#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "def.h"


int  read_parameter
//=======================================================
//
//  READ PARAMETER FILE
//
//
//
(
    process   *prc,
    domain    *cdo,
    lattice   *ltc,
    cell      *cel,
    Wall      *wall,
    char      *filename
)
//-------------------------------------------------------
{
    FILE      *fp;

    fp = fopen(filename,"r"); if (fp == NULL) error(2);
    printf("Read file is %s\n",filename);

    fscanf(fp,"%lf %lf %lf %lf %lf    \n",&ltc->Ca,&ltc->Re,&ltc->KL,&ltc->Df,&ltc->Cr);
    fscanf(fp,"%lf %lf %lf            \n",&cdo->D,&cdo->U,&cdo->L);
    fscanf(fp,"%lf %lf %lf            \n",&cdo->LX,&cdo->LY,&cdo->LZ);
    fscanf(fp,"%d %d %d %d %d         \n",&cel->n,&cel->n_irbc,&cel->n_w,&cel->n_c,&cel->n_p);
    fscanf(fp,"%lf %lf %lf %lf %lf    \n",&ltc->CKL,&cel->CR,&cel->CW,&cel->CC,&cel->CP);
    fscanf(fp,"%lf %lf %lf %lf        \n",&cel->Gs,&cel->RGw,&cel->RGc,&cel->RGp);
    fscanf(fp,"%lf %lf %lf            \n",&cdo->dt,&cdo->dx,&cdo->tau);
    fscanf(fp,"%d %d %d %d %d         \n",&cdo->nx,&cdo->ny,&cdo->nz,&cdo->n,&wall->n_wall);
    fscanf(fp,"%d %d %d %d %d         \n",&prc->IterLimit,&prc->Fout,&prc->STout,&prc->FlowS,&prc->Num_Flag);
    fscanf(fp,"%lf %lf                \n",&cdo->lam1,&cdo->lam2);
    fscanf(fp,"%lf %d                 \n",&cdo->freq,&cdo->flag); // flag = 0 -> cos || flag = 1 -> sin
    fclose(fp);

    return 0;
}

int  read_wall_bin
//=======================================================
//
//  READ BINARY DATAFILE OF WALL
//
//
//
(
    Wall   *wall,
    char   *filename
)
//-------------------------------------------------------
{
    FILE       *fp;

    wall->x = (double*)malloc(sizeof(double)*wall->n_wall*3);

    fp = fopen(filename,"rb"); if(fp == NULL) return -1;
    fread(wall->x,sizeof(double),wall->n_wall*3,fp);
    fclose(fp);
    printf("Read file is %s\n",filename);

    return 0;
}


int  read_fluid_bin
//=======================================================
//
//  READ BINARY DATAFILE OF CELL
//
//
//
(
    int       icnt,
    process   *prc,
    domain    *cdo,
    lattice   *ltc,
    cell      *cel,
    char      *filename
)
//-------------------------------------------------------
{
    FILE       *fp;
    static int file_start = 0;

    cdo->comp = 0;

    fp = fopen(filename,"rb");
    if (fp == NULL) return -1;
    else file_start++;

    fread(&prc->iter,sizeof(int),1,fp);
  /*  
    fread(&cdo->dt,  sizeof(double),1,fp);
    fread(&cdo->dx,  sizeof(double),1,fp);
    fread(&cdo->nx,  sizeof(int   ),1,fp);
    fread(&cdo->ny,  sizeof(int   ),1,fp);
    fread(&cdo->nz,  sizeof(int   ),1,fp);
    fread(&cdo->n,   sizeof(int   ),1,fp);
  */  

    if (file_start == 1) {
      ltc->d   = (double *)malloc(sizeof(double)*cdo->n);
      ltc->u   = (double *)malloc(sizeof(double)*cdo->n);
      ltc->v   = (double *)malloc(sizeof(double)*cdo->n);
      ltc->w   = (double *)malloc(sizeof(double)*cdo->n);
      ltc->vf  = (double *)malloc(sizeof(double)*cdo->n);
      ltc->vf2 = (double *)malloc(sizeof(double)*cdo->n);
      ltc->bc  = (int    *)malloc(sizeof(int   )*cdo->n);
    }

    fread(ltc->d,  sizeof(double),cdo->n,fp);
    fread(ltc->u,  sizeof(double),cdo->n,fp);
    fread(ltc->v,  sizeof(double),cdo->n,fp);
    fread(ltc->w,  sizeof(double),cdo->n,fp);
    fread(ltc->vf, sizeof(double),cdo->n,fp);
    fread(ltc->vf2,sizeof(double),cdo->n,fp);
    fread(ltc->bc, sizeof(int   ),cdo->n,fp);
    fclose(fp);

    printf("Read file is %s\n",filename);

    return 0;
}

int  read_cell_bin
//=======================================================
//
//  READ BINARY DATAFILE OF CELL
//  2012/06/05
//
//
//
(
    int       icnt,
    cell      *cel,
    Wall      *wall,
    char      *filename
)
//-------------------------------------------------------
{
    FILE       *fp;
    static int file_start = 0;

    fp = fopen(filename,"rb");
    if (fp == NULL) return -1;
    else            file_start++;

// Read data about cells
   if (cel->n > 0) {
     //fread(&cel->flag,   sizeof(int),1,fp); // For platelets or MPs at D = 8 um
     fread(&cel->n,      sizeof(int),1,fp);
     fread(&cel->vertex, sizeof(int),1,fp);
     fread(&cel->element,sizeof(int),1,fp);

     if (file_start == 1) {
       cel->x   = (double *)malloc(sizeof(double)*cel->n*cel->vertex   );
       cel->y   = (double *)malloc(sizeof(double)*cel->n*cel->vertex   );
       cel->z   = (double *)malloc(sizeof(double)*cel->n*cel->vertex   );
       cel->u   = (double *)malloc(sizeof(double)*cel->n*cel->vertex   );
       cel->v   = (double *)malloc(sizeof(double)*cel->n*cel->vertex   );
       cel->w   = (double *)malloc(sizeof(double)*cel->n*cel->vertex   );
       cel->q   = (double *)malloc(sizeof(double)*cel->n*cel->vertex*3 );
       cel->t   = (double *)malloc(sizeof(double)*cel->n*cel->element*2);
       cel->ele = (int    *)malloc(sizeof(int   )*cel->n*cel->element*3);
       #if (REPULSION == true)
       cel->r_ep   = (double *)malloc(sizeof(double)*cel->n*cel->vertex  );
       cel->is_rep = (bool   *)malloc(sizeof(bool  )*cel->n*cel->vertex  );
       cel->f_rep  = (double *)malloc(sizeof(double)*cel->n*cel->vertex*3);
       #endif
     }

     fread(cel->x,   sizeof(double),cel->n*cel->vertex,   fp);
     fread(cel->y,   sizeof(double),cel->n*cel->vertex,   fp);
     fread(cel->z,   sizeof(double),cel->n*cel->vertex,   fp);
     fread(cel->u,   sizeof(double),cel->n*cel->vertex,   fp);
     fread(cel->v,   sizeof(double),cel->n*cel->vertex,   fp);
     fread(cel->w,   sizeof(double),cel->n*cel->vertex,   fp);
     fread(cel->q,   sizeof(double),cel->n*cel->vertex*3, fp);
     fread(cel->t,   sizeof(double),cel->n*cel->element*2,fp);
     fread(cel->ele, sizeof(int   ),cel->n*cel->element*3,fp);
     #if (REPULSION == true)
     fread(cel->r_ep,   sizeof(double),cel->n*cel->vertex,   fp);
     fread(cel->is_rep, sizeof(bool  ),cel->n*cel->vertex,   fp);
     fread(cel->f_rep,  sizeof(double),cel->n*cel->vertex*3, fp);
     #endif

   } else { cel->vertex = 0; cel->element = 0; }

// Read data about nucleus
   if (cel->n_w > 0) {
     fread(&cel->n_w,      sizeof(int),1,fp);
     fread(&cel->vertex_w, sizeof(int),1,fp);
     fread(&cel->element_w,sizeof(int),1,fp);

     if (file_start == 1) {
       cel->x_w   = (double *)malloc(sizeof(double)*cel->n_w*cel->vertex_w   );
       cel->y_w   = (double *)malloc(sizeof(double)*cel->n_w*cel->vertex_w   );
       cel->z_w   = (double *)malloc(sizeof(double)*cel->n_w*cel->vertex_w   );
       cel->u_w   = (double *)malloc(sizeof(double)*cel->n_w*cel->vertex_w   );
       cel->v_w   = (double *)malloc(sizeof(double)*cel->n_w*cel->vertex_w   );
       cel->w_w   = (double *)malloc(sizeof(double)*cel->n_w*cel->vertex_w   );
       cel->q_w   = (double *)malloc(sizeof(double)*cel->n_w*cel->vertex_w*3 );
       cel->t_w   = (double *)malloc(sizeof(double)*cel->n_w*cel->element_w*2);
       cel->ele_w = (int    *)malloc(sizeof(int   )*cel->n_w*cel->element_w*3);
       #if (REPULSION == true)
       cel->r_ep_w   = (double *)malloc(sizeof(double)*cel->n_w*cel->vertex_w  );
       cel->is_rep_w = (bool   *)malloc(sizeof(bool  )*cel->n_w*cel->vertex_w  );
       cel->f_rep_w  = (double *)malloc(sizeof(double)*cel->n_w*cel->vertex_w*3);
       #endif
     }

     fread(cel->x_w,   sizeof(double),cel->n_w*cel->vertex_w,   fp);
     fread(cel->y_w,   sizeof(double),cel->n_w*cel->vertex_w,   fp);
     fread(cel->z_w,   sizeof(double),cel->n_w*cel->vertex_w,   fp);
     fread(cel->u_w,   sizeof(double),cel->n_w*cel->vertex_w,   fp);
     fread(cel->v_w,   sizeof(double),cel->n_w*cel->vertex_w,   fp);
     fread(cel->w_w,   sizeof(double),cel->n_w*cel->vertex_w,   fp);
     fread(cel->q_w,   sizeof(double),cel->n_w*cel->vertex_w*3, fp);
     fread(cel->t_w,   sizeof(double),cel->n_w*cel->element_w*2,fp);
     fread(cel->ele_w, sizeof(int   ),cel->n_w*cel->element_w*3,fp);
     #if (REPULSION == true)
     fread(cel->r_ep_w,   sizeof(double),cel->n_w*cel->vertex_w,   fp);
     fread(cel->is_rep_w, sizeof(bool  ),cel->n_w*cel->vertex_w,   fp);
     fread(cel->f_rep_w,  sizeof(double),cel->n_w*cel->vertex_w*3, fp);
     #endif

   } else { cel->vertex_w = 0; cel->element_w = 0; }
   
// Read data about adhesion
   if (file_start == 1) {
     if (cel->n == 1 && cel->n_w == 0) {
       cel->adnum = cel->n;
       cel->advertex = cel->vertex;
     } else {
       cel->adnum = cel->n_w;
       cel->advertex = cel->vertex_w;
     }

     #if (CYTOADHESION == true)
     cel->f_adh            = (double *)malloc(sizeof(double)*cel->adnum*cel->advertex*3                );
     cel->is_knob          = (bool   *)malloc(sizeof(bool  )*cel->adnum*cel->advertex                  );
     cel->attachment_point = (int    *)malloc(sizeof(int   )*cel->adnum*cel->advertex*cel->nspring_knob);
     wall->is_attached     = (int    *)malloc(sizeof(int   )*wall->n_wall                              );
     #endif
     #if (POLYMERIZATION == true)
     cel->point = (bool   *)malloc(sizeof(bool  )*cel->adnum*cel->advertex  );
     cel->vg    = (double *)malloc(sizeof(double)*cel->adnum*3              );
     cel->vec   = (double *)malloc(sizeof(double)*cel->adnum*cel->advertex*3);
//     cel->ddlp  = (bool   *)malloc(sizeof(bool  )*cel->adnum*cel->advertex  );
     cel->ddlx  = (double *)malloc(sizeof(double)*cel->adnum*cel->advertex*3);
     #endif
   }

   #if (CYTOADHESION == true)
   fread(cel->f_adh,            sizeof(double),cel->adnum*cel->advertex*3,                 fp);
   fread(cel->is_knob,          sizeof(bool  ),cel->adnum*cel->advertex,                   fp);
   fread(cel->attachment_point, sizeof(int   ),cel->adnum*cel->advertex*cel->nspring_knob, fp);
   fread(wall->is_attached,     sizeof(int   ),wall->n_wall,                               fp);
   #endif
   #if (POLYMERIZATION == true)
   fread(cel->point, sizeof(bool  ),cel->adnum*cel->advertex,   fp);
   fread(cel->vg,    sizeof(double),cel->adnum*3,               fp);
   fread(cel->vec,   sizeof(double),cel->adnum*cel->advertex*3, fp);
//   fread(cel->ddlp,  sizeof(bool  ),cel->adnum*cel->advertex,   fp);
   fread(cel->ddlx,  sizeof(double),cel->adnum*cel->advertex*3, fp);
   #endif

   fclose(fp);
   printf("Read file is %s\n",filename);

   cel->vertex_p = 0; cel->element_p = 0;
   cel->vertex_c = 0; cel->element_c = 0;

   return 0;
}
