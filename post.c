#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "def.h"
//#include "f2c.h"
//#include "clapack.h"

int  BinToDat_Cell
//=======================================================
//  1. (5th Sep. 2022)
//  CONVERT BIN DATA TO DAT DATA
//
//
(
    int       icnt,
    domain    *cdo,
    lattice   *ltc,
    cell      *cel,
    Wall      *wall,
    char      *fpoint,
    char      *fname
)
//-------------------------------------------------------
{
// Output cell file ---
    FILE      *fp;
    char      filename[256];
    int       ic, j, e,
              file_check = 0, c = 3, b = 5;

    sprintf(filename,"%s%s/result/cell/cell%04d.bin",fpoint,fname,icnt);
    file_check = read_cell_bin(icnt,cel,wall,filename);
    if (file_check == -1) return 1; 

    sprintf(filename,"./result/%s/cell%04d.vtk",fname,icnt);
    fp = fopen(filename,"w"); if(fp == NULL) error(1);

    fprintf(fp,"# vtk DataFile Version 2.0 \n");
    fprintf(fp,"# vtk output \n");
    fprintf(fp,"ASCII \n");
    fprintf(fp,"DATASET UNSTRUCTURED_GRID \n");
    fprintf(fp,"POINTS %04d float\n",cel->vertex*cel->n);

    for (ic = 0; ic < cel->n; ic++) {
      for (j = 0; j < cel->vertex; j++) {
        fprintf(fp,"%f %f %f \n",
        cel->x[ic*cel->vertex + j],
        cel->y[ic*cel->vertex + j],
        cel->z[ic*cel->vertex + j]);
      }
    }

    fprintf(fp,"CELLS %04d %04d\n ",cel->n*cel->element,cel->n*cel->element*4);
    for(ic = 0; ic < cel->n; ic++) {
      for(e = 0;e < cel->element;e++){
        fprintf(fp,"%d %10d %10d %10d\n",
        c,
        cel->ele[ic*cel->element*TRI + e*TRI  ],
        cel->ele[ic*cel->element*TRI + e*TRI+1],
        cel->ele[ic*cel->element*TRI + e*TRI+2]);
      }
    }

    fprintf(fp,"CELL_TYPES %04d\n",cel->n*cel->element);
    for (ic = 0; ic < cel->n*cel->element; ic++) fprintf(fp,"%d \n",b);
    
    fclose(fp);

    return 0;
}

int Velocity
//============================================================================== 
(
// 2. (5th Sep. 2022)
// CENTROID VELOCITY 
//
//
    int       icnt,
    process   *prc,
    domain    *cdo,
    lattice   *ltc,
    cell      *cel,
    Wall      *wall,
    char      *fpoint,
    char      *fname
)
//------------------------------------------------------------------------------   
{
    FILE          *fp;
    char          filename[256];
    int           file_check = 0;
    double        a, r, rt, DT, T;
    
    a  = DX*ltc->CKL;                            // Radius of a RBC               [m]
    // SHEAR FLOW
    ///r  = ltc->Ca*cel->Gs/(MU*a);                // Shear rate                    [1/s]
    // CHANNEL FLOW
    r  = 2.0*ltc->Ca*cel->Gs/(MU*a);             // Shear rate                    [1/s]
    rt = cdo->dt*(double)icnt*(double)prc->Fout; // Non-dimensional time          [-]
    DT = cdo->dt/r;                              // Dimensional delta time        [s]
    T  = DT*(double)icnt*(double)prc->Fout;      // Dimensional time              [s]

    // Reading cell file
    sprintf(filename,"%s%s/result/cell/cell%04d.bin",fpoint,fname,icnt);
    file_check = read_cell_bin(icnt,cel,wall,filename);
    if (file_check == -1) return 1; 

    // Set variables
    if (icnt == 0) {
      cel->rc = (double *)malloc(sizeof(double)*cel->n*3); if (cel->rc == NULL) error(1);
      cel->vc = (double *)malloc(sizeof(double)*cel->n*3); if (cel->vc == NULL) error(1);
      cel->gp = (double *)malloc(sizeof(double)*cel->n*3); if (cel->gp == NULL) error(1);
    } // --- End of initializing ---

    // Calculating the centroid & velocity of each RBC
    Centroid_Capsules(cel,cdo);

    //Radial position
    //rc_xy = sqrt(cel->rc[0]*cel->rc[0] + cel->rc[1]*cel->rc[1]);

    // File output
    sprintf(filename,"./centroid_velocity/%s.dat",fname);
    if (icnt == 0) fp = fopen(filename,"w");
    else           fp = fopen(filename,"a");
    if (fp == NULL) error(2);
    // rt    : non-dilensional time[-] 
    // T     : dilensional time[s] 
    // rc[2] : capsule centoid coordinate along flow-direction (z-direction)[-] 
    // vc[2] : capsule centoid velocity along flow-direction (z-direction)[-] 
    // sin() : input velocity at time [-]
    fprintf(fp,"%15.9e %15.9e %15.9e %15.9e %15.9e\n",rt,T,cel->rc[2],cel->vc[2],sin(M_PI*rt/cdo->freq));
    fclose(fp);

    return 0;
}
