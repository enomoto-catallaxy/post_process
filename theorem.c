#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "def.h"

void  Add_Nlink
//==========================================================
//
//  SET LINKING DATA
//
//
(
    int       na,            // first node
    int       nb,            // second node
    int       *nln,          // number of linking node
    int       *ln            // linking node id number
)
//----------------------------------------------------------
{
    int       flag = 0, i;

    for(i = 0; i < nln[na]; i++){
        if(nb == ln[7*na + i]){
           flag = 1;
           continue;
        }
    }
    if(flag == 0){
       ln[7*na + nln[na]] = nb;
       nln[na]++;
    }

    return;
}


void  Add_Elink
//==========================================================
//
//  SET LINKING DATA
//
//
(
    int       ni,            // id number of node
    int       ei,            // id number of element
    int       *nle,          // number of linking element
    int       *le            // linking element id number
)
//----------------------------------------------------------
{
    le[6*ni + nle[ni]] = ei;
    nle[ni]++;

    return;
}


int  Checklink
//==========================================================
//
//  CHECK WHETHER LINKED OR NOT
//
//
(
    int       na,            // first node
    int       nb,            // second node
    int       *nln,          // number of linking node
    int       *ln            // linking node id number
)
//----------------------------------------------------------
{
    int       i;
    int       flag = 0;

    for(i = 0; i < nln[na]; i++){
        if(nb == ln[7*na + i]) flag = 1;
    }

    return flag;
}


void  Sort_Link_Spring
//==========================================================
//
//  SORT LINKING DATA
//
//
(
    double    *x,            // position coordinate at n step
    int       *nln,          // number of linking node
    int       *ln,           // linking node id number
    int       vertex,
    int       ic
)
//----------------------------------------------------------
{
    int       j, k, l, n, na, nb, jid;
    int       lns[7];
    double    va[3], vb[3], c[3], inpro;
    double    gx = 0.0, gy = 0.0, gz = 0.0;


    for(j = 0;j < vertex;j++){
        jid = ic*vertex + j;
        gx += x[3*jid + 0];
        gy += x[3*jid + 1];
        gz += x[3*jid + 2];
    }
    gx /= (double)vertex;
    gy /= (double)vertex;
    gz /= (double)vertex;

    for(j = 0;j < vertex;j++){
        jid = ic*vertex + j;
        n = 0;
        lns[0] = ln[7*jid + 0];

        for(k = 1;k < nln[jid];k++){
            na = lns[k - 1];

            va[0] = x[3*na + 0] - x[3*jid + 0];
            va[1] = x[3*na + 1] - x[3*jid + 1];
            va[2] = x[3*na + 2] - x[3*jid + 2];

            for(l = 0;l < nln[jid];l++){
                nb = ln[7*jid + l];
                if(na == nb) continue;

                if(Checklink(na,nb,nln,ln) == 1){
                    vb[0] = x[3*nb + 0] - x[3*jid + 0];
                    vb[1] = x[3*nb + 1] - x[3*jid + 1];
                    vb[2] = x[3*nb + 2] - x[3*jid + 2];

                    c[0] = va[1]*vb[2] - va[2]*vb[1];
                    c[1] = va[2]*vb[0] - va[0]*vb[2];
                    c[2] = va[0]*vb[1] - va[1]*vb[0];

                    inpro = c[0]*(x[3*jid + 0] - gx)
                          + c[1]*(x[3*jid + 1] - gy)
                          + c[2]*(x[3*jid + 2] - gz);

                    if(inpro > 0){
                        lns[k] = nb;
                        n++;
                        continue;
                    }
                }
            }
        }

        if(n + 1 != nln[jid]) error(4);
        else                  for(k = 0;k < nln[jid];k++) ln[7*jid + k] = lns[k];
    }

    return;
}


void  Sort_Link_Order
//==========================================================
//
//  SORT LINKING DATA
//
//
(
    int       *nln,          // number of linking node
    int       *ln,           // linking node id number
    int       vertex,
    int       ic
)
//----------------------------------------------------------
{
    int       j, k, l, jid;
    int       temp;


    for(j = 0; j < vertex; j++){
        jid = ic*vertex + j;
        ln[7*jid + nln[jid]] = jid;

        for(k = 0; k < nln[jid] + 1; k++){
            for(l = nln[jid]; l > k; l--){
                if(ln[7*jid + l - 1] > ln[7*jid + l]){
                    temp = ln[7*jid + l];
                    ln[7*jid + l] = ln[7*jid + l - 1];
                    ln[7*jid + l - 1] = temp;
                }
            }
        }
    }

    return;
}

void  Group_Node_Link
//=======================================================
//  12th Jun 2014
//  GROUPING ADJECENT NODES
//
//
//
(
    int       nc,
    int       nv,
    int       ne,
    int       *e,
    double    *x,
    double    *y,
    double    *z,
    cell      *cel,
    domain    *cdo
)
//-------------------------------------------------------
{
    int        i, j, k, ea, eb, ec, jid, eid, n_edge, division, n_th, checker;
    static int flag = 0;

    if (flag == 0) flag ++; 
    else return;

// Parameters ---
    n_edge   = nv + ne - 2;
    division = (int)ceil((double)nv*(double)nc/(double)DIV);
    n_th     = division*DIV;

// Allocate variables ---
    cel->xyz = (double*)malloc(sizeof(double)*n_th*3); if (cel->xyz == NULL) error(0);
    cel->nln = (int   *)malloc(sizeof(int   )*n_th  ); if (cel->nln == NULL) error(0);
    cel->ln  = (int   *)malloc(sizeof(int   )*n_th*7); if (cel->ln  == NULL) error(0);
    cel->lo  = (int   *)malloc(sizeof(int   )*n_th*7); if (cel->lo  == NULL) error(0);
    cel->nle = (int   *)malloc(sizeof(int   )*n_th  ); if (cel->nle == NULL) error(0);
    cel->le  = (int   *)malloc(sizeof(int   )*n_th*6); if (cel->le  == NULL) error(0);

    for (j = 0; j < n_th*3; j++) cel->xyz[j] = 0.0e+00;
    for (j = 0; j < n_th;   j++) cel->nln[j] = 0;
    for (j = 0; j < n_th*7; j++) cel->ln[j]  = 0;
    for (j = 0; j < n_th*7; j++) cel->lo[j]  = 0;
    for (j = 0; j < n_th;   j++) cel->nle[j] = 0;
    for (j = 0; j < n_th*6; j++) cel->le[j]  = 0;

    for (i = 0; i < nc; i++) {
      for (j = 0; j < nv; j++) {
        jid = i*nv + j;
        cel->xyz[jid*3 + 0] = x[jid];
        cel->xyz[jid*3 + 1] = y[jid];
        cel->xyz[jid*3 + 2] = z[jid];
      }
    }

// Linking information ---
    for (i = 0; i < nc; i++) {
      for (j = 0; j < ne; j++) {
        eid = i*ne + j;
        ea  = e[3*eid + 0];
        eb  = e[3*eid + 1];
        ec  = e[3*eid + 2];

        Add_Nlink(ea,eb,cel->nln,cel->ln);
        Add_Nlink(ea,ec,cel->nln,cel->ln);
        Add_Nlink(eb,ec,cel->nln,cel->ln);
        Add_Nlink(eb,ea,cel->nln,cel->ln);
        Add_Nlink(ec,ea,cel->nln,cel->ln);
        Add_Nlink(ec,eb,cel->nln,cel->ln);

        Add_Elink(ea,eid,cel->nle,cel->le);
        Add_Elink(eb,eid,cel->nle,cel->le);
        Add_Elink(ec,eid,cel->nle,cel->le);
      }

      checker = 0;
      for (j = 0; j < nv; j++) {
        jid = i*nv + j;
        checker += cel->nln[jid];
      }
      if (checker != 2*n_edge) error(4);

      for (j = 0; j < nv; j++) {
        jid = i*nv + j;
        for (k = 0; k < cel->nln[jid]; k++) {
          cel->lo[7*jid + k] = cel->ln[7*jid + k];
        }
      }
      Sort_Link_Spring(cel->xyz,cel->nln,cel->ln,nv,i);
      Sort_Link_Order(cel->nln,cel->lo,nv,i);
    }

    return;
}

void  Storage_Periodic_Process
//=======================================================
//  27th May 2014
//  PROCUDURE FOR PERIODIC
//  
//
//
(
    int       nc,
    int       nv,
    int       ne,
    int       *e,
    double    *x,
    double    *y,
    double    *z,
    double    *u,
    double    *v,
    double    *w,
    cell      *cel,
    domain    *cdo
)
//-------------------------------------------------------
{
    int       i, j, k, j1, j2, flag = 0;
    double    zz;

    // Storage variables for each capsule
    if (cel->n_w > 0) {
      for (i = 0; i < nc*nv; i++) {
        x[i] = cel->x_w[i]; y[i] = cel->y_w[i]; z[i] = cel->z_w[i];
        u[i] = cel->u_w[i]; v[i] = cel->v_w[i]; w[i] = cel->w_w[i];
      }
      for (j = 0; j < nc*ne; j++) {
        e[TRI*j+0] = cel->ele_w[TRI*j+0];
        e[TRI*j+1] = cel->ele_w[TRI*j+1];
        e[TRI*j+2] = cel->ele_w[TRI*j+2];
      }

    } else if (cel->n_p > 0) {
      for (i = 0; i < nc*nv; i++) {
        x[i] = cel->x_p[i]; y[i] = cel->y_p[i]; z[i] = cel->z_p[i];
        u[i] = cel->u_p[i]; v[i] = cel->v_p[i]; w[i] = cel->w_p[i];
      }
      for (j = 0; j < nc*ne; j++) {
        e[TRI*j+0] = cel->ele_p[TRI*j+0];
        e[TRI*j+1] = cel->ele_p[TRI*j+1];
        e[TRI*j+2] = cel->ele_p[TRI*j+2];
      }

    } else if (cel->n_c > 0) {
      for (i = 0; i < nc*nv; i++) {
        x[i] = cel->x_c[i]; y[i] = cel->y_c[i]; z[i] = cel->z_c[i];
        u[i] = cel->u_c[i]; v[i] = cel->v_c[i]; w[i] = cel->w_c[i];
      }
      for (j = 0; j < nc*ne; j++) {
        e[TRI*j+0] = cel->ele_c[TRI*j+0];
        e[TRI*j+1] = cel->ele_c[TRI*j+1];
        e[TRI*j+2] = cel->ele_c[TRI*j+2];
      }
    }
//    return;

    // Calculation : Centroid of a capsule (from vertex)
    for (i = 0; i < nc; i++) {
      flag = 0;
      for (j = 0; j < ne; j++) {
        for (k = 0; k < TRI; k++) {
          j1 = e[i*ne*TRI + j*TRI+k];
          j2 = e[i*ne*TRI + j*TRI+(k+1)%TRI];
          zz = z[j2] - z[j1];
          if      (zz >  cdo->LZ*0.5){ z[j1] += cdo->LZ; flag++; }
          else if (zz < -cdo->LZ*0.5){ z[j2] += cdo->LZ; flag++; }
        }
      }
      if (flag > 0) {
        for (j = 0; j < nv; j++) {
            if(z[i*nv+j] < 0.0e+10) z[i*nv+j] += cdo->LZ;
        }
      }
    }

    return;
}

void  Centroid_Capsules
//=======================================================
//
//  CALCULATE THE CENTROID OF RBCS
//  2016.05.04
//
//
(
    cell      *cel,
    domain    *cdo
)
//-------------------------------------------------------
{
    int       e, i, j, k, j1, j2,  ea, eb, ec;
    bool      xflag, zflag;
    double    x0,  y0,  z0,  x1,  y1,  z1,  x2,  y2,  z2,
              x01, y01, z01, x02, y02, z02, vex, vey, vez, 
              uvx, uvy, uvz, rx, ry, rz, uu, uv, uw, ds, xx, zz,
              volume, jr[3];

    for (i = 0; i < cel->n*3; i++) {
      cel->rc[i] = 0.0e+10;
      cel->vc[i] = 0.0e+10;
      cel->gp[i] = 0.0e+10;
    }

// Procedure for periodic boundary
    for (i = 0; i < cel->n; i++) {
      xflag = false;
      zflag = false;

      for (e = 0; e < cel->element; e++) {
        for (k = 0; k < TRI; k++) {
          j1 = cel->ele[i*cel->element*TRI + e*TRI+k];
          j2 = cel->ele[i*cel->element*TRI + e*TRI+(k+1)%TRI];
          xx = cel->x[j2] - cel->x[j1];
          zz = cel->z[j2] - cel->z[j1];
          if      (xx >  cdo->LX*0.5) { cel->x[j1] += cdo->LX; xflag = true; }
          else if (xx < -cdo->LX*0.5) { cel->x[j2] += cdo->LX; xflag = true; }
          if      (zz >  cdo->LZ*0.5) { cel->z[j1] += cdo->LZ; zflag = true; }
          else if (zz < -cdo->LZ*0.5) { cel->z[j2] += cdo->LZ; zflag = true; }
        }
      }

      if (xflag == true && zflag == true) {
        for (j = 0; j < cel->vertex; j++) {
          if (cel->x[i*cel->vertex + j] < 0.0e+10) cel->x[i*cel->vertex + j] += cdo->LX;
          if (cel->z[i*cel->vertex + j] < 0.0e+10) cel->z[i*cel->vertex + j] += cdo->LZ;
          cel->gp[i*TRI+0] += cel->x[i*cel->vertex + j];
          cel->gp[i*TRI+1] += cel->y[i*cel->vertex + j];
          cel->gp[i*TRI+2] += cel->z[i*cel->vertex + j];
        }
      } else if (xflag == true && zflag == false) {
        for (j = 0; j < cel->vertex; j++) {
          if (cel->x[i*cel->vertex + j] < 0.0e+10) cel->x[i*cel->vertex + j] += cdo->LX;
          cel->gp[i*TRI+0] += cel->x[i*cel->vertex + j];
          cel->gp[i*TRI+1] += cel->y[i*cel->vertex + j];
          cel->gp[i*TRI+2] += cel->z[i*cel->vertex + j];
        }
      } else if (xflag == false && zflag == true) {
        for (j = 0; j < cel->vertex; j++) {
          if (cel->z[i*cel->vertex + j] < 0.0e+10) cel->z[i*cel->vertex + j] += cdo->LZ;
          cel->gp[i*TRI+0] += cel->x[i*cel->vertex + j];
          cel->gp[i*TRI+1] += cel->y[i*cel->vertex + j];
          cel->gp[i*TRI+2] += cel->z[i*cel->vertex + j];
        }
      } else {
        for (j = 0; j < cel->vertex; j++) { 
          cel->gp[i*TRI+0] += cel->x[i*cel->vertex + j];
          cel->gp[i*TRI+1] += cel->y[i*cel->vertex + j];
          cel->gp[i*TRI+2] += cel->z[i*cel->vertex + j];
        }
      }

      cel->gp[i*TRI+0] /= (double)cel->vertex;
      cel->gp[i*TRI+1] /= (double)cel->vertex;
      cel->gp[i*TRI+2] /= (double)cel->vertex;
    }

// Calculation : Centoid position
    for (i = 0; i < cel->n; i++) {
      volume = 0.0e+10;

      for (e = 0; e < cel->element; e++) {
        ea = cel->ele[i*cel->element*TRI + e*TRI + 0];
        eb = cel->ele[i*cel->element*TRI + e*TRI + 1];
        ec = cel->ele[i*cel->element*TRI + e*TRI + 2];

        x0 = cel->x[ea]; x1 = cel->x[eb]; x2 = cel->x[ec];
        y0 = cel->y[ea]; y1 = cel->y[eb]; y2 = cel->y[ec];
        z0 = cel->z[ea]; z1 = cel->z[eb]; z2 = cel->z[ec];

        x01 = x1 - x0;   x02 = x2 - x0;
        y01 = y1 - y0;   y02 = y2 - y0;
        z01 = z1 - z0;   z02 = z2 - z0;

        vex = y01*z02 - z01*y02;
        vey = z01*x02 - x01*z02;
        vez = x01*y02 - y01*x02;
       
        ds = 0.5*sqrt(vex*vex + vey*vey + vez*vez);

        // Normal unit vector of a triangle-element
        uvx = vex/sqrt(vex*vex + vey*vey + vez*vez);
        uvy = vey/sqrt(vex*vex + vey*vey + vez*vez);
        uvz = vez/sqrt(vex*vex + vey*vey + vez*vez);

        // Centroid vector of a triangle-element
        rx = (x0 + x1 + x2)/3.0;
        ry = (y0 + y1 + y2)/3.0;
        rz = (z0 + z1 + z2)/3.0;

        volume += (uvx*rx + uvy*ry + uvz*rz)*ds;

        cel->rc[i*TRI+0] += (rx*rx + ry*ry + rz*rz)*uvx*ds;
        cel->rc[i*TRI+1] += (rx*rx + ry*ry + rz*rz)*uvy*ds;
        cel->rc[i*TRI+2] += (rx*rx + ry*ry + rz*rz)*uvz*ds;
      }
      // Volume of a capsule
      volume /= 3.0;

      // Centroid of a casule
      cel->rc[i*TRI+0] /= (2.0*volume);
      cel->rc[i*TRI+1] /= (2.0*volume);
      cel->rc[i*TRI+2] /= (2.0*volume);
    }

// Calculation : Volume average velocity
    for (i = 0; i < cel->n; i++) {
      volume = 0.0e+10;

      for (e = 0; e < cel->element; e++) {
        ea = cel->ele[i*cel->element*TRI + e*TRI + 0];
        eb = cel->ele[i*cel->element*TRI + e*TRI + 1];
        ec = cel->ele[i*cel->element*TRI + e*TRI + 2];

        x0 = cel->x[ea]; x1 = cel->x[eb]; x2 = cel->x[ec];
        y0 = cel->y[ea]; y1 = cel->y[eb]; y2 = cel->y[ec];
        z0 = cel->z[ea]; z1 = cel->z[eb]; z2 = cel->z[ec];

        x01 = x1 - x0;   x02 = x2 - x0;
        y01 = y1 - y0;   y02 = y2 - y0;
        z01 = z1 - z0;   z02 = z2 - z0;

        vex = y01*z02 - z01*y02;
        vey = z01*x02 - x01*z02;
        vez = x01*y02 - y01*x02;
       
        ds = 0.5*sqrt(vex*vex + vey*vey + vez*vez);

        // Normal unit vector of a triangle-element
        uvx = vex/sqrt(vex*vex + vey*vey + vez*vez);
        uvy = vey/sqrt(vex*vex + vey*vey + vez*vez);
        uvz = vez/sqrt(vex*vex + vey*vey + vez*vez);

        // Centroid vector of a triangle-element
        rx = (x0 + x1 + x2)/3.0;
        ry = (y0 + y1 + y2)/3.0;
        rz = (z0 + z1 + z2)/3.0;

        jr[0] = rx - cel->gp[i*TRI+0];
        jr[1] = ry - cel->gp[i*TRI+1];
        jr[2] = rz - cel->gp[i*TRI+2];

        volume += (uvx*jr[0] + uvy*jr[1] + uvz*jr[2])*ds;

        // Centroid velocity of a triangle-element
        uu = (cel->u[ea] + cel->u[eb] + cel->u[ec])/3.0;
        uv = (cel->v[ea] + cel->v[eb] + cel->v[ec])/3.0;
        uw = (cel->w[ea] + cel->w[eb] + cel->w[ec])/3.0;

        cel->vc[i*TRI+0] += (uvx*uu + uvy*uv +uvz*uw)*jr[0]*ds;
        cel->vc[i*TRI+1] += (uvx*uu + uvy*uv +uvz*uw)*jr[1]*ds;
        cel->vc[i*TRI+2] += (uvx*uu + uvy*uv +uvz*uw)*jr[2]*ds;
      }
      // Volume of a capsule
      volume /= 3.0;

      // Colume-average velocity of a capsule
      cel->vc[i*TRI+0] /= volume;
      cel->vc[i*TRI+1] /= volume;
      cel->vc[i*TRI+2] /= volume;
    }

    return;
}

void Bending
//============================================================================== 
(
// 23th Jun 2018
// BENDING RESISTANCE 
//
//
    int       icnt,
    double    *fx,
    double    *fy,
    double    *fz,
    int       *ec,
    lattice   *ltc,
    cell      *cel,
    domain    *cdo
)
//------------------------------------------------------------------------------   
{
    int       e, i, ic, k, jid, eid, e1, e2,
              i1, i2, i3, j1, j2, j3,
	      k1 = 0, k2 = 0, k3 = 0, k4 = 0, ecn;
    double    xx, zz, a, br, //a, r, rt, DT, T, S = 0.0,
	      //ave_q = 0.0, total_q = 0.0,
	      // Parameter for bending
	      cost0, sint0, cost, cost2, sint, b11, b12, b22, beta, eps = 1.0e-6,
              a1x, a2x, a3x, a4x, a21x, a31x, a34x, a24x, a23x, t1x, t2x, xix, zex, fxb,
              a1y, a2y, a3y, a4y, a21y, a31y, a34y, a24y, a23y, t1y, t2y, xiy, zey, fyb,
              a1z, a2z, a3z, a4z, a21z, a31z, a34z, a24z, a23z, t1z, t2z, xiz, zez, fzb;
    
    a  = DX*ltc->CKL;      // Radius of a RBC  [m]
    br = KC/(cel->Gs*a*a); // Bending regidity [-]

    // Set initial angles
    if (icnt == 0) {
      // Grouping linked node (5 or 6)
      Group_Node_Link (cel->n,cel->vertex,cel->element,cel->ele,cel->x,cel->y,cel->z,cel,cdo);

      // Element connectivity ---
      for(ic = 0; ic < cel->n; ic++) {
        for(e = 0; e < cel->element; e++) {
          ecn = 0;
          e1 = ic*cel->element + e;
          i1 = cel->ele[e1*TRI  ];
          i2 = cel->ele[e1*TRI+1];
          i3 = cel->ele[e1*TRI+2];

          for (k = 0; k < cel->element; k++) {
            if (e == k) continue;
            e2 = ic*cel->element + k;
            j1 = cel->ele[e2*TRI  ];
            j2 = cel->ele[e2*TRI+1];
            j3 = cel->ele[e2*TRI+2];

            if ((i1 == j1 && i2 == j2)
                || (i1 == j1 && i2 == j3)
                || (i1 == j2 && i2 == j3)
                || (i1 == j2 && i2 == j1)
                || (i1 == j3 && i2 == j1)
                || (i1 == j3 && i2 == j2)
                || (i2 == j1 && i3 == j2)
                || (i2 == j1 && i3 == j3)
                || (i2 == j2 && i3 == j3)
                || (i2 == j2 && i3 == j1)
                || (i2 == j3 && i3 == j1)
                || (i2 == j3 && i3 == j2)
                || (i3 == j1 && i1 == j2)
                || (i3 == j1 && i1 == j3)
                || (i3 == j2 && i1 == j3)
                || (i3 == j2 && i1 == j1)
                || (i3 == j3 && i1 == j1)
                || (i3 == j3 && i1 == j2)){
              ec[e1*TRI+ecn] = e2;
              ecn++;
            }
            else continue;
          }
        }
      }
    }

    // Initializing
    for (jid = 0; jid < cel->n*cel->vertex; jid++) {
      fx[jid] = 0.0e+00; 
      fy[jid] = 0.0e+00;
      fz[jid] = 0.0e+00; 
    }
 
    // Procedure for periodic boundary
    for (i = 0; i < cel->n; i++) {
      for (e = 0; e < cel->element; e++) {
        for (k = 0; k < TRI; k++) {
          j1 = cel->ele[i*cel->element*TRI + e*TRI+k];
          j2 = cel->ele[i*cel->element*TRI + e*TRI+(k+1)%TRI];
          xx = cel->x[j2] - cel->x[j1];
          zz = cel->z[j2] - cel->z[j1];
          if      (xx >  cdo->LX*0.5) { cel->x[j1] += cdo->LX; }
          else if (xx < -cdo->LX*0.5) { cel->x[j2] += cdo->LX; }
          if      (zz >  cdo->LZ*0.5) { cel->z[j1] += cdo->LZ; }
          else if (zz < -cdo->LZ*0.5) { cel->z[j2] += cdo->LZ; }
        }
      }
    }

    // Calculating bending rigidity
    for (eid = 0; eid < cel->n*cel->element; eid++) {
      e1 = eid;
      i1 = cel->ele[e1*TRI  ];
      i2 = cel->ele[e1*TRI+1];
      i3 = cel->ele[e1*TRI+2];

      for (k = 0;k < TRI;k++) {
        e2 = ec[e1*TRI+k];
        j1 = cel->ele[e2*TRI  ];
        j2 = cel->ele[e2*TRI+1];
        j3 = cel->ele[e2*TRI+2];

        if      ((i1 == j1 && i2 == j2) || (i1 == j2 && i2 == j1)) { k1 = i3; k2 = i1; k3 = i2; k4 = j3; }
        else if ((i1 == j2 && i2 == j3) || (i1 == j3 && i2 == j2)) { k1 = i3; k2 = i1; k3 = i2; k4 = j1; }
        else if ((i1 == j3 && i2 == j1) || (i1 == j1 && i2 == j3)) { k1 = i3; k2 = i1; k3 = i2; k4 = j2; }
        else if ((i2 == j1 && i3 == j2) || (i2 == j2 && i3 == j1)) { k1 = i1; k2 = i2; k3 = i3; k4 = j3; }
        else if ((i2 == j2 && i3 == j3) || (i2 == j3 && i3 == j2)) { k1 = i1; k2 = i2; k3 = i3; k4 = j1; }
        else if ((i2 == j3 && i3 == j1) || (i2 == j1 && i3 == j3)) { k1 = i1; k2 = i2; k3 = i3; k4 = j2; }
        else if ((i3 == j1 && i1 == j2) || (i3 == j2 && i1 == j1)) { k1 = i2; k2 = i3; k3 = i1; k4 = j3; }
        else if ((i3 == j2 && i1 == j3) || (i3 == j3 && i1 == j2)) { k1 = i2; k2 = i3; k3 = i1; k4 = j1; }
        else if ((i3 == j3 && i1 == j1) || (i3 == j1 && i1 == j3)) { k1 = i2; k2 = i3; k3 = i1; k4 = j2; }

        a1x = cel->x[k1]; a2x = cel->x[k2]; a3x = cel->x[k3]; a4x = cel->x[k4];
        a1y = cel->y[k1]; a2y = cel->y[k2]; a3y = cel->y[k3]; a4y = cel->y[k4];
        a1z = cel->z[k1]; a2z = cel->z[k2]; a3z = cel->z[k3]; a4z = cel->z[k4];

        if (a2x - a1x > cdo->LX*0.5) a2x -= cdo->LX; else if (a2x - a1x < -cdo->LX*0.5) a2x += cdo->LX;
        if (a3x - a1x > cdo->LX*0.5) a3x -= cdo->LX; else if (a3x - a1x < -cdo->LX*0.5) a3x += cdo->LX;
        if (a4x - a1x > cdo->LX*0.5) a4x -= cdo->LX; else if (a4x - a1x < -cdo->LX*0.5) a4x += cdo->LX;
        if (a2z - a1z > cdo->LZ*0.5) a2z -= cdo->LZ; else if (a2z - a1z < -cdo->LZ*0.5) a2z += cdo->LZ;
        if (a3z - a1z > cdo->LZ*0.5) a3z -= cdo->LZ; else if (a3z - a1z < -cdo->LZ*0.5) a3z += cdo->LZ;
        if (a4z - a1z > cdo->LZ*0.5) a4z -= cdo->LZ; else if (a4z - a1z < -cdo->LZ*0.5) a4z += cdo->LZ;

        a21x = a2x - a1x; a31x = a3x - a1x; a34x = a3x - a4x; a24x = a2x - a4x; a23x = a2x - a3x;
        a21y = a2y - a1y; a31y = a3y - a1y; a34y = a3y - a4y; a24y = a2y - a4y; a23y = a2y - a3y;
        a21z = a2z - a1z; a31z = a3z - a1z; a34z = a3z - a4z; a24z = a2z - a4z; a23z = a2z - a3z;

        xix = a21y*a31z - a21z*a31y; zex = a34y*a24z - a34z*a24y;
        xiy = a21z*a31x - a21x*a31z; zey = a34z*a24x - a34x*a24z;
        xiz = a21x*a31y - a21y*a31x; zez = a34x*a24y - a34y*a24x;

        t1x = (a1x + a2x + a3x)/3.0; t2x = (a4x + a2x + a3x)/3.0;
        t1y = (a1y + a2y + a3y)/3.0; t2y = (a4y + a2y + a3y)/3.0;
        t1z = (a1z + a2z + a3z)/3.0; t2z = (a4z + a2z + a3z)/3.0;

        cost0 = 1.0;
        sint0 = 0.0;

        cost  = (xix*zex + xiy*zey + xiz*zez) 
               /sqrt(xix*xix + xiy*xiy + xiz*xiz)
               /sqrt(zex*zex + zey*zey + zez*zez);
        cost2 = 1.0 - cost*cost; if (cost2 < eps) continue;

        if (((xix-zex)*(t1x-t2x) + (xiy-zey)*(t1y-t2y) + (xiz-zez)*(t1z-t2z)) >= 0.0)
        { sint =  sqrt(cost2); beta =  2.0/sqrt(3.0)*br*(sint*cost0 - cost*sint0)/sqrt(cost2); }
        else
        { sint = -sqrt(cost2); beta = -2.0/sqrt(3.0)*br*(sint*cost0 - cost*sint0)/sqrt(cost2); }

        b11 = -beta*cost/(xix*xix + xiy*xiy + xiz*xiz);
        b12 =  beta/sqrt(xix*xix + xiy*xiy + xiz*xiz)/sqrt(zex*zex + zey*zey + zez*zez);
        b22 = -beta*cost/(zex*zex + zey*zey + zez*zez);

        fxb = b11*(xiy*-a23z - xiz*-a23y) + b12*(zey*-a23z - zez*-a23y);
        fyb = b11*(xiz*-a23x - xix*-a23z) + b12*(zez*-a23x - zex*-a23z);
        fzb = b11*(xix*-a23y - xiy*-a23x) + b12*(zex*-a23y - zey*-a23x);

        fx[k1] += fxb;
        fy[k1] += fyb;
        fz[k1] += fzb;

        fxb = b11*(xiy*-a31z - xiz*-a31y) + b22*(zey*a34z - zez*a34y)
            + b12*((xiy*a34z - xiz*a34y) + (zey*-a31z - zez*-a31y));
        fyb = b11*(xiz*-a31x - xix*-a31z) + b22*(zez*a34x - zex*a34z)
            + b12*((xiz*a34x - xix*a34z) + (zez*-a31x - zex*-a31z));
        fzb = b11*(xix*-a31y - xiy*-a31x) + b22*(zex*a34y - zey*a34x)
            + b12*((xix*a34y - xiy*a34x) + (zex*-a31y - zey*-a31x));

        fx[k2] += 0.5*fxb;
        fy[k2] += 0.5*fyb;
        fz[k2] += 0.5*fzb;

        fxb = b11*(xiy*a21z - xiz*a21y) + b22*(zey*-a24z - zez*-a24y)
            + b12*((xiy*-a24z - xiz*-a24y) + (zey*a21z - zez*a21y));
        fyb = b11*(xiz*a21x - xix*a21z) + b22*(zez*-a24x - zex*-a24z)
            + b12*((xiz*-a24x - xix*-a24z) + (zez*a21x - zex*a21z));
        fzb = b11*(xix*a21y - xiy*a21x) + b22*(zex*-a24y - zey*-a24x)
            + b12*((xix*-a24y - xiy*-a24x) + (zex*a21y - zey*a21x));

        fx[k3] += 0.5*fxb;
        fy[k3] += 0.5*fyb;
        fz[k3] += 0.5*fzb;
      }
    }

    return;
}

void ShearStress
//============================================================================== 
(
// 24th Des 2018
// SHEAR RATE/STRESS FOR SINGLE CELL
//
//
    double    *gex,
    double    *gey,
    double    *gez,
    lattice   *ltc,
    cell      *cel,
    domain    *cdo
)
//------------------------------------------------------------------------------   
{

    int       ic, i, j, k, e, ea, eb, ec, ix, iy, iz,
	      i___, im__, i_m_, i__m, i_mm, im_m, imm_, immm,
	      nx = cdo->nx, ny = cdo->ny, nz = cdo->nz;
    double    x0, y0, z0, x1, y1, z1, x2, y2, z2, x01, y01, z01, x02, y02, z02,
              vx, vy, vz, vn, vnx, vny, vnz, xg, yg, zg,xtr, ytr, ztr,
              s, t, r, u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z, shrx, shry, shrz,
	      dx = cdo->dx;

// Initializing
    cel->ave_sh = 0.0e+00,
    cel->sd_sh  = 0.0e+00,
    cel->max_sh =-1.0e+10,
    cel->min_sh = 1.0e+10;

// Centroid of an element (gex, gey, gez)
    for (ic = 0; ic < cel->n; ic++) {
      for (e = 0; e < cel->element; e++) {
        ea = cel->ele[ic*cel->element*TRI + e*TRI + 0];
        eb = cel->ele[ic*cel->element*TRI + e*TRI + 1];
        ec = cel->ele[ic*cel->element*TRI + e*TRI + 2];
        
        x0 = cel->x[ea];  x1 = cel->x[eb];  x2 = cel->x[ec];
        y0 = cel->y[ea];  y1 = cel->y[eb];  y2 = cel->y[ec];
        z0 = cel->z[ea];  z1 = cel->z[eb];  z2 = cel->z[ec];     

	x01 = x1 - x0;  x02 = x2 - x0;
	y01 = y1 - y0;  y02 = y2 - y0;
   	z01 = z1 - z0;  z02 = z2 - z0;

	// Periodic process
        if (z01 > cdo->LZ*0.5 && z02 > cdo->LZ*0.5) {z0 += cdo->LZ;}
        if (z01 > cdo->LZ*0.5 && z02 < 0.0e+00)     {z0 += cdo->LZ; z2 += cdo->LZ;}
        if (z01 < 0.0e+00     && z02 > cdo->LZ*0.5) {z0 += cdo->LZ; z1 += cdo->LZ;}
      
        gex[ic*cel->element + e] = (x0 + x1 + x2)/3.0;
        gey[ic*cel->element + e] = (y0 + y1 + y2)/3.0;
        gez[ic*cel->element + e] = (z0 + z1 + z2)/3.0; 

	// Periodic process
        if (gez[ic*cel->element + e] > cdo->LZ*0.5) gez[ic*cel->element + e] -= cdo->LZ;
      }
    }

// Calculating shear rate on the membrane
    for (ic = 0; ic < cel->n; ic++) {
      for (e = 0; e < cel->element; e++) {
        // Calculation of normal vector
	ea = cel->ele[ic*cel->element*TRI + e*TRI + 0];
	eb = cel->ele[ic*cel->element*TRI + e*TRI + 1];
	ec = cel->ele[ic*cel->element*TRI + e*TRI + 2];

	x0 = cel->x[ea];  x1 = cel->x[eb];  x2 = cel->x[ec];
        y0 = cel->y[ea];  y1 = cel->y[eb];  y2 = cel->y[ec];
	z0 = cel->z[ea];  z1 = cel->z[eb];  z2 = cel->z[ec];     

	x01 = x1 - x0;  x02 = x2 - x0;
	y01 = y1 - y0;  y02 = y2 - y0;
   	z01 = z1 - z0;  z02 = z2 - z0;
 
     	vx = (y01*z02 - z01*y02);
     	vy = (z01*x02 - x01*z02);
     	vz = (x01*y02 - y01*x02); 

	vn  = sqrt(vx*vx + vy*vy + vz*vz); 
        vnx = vx/vn;
        vny = vy/vn;
        vnz = vz/vn;

        // Adjacent fluid latiices (Euler points) to an element ---------
        xg = gex[ic*cel->element + e];
        yg = gey[ic*cel->element + e];
        zg = gez[ic*cel->element + e];

        xtr = xg + (double)(nx - 1)*dx/2.0;
        ytr = yg + (double)(ny - 1)*dx/2.0;
        ztr = zg + (double)(nz - 1)*dx/2.0;

        ix = (int)ceil(xtr/dx); if (ix > nx-5) ix -= (nx-5); 
        iy = (int)ceil(ytr/dx); if (iy > ny-5) iy -= (ny-5);
        iz = (int)ceil(ztr/dx); if (iz > nz-5) iz -= (nz-5);

        i  = ix + nx*iy + nx*ny*iz;

	// Level-set (distance) function
	s = xg/dx - (double)(ix - 1);
	t = yg/dx - (double)(iy - 1);
	r = zg/dx - (double)(iz - 1);

        // Calculating velocity gradient (u_, v_, w_)
        i___ = i;
	
       	if (ix <= 2) im__ = i + (nx-5);
	else         im__ = i - 1;

	if (iy <= 2) i_m_ = i + nx*(ny-5);
	else         i_m_ = i - nx;

	if (iz <= 2) i__m = i + nx*ny*(nz-5);
	else         i__m = i - nx*ny;

	if (ix <= 2 && iy <= 2) imm_ = i + (nx-5) + nx*(ny-5);
	else if (ix <= 2)       imm_ = i + (nx-5) - nx;
	else if (iy <= 2)       imm_ = i - 1 + nx*(ny-5);
	else                    imm_ = i - 1 - nx;

	if (ix <= 2 && iz <= 2) im_m = i + (nx-5) + nx*ny*(nz-5);
	else if (ix <= 2)       im_m = i + (nx-5) - nx*ny;
	else if (iz <= 2)       im_m = i - 1 + nx*ny*(nz-5);
	else                    im_m = i - 1 - nx*ny;

        if (iy <= 2 && iz <= 2) i_mm = i + nx*(ny-5) + nx*ny*(nz-5);
	else if (iy <= 2)       i_mm = i + nx*(ny-5) - nx*ny;
	else if (iz <= 2)       i_mm = i - nx + nx*ny*(nz-5);
	else	                i_mm = i - nx - nx*ny;

	if (ix <= 2 && iy <= 2 && iz <= 2) immm = i + (nx-5) + nx*(ny-5) + nx*ny*(nz-5);
	else if (ix <= 2 && iy <= 2)       immm = i + (nx-5) + nx*(ny-5) - nx*ny;
	else if (ix <= 2 && iz <= 2)       immm = i + (nx-5) - nx + nx*ny*(nz-5);
	else if (iy <= 2 && iz <= 2)       immm = i - 1 + nx*(ny-5) + nx*ny*(nz-5);
	else if (ix <= 2)                  immm = i + (nx-5) - nx - nx*ny;
	else if (iy <= 2)                  immm = i - 1 + nx*(ny-5) - nx*ny;
	else if (iz <= 2)                  immm = i - 1 - nx + nx*ny*(nz-5);
	else 	                           immm = i - 1 - nx - nx*ny;

	u_x = (1.0-r)*( -(1.0-t)*ltc->u[immm] + (1.0-t)*ltc->u[i_mm] - t*ltc->u[im_m] + t*ltc->u[i__m] ) 
	          + r*( -(1.0-t)*ltc->u[imm_] + (1.0-t)*ltc->u[i_m_] - t*ltc->u[im__] + t*ltc->u[i___] );

	u_y = (1.0-r)*( -(1.0-s)*ltc->u[immm] - s*ltc->u[i_mm] + (1.0-s)*ltc->u[im_m] + s*ltc->u[i__m] ) 
	          + r*( -(1.0-s)*ltc->u[imm_] - s*ltc->u[i_m_] + (1.0-s)*ltc->u[im__] + s*ltc->u[i___] );

	u_z = -( (1.0-s)*(1.0-t)*ltc->u[immm] + s*(1.0-t)*ltc->u[i_mm] + (1.0-s)*t*ltc->u[im_m] + s*t*ltc->u[i__m])
	      +( (1.0-s)*(1.0-t)*ltc->u[imm_] + s*(1.0-t)*ltc->u[i_m_] + (1.0-s)*t*ltc->u[im__] + s*t*ltc->u[i___] );

	v_x = (1.0-r)*( -(1.0-t)*ltc->v[immm] + (1.0-t)*ltc->v[i_mm] - t*ltc->v[im_m] + t*ltc->v[i__m] ) 
	          + r*( -(1.0-t)*ltc->v[imm_] + (1.0-t)*ltc->v[i_m_] - t*ltc->v[im__] + t*ltc->v[i___] );

	v_y = (1.0-r)*( -(1.0-s)*ltc->v[immm] - s*ltc->v[i_mm] + (1.0-s)*ltc->v[im_m] + s*ltc->v[i__m] ) 
	          + r*( -(1.0-s)*ltc->v[imm_] - s*ltc->v[i_m_] + (1.0-s)*ltc->v[im__] + s*ltc->v[i___] );

	v_z = -( (1.0-s)*(1.0-t)*ltc->v[immm] + s*(1.0-t)*ltc->v[i_mm] + (1.0-s)*t*ltc->v[im_m] + s*t*ltc->v[i__m] )
	      +( (1.0-s)*(1.0-t)*ltc->v[imm_] + s*(1.0-t)*ltc->v[i_m_] + (1.0-s)*t*ltc->v[im__] + s*t*ltc->v[i___] );

	w_x = (1.0-r)*( -(1.0-t)*ltc->w[immm] + (1.0-t)*ltc->w[i_mm] - t*ltc->w[im_m] + t*ltc->w[i__m] ) 
	          + r*( -(1.0-t)*ltc->w[imm_] + (1.0-t)*ltc->w[i_m_] - t*ltc->w[im__] + t*ltc->w[i___] );

	w_y = (1.0-r)*( -(1.0-s)*ltc->w[immm] - s*ltc->w[i_mm] + (1.0-s)*ltc->w[im_m] + s*ltc->w[i__m] ) 
	          + r*( -(1.0-s)*ltc->w[imm_] - s*ltc->w[i_m_] + (1.0-s)*ltc->w[im__] + s*ltc->w[i___] );

	w_z = -( (1.0-s)*(1.0-t)*ltc->w[immm] + s*(1.0-t)*ltc->w[i_mm] + (1.0-s)*t*ltc->w[im_m] + s*t*ltc->w[i__m] )
	      +( (1.0-s)*(1.0-t)*ltc->w[imm_] + s*(1.0-t)*ltc->w[i_m_] + (1.0-s)*t*ltc->w[im__] + s*t*ltc->w[i___] );

	u_x = u_x/dx;  u_y = u_y/dx;  u_z = u_z/dx;
	v_x = v_x/dx;  v_y = v_y/dx;  v_z = v_z/dx;
	w_x = w_x/dx;  w_y = w_y/dx;  w_z = w_z/dx;
    
        // Calculating shear rates   
	shrx = ( vnx*(u_x + u_x) + vny*(u_y + v_x) + vnz*(u_z + w_x) )*(1.0-vnx*vnx)
	      +( vnx*(v_x + u_y) + vny*(v_y + v_y) + vnz*(v_z + w_y) )*(-vny*vnx)
	      +( vnx*(w_x + u_z) + vny*(w_y + v_z) + vnz*(w_z + w_z) )*(-vnz*vnx);

	shry = ( vnx*(u_x + u_x) + vny*(u_y + v_x) + vnz*(u_z + w_x) )*(-vnx*vny)
	      +( vnx*(v_x + u_y) + vny*(v_y + v_y) + vnz*(v_z + w_y) )*(1.0-vny*vny)
	      +( vnx*(w_x + u_z) + vny*(w_y + v_z) + vnz*(w_z + w_z) )*(-vnz*vny);

	shrz = ( vnx*(u_x + u_x) + vny*(u_y + v_x) + vnz*(u_z + w_x) )*(-vnx*vnz)
	      +( vnx*(v_x + u_y) + vny*(v_y + v_y) + vnz*(v_z + w_y) )*(-vny*vnz)
	      +( vnx*(w_x + u_z) + vny*(w_y + v_z) + vnz*(w_z + w_z) )*(1.0-vnz*vnz);

	cel->sh[ic*cel->element + e] = sqrt(shrx*shrx + shry*shry + shrz*shrz);
	cel->ave_sh                 += sqrt(shrx*shrx + shry*shry + shrz*shrz);

	if (cel->max_sh < cel->sh[ic*cel->element + e]) cel->max_sh = cel->sh[ic*cel->element + e]; 
	if (cel->min_sh > cel->sh[ic*cel->element + e]) cel->min_sh = cel->sh[ic*cel->element + e]; 
	
	for (k = 0; k < TRI; k++) {
	  if      (k == 0) j = ea;
	  else if (k == 1) j = eb;
	  else if (k == 2) j = ec;
	  cel->sh_node[ic*cel->vertex*3 + j*3 + 0] += shrx/3.0;
	  cel->sh_node[ic*cel->vertex*3 + j*3 + 1] += shry/3.0;
	  cel->sh_node[ic*cel->vertex*3 + j*3 + 2] += shrz/3.0;
	}

      } // element
    } // n       

    cel->ave_sh /= (double)(cel->n*cel->element);

    for (ic = 0; ic < cel->n; ic++) {
      for (e = 0; e < cel->element; e++) {
        cel->sd_sh += (cel->sh[ic*cel->element + e] - cel->ave_sh)*(cel->sh[ic*cel->element + e] - cel->ave_sh);
      }
    }

    cel->sd_sh = sqrt(cel->sd_sh);

    return;
}
