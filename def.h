
// Definition of parameters ---
#define MAX(x,y)       ((x) > (y) ? (x) : (y)) 
#define MIN(x,y)       ((x) < (y) ? (x) : (y)) 

#define CYTOADHESION             false
#define REPULSION                false
#define CONSTRICTED_CHANNELFLOW  false
#define POLYMERIZATION           false

#define DEV            0
#define Q              19
#define TRI            3
#define HEX            6
#define DOB            2
#define DIV            256
#define N              9  // Number of matrix components
#define TH            10  // Number of neighboring lattices for interporation

#define COMP           0
#define BUFF           1
#define INLET          2
#define OUTLET         3
#define WALL           4
#define MU             1.2e-03 // Viscosity of plasma [Pa*s]
#define KC             5.0e-19 // [J]
#define DP             1.025e+03 // Density of plasma [kg/(m*m*m)]

#define DX             2.500e-07 // Lattice space (=0.250[um]) [m]
//#define DX             1.250e-07 // Lattice space (=0.125[um]) [m]

#define DKNOB          1.25e+00  // [Microvillus/um^2]
#define KBT            (310.0*1.38e-23) // [J]
#define C              1.0e+02  // Area dilation modulus

// Information of computation process ---
typedef struct
{
    int       iter,                    // Iteration
              iter_limit,              // Limit of iteration
              file_interval,           // Iteration interval for file output
              fnum_ad,
              fnum_iter,
              IterLimit,
              Fout,
              FoutAd,
              STout,
              Num_Flag,
              FlowS,
              END_FILE;
} process;


// Information of computational domain ---
typedef struct
{
    double    dt,                      // Time intarbal
              dx,                      // Size of lattice
              sdT,                     // Dimensional time for smaller separation
              sr,                      // Sheare rate
              tau,                     // Relaxation time
              freq,                    // fequency
	      lam1,                    // Viscosity ratio of inner fluid (RBC) & outer fluid (plasma)
	      lam2,                    // Viscosity ratio of inner fluid (WBC) & outer fluid (plasma)
              xmin,                    // Domain size
              xmax,
              ymin,
              ymax,
              zmin,
              zmax,
              D,
              U,
              L,
              LX,
              LY,
              LZ;
    int       nx,                      // A number of nodes
              ny,
              nz,
              nw,                      // A number of all wall nodes
              n,                       // A number of all nodes
	      flag,
              comp;                    // A number of computational domain
} domain;


// Information of lattice Boltzmann nodes ---
typedef struct
{
    double    *buV,
              *avel,                   // Average blood velocity z-direction
              *d,                      // Density
              *u,                      // Velocity
              *v,
              *w,
              *vf,                     // VOF of a red blood cell
              *vf2,                    // VOF of other cells
              Ca,                      // Capillary number
              Re,                      // Reynolds number
              Df,                      // Diffusion number
              Cr,                      // Courant number
              KL,                      // Resolution of a lattice (capusle)
              CKL;                     // Resolution of a capsule
    int       *bc;                     // Domain definition
// Variables about adhesion
    double    *xw,                     // Corrdinates of wall
              *yw,
              *zw;
} lattice;


// Information of cell ---
typedef struct
{
// Variables for RBC & WBC & Platelet & Cancer cell
    double    CR,               // Radius of a RBC
              Gs,               // Surface shear elastic modulus
              *xr,              // Reference coordinates
              *x, *xD,          // Coordinates
	      *y, *yD, 
	      *z, *zD,    
              *u,               // Velocity
	      *v,
	      *w,       
              *q,               // Load
              *f,               // Force acting on element center
              *t,               // Tension
              *maxr, *maxrD,    // Maximum length between two nodes on a cell
              *minr, *minrD;    // Minimum length between two nodes on a cell
    int       flag,             // Flag for inflation process
              n,                // A number of cells
              n_r,              // A number of RBCs
              n_irbc,           // A number of Infected RBC
              vertex,           // A number of vertexes
              element,          // A number of triangular elements
              *ele, *eleD,      // RBC's vertex number consisting of triangular element
              model;
// Variables for WBCs
    double    CW,               // Radius of a WBC
              RGw,              // Ratio of Gs
              *x_w, *y_w, *z_w, // Coordinates
              *u_w, *v_w, *w_w, // Velocity
              *q_w,             // Load
              *f_w,             // Force acting on element center
              *t_w,             // Tension
              *rest_Dz_w;       // Z-diameter of resting shape
    int       n_w,              // A number of WBCs
              vertex_w,         // A number of vertexes
              element_w,        // A number of triangular elements
              *ele_w,           // WBC's vertex number consisting of triangular element
              model_w;
// Variables for platelet
    double    CP,               // Radius of a platelet
              RGp,              // Ratio of Gs
              *x_p, *y_p, *z_p, // Coordinates
              *u_p, *v_p, *w_p, // Velocity
              *q_p,             // Load
              *f_p,             // Force acting on element center
              *t_p,             // Tension
              *rest_Dz_p,       // Z-diameter of resting shape
              *dl,              // Averaged distance between the centroid & wall
              *dr;              // Averaged distance between the membrane & wall
    int       n_p,              // A number of platelets
              vertex_p,         // A number of vertexes
              element_p,        // A number of triangular elements
              *ele_p,           // Platelet's vertex number consisting of triangular element
              model_p;
// Variables for carcinoma
    double    CC,               // Radius of a cancer cell
              RGc,              // Ratio of Gs
              *x_c, *y_c, *z_c, // Coordinates
              *u_c, *v_c, *w_c, // Velocity
              *q_c,             // Load
              *f_c,             // Force acting on element center
              *t_c,             // Tension
              *rest_Dz_c;       // Z-diameter of resting shape
    int       n_c,              // A number of cancer cells
              vertex_c,         // A number of vertexes
              element_c,        // A number of triangular elements
              *ele_c,           // Cancer cell's vertex number consisting of triangular element
              model_c;
// Variables for capsule centroid
    double    *xc,              // Centroid position on axial direction
              *nxc,             // Non-Centroid position on axial direction
              *xd,              // Dstance between membrane and wall surface
              *xv,              // Radial velocity
              *rc,              // RBCs centroid position
              *vc,              // RBCs centroid velocity
              *gp,              // RBCs centroid position calculated by vertex
              *nont,            // Nondimensional time
              *ave_cent,        // RBCs average velocity z-direction
              *vave,            // Relative average velocity z-direction
              *vmax,            // Relative maximum velocity z-direction
              *difv,            // Difference between cell velocity and RBCs average velocity
              *ave2,            // Number of average for RBC & time aberage
              *rbc_bv,          // Velocity of averaged RBCs relative to blood velocity
              *wbc_bv,          // Velocity of WBCs relative to blood velocity
              *p_bv,            // Velocity of average platelets relative to blood velocity
              *rwc_bv,          // Velocity of (RBC - WBC)/Blood
              *rpc_bv,          // Velocity of (RBC - Platelets)/Blood
              *wrc_bv,          // Velocity of (WBC - RBC)/Blood
              *prc_bv,          // Velocity of (Platelets - RBC)/Blood
              *nwbc_bv,         // non-Velocity of WBC/Blood
              *np_bv,           // non-Velocity of average platelets/Blood
              *nrbc_bv;         // non-Velocity of RBC/Blood
// Variables for adhesion model
    double    *fadx,            // Adhesive force acting on the node
              *fady,
              *fadz,
              *detach,          // Number of detach length
              *touch,           // Number of touch length
              npua,             // Number of knob on an unit area
              lim,              // Limitation of spring length [mesh]
              sl;               // Natural length of spring [mesh]
    int       *flagad,          // Adhesion node point
              *knob,            // ID of knob
              knobn,            // Number of knobs
              nure,             // Number of spring at a ligad position
              max_knob,         // Maximum number of knob
              max_det_cat;      // maximum number of detachment & catch point
// Variables for sorting node
   double     *xyz;             // Coordinate for x, y, z-direction
   int        *nln,             // Number of linking node
              *ln,              // Linking node id number
              *lo,              // Order of linking node
              *nle,             // Number of linking element
              *le;              // Linking element id number
// Variables about adhesion
    double    *f_adh,                  // adhesive force acting on the node
              *k_on,                   // association rate [-]
              *k_off,                  // dissociation rate [-]
              *r_ij,                   // distance between membrane & wall [-]
              *r_wall,                 // distance between membrane & wall [-]
              *avev,                   // dissociation rate []
              k_on0,                   // stocastic parameter for adhesion [1/s]
              k_off0,                  // stocastic parameter for adhesion [1/s]
              p_k_on0,                 // stocastic parameter for polymerization [1/(muM*s)]
              p_k_off0,                // stocastic parameter for adhesion [1/s]
              G_EC,
              G_EL,
              kad,                     // spring constant [N/m]
              chib,                    // reactive compliance [m]
              limit,                   // limitation of spring length [mesh]
              l0,                      // natural length of spring length [mesh]
              Drecep,                  // receptor density [/m^2]
              Nknob;                   // number of knob (microvilli) [/cell]
    int       adnum,
	      advertex,
	      *attachment_point,       // adhesion wall point
              nspring_knob;
    bool      *is_knob;
// Variables about polymerization
    double    *f_pol,
              *vg,
              *vec,
	      *ddlx;
    bool      *point,
	      *ddlp;
// Variables about repulsion force
    bool      *is_rep, *is_rep_w;
    double    *r_ep,   *r_ep_w,
              *f_rep,  *f_rep_w;
// Variables about repulsion force
    double    *sh,
	      *sh_node,
	      ave_sh,
	      sd_sh,
	      max_sh,
	      min_sh;
} cell;

typedef struct
{
    int     nz,
            ntheta,
            n_wall,
            max_neighwall,
           *is_attached;
    double  radius,
            dz,
            dtheta,
            *x;
} Wall;

// err.c ---
void  error
//=======================================================
//
//  PRINT ERROR MESSAGE AND EXIT THE PROGRAM
//
//
//
(
    int       error_type
);
//-------------------------------------------------------


// file.c ---
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
);
//-------------------------------------------------------

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
);
//-------------------------------------------------------

int  read_fluid_bin
//=======================================================
//
//  READ BINARY DATAFILE OF CELL
//
//
//
(
    int       icnt,
//    int       post_process,
    process   *prc,
    domain    *cdo,
    lattice   *ltc,
    cell      *cel,
    char      *filename
);
//-------------------------------------------------------

int  read_cell_bin
//=======================================================
//
//  READ TEXT DATAFILE OF CELL
//
//
(
    int       icnt,
    cell      *cel,
    Wall      *wall,
    char      *filename
);
//-------------------------------------------------------

// theorem.c ---
void  Add_Nlink (int na, int nb, int *nln, int *ln);
void  Add_Elink (int ni, int ei, int *nle, int *le);
int   Checklink (int na, int nb, int *nln, int *ln);
void  Sort_Link_Spring (double *x, int *nln, int *ln, int vertex, int ic);
void  Sort_Link_Order  (int *nln, int *ln, int vertex, int ic);

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
    int       *ele,
    double    *x,
    double    *y,
    double    *z,
    cell      *cel,
    domain    *cdo
);
//-------------------------------------------------------

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
    int       *ele,
    double    *x,
    double    *y,
    double    *z,
    double    *u,
    double    *v,
    double    *w,
    cell      *cel,
    domain    *cdo
);
//-------------------------------------------------------

void  Centroid_Capsules
//=======================================================
//
//  CALCULATE THE CENTROID OF RBCS
//
//
//
(
    cell      *cel,
    domain    *cdo
);
//-------------------------------------------------------

void Bending
//=======================================================
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
);
//-------------------------------------------------------

// post.c ---
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
);
//-------------------------------------------------------

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
);
//------------------------------------------------------------------------------   
