#ifndef LISFLOOD_H
#define LISFLOOD_H

/*
  #####################################################################################
  LISFLOOD-FP flood inundation model
  #####################################################################################
  extended with Basic Model Interface (BMI) compatibility
  #####################################################################################

  © copyright Bristol University Hydrology Research Group 2008

  webpage -	http://www.ggy.bris.ac.uk/research/hydrology/models/lisflood
  contact -	Professor Paul Bates, email: paul.bates@Bristol.ac.uk,
  Tel: +44-117-928-9108, Fax: +44-117-928-7878
  
  for inquiries and questions regarding BMI compatibility, please contact
  Jannis Hoch, M.Sc. Utrecht University, Faculty of Geosciences, Department Physical Geography): j.m.hoch@uu.nl

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
//#include <omp.h>
#include <vector> // CCS
#include <iostream> // CCS
//#include <netcdf.h> // JCN

// define basic constants
#define ON 1
#define OFF 0
#define CHKINTERVAL 1.0 // default checkpoint interval in hours runtime
#define NULLVAL -9999.0 // MT: define ascii file NULL value as constant

// limits
#define MAXPI 20000

// MT: disable visual C++ depreciation warnings so we can see other warnings
#pragma warning( disable : 4996)

using namespace std; // CCS

/*
*****************************************************************************

Define the structures
---------------------
Fnames - Contains all the filenames from the .par file
States - Contains all the state parameters for the simulation
Pars - Contains the parameter values specified in the .par file
Solver - Defines the solution settings
Arrays - Defines the global arrays
ChannelSegmentType - Defines the channel
Stage - Variables for outputting stage information at specified locations
Files - General output file pointers
BoundCs - Boundary conditions

*****************************************************************************
*/

/*! \struct Arrays
  Stores the pointers to arrays required globally in the computation. Defined as 1D
  vectors but stores 2D data determined by the array subscripts.
*/
struct Arrays{
  /*! DEM, Water height, Flow in x-direction and Flow in y-direction */
  double *DEM; // Digital elevation model
  double *H,*Qx,*Qy,*Qxold,*Qyold,*U,*V,*Hflowx,*Hflowy;
  /* TRENT additions  */
  double *HU,*HV,*RSHU,*LSHU,*RSHV,*LSHV,*BSHU,*TSHU,*BSHV,*TSHV,*FHx,*FHUx,*FHVx,*FHy,*FHUy,*FHVy;
  /* Flow direction map for Rainfall*/
  int *FlowDir; // CCS: added to hold DEM flow direction map for routing shallow rainfall flow 13/03/2012
  double *Route_dH, *RouteInt; // CCS: added to record routing scheme dH and interval

  /* ---------------- */
  double *maxH,*maxHtm,*initHtm,*totalHtm;
  double *Manningsn, *SGCManningsn;
  double *paerial, *pbound;
  double *Weir_hc,*Weir_Cd,*Weir_m,*Weir_w;
  int *Weir_Identx,*Weir_Identy,*Weir_Fixdir,*Weir_Typ;
  double *evap, *rain;
  int *ChanMask,*SegMask;
  double *TRecx,*TRecy; // MT: add to record TStep
  double *LimQx,*LimQy; // MT: add to record Qlimits
  double *Vx, *Vy, *maxVx, *maxVy, *Vc; // JCN: added to record velocity
  double *maxVc, *maxVcH, *maxHaz; // JCN added to calculate hazard
  double *SGCwidth, *SGCz, *QxSGold, *QySGold, *SGCbfH, *SGCVol, *SGCdVol, *SGCbfV, *SGCc, *SGCFlowWidth, *SGCdx, *SGCcat_area; // JCN added to store widths and depths
  double *SGCQin; //JMH
  double *dx, *dy, *dA; // CCS added for lat long data
  int  *SGCgroup;
} ;

//-------------------------------------------
// Files
struct Files{
  FILE *mass_fp;
  FILE *stage_fp;
  FILE *vel_fp;
  FILE *gau_fp;
} ;

//-------------------------------------------
// Fnames
struct Fnames{
  char demfilename[80];
  char startfilename[80];
  char chanfilename[80];
  //char Qfilename[80]; Not used! JCN
  char resrootname[255];
  char dirrootname[255];
  char qfilename[80];
  char nfilename[80];
  char SGCnfilename[80];
  char porfilename[80];
  char rivername[80];
  char bcifilename[80];
  char bdyfilename[80];
  char weirfilename[80];
  char opfilename[80];
  char stagefilename[80];
  char ascheaderfilename[80];
  char multiriverfilename[80]; // CCS
  char checkpointfilename[255]; // used to write a checkpoint file (note it is placed in the input file dir not the results dir)
  char loadCheckpointFilename[255]; // explicit checkpoint file to start run (specify with -loadcheck option, defaults to checkpointfilename)
  char evapfilename[255];
  char rainfilename[255];
  char logfilename[255];
  char SGCwidthfilename[255]; // JN sub grid channel widths
  char SGCbankfilename[255]; // JN sub grid channel bank elevations
  char SGCbedfilename[250];  // JN sub grid channel bed elevation
  char SGCcat_areafilename[250]; // JN sub grid channel accumulation area
  char SGCchangroupfilename[250];
  char SGCchanpramsfilename[250];
  char gaugefilename[255];
} ;

//-------------------------------------------
// Boundary Conditions
struct BoundCs{
  int   *xpi;
  int   *ypi;
  int   *PS_Ident;
  int   numPS;
  double *PS_Val;
  double *PS_qold;
  double *PS_qSGold;
  char  *PS_Name;
  int   *BC_Ident;
  double *BC_Val;
  double **BCVarlist;
  char  *BC_Name;
  double Qpoint;
  double Qin;
  double Qout;
  double QChanOut;
  double VolInMT; // added by JCN stores volume in over mass inteval
  double VolOutMT; // added by JCN stores volume out over mass inteval
} ;

//-------------------------------------------
// Stage
struct Stage{
  int Nstages, Ngauges;
  double *stage_loc_x,*stage_loc_y;
  double *gauge_loc_x, *gauge_loc_y, *gauge_dist;
  int *stage_grid_x,*stage_grid_y, *vel_grid_xy, *stage_check;
  int *gauge_grid_xy, *gauge_grid_x, *gauge_grid_y;
  int *gauge_dir, *gauge_cells;
} ;

// SGC parameters
struct SGCprams{
  int NSGCprams;
  int *SGCchantype;
  double SGCbetahmin;
  double *SGCp, *SGCr, *SGCs, *SGCn, *SGCm, *SGCa;
  double *SGCgamma, *SGCbeta1, *SGCbeta2, *SGCbeta3, *SGCbeta4, *SGCbeta5;
} ;

//-------------------------------------------
// Simulation States
struct States{
  int ChannelPresent;
  int TribsPresent;
  int NCFS;
  int save_depth;
  int save_elev;
  int out_dir;
  int single_op;
  int multi_op;
  int calc_area;
  int calc_meandepth;
  int calc_volume;
  int save_stages;
  int adaptive_ts;
  int acceleration; // PB: Flag to switch to acceleration formulation
  int qlim; // TJF: Flag for qlim version
  int debugmode;
  int save_Qs;
  int calc_infiltration;
  int call_gzip;
  int alt_ascheader;
  int checkpoint;
  int checkfile;
  int calc_evap;
  int rainfall; // TJF: added for time varying, spatially uniform rainfall
  int routing; // CCS: added for routing routine.
  int reset_timeinit;
  int profileoutput;
  int porosity;
  int weirs;
  int save_Ts;   // MT: added flag to output adaptive timestep
  int save_QLs;  // MT: added flag to output Qlimits
  int diffusive; // MT: added flag to indicate wish to use diffusive channel solver instead of default kinematic
  int startq;    // MT: added flag to indicate wish to use start flow to calculate initial water depths throughout channel
  int logfile;   // MT: added flag to record logfile
  int startfile; // MT: added flag to note use of startfile
  int start_ch_h; // MT: added flag to note use of starting H in channel
  int comp_out; // TJF: added to make computational output information optional
  int chainagecalc; // MT: added so user can switch off grid independent river chainage calculation
  int mint_hk; // JN: added to request maxH, maxHtm totalHtm and initHtm be calulated at the mass interval
  int Roe; // JN/IV: added to use Roe solver
  int killsim; // MDW: added to flag kill of simulation after specified run time
  int dhoverw; // TJF: added as a switch for dhlin (ON - dhlin set by command line/parfile; OFF - dhlin prescribed by gradient 0.0002 Cunge et al. 1980)
  int drychecking; //JN Option to turn DryCheck off
  int voutput; // exports velocity esimates based on Q's of Roe velocity (JCN)
  int steadycheck; // MDW: added flag to check for model steady state
  int hazard; // JN additional module for calculating hazards
  int startq2d; // JN: initalises inertial model with uniform flow for Qold
  int Roe_slow; // JN: ghost cell version of Roe solver
  int multiplerivers; // CCS multiple river switch
  int SGC; // JN sub gird channels r
  int SGCbed; // JN sub grid bed elevation file to override hydraulic geometry
  int SGCcat_area; // JN sub grid channel accumulated area to override hydraulic geometry based on width
  int SGCchangroup; // turns on distributed channe groups
  int SGCchanprams; // parameters for distributed channel groups
  int binary_out; // JN binary raster output
  int gsection; // JN virtual gauge sections
  int binarystartfile; // JN load a binary start file
  int startelev; // used to use an elevation file for the startfile
  int latlong; // CCS: added for lat-long coordinate systems
  int SGCbfh_mode; // JCN switches model to use parameter p as bank full depth
  int SGCA_mode; // JCN switches model to use parameter p as bank full Area
  int dist_routing; // JCN turnes on spatially distributed routing velocity
  int SGCvoutput; // JCN Turns on sub-grid channel velocity output
} ;

//-------------------------------------------
// Model Parameters
struct Pars{
  int xsz,ysz;
  double dx, dx_sqrt;
  double dy,dA;
  double FPn;
  double tlx,tly,blx,bly;
  double SaveInt, MassInt;
  double SaveTotal, MassTotal;
  int SaveNo;
  int op_multinum;
  double *op_multisteps;
  int *op_multiswitch;
  double op;
  double InfilRate;
  double InfilLoss, EvapLoss, RainLoss; // previous mass interval loss
  double InfilTotalLoss, EvapTotalLoss, RainTotalLoss; // cumulative loss
  double checkfreq, nextcheck;
  double reset_timeinit_time;
  double maxelev,zlev; // Water depth dependent porosity
  int zsz; // Water depth dependent porosity
  int Por_Ident;
  double dAPor;
  char **ascheader;
  double ch_start_h; // starting depth of channel flow. default to 2m or read from par file.
  double killsim_time; // time to kill simulation
  double steadyQdiff,steadyQtol,steadyInt,steadyTotal; // used for checking steady-state
  double SGC_p; // sub grid channel width depth exponent
  double SGC_r; // sub grid channel width depth mutiplier
  int SGCchan_type; // JCN bank slop for trapazoidal channel
  double SGC_s, SGC_2, SGC_n; // JCN trapazodal channel slope dependent constant
  double *SGCprams; // pointer to table of SGC parameters
  double Routing_Speed, RouteInt; // CCS variables controlling routing speed in rainfall routing routine
  double RouteSfThresh; // CCS water surface slope at which routing scheme takes over from shallow water eqn is SGC mode (when routing==ON).
  double SGC_m, SGC_a; // allows a meander coefficient to be set for the sub-grid model, default 1, allows channel upstream area to be set, defaul -1;
  double min_dx,min_dy,min_dx_dy; // CCS added to hold minimum values of dx and dy when using lat-long projected grids.
} ;

// Solver settings
struct Solver{
  double t;
  double g;
  double divg;
  double cfl;
  int    ts_multiple; // channel timestep multiple for running 1D decoupled from 2D
  long   Nit,itCount;
  double Sim_Time;
  double InitTstep; // Maximum timestep
  double Tstep;  // Adapting timestep
  double MinTstep;  // Stores minimum timestep during simulation
  double SolverAccuracy;
  int dynsw; // Switch for full dynamic steady state or diffusive steady state
  double Impfactor;
  double Hds;
  double vol1,vol2;
  double Qerror;
  double Verror;
  double FArea; // Store flooded area
  double DepthThresh,MomentumThresh,MaxHflow;
  double dhlin;
  double htol;
  double Qlimfact; // MT added to allow user to relax Qlimit
  double itrn_time;
  double itrn_time_now;
  double SGCtmpTstep; // JCN added to enable time step calculatin in UpdateH for SGC method
  time_t time_start;
  time_t time_finish;
  time_t time_check;
  double theta; //GAMA added for q-centred scheme
  int fricSolver2D; //GAMA: Solves the friction term using the vectorial (2D) scheme
} ;

//-------------------------------------------
// ChannelSegmentType
struct ChannelSegmentType{
  double *Chandx;
  double *Shalf;
  double *Chainage;
  double *ChanQ; // only for recording Q for output in profile
  double *A;
  double *NewA;
  double *ChanWidth;
  double *ChanN;
  int   *ChanX;
  int 	*ChanY;
  double *Q_Val;
  double **QVarlist;
  int   *Q_Ident;
  char  *Q_Name;
  double *BankZ;
  int chsz;
  int Next_Segment;
  int Next_Segment_Loc;
  int N_Channel_Segments;
  double JunctionH; // allows recording of H data for dummy junction point of tributary, without overwriting main channel info
  double JunctionDEM; // allows recording of DEM data for dummy junction point of tributary, without overwriting main channel info
} ;

//-------------------------------------------
/* QID7_Store // CCS for temp storage of trib boundary condition info when (Q_Ident_tmp[i]==7) in LoadRiver. A vector containing these
   structures is built in LoadRiver function and used in UpdateChannelsVector function. */
struct QID7_Store{
  int trib;
  int Next_Segment_Loc;
  int chseg;
  int RiverID;
};
//-------------------------------------------


/*
*****************************************************************************

Define the function prototypes
---------------------
Prototypes split into approximate groups that correspond to file locations
for easier editing.

*****************************************************************************
*/

// Input prototypes - input.cpp
void ReadParamFile(const char *, Fnames *, States *, Pars *, Solver*, int*);
void LoadDEM(Fnames *,States *,Pars *,Arrays *,int*);
void LoadManningsn(Fnames *,Pars *,Arrays *,int*);
void LoadSGCManningsn(Fnames *,Pars *,Arrays *,int*);
void LoadRiverNetwork(Fnames *,States *, Pars *, vector<ChannelSegmentType> *, Arrays *,vector<QID7_Store> *,vector<int> *, int*); // CCS
void LoadRiver(Fnames *,States *, Pars *,vector<ChannelSegmentType> *,Arrays *,vector<QID7_Store> *, vector<int> *, int*);
void UpdateChannelsVector(States *, ChannelSegmentType *, vector<QID7_Store> *, QID7_Store *,int *); // CCS
void LoadStart(Fnames *,States *,Pars *,Arrays *,SGCprams *, int *);
void LoadStart(Fnames *,Pars *,Arrays *, int *);
void LoadBCs(Fnames *,States *,Pars *,BoundCs *,Arrays *,int *);
void LoadBCVar(Fnames *,States *,Pars *,BoundCs *,ChannelSegmentType *,Arrays *,vector<ChannelSegmentType> *, int *);
void LoadWeir(Fnames *,States *,Pars *,Arrays *,int*);
void LoadStages(Fnames *,States *,Pars *,Stage *,int*);
void LoadGauges(Fnames *,States *,Pars *,Stage *,int*);
void LoadPor(Fnames *,States *,Pars *,Arrays *, int *);
void LoadEvap(Fnames *,Arrays *, int *);
void LoadRain(Fnames *,Arrays *, int *);
void LoadSGC(Fnames *,Pars *,Arrays *,States *, SGCprams *, int *);
void LoadBinaryStart(Fnames *,States *,Pars *,Arrays *, int *);
void LoadSGCChanPrams(Fnames *,States *,Pars *,SGCprams *, int *);

// LISFLOOD Solution prototypes - iterateq.cpp
void IterateQ(); // Fnames *, Files *,States *, Pars *, Solver*, BoundCs *, Stage *, ChannelSegmentType *, Arrays *, vector<int> *, int *, vector<ChannelSegmentType> *, int *);
void UpdateH(States *, Pars *, Solver *,BoundCs *,ChannelSegmentType *,Arrays *);

// Floodplain prototypes - fp_flow.cpp
void FloodplainQ(States *,Pars *,Solver *,Arrays *,SGCprams *);
double CalcFPQx(int i,int j,States *,Pars *, Solver *, Arrays *, double * TSptr);
double CalcFPQy(int i,int j,States *,Pars *, Solver *, Arrays *, double * TSptr);
int MaskTest(int m1,int m2);
int MaskTestAcc(int m1);

// Channel prototypes - ch_flow.cpp
void SetChannelStartH(States *Statesptr, Pars *Parptr, Arrays *Arrptr,ChannelSegmentType *ChannelSegments, vector<int> *, int *);
void SetChannelStartHfromQ(States *Statesptr, Pars *Parptr, Arrays *Arrptr,ChannelSegmentType *ChannelSegments,Solver *, vector<int> *, int *);
void CalcChannelStartQ(States *Statesptr, Pars *Parptr, Arrays *Arrptr,ChannelSegmentType *ChannelSegments, vector<int> *, int *);
void ChannelQ(double deltaT, States *, Pars *,Solver *,BoundCs *,ChannelSegmentType *,Arrays *, vector<int> *, int *);
double CalcA(double n,double s,double w,double Q);
double BankQ(int chani,ChannelSegmentType *, Pars *,Arrays *);
double ChannelVol(States *,Pars *, ChannelSegmentType *, Arrays *);
double CalcQ(double n,double s,double w,double h);
double Newton_Raphson(double Ai,double dx,double a0,double a1,double c,Solver *);

// Diffusive channel solver specific functions
void ChannelQ_Diff(double deltaT, States *, Pars *,Solver *,BoundCs *,ChannelSegmentType *,Arrays *, vector<int> *, int *);
void bandec(double **a, int n, int m1, int m2, double **al, int indx[], double &d);
void banbks(double **a, int n, int m1, int m2, double **al, int indx[], double b[]);
void SWAP(double &a,double &b);
void calcF(double *x,double *xn,double *f, double dt, ChannelSegmentType *csp, Pars *Parptr, Arrays *Arrptr, double Qin, int chseg, double WSout, int HoutFREE,  Solver *Solverptr, int low);
void calcJ(double *x,double *xn,double **J, double dt, ChannelSegmentType *csp, int chseg, int HoutFREE);
double norm(double *x,int n);
double CalcEnergySlope(double n,double w,double h,double Q);
//void precond(double **a, int n);

// Acceleration floodplain solver
double CalcFPQxAcc(int i,int j, States *,Pars *, Solver *, Arrays *);
double CalcFPQyAcc(int i,int j, States *,Pars *, Solver *, Arrays *);
double CalcMaxH(Pars *, Arrays *);
void CalcT(Pars *, Solver *, Arrays *);
void UpdateQs(Pars *, Arrays *);

// Sub gid channel prototypes
void SGC_FloodplainQ(States *,Pars *,Solver *,Arrays *,SGCprams *);
void CalcSGCz(Fnames *,States *,Pars *, Arrays *,SGCprams *, int *);
double CalcFPQxSGC(int i,int j, States *,Pars *, Solver *, Arrays *,SGCprams *);
double CalcFPQySGC(int i,int j, States *,Pars *, Solver *, Arrays *,SGCprams *);
void SGC_UpdateH(States *, Pars *, Solver *,BoundCs *,ChannelSegmentType *,Arrays *,SGCprams *);
void SGC_BCs(States *,Pars *, Solver *,BoundCs *, ChannelSegmentType *, Arrays *,SGCprams *);
void SGC_hotstart(States *,Pars *,Solver *,Arrays *);
void SGC_Evaporation(Pars *, Solver *, Arrays *,SGCprams *);
void SGC_Rainfall(Pars *, Solver *, Arrays *); // CCS May 2013
void SGC_Routing(States *,Pars *, Solver *, Arrays *);// CCS May 2013
void CalcSGC_A  (int, double, double, double *, double *, SGCprams *);
double CalcSGC_R  (int, double, double, double, double, double, SGCprams *);
double CalcSGC_UpH(int, double, double, double);
double CalcSGC_UpV(int, double, double, double);
void SGC_wp_prams (SGCprams *);
void CalcSGC_pointFREE(double, double, double, double, double, double, double, double, double, double, int, int, double *, double *, double *, SGCprams *);

// Boundary prototypes - boundary.cpp
void BCs(States *,Pars *, Solver *,BoundCs *, ChannelSegmentType *, Arrays *);
void BoundaryFlux(States *,Pars *,Solver *,BoundCs *,ChannelSegmentType *,Arrays *, vector<ChannelSegmentType> *);
double InterpBC(double *varlist,double t);
double RoeBCy(int edge, int p0, int p1, int pq0, double z0, double z1, double hl, double hr, double hul, double hur, double hvl, double hvr,States *Statesptr,Pars *Parptr,Solver *Solverptr,Arrays *Arrptr);
double RoeBCx(int edge, int p0, int p1, int pq0, double z0, double z1, double hl, double hr, double hul, double hur, double hvl, double hvr,States *Statesptr,Pars *Parptr,Solver *Solverptr,Arrays *Arrptr);

// Optional addon protoypes
// chkpnt.cpp
void ReadCheckpoint(Fnames *, States *, Pars *, Solver *, BoundCs *, ChannelSegmentType *, Arrays *, int *);
void WriteCheckpoint(Fnames *, States *, Pars *, Solver *, BoundCs *, ChannelSegmentType *, Arrays *,  int *);
// infevap.cpp
void FPInfiltration(Pars *, Solver *, Arrays *);
void Evaporation(Pars *, Solver *, Arrays *);
void Rainfall(Pars *, Solver *, Arrays *);
void FlowDirDEM(Pars *, Arrays *, States *,BoundCs *); // Calculate routing intervals and flow directions from DEM for rainfall component CCS 13/03/2012
void Routing(States *, Pars *, Solver *, Arrays *); // Route shallow flows from rainfall CCS 14/03/2012
// por_flow.cpp
double CalcFPQxPor(int i,int j, States *,Pars *, Solver *, Arrays *);
double CalcFPQyPor(int i,int j, States *,Pars *, Solver *, Arrays *);
double PorArea(int t,int j,Pars *, Arrays *);
// weir_flow.cpp
double CalcWeirQx(int i,int j,Pars *, Arrays *, Solver *, States *, SGCprams *);
double CalcWeirQy(int i,int j,Pars *, Arrays *, Solver *, States *, SGCprams *);

// Utility prototypes - util.cpp
double DomainVol(States *, Pars *, ChannelSegmentType *, Arrays *, vector<ChannelSegmentType> *);
void SmoothBanks(Pars *,Solver *, ChannelSegmentType *, Arrays *, vector<ChannelSegmentType> *, int *);
double getmax(double a,double b);
double getmin(double a,double b);
void DryCheck(Pars *, Solver *,Arrays *);
int signR(double a);
void UpdateV(States *, Pars *, Solver *,BoundCs *,ChannelSegmentType *,Arrays *);
void InitFloodplainQ(States *,Pars *,Solver *,Arrays *);
double CalcVirtualGauge(int i, Pars *, Arrays *, Stage *);
void CalcArrayDims(States *, Pars *, Arrays *);

// Ouput prototypes - output.cpp
void fileoutput(Fnames *, States *, Pars *, Arrays *);
void write_regular_output( Fnames *, Solver *, States *, Pars *, Arrays *, SGCprams *);

// MT new general purpose functions
void write_ascfile(char *root, int SaveNumber, char *extension, double *data, double *dem, int outflag, States *, Pars *); // general purpose ascii write routine
void write_binrasterfile(char *root, int SaveNumber, char *extension, double *data, double *dem, int outflag, States *, Pars *); // general purpose ascii write routine
void write_ascfile_SGCf(char *root, int SaveNumber, char *extension, double *data, double *SGCbfH, States *,Pars *); // specific routine for SGC floodplain depth export ascii
void write_binrasterfile_SGCf(char *root, int SaveNumber, char *extension, double *data, double *SGCbfH, States *,Pars *); // specific routine for SGC floodplain depth export binary
void write_profile(char *root, int SaveNumber, char *extension, States *Statesptr,ChannelSegmentType *ChannelSegments, Arrays *Arrptr,Pars *Parptr,  vector<int> *, int *); // write river channel profiles
void debugfileoutput(Fnames *, States *, Pars *, Arrays *); // Debug option file output (currently the modified DEM and the channel and trib masks
void printversion(); // output program version header
int fexist( char *filename ); // check if file exists

// TRENT functions
double maximum(double a,double b,double c);
double CalcFPQxRoe(int i,int j, States *,Pars *, Solver *, Arrays *);
double CalcFPQyRoe(int i,int j, States *,Pars *, Solver *, Arrays *);
void UpdateQsRoe(Pars *, Solver *, Arrays *);
//void UpdateHRoe(States *, Pars *, Solver *,BoundCs *,ChannelSegmentType *,Arrays *);
void CalcTRoe(Pars *, Solver *, Arrays *);


#endif /* LISFLOOD_H */
