#include "initialize.h"
#include "finalize.h"

using namespace std; // CCS

// Global variables that we initiate and declare.
SGCprams *SGCptr;
States *Statesptr;
Pars *Parptr;
Fnames *Fnameptr;
Solver *Solverptr;
BoundCs *BCptr;
Stage *Stageptr;
ChannelSegmentType *CSTypePtr;
Arrays *Arrptr;
vector<int> *RiversIndexVecPtr;
int *RiversIndexPtr;
vector<ChannelSegmentType> *ChannelSegmentsVecPtr;
int *verbose;
char t1[80];
char tmp_sys_com[255]; // temporary string to hold system command

// Fnames ParFp;
Files Fps;
// Solver ParSolver;

double Previous_t;      // previous time channel was calculated

int steadyCount;
int tstep_counter;   // start at -1 so that in first run through we calculate river

double tstep_channel; // channel timestep


extern "C" int init(int argc, const char *argv[]) {

  int result = 0;

  int i,chseg;
  FILE *tmp_fp;
  char strtmp[256];
  double tmp;

  //Instances of Vectors
  ChannelSegmentsVecPtr= new vector<ChannelSegmentType>(); // CCS

  vector<QID7_Store> QID7; //#CCS A temporary store for some ChannelSegments variables that cannot be written to the correct location during LoadRiver().
  vector<QID7_Store> *QID7_Vec_Ptr; // CCS
  QID7_Vec_Ptr=&QID7; // CCS

  RiversIndexVecPtr= new vector<int>(); // CCS

  // DEFINE & DECLARE: Pointers to Structures
  Arrptr = new Arrays();
  Fnameptr = new Fnames();
  Statesptr= new States();
  Parptr = new Pars();
  Solverptr = new Solver();
  BCptr = new BoundCs();
  Stageptr = new Stage();
  SGCptr = new SGCprams();

  // Define initial value for common simulation states (eg. verbose)
  // Pointers for common simulation states
  verbose = new int;
  *verbose = OFF;

  // Define initial value for parameters
  Parptr->dx=10.0;
  Parptr->dy=10.0;
  Parptr->dA=100.0;
  Parptr->tlx=0.0;
  Parptr->tly=0.0;
  Parptr->blx=0.0;
  Parptr->bly=0.0;
  Parptr->FPn=0.06;
  Parptr->SaveInt=1000.0;
  Parptr->SaveTotal=0.0;
  Parptr->MassInt=100.0;
  Parptr->MassTotal=0.0;
  Parptr->SaveNo=0;
  Parptr->op=100.0;
  Parptr->InfilLoss=0.0;
  Parptr->EvapLoss=0.0;
  Parptr->RainLoss=0.0;
  Parptr->InfilTotalLoss=0.0;
  Parptr->EvapTotalLoss=0.0;
  Parptr->RainTotalLoss=0.0;
  Parptr->checkfreq=CHKINTERVAL;  // set default checkpointing file output interval
  Parptr->nextcheck=0.0;
  Parptr->reset_timeinit_time=0;
  Parptr->op_multinum = 0; // default to zero or can cause problems with checkpointing if multipass not used
  Parptr->ch_start_h = 2.0; // default start water depth for channel
  Parptr->steadyQtol = 0.0005; // tolerance for steady-state definition
  Parptr->steadyInt=1800.0; // interval at which to assess steady-state
  Parptr->steadyTotal=0.0;
  Parptr->SGC_p=0.78; // default for sub grid channel exponent
  Parptr->SGC_r=0.12; // default for sub grid channel multiplier (British rivers average Hey and Thorne (1986))
  Parptr->SGCchan_type=1; // defines the type of channel used by the sub_grid model, default is rectangular channel.
  Parptr->SGC_s=2.0; // sub-grid channel parameter used for some of the channel types, parabolic channel default.
  Parptr->SGC_2=0.0; // sub-grid channel parameter used for some of the channel types, meaningless default.
  Parptr->SGC_n=0.035; // sub-grid channel parameter used for some of the channel types, meaningless default.
  Parptr->Routing_Speed=0.1; // CCS default routing speed for shallow rainfall flows 0.1 m/s
  Parptr->RouteInt=0.0; // CCS will be reasigned when FlowDirDEM function is called
  Parptr->RouteSfThresh=0.1; // CCS water surface slope at which routing takes over from shallow water equations when routing is enabled.
  Parptr->SGC_m=1;  // JCN meander coefficient for sub-grid model.
  Parptr->SGC_a=-1; // JCN upstream area for sub-grid model.
  Parptr->min_dx=10.0; // CCS Holds min_dx value (needed for variable dimension lat-long grids)
  Parptr->min_dy=10.0; // CCS Holds min_dy value (needed for variable dimension lat-long grids)
  Parptr->min_dx_dy=10.0; // CCS Holds min of min_dx and min_dy values (needed for variable dimension lat-long grids)


  // Define initial values for boundary conditions
  BCptr->Qin=0.0;
  BCptr->Qout=0.0;
  BCptr->VolInMT=0.0;
  BCptr->VolOutMT=0.0;

  // Define initial values for arrays
  Arrptr->Manningsn=NULL;
  Arrptr->SGCManningsn=NULL;

  // Define initial values for solver settings
  Solverptr->Sim_Time=3600.0;
  Solverptr->InitTstep=10.0;		// Maximum timestep
  Solverptr->Nit=360;
  Solverptr->itCount=0;
  Solverptr->t=0.0;
  Solverptr->g=9.8065500000000;
  Solverptr->divg=(1/(2*Solverptr->g));
  Solverptr->cfl=0.7;
  Solverptr->SolverAccuracy=1e-4;
  Solverptr->dynsw=0; // Switch for full dynamic steady state (1) or diffusive steady state (0)
  Solverptr->DepthThresh=1e-3;
  Solverptr->MomentumThresh=1e-2;
  Solverptr->MaxHflow=10.0;
  Solverptr->Hds=0.0;
  Solverptr->Qerror=0.0;
  Solverptr->Verror=0.0;
  Solverptr->dhlin=0.01;
  Solverptr->htol=1.0;
  Solverptr->Qlimfact=1.0;
  Solverptr->itrn_time=0.0;
  Solverptr->itrn_time_now=0.0;
  Solverptr->ts_multiple=1;  // default to x1 timestep decouple multiple
  Solverptr->SGCtmpTstep =1; // JCN any number
  Solverptr->theta=1.0; // GAMA (for q-centred numerical scheme), 1.0= semi-implicit version (Bates et al 2010);
  Solverptr->fricSolver2D=ON; //GAMA: uses the 2D friction scheme as default

  // Define default values for SimStates instance of States
  Statesptr->diffusive=OFF;	// CCS added default state
  Statesptr->ChannelPresent=OFF;
  Statesptr->TribsPresent=ON;
  Statesptr->NCFS=ON;
  Statesptr->save_depth=ON;
  Statesptr->save_elev=ON;
  Statesptr->out_dir=OFF;
  Statesptr->single_op=OFF;
  Statesptr->multi_op=OFF;
  Statesptr->calc_area=OFF;
  Statesptr->calc_meandepth=OFF;
  Statesptr->calc_volume=OFF;
  Statesptr->save_stages=OFF;
  Statesptr->adaptive_ts=ON;
  Statesptr->qlim=OFF; //TJF: Switch for qlim version, default is OFF
  Statesptr->acceleration=OFF; //PB: Switch for acceleration version, default is OFF
  Statesptr->debugmode=OFF;
  Statesptr->save_Qs=OFF;
  Statesptr->calc_infiltration=OFF;
  Statesptr->call_gzip=OFF;
  Statesptr->alt_ascheader=OFF;
  Statesptr->checkpoint=OFF;
  Statesptr->checkfile=OFF;
  Statesptr->calc_evap=OFF;
  Statesptr->routing=OFF; //CCS: Switch for rainfall routing routine
  Statesptr->rainfall=OFF;
  Statesptr->reset_timeinit=OFF;
  Statesptr->profileoutput=OFF;
  Statesptr->porosity=OFF;
  Statesptr->weirs=OFF;
  Statesptr->save_Ts=OFF;
  Statesptr->save_QLs=OFF;
  Statesptr->startq=OFF;
  Statesptr->logfile=OFF;
  Statesptr->startfile=OFF;
  Statesptr->start_ch_h=OFF;
  Statesptr->comp_out=OFF;
  Statesptr->chainagecalc=ON;
  Statesptr->mint_hk=OFF;
  Statesptr->Roe=OFF;
  Statesptr->killsim=OFF;
  Statesptr->dhoverw=OFF;
  Statesptr->drychecking=ON;
  Statesptr->voutput=OFF;
  Statesptr->steadycheck=OFF;
  Statesptr->hazard=OFF;
  Statesptr->startq2d=OFF;
  Statesptr->Roe_slow=OFF;
  Statesptr->multiplerivers=OFF;
  Statesptr->SGC=OFF;
  Statesptr->SGCbed=OFF;
  Statesptr->SGCcat_area=OFF;
  Statesptr->SGCchangroup=OFF;
  Statesptr->SGCchanprams=OFF;
  Statesptr->SGCbfh_mode=OFF;
  Statesptr->SGCA_mode=OFF;
  Statesptr->binary_out=OFF;
  Statesptr->gsection=OFF;
  Statesptr->binarystartfile=OFF;
  Statesptr->startelev=OFF;
  Statesptr->latlong=OFF;
  Statesptr->dist_routing=OFF;
  Statesptr->SGCvoutput=OFF; // switch for sub-grid channel velocity output

  SGCptr->NSGCprams=0;
  SGCptr->SGCbetahmin=0.2;

  /*default resrootname*/
  strcpy(Fnameptr->resrootname,"out");

  //print version
  printversion();

  // if user only wants to know the version then exit
  for(i=1;i<argc;i++) if(!strcmp(argv[i],"-version")) return(0);

  // check for verbose and debug modes
  for(i=1;i<argc;i++) if(!strcmp(argv[i],"-v")) *verbose=ON;
  for(i=1;i<argc;i++) if(!strcmp(argv[i],"-debug")) Statesptr->debugmode=ON;

  // Assume the last command line argument is the parameter file
  ReadParamFile(argv[argc-1], Fnameptr, Statesptr, Parptr, Solverptr, verbose);

  // allow output folder to be determined by commandline
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-dir"))
                        {
                          sscanf(argv[i+1],"%s",Fnameptr->dirrootname);
                          Statesptr->out_dir=ON;
                          if(*verbose==ON) printf("Output folder set by command line: %s\n",Fnameptr->dirrootname);
                        }
  // TF: allow results root to be determined by commandline
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-resroot"))
                        {
                          sscanf(argv[i+1],"%s",Fnameptr->resrootname);
                          if(*verbose==ON) printf("Results root name reset by command line: %s\n",Fnameptr->resrootname);
                        }

  // PB: switch to acceleration version
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-acceleration"))
                        {
                          Statesptr->acceleration=ON;
                          Statesptr->adaptive_ts=OFF;
                          Statesptr->qlim=OFF;
                          if(*verbose==ON) printf("\nUsing acceleration formulation for floodplain flow\n");
                        }

  // use output folder if requested in parameter file or command line
  if(Statesptr->out_dir==ON)
  {
    if (fexist(Fnameptr->dirrootname)==0) // check if it doesn't exist
    {
      //create output folder
      sprintf(tmp_sys_com,"%s%s","mkdir ",Fnameptr->dirrootname);
      system(tmp_sys_com);
    }
    //update the resroot to include the folder information
    strcpy(strtmp,Fnameptr->resrootname);

#ifdef unix
    sprintf(Fnameptr->resrootname,"%s/%s",Fnameptr->dirrootname,strtmp);
#elif __APPLE__
    sprintf(Fnameptr->resrootname,"%s/%s",Fnameptr->dirrootname,strtmp);
#else
    sprintf(Fnameptr->resrootname,"%s\\%s",Fnameptr->dirrootname,strtmp);
#endif

  }

  // (MT) redirect all console output to logfile if requested
  for(i=1;i<argc;i++) if(!strcmp(argv[i],"-log"))
                      {
                        sprintf(Fnameptr->logfilename,"%s%s",Fnameptr->resrootname,".log");  //Default log filename
                        printf("Redirecting all console output to %s\n\n",Fnameptr->logfilename);
                        printf("Lisflood is running ......\n");
                        freopen( Fnameptr->logfilename, "w", stdout ); // redirect stdout to log file
                        setvbuf(stdout,NULL,_IONBF,0); // set buffer to zero so log file is always up to date
                        printversion(); // output version here as well (so we get on screen and in file)
                        Statesptr->logfile=ON;
                      }

  for(i=1;i<argc;i++) if(!strcmp(argv[i],"-checkpoint")) Statesptr->checkpoint=ON;

  // A different sim_time if requested
  for(i=1;i<argc;i++) if(!strcmp(argv[i],"-simtime")) sscanf(argv[i+1],"%lf",&Solverptr->Sim_Time);

  sprintf(Fnameptr->checkpointfilename,"%s%s",Fnameptr->resrootname,".chkpnt");  //Default checkpoint filename

  for(i=1;i<argc;i++)	if(!strcmp(argv[i],"-loadcheck"))
                      {
                        strcpy(Fnameptr->loadCheckpointFilename,argv[i+1]);
                        Statesptr->checkpoint=ON;
                        Statesptr->checkfile=ON;
                      }
  if(Statesptr->checkpoint==ON && *verbose==ON) printf("Running in checkpointing mode: frequency %lf hours\n",Parptr->checkfreq);

  //override parameter file settings if given on command line
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-nfp")) {
      sscanf(argv[i+1],"%lf",&Parptr->FPn);
      if(*verbose==ON) printf("Floodplain friction reset by command line: %lf\n",Parptr->FPn);
    }
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-weir")) {
      sscanf(argv[i+1],"%s",Fnameptr->weirfilename);
      if(*verbose==ON) printf("Weir file set by command line: %s\n",Fnameptr->weirfilename);
    }
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-inf")) {
      sscanf(argv[i+1],"%lf",&Parptr->InfilRate);
      Statesptr->calc_infiltration=ON;
    }

  // MDW: set a kill time for the simulation
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-kill")) {
      sscanf(argv[i+1],"%lf",&Parptr->killsim_time);
      if(*verbose==ON) printf("\n***** WARNING: Simulation will be killed after %.2lf hours *****\n", Parptr->killsim_time);
      Parptr->killsim_time*=3600; //convert hours to seconds
      Statesptr->killsim=ON;
    }

  // Reset CFL value for acceleration version on the command line (TJF)
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-cfl")){
      if(Statesptr->acceleration==OFF && *verbose==ON) printf("WARNING: CFL value changed on command line but acceleration version off\n");
      sscanf(argv[i+1],"%lf",&Solverptr->cfl);
      if(*verbose==ON) printf("CFL value reset by command line: %lf\n",Solverptr->cfl);
    }

  // Reset theta value for acceleration version on the command line (GAMA)
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-theta")){
      if(Statesptr->acceleration==OFF && *verbose==ON) printf("WARNING: theta value changed on command line but acceleration version off\n");
      sscanf(argv[i+1],"%lf",&Solverptr->theta);
      if(*verbose==ON) printf("theta value reset by command line: %lf\n",Solverptr->theta);
    }

  // Reset dhlin value for adaptive version on the command line (TJF)
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-dhlin")){
      sscanf(argv[i+1],"%lf",&Solverptr->dhlin);
      Statesptr->dhoverw=ON;
      if(*verbose==ON) printf("dhlin value reset by command line: %lf\n",Solverptr->dhlin);
    }

  // Turn on steady-state checking, and check for a specified tolerance (MDW)
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-steady")) Statesptr->steadycheck=ON;
  if(Statesptr->steadycheck==ON) {
    for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-steadytol")) sscanf(argv[i+1],"%lf",&Parptr->steadyQtol); // optional tolerance
    if(*verbose==ON) printf("\nWARNING: simulation will stop on steady-state (tolerance: %.6lf), or after %.1lfs.\n",Parptr->steadyQtol,Solverptr->Sim_Time);
  }

  //code to load in alternative ASCII header for output files
  if(Statesptr->alt_ascheader==ON){
    Parptr->ascheader=new char*[6];//6 lines in the file
    tmp_fp=fopen(Fnameptr->ascheaderfilename,"r");
    for(i=0;i<6;i++) {
      Parptr->ascheader[i]=new char[256];//255 characters per line
      fgets(Parptr->ascheader[i],255,tmp_fp);
    }
    if(*verbose==ON) printf("Using alternative ASCII header for output\n");
    fclose(tmp_fp);
  }
  // Turn on full dynamic channel steady state initial solution
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-dynsw")){
      Solverptr->dynsw=1;
      if(*verbose==ON) printf("Full dynamic steady state channel solver turned on for initial solution\n");
    }

  // get system time and echo for user
  if(*verbose==ON){
    time_t ts = time(0);
    tm timeS = *localtime(&ts);
    printf("\nStart Date: %d/%d/%d \n",timeS.tm_mday,timeS.tm_mon + 1,timeS.tm_year + 1900);
    printf("Start Time: %d:%d:%d \n\n",timeS.tm_hour, timeS.tm_min, timeS.tm_sec);
  }

  LoadDEM(Fnameptr,Statesptr,Parptr,Arrptr,verbose);

  CalcArrayDims(Statesptr,Parptr,Arrptr); // CCS populates dx, dy and dA arrays (calcs correct dimensions if using lat long grid)

  // dhlin value calculated "on the fly" as a function of dx and gradient (0.0002) from Cunge et al. 1980
  if(Statesptr->dhoverw==OFF) Solverptr->dhlin=Parptr->dx*0.0002;

  LoadRiverNetwork(Fnameptr,Statesptr,Parptr, ChannelSegmentsVecPtr,Arrptr,QID7_Vec_Ptr,RiversIndexVecPtr,verbose); // CCS
  //if(Statesptr->ChannelPresent==OFF) ChannelSegments.resize(1); // temp fix to prevent visual studio debuger exiting on the next line (JCN); here some ISSUES with BMI

  CSTypePtr = &(ChannelSegmentsVecPtr->front()); // CCS has to be defined after LoadRiverNetwork has completed.
  RiversIndexPtr=&(RiversIndexVecPtr->front());  // CCS has to be defined after LoadRiverNetwork has completed.

  if(QID7.size()!=0) // CCS If there are any tribs then we need to copy the terms from the temp store to the correct place.
  {
    QID7_Store *QID7Ptr=&QID7[0]; // CCS
    UpdateChannelsVector(Statesptr,CSTypePtr,QID7_Vec_Ptr,QID7Ptr,RiversIndexPtr); // CCS
  }

  //override river file friction if specified on command line
  for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-nch")) {
      sscanf(argv[i+1],"%lf",&tmp);
      if(*verbose==ON) printf("Channel friction reset by command line: %lf\n\n",tmp);
      for(chseg=0;chseg<CSTypePtr->N_Channel_Segments;chseg++) for(i=0;i<CSTypePtr[chseg].chsz;i++) CSTypePtr[chseg].ChanN[i]=tmp;
    }
  if(Statesptr->ChannelPresent==ON) SmoothBanks(Parptr, Solverptr, CSTypePtr, Arrptr, ChannelSegmentsVecPtr, verbose);

  if(Statesptr->SGC==ON) LoadSGC(Fnameptr,Parptr,Arrptr,Statesptr,SGCptr,verbose); // load sub grid channels
  if(Statesptr->SGC==ON && Statesptr->SGCchanprams==ON) LoadSGCChanPrams(Fnameptr,Statesptr,Parptr,SGCptr, verbose); // This loads the parameters for the SGC group information
  if(Statesptr->SGC==ON) CalcSGCz(Fnameptr,Statesptr,Parptr,Arrptr,SGCptr,verbose);

  if(Statesptr->startfile==ON) LoadStart(Fnameptr,Statesptr,Parptr,Arrptr, SGCptr,verbose);
  if(Statesptr->binarystartfile==ON) LoadBinaryStart(Fnameptr,Statesptr,Parptr,Arrptr,verbose);

  LoadBCs(Fnameptr,Statesptr,Parptr,BCptr,Arrptr,verbose);
  LoadBCVar(Fnameptr,Statesptr,Parptr,BCptr,CSTypePtr,Arrptr,ChannelSegmentsVecPtr,verbose);
  LoadManningsn(Fnameptr,Parptr,Arrptr,verbose);
  LoadSGCManningsn(Fnameptr,Parptr,Arrptr,verbose);
  LoadPor(Fnameptr,Statesptr,Parptr,Arrptr,verbose);
  LoadWeir(Fnameptr,Statesptr,Parptr,Arrptr,verbose);
  if(Statesptr->calc_evap==ON) LoadEvap(Fnameptr,Arrptr,verbose);
  if(Statesptr->rainfall==ON) LoadRain(Fnameptr,Arrptr,verbose);
  if(Statesptr->save_stages==ON) LoadStages(Fnameptr, Statesptr, Parptr, Stageptr, verbose);
  if(Statesptr->gsection==ON) LoadGauges(Fnameptr, Statesptr, Parptr, Stageptr, verbose);

  if (Statesptr->routing==ON) // Call FlowDirDEM to generate flow direction map from DEM before main loop CCS
  {
    FlowDirDEM(Parptr, Arrptr, Statesptr, BCptr);
    if(*verbose==ON) printf("Flow direction map generated from DEM\n\n");
  }

  // apply different starting methods for channel
  if(Statesptr->ChannelPresent==ON)
  {
    // calc initial steady state flows down channel
    CalcChannelStartQ(Statesptr, Parptr, Arrptr, CSTypePtr, RiversIndexVecPtr, RiversIndexPtr);

    if(Statesptr->startfile == ON)
    {
      // start file is specified. Do nothing, as starting H values for channel already read in from the startfile.
    }
    else if(Statesptr->startq == ON)
    {
      // Kinematic: Uses the kinematic initial solution to calculate H from Q
      // Diffusive: Uses diffusive steady state initial solution (default) or can use full dynamic steady state
      // initial if turned on using -dynsw on command line or "ch_dynamic" in the parameter file

      // use the flows to calculate a starting H
      SetChannelStartHfromQ(Statesptr, Parptr, Arrptr, CSTypePtr,Solverptr,RiversIndexVecPtr, RiversIndexPtr);
    }
    else
    {
      // set channel start H to default or user defined H
      SetChannelStartH(Statesptr, Parptr, Arrptr, CSTypePtr, RiversIndexVecPtr, RiversIndexPtr);
    }
  }
  // apply hot starting methods to SGC model
  if(Statesptr->startq == ON && Statesptr->SGC == ON)
  {
    SGC_hotstart(Statesptr,Parptr,Solverptr,Arrptr);
    if(*verbose==ON) printf("\nStartq for SGC model implemented\n");
  }

  if(*verbose==ON) if(Statesptr->calc_infiltration==ON) printf("Floodplain infiltration set at: %.10f ms-1\n\n",Parptr->InfilRate);

  //get multiple overpass timings from file
  if(Statesptr->multi_op==ON) {
    tmp_fp=fopen(Fnameptr->opfilename,"r");
    if(tmp_fp!=NULL)
    {
      fscanf(tmp_fp,"%i",&Parptr->op_multinum);
      if(*verbose==ON) printf("\nMultiple overpass files to be output: %d\n",Parptr->op_multinum);
      Parptr->op_multisteps=new double[Parptr->op_multinum];
      Parptr->op_multiswitch=new int[Parptr->op_multinum];
      for(i=0;i<Parptr->op_multinum;i++){
        if(fscanf(tmp_fp,"%lf",&Parptr->op_multisteps[i])!=1) // read in value and check if one value read in successfully
        {
          printf("\nWARNING: overpass file read error at line %i\n",i+1);
          Parptr->op_multinum=i; // reset to number of values actually read in
          break;
        }
        Parptr->op_multiswitch[i]=0;
        if(*verbose==ON) printf("Overpass %d at %lf seconds\n",i,Parptr->op_multisteps[i]);
      }
      fclose(tmp_fp);
    } else {
      Statesptr->multi_op=OFF;
      if(*verbose==ON) printf("\nUnable to open multiple overpass output file: %s\n",Fnameptr->opfilename);
    }
  }

  //Load checkpointed data if this job has been restarted
  if(Statesptr->checkpoint==ON) {
    ReadCheckpoint(Fnameptr,Statesptr,Parptr,Solverptr,BCptr,CSTypePtr,Arrptr, verbose);
    if(*verbose==ON) printf(" - checkpoint output file: %s\n",Fnameptr->checkpointfilename);
  }

  //mass balance
  sprintf(t1,"%s%s",Fnameptr->resrootname,".mass");
  if(Statesptr->checkpoint==ON && Solverptr->t>0) { //if this is a checkpointed job, we only need to amend the .mass file
    Fps.mass_fp=fopen(t1,"a");
  } else {
    Fps.mass_fp=fopen(t1,"w");
  }
  if(Fps.mass_fp!=NULL)
  {
    if(Solverptr->t==0) fprintf(Fps.mass_fp,"Time         Tstep      MinTstep   NumTsteps    Area         Vol         Qin         Hds        Qout          Qerror       Verror       Rain-Inf+Evap\n");
    else
    {
      // make a note in the mass file that this is a restart point - user can then edit the overlap out if they want a continuous mass file record.
      fprintf(Fps.mass_fp,"####################################################### Checkpoint restart ########################################################\n");
      fprintf(Fps.mass_fp,"Time         Tstep      MinTstep   NumTsteps    Area         Vol         Qin         Hds        Qout          Qerror       Verror       Rain-Inf+Evap\n");
      fflush(Fps.mass_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
    }
  }
  else
  {
    if(*verbose==ON)
    {
      printf("Unable to open mass balance file: %s",t1);
      result = 1;
      return result;
    }
  }

  //stage output file
  if(Statesptr->save_stages==ON) {
    sprintf(t1,"%s%s",Fnameptr->resrootname,".stage");
    if(Statesptr->checkpoint==ON && Solverptr->t>0) { //if this is a checkpointed job, we only need to amend the .stage file
      Fps.stage_fp=fopen(t1,"a");
    } else {
      Fps.stage_fp=fopen(t1,"w");
    }
    if(Fps.stage_fp!=NULL)
    {
      if (Solverptr->t==0.0 || Statesptr->checkpoint==OFF)
      {
        fprintf(Fps.stage_fp,"Stage output, depth (m). Stage locations from: %s\n\n",Fnameptr->stagefilename);
        fprintf(Fps.stage_fp,"Stage information (stage,x,y,elev):\n");
        for(i=0;i<Stageptr->Nstages;i++)
        {
          if (Statesptr->SGC==ON && Arrptr->SGCwidth[Stageptr->stage_grid_x[i]+Stageptr->stage_grid_y[i]*Parptr->xsz] > 0) // if a SUB GRID channel is present export the channel bed elevation)
          {
            if(Stageptr->stage_check[i]==1) fprintf(Fps.stage_fp,"%d\t%.4f\t%.4f\t%.4f\n",i+1,Stageptr->stage_loc_x[i],Stageptr->stage_loc_y[i],Arrptr->SGCz[Stageptr->stage_grid_x[i]+Stageptr->stage_grid_y[i]*Parptr->xsz]);
            else fprintf(Fps.stage_fp,"%d\t%.4f\t%.4f\tn/a\n",i+1,Stageptr->stage_loc_x[i],Stageptr->stage_loc_y[i]);
          }
          else
          {
            if(Stageptr->stage_check[i]==1) fprintf(Fps.stage_fp,"%d\t%.4f\t%.4f\t%.4f\n",i+1,Stageptr->stage_loc_x[i],Stageptr->stage_loc_y[i],Arrptr->DEM[Stageptr->stage_grid_x[i]+Stageptr->stage_grid_y[i]*Parptr->xsz]);
            else fprintf(Fps.stage_fp,"%d\t%.4f\t%.4f\tn/a\n",i+1,Stageptr->stage_loc_x[i],Stageptr->stage_loc_y[i]);
          }
        }
        fprintf(Fps.stage_fp,"\nOutput, depths:\n");
        fprintf(Fps.stage_fp,"Time; stages 1 to %d\n",Stageptr->Nstages);
      }
      else
      {
        fprintf(Fps.stage_fp,"####################################################### Checkpoint restart ########################################################\n");
        fflush(Fps.stage_fp);
      }
    }
    else
    {
      if(*verbose==ON) printf("Unable to open stage output file: %s",t1);
      Statesptr->save_stages=OFF;
    }

  }
  //velocity output file
  if(Statesptr->save_stages==ON && Statesptr->voutput==ON)
  {
    sprintf(t1,"%s%s",Fnameptr->resrootname,".velocity");
    if(Statesptr->checkpoint==ON && Solverptr->t>0) { //if this is a checkpointed job, we only need to amend the .stage file
      Fps.vel_fp=fopen(t1,"a");
    } else {
      Fps.vel_fp=fopen(t1,"w");
    }
    if(Fps.vel_fp!=NULL){
      if(Solverptr->t==0) {
        fprintf(Fps.vel_fp,"Velocity output, velocity (ms-1). Velocity locations from: %s\n\n",Fnameptr->stagefilename);
        fprintf(Fps.vel_fp,"Stage information (stage,x,y,elev):\n");
        for(i=0;i<Stageptr->Nstages;i++){
          if(Stageptr->stage_check[i]==1) fprintf(Fps.vel_fp,"%d\t%.4f\t%.4f\t%.4f\n",i+1,Stageptr->stage_loc_x[i],Stageptr->stage_loc_y[i],Arrptr->DEM[Stageptr->stage_grid_x[i]+Stageptr->stage_grid_y[i]*Parptr->xsz]);
          else fprintf(Fps.vel_fp,"%d\t%.4f\t%.4f\tn/a\n",i+1,Stageptr->stage_loc_x[i],Stageptr->stage_loc_y[i]);
        }
        fprintf(Fps.vel_fp,"\nOutput, depths:\n");
        fprintf(Fps.vel_fp,"Time; velocities 1 to %d\n",Stageptr->Nstages);
      }
      else {
        fprintf(Fps.vel_fp,"####################################################### Checkpoint restart ########################################################\n");
        fflush(Fps.vel_fp);
      }
    } else {
      if(*verbose==ON) printf("Unable to open velocity output file: %s",t1);
      Statesptr->save_stages=OFF;
    }
  }

  //velocity output file
  if(Statesptr->gsection==ON)
  {
    sprintf(t1,"%s%s",Fnameptr->resrootname,".discharge");
    if(Statesptr->checkpoint==ON && Solverptr->t>0) { //if this is a checkpointed job, we only need to amend the .stage file
      Fps.gau_fp=fopen(t1,"a");
    } else {
      Fps.gau_fp=fopen(t1,"w");
    }
    if(Fps.gau_fp!=NULL){
      if(Solverptr->t==0) {
        fprintf(Fps.gau_fp,"Discharge output, discharge (m3s-1). Discharge locations from: %s\n\n",Fnameptr->gaugefilename);
        fprintf(Fps.gau_fp,"Time; discharge 1 to %d\n",Stageptr->Ngauges);
      }
      else {
        fprintf(Fps.gau_fp,"####################################################### Checkpoint restart ########################################################\n");
        fflush(Fps.gau_fp);
      }
    } else {
      if(*verbose==ON) printf("Unable to open discharge output file: %s",t1);
      Statesptr->gsection=OFF;
    }
  }

  //find out if we are going to compress output on the fly
  for(i=1;i<argc;i++) {
    if(!strcmp(argv[i],"-gzip")) {
      Statesptr->call_gzip=ON;
      if(*verbose==ON) printf("\nOutput will be compressed using Gzip\n");
    }
  }

  // output debug files (used DEM, channel mask seg mask) if required
  if(Statesptr->debugmode==ON) debugfileoutput(Fnameptr,Statesptr,Parptr,Arrptr);

  //start simulation
  time(&Solverptr->time_start);

  return result;
}

extern "C" void init_iterateq()
{
  Solverptr->vol1=DomainVol(Statesptr,Parptr,&(ChannelSegmentsVecPtr->front()),Arrptr,ChannelSegmentsVecPtr);

  if(*verbose==ON)
  {
    printf("\nStarting time steps: ");
    fflush(stdout);
  }
  Solverptr->itrn_time_now=Solverptr->itrn_time;

  // Populating Tstep variables prior to start of simulation
  if(Statesptr->adaptive_ts==ON)
  {
    if(Solverptr->t==0)
    {
      Solverptr->Tstep=Solverptr->InitTstep;
      Solverptr->MinTstep=Solverptr->InitTstep;
    }
    if(*verbose==ON) printf("adaptive mode\n\n");
    fflush(stdout);
  }
  else if(Statesptr->acceleration==ON)
  {
    if(Solverptr->t==0)
    {
      Solverptr->Tstep=Solverptr->InitTstep;
      Solverptr->MinTstep=Solverptr->InitTstep;
    }
    if(*verbose==ON) printf("acceleration mode\n\n");
    fflush(stdout);
  }
  else if(Statesptr->Roe==ON)
  {
    if(Solverptr->t==0)
    {
      Solverptr->Tstep=Solverptr->InitTstep;
      Solverptr->MinTstep=Solverptr->InitTstep;
    }
    if(*verbose==ON) printf("Roe mode\n\n");
    fflush(stdout);
  }
  else
  {
    if(Solverptr->t==0)
    {
      Solverptr->Tstep=Solverptr->InitTstep;
      Solverptr->MinTstep=Solverptr->InitTstep;
    }
    if(*verbose==ON) printf("non-adaptive mode\n\n");
    fflush(stdout);
  }
  if (Statesptr->SGC==ON)
  {
    // because the SGC model calculates the time step in UpdateH rather than during calcFPflow it needs to initalise
    // SGCtmpTstep, which would usually be calculated in UpdateH
    Solverptr->Tstep=Solverptr->InitTstep;
    CalcT(Parptr, Solverptr, Arrptr);
    Solverptr->SGCtmpTstep = Solverptr->Tstep;
    if(*verbose==ON) printf("SGC mode\n\n");
    fflush(stdout);
  }
  // set previous time to one timestep backwards so that first river calcs uses Solverptr->Tstep for 1st iteration
  // this is because of the way the timestep is calculated as the time difference from the last time the river was run
  //Previous_t = -Solverptr->Tstep;
  Previous_t = Solverptr->t-Solverptr->Tstep;

  steadyCount=0;
  tstep_counter=-1;   // start at -1 so that in first run through we calculate river

  tstep_channel=0; // channel timestep

}


//---------------------------------------------------------------------------
void printversion()
// printout header with program and version number
{
  printf("***************************\n");
  printf(" LISFLOOD-FP version %d.%d.%d\n",VersionMajor,VersionMinor,VersionInc);
  printf("***************************\n\n");
  printf(" This is a modified version 5.9 including Basic Model Interface (BMI) compatibility\n\n");
  printf("***************************\n\n");
  printf(" The model needs therefore be executed accordingly!\n\n");
  printf("***************************\n\n");
}
