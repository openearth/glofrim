/*
*****************************************************************************
LOAD PARAMETER FILE PARAMETERS
---------------------

Load the filenames and options from the parameter file (.par) and defines the
simulation states. Detailed explanation of model parameters, files and simulation
states found in the manual. TJF.

*****************************************************************************
*/

#include "lisflood.h"

//-----------------------------------------------------------------------------
// LOAD COMMAND LINE PARAMETERS INTO GLOBAL VARIABLES
void ReadParamFile(const char *fname, Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, int *verbose)
{
  char buffer[80];
  FILE *par_fp;
  par_fp=fopen(fname,"r");

  if(par_fp==NULL)
  {
    fprintf(stderr,"Parameter file not found. Aborting.\n");
    exit(0);
  }

  if(*verbose==ON) printf("Loading parameters... \n");

  while(!feof(par_fp))
  {
    fscanf(par_fp,"%s",buffer); // get next line of parameter file

    if(!strcmp(buffer,"DEMfile")) fscanf(par_fp,"%s",Fnameptr->demfilename);
    if(!strcmp(buffer,"startfile"))
    {
      fscanf(par_fp,"%s",Fnameptr->startfilename);
      Statesptr->startfile = ON;
	  Statesptr->binarystartfile = OFF;
    }
	if(!strcmp(buffer,"binarystartfile"))
    {
      fscanf(par_fp,"%s",Fnameptr->startfilename);
      Statesptr->startfile = OFF;
	  Statesptr->binarystartfile = ON;
	  if(*verbose==ON) printf("Using binary startfile\n");
    }
	if(!strcmp(buffer,"startelev"))
    {
		Statesptr->startelev = ON;
		if(*verbose==ON) printf("Using water surface elevations for restart\n");
	}
    if(!strcmp(buffer,"resroot")) fscanf(par_fp,"%s",Fnameptr->resrootname);
    if(!strcmp(buffer,"dirroot"))
    {
      fscanf(par_fp,"%s",Fnameptr->dirrootname);
      Statesptr->out_dir = ON;
    }
    if(!strcmp(buffer,"fpfric")) fscanf(par_fp,"%lf",&Parptr->FPn);
    if(!strcmp(buffer,"startq")) Statesptr->startq=ON;
    if(!strcmp(buffer,"ch_start_h"))
    {
      fscanf(par_fp,"%lf",&Parptr->ch_start_h);
      Statesptr->start_ch_h = ON;
    }

    if(!strcmp(buffer,"sim_time")) fscanf(par_fp,"%lf",&Solverptr->Sim_Time);

    if(!strcmp(buffer,"initial_tstep")) fscanf(par_fp,"%lf",&Solverptr->InitTstep);
    if(!strcmp(buffer,"massint")) fscanf(par_fp,"%lf",&Parptr->MassInt);
    if(!strcmp(buffer,"saveint")) fscanf(par_fp,"%lf",&Parptr->SaveInt);
    if(!strcmp(buffer,"htol")) fscanf(par_fp,"%lf",&Solverptr->htol);
    if(!strcmp(buffer,"qlimfact"))
    {
       fscanf(par_fp,"%lf",&Solverptr->Qlimfact);
       if(*verbose==ON) printf("Relaxing Q limit by factor x%f \n\n",Solverptr->Qlimfact);
    }

    if(!strcmp(buffer,"ts_multiple"))
    {
       fscanf(par_fp,"%d",&Solverptr->ts_multiple);
       if(*verbose==ON) printf("Running channel at x%d multiples of floodplain timestep\n\n",Solverptr->ts_multiple);
    }

    if(!strcmp(buffer,"overpass"))
    {
      fscanf(par_fp,"%lf",&Parptr->op);
      Statesptr->single_op=ON;
    }
    if(!strcmp(buffer,"overpassfile"))
    {
      fscanf(par_fp,"%s",Fnameptr->opfilename);
      Statesptr->multi_op=ON;
    }
    if(!strcmp(buffer,"Qfile")) fscanf(par_fp,"%s",Fnameptr->qfilename);
    if(!strcmp(buffer,"manningfile")) fscanf(par_fp,"%s",Fnameptr->nfilename);
	if(!strcmp(buffer,"SGCmanningfile"))
	{
		fscanf(par_fp,"%s",Fnameptr->SGCnfilename);
		if(*verbose==ON) printf("SGC Manning file %s\n",Fnameptr->SGCnfilename);
	}
    if(!strcmp(buffer,"riverfile")) fscanf(par_fp,"%s",Fnameptr->rivername);
    if(!strcmp(buffer,"bcifile")) fscanf(par_fp,"%s",Fnameptr->bcifilename);
    if(!strcmp(buffer,"bdyfile")) fscanf(par_fp,"%s",Fnameptr->bdyfilename);
    if(!strcmp(buffer,"weirfile")) fscanf(par_fp,"%s",Fnameptr->weirfilename);

	if(!strcmp(buffer,"ch_dynamic")) Solverptr->dynsw=1; // Full dynamic steady state channel solver used for initial condition

    if(!strcmp(buffer,"porfile"))
    {
      fscanf(par_fp,"%s",Fnameptr->porfilename);
      Statesptr->porosity=ON;
    }

    if(!strcmp(buffer,"depthoff")) Statesptr->save_depth=OFF;
    if(!strcmp(buffer,"elevoff")) Statesptr->save_elev=OFF;

	if(!strcmp(buffer,"adaptoff"))
	{
		Statesptr->adaptive_ts=OFF;
		Statesptr->qlim=ON;
		if(*verbose==ON) printf("Using Qlim formulation for floodplain flow\n");

	}
	if(!strcmp(buffer,"acceleration"))
	{
		Statesptr->acceleration=ON;
		Statesptr->adaptive_ts=OFF;
		Statesptr->qlim=OFF;
		if(*verbose==ON) printf("Using acceleration formulation for floodplain flow\n");
	}

	// Don't combine adaptoff and acceleration keywords. Program aborts.
	if((!strcmp(buffer,"acceleration") && Statesptr->qlim==ON) || (!strcmp(buffer,"adaptoff") && Statesptr->acceleration==ON))
	{
		if(*verbose==ON) printf("\nIncompatible timestep functions chosen (adaptoff and acceleration) - Check parameter file. Aborting\n");
		exit (1);
	}
    if(!strcmp(buffer,"chainageoff")) Statesptr->chainagecalc=OFF;
    if(!strcmp(buffer,"qoutput")) Statesptr->save_Qs=ON;
    if(!strcmp(buffer,"voutput")) Statesptr->voutput=ON;
	if(!strcmp(buffer,"hazard"))
	{
		Statesptr->voutput=ON;
		Statesptr->hazard=ON;
		if(*verbose==ON) printf("\nHazard mode on\n");
	}
    if(!strcmp(buffer,"profiles")) Statesptr->profileoutput=ON;
    if(!strcmp(buffer,"debug")) Statesptr->debugmode=ON;
    if(!strcmp(buffer,"comp_out")) Statesptr->comp_out=ON;
    if(!strcmp(buffer,"binary_out"))
    {
	    Statesptr->binary_out=ON;
	    printf("Using binary raster output\n");
    }

    // JN: added to request maxH, maxHtm totalHtm and initHtm be calulated at the mass interval
	if(!strcmp(buffer,"mint_hk")) Statesptr->mint_hk=ON;
	// JN/IV: added to turn on Roe solever, don't combine with adaptoff
	if(!strcmp(buffer,"Roe"))
	{
		Statesptr->Roe=ON;
		Statesptr->Roe_slow=OFF;
		Statesptr->adaptive_ts=OFF;
		Statesptr->acceleration=OFF;
		Statesptr->qlim=OFF;
		if(*verbose==ON) printf("\nUsing Roe formulation for floodplain flow\n");
	}
		if(!strcmp(buffer,"Roe_slow"))
	{
		Statesptr->Roe_slow=ON; // uses ghost cell version of Roe solver.
		Statesptr->Roe=ON;
		Statesptr->adaptive_ts=OFF;
		Statesptr->acceleration=OFF;
		Statesptr->qlim=OFF;
		if(*verbose==ON) printf("\nUsing Roe formulation for floodplain flow\n");
	}
	if((!strcmp(buffer,"Roe") && Statesptr->qlim==ON) || (!strcmp(buffer,"adaptoff") && Statesptr->Roe==ON) || (!strcmp(buffer,"acceleration") && Statesptr->Roe==ON))
	{
		if(*verbose==ON) printf("\nIncompatible timestep functions chosen (adaptoff/qlim/acceleration and Roe) - Check parameter file. Aborting\n");
		exit (1);
	}

    // MT: added to output adaptive timestep & qlimits & diffusive channel solver
    if(!strcmp(buffer,"toutput")) Statesptr->save_Ts=ON;
    if(!strcmp(buffer,"qloutput")) Statesptr->save_QLs=ON;
    if(!strcmp(buffer,"diffusive"))
    {
      Statesptr->diffusive=ON;
      if(*verbose==ON) printf("\nUsing diffusive solver for channel flow\n");
    }

    if(!strcmp(buffer,"ascheader"))
    {
      fscanf(par_fp,"%s",Fnameptr->ascheaderfilename);
      Statesptr->alt_ascheader=ON;
    }

    if(!strcmp(buffer,"calcarea")) Statesptr->calc_area=ON;
    if(!strcmp(buffer,"calcmeandepth")) Statesptr->calc_meandepth=ON;
    if(!strcmp(buffer,"calcvolume")) Statesptr->calc_volume=ON;

    if(!strcmp(buffer,"stagefile"))
    {
      fscanf(par_fp,"%s",Fnameptr->stagefilename);
      Statesptr->save_stages=ON;
    }
	if(!strcmp(buffer,"gaugefile"))
    {
      fscanf(par_fp,"%s",Fnameptr->gaugefilename);
      Statesptr->gsection=ON;
    }
    if(!strcmp(buffer,"infiltration"))
    {
      Statesptr->calc_infiltration=ON;
      fscanf(par_fp,"%lf",&Parptr->InfilRate);
    }

    if(!strcmp(buffer,"evaporation"))
    {
      Statesptr->calc_evap=ON;
      fscanf(par_fp,"%s",Fnameptr->evapfilename);
    }

	if(!strcmp(buffer,"rainfall")) // Enable rainfall CCS
    {
      Statesptr->rainfall=ON;
      fscanf(par_fp,"%s",Fnameptr->rainfilename);
    }

	if(!strcmp(buffer,"routing")) // Enable routing scheme CCS
    {
      Statesptr->routing=ON;
    }

	if(!strcmp(buffer,"rainfallrouting")) // Legacy routing keyword CCS
    {
      Statesptr->rainfall=ON;
	  Statesptr->routing=ON;
	  fscanf(par_fp,"%s",Fnameptr->rainfilename);
    }

	if(!strcmp(buffer,"routingspeed")) // Set speed at which shallow water (<depththresh) is routed across DEM CCS 14/03/2012
    {
	  Statesptr->routing=ON;
	  fscanf(par_fp,"%lf",&Parptr->Routing_Speed);
    }

	if(!strcmp(buffer,"routesfthresh")) // Set slope where routing replaces shallow water eqn  CCS 14/03/2012
    {
		Statesptr->routing=ON;
		fscanf(par_fp,"%lf",&Parptr->RouteSfThresh);
    }

    if(!strcmp(buffer,"checkpoint"))
    {
      Statesptr->checkpoint=ON;
      fscanf(par_fp,"%lf",&Parptr->checkfreq);
      if(Parptr->checkfreq<=0.0) Parptr->checkfreq=CHKINTERVAL; // if interval is less than zero set to default
    }
    if(!strcmp(buffer,"loadcheck"))
    {
      Statesptr->checkpoint=ON;
      Statesptr->checkfile=ON;
      fscanf(par_fp,"%s",Fnameptr->loadCheckpointFilename);
    }
    if(!strcmp(buffer,"resettimeinit"))
    {
      Statesptr->reset_timeinit=ON;
      fscanf(par_fp,"%lf",&Parptr->reset_timeinit_time);
      if(*verbose==ON) printf("\n Time of initial inundation will be reset at %lf seconds\n",Parptr->reset_timeinit_time);
    }
	// Reset the cfl condition for acceleration version - reduce to increase stability
	if(!strcmp(buffer,"cfl")) {
		fscanf(par_fp,"%lf",&Solverptr->cfl);
		if(*verbose==ON) printf("cfl changed to %f \n",Solverptr->cfl);
	}

	// Reset the theta parameter for acceleration version - reduce to increase numerical diffusion
	if(!strcmp(buffer,"theta")) {
		fscanf(par_fp,"%lf",&Solverptr->theta);
		if(*verbose==ON) printf("theta changed to %f \n",Solverptr->theta);
	}

	// Uses the 1D (old) version of the friction term
	if(!strcmp(buffer,"1Dfriction")) {
		Solverptr->fricSolver2D=OFF;
		if(*verbose==ON) printf("Using the 1D version of the friction term\n");
	}

	// Reset the dhlin condition for adaptive version - reduce to increase stability/ increase to reduce computation time
	// Overwrites the standard value which is set as gradient = 0.0002 (i.e. dhlin is a function of dx)
	if(!strcmp(buffer,"dhlin")) {
		Statesptr->dhoverw=ON;
		fscanf(par_fp,"%lf",&Solverptr->dhlin);
	}
	// Options to change depth and momentum thresholds defaults used if not set
	if(!strcmp(buffer,"depththresh"))
	{
		fscanf(par_fp,"%lf",&Solverptr->DepthThresh);
		if(*verbose==ON) printf("Depth threshold changed to %f \n",Solverptr->DepthThresh);
	}
    if(!strcmp(buffer,"momentumthresh"))
	{
		fscanf(par_fp,"%lf",&Solverptr->MomentumThresh);
		if(*verbose==ON) printf("Momentum threshold changed to %f \n",Solverptr->MomentumThresh);
	}
	if(!strcmp(buffer,"drycheckoff")) { // turns DryCheck off use at own risk!
		Statesptr->drychecking=OFF;
		if(*verbose==ON) printf("DryCheck is off\n");
	}
	if(!strcmp(buffer,"multiriverfile")) // CCS read in multiriver index file
	{
		fscanf(par_fp,"%s",Fnameptr->multiriverfilename);
		Statesptr->multiplerivers=ON;
		if(*verbose==ON) printf("\nMultiple river mode selected\n");
	}
	if(!strcmp(buffer,"SGCp"))
	{
		Statesptr->SGC=ON;
		fscanf(par_fp,"%lf",&Parptr->SGC_p);
		if(*verbose==ON) printf("SGC exponent changed to %.5f \n",Parptr->SGC_p);
	}
	if(!strcmp(buffer,"SGCr"))
	{
		Statesptr->SGC=ON;
		fscanf(par_fp,"%lf",&Parptr->SGC_r);
		if(*verbose==ON) printf("SGC multiplier changed to %.5f \n",Parptr->SGC_r);
	}
	if(!strcmp(buffer,"SGCm"))
	{
		fscanf(par_fp,"%lf",&Parptr->SGC_m);
		if(*verbose==ON) printf("SGC meander coefficient changed to %.5f \n",Parptr->SGC_m);
	}
	if(!strcmp(buffer,"SGCchan"))
	{
		fscanf(par_fp,"%i",&Parptr->SGCchan_type);
		if(*verbose==ON) printf("SGC channel changed to %i \n",Parptr->SGCchan_type);
	}
	if(!strcmp(buffer,"SGCs"))
	{
		fscanf(par_fp,"%lf",&Parptr->SGC_s);
		if(*verbose==ON) printf("SGCs changed to %.5f \n",Parptr->SGC_s);
	}
	if(!strcmp(buffer,"SGC2"))
	{
		fscanf(par_fp,"%lf",&Parptr->SGC_2);
		if(*verbose==ON) printf("SGC2 changed to %.5f \n",Parptr->SGC_2);
	}
	if(!strcmp(buffer,"SGCn"))
	{
		fscanf(par_fp,"%lf",&Parptr->SGC_n);
		if(*verbose==ON) printf("SGCn changed to %.5f \n",Parptr->SGC_n);
	}
	if(!strcmp(buffer,"SGCa"))
	{
		fscanf(par_fp,"%lf",&Parptr->SGC_a);
		if(*verbose==ON) printf("SGCa changed to %.5f \n",Parptr->SGC_a);
	}
	if(!strcmp(buffer,"SGCbfh_mode"))
	{
		Statesptr->SGCbfh_mode=ON; // Turns on mode where parameter p is used as the channel bank full depth
		if(*verbose==ON) printf("Using sub-grid parameter p as bank full depth\n");
	}
	if(!strcmp(buffer,"SGCA_mode"))
	{
		Statesptr->SGCbfh_mode=ON; // Turns on mode where parameter p is used as the channel bank full area
		if(*verbose==ON) printf("Using sub-grid parameter p as bank full area\n");
	}
	if(!strcmp(buffer,"SGCvoutput"))
	{
		Statesptr->SGCvoutput=ON; // Turns on mode where parameter p is used as the channel bank full area
		if(*verbose==ON) printf("Output sub-grid channel velocity\n");
	}
	if(!strcmp(buffer,"SGCwidth"))
	{
		Statesptr->SGC=ON;
		fscanf(par_fp,"%s",Fnameptr->SGCwidthfilename);
		Statesptr->acceleration=ON; // sug-grid channels only work with acceleration model so ensure its in this model.
		Statesptr->adaptive_ts=OFF;
		Statesptr->qlim=OFF;
		if(*verbose==ON) printf("Using sub-grid channels and acceleration formulation\n");
	}
	if(!strcmp(buffer,"SGCbed"))
	{
		Statesptr->SGCbed=ON;
		fscanf(par_fp,"%s",Fnameptr->SGCbedfilename);
		if(*verbose==ON) printf("SGC bed elevation read\n");
	}
	if(!strcmp(buffer,"SGCcat_area"))
	{
		Statesptr->SGCcat_area=ON;
		fscanf(par_fp,"%s",Fnameptr->SGCcat_areafilename);
		if(*verbose==ON) printf("SGC catchment area read\n");
	}
	if(!strcmp(buffer,"SGCchangroup"))
	{
		Statesptr->SGCchangroup=ON;
		fscanf(par_fp,"%s",Fnameptr->SGCchangroupfilename);
		if(*verbose==ON) printf("SGC channel type\n");
	}
	if(!strcmp(buffer,"SGCchanprams"))
	{
		Statesptr->SGCchanprams=ON;
		fscanf(par_fp,"%s",Fnameptr->SGCchanpramsfilename);
		if(*verbose==ON) printf("SGC channel parameters\n");
	}
	if(!strcmp(buffer,"SGCbank")) fscanf(par_fp,"%s",Fnameptr->SGCbankfilename);
	if(!strcmp(buffer,"tstart"))
	{
		 fscanf(par_fp,"%lf",&Solverptr->t);
		 Parptr->SaveTotal = Solverptr->t;
		 if(*verbose==ON) printf("Simulation start time changed to %f \n",Solverptr->t);
	}
	if(!strcmp(buffer,"gravity"))
	{
		 fscanf(par_fp,"%lf",&Solverptr->g);
		 if(*verbose==ON) printf("Simulation gravity changed to %f \n",Solverptr->g);
	}
	if(!strcmp(buffer,"latlong"))
	{
		Statesptr->latlong=ON;
		if(*verbose==ON) printf("\nLat-Long coordinate system on.\n");
	}
	if(!strcmp(buffer,"dist_routing"))
	{
		Statesptr->dist_routing=ON;
		if(*verbose==ON) printf("\nUsing slope dependent routing velocity.\n");
	}
  }

  fclose(par_fp);

  if(*verbose==ON) printf("Reading parameters done.\n\n");

  // do some basic checks to warn user of option conflicts
  if (Statesptr->start_ch_h == ON && Statesptr->startq == ON && *verbose==ON) printf("\nWARNING: startq option overides ch_start_h values\n\n");
  if (Statesptr->startfile == ON && Statesptr->startq == ON && *verbose==ON) printf("\nWARNING: startfile H values overide startq values\n\n");
  if (Statesptr->startfile == ON && Statesptr->start_ch_h == ON && *verbose==ON) printf("\nWARNING: startfile H values overide ch_start_h values\n\n");
  if (Statesptr->routing == ON && Statesptr->acceleration == OFF && Statesptr->SGC == OFF && *verbose==ON) // CCS disable routing if not being used with inertial or SG model:
  {
	  Statesptr->routing = OFF;
	  printf("\nWARNING: Routing must be used with inertial or subgrid models. Routing disabled.\n\n");
  }
  if (Statesptr->latlong == ON && Statesptr->SGC == OFF && *verbose==ON) // CCS abort if trying to use latlong without SG:
  {
	  printf("\nWARNING: Latlong must be used with subgrid model. Aborting...\n\n");
	  exit (1);
  }

  return;
}
