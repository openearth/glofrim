#include "global.h"
void iterateq_step() {

  int i,j;
  size_t ptr;

  Files *Fptr = &Fps;
  Stage *Locptr = Stageptr;

  ChannelSegmentType *ChannelSegments;
  ChannelSegments = &(ChannelSegmentsVecPtr->front());

// ChannelQ routine is run using Tstep calculated for previous iteration
  // Use default kinematic solver unless diffusive flag set
  if(Statesptr->ChannelPresent==ON)
  {
    if(tstep_counter == -1 || tstep_counter == Solverptr->ts_multiple)
      // do river calc if start of run (-1) or if timestep counter is equal to ts_multiple
    {
      tstep_channel = Solverptr->t - Previous_t; // calculate river timestep (deltaT since last river calc)
      if(Statesptr->diffusive==ON) ChannelQ_Diff(tstep_channel,Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr, RiversIndexVecPtr, RiversIndexPtr);
      else ChannelQ(tstep_channel,Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr,RiversIndexVecPtr, RiversIndexPtr);
      tstep_counter = 0;  // set timestep counter to zero
      Previous_t = Solverptr->t; // record previous timestep
    }
    tstep_counter ++; // increment timestep counter
  }
  // Chose between a sub-grid or conventional floodplain model
  if (Statesptr->SGC==ON)
  {
    // sub grid floodplain models
    SGC_FloodplainQ(Statesptr,Parptr,Solverptr,Arrptr, SGCptr);
    SGC_BCs(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr,SGCptr);
    // Infiltration, evaporation and rainfall routines after time step update (TJF)11
	if(Statesptr->calc_evap==ON) SGC_Evaporation(Parptr,Solverptr,Arrptr,SGCptr);
	if(Statesptr->rainfall==ON && Statesptr->routing==OFF) SGC_Rainfall(Parptr,Solverptr,Arrptr); // CCS rainfall with routing scheme disabled
	if(Statesptr->routing==ON) SGC_Routing(Statesptr,Parptr,Solverptr,Arrptr);	// CCS Routing scheme (controls rainfall if enabled; can be used without rainfall)
	if(Statesptr->hazard==ON) UpdateV(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr);
	SGC_UpdateH(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr,SGCptr);
	BoundaryFlux(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr, ChannelSegmentsVecPtr);
    // NOTES: Update Q's handeled within SGC flux equations (SGC_FloodplainQ etc.) time step calculateion intergrated with UpdateH
  }
  else // conventional floodplain models
  {
    // Time step is reset to initial time step at the start of each FloodplainQ calculation
    if(Solverptr->t>0.0) Solverptr->Tstep=Solverptr->InitTstep;

    FloodplainQ(Statesptr,Parptr,Solverptr,Arrptr,SGCptr);

    BCs(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr);
		if(Statesptr->drychecking==ON) DryCheck(Parptr,Solverptr,Arrptr);

		// Infiltration, evaporation and rainfall routines after time step update (TJF)11
		if(Statesptr->calc_infiltration==ON) FPInfiltration(Parptr,Solverptr,Arrptr);
		if(Statesptr->calc_evap==ON) Evaporation(Parptr,Solverptr,Arrptr);
		if(Statesptr->rainfall==ON && Statesptr->routing==OFF) Rainfall(Parptr,Solverptr,Arrptr); // CCS rainfall with routing scheme disabled
		if(Statesptr->routing==ON) Routing(Statesptr,Parptr,Solverptr,Arrptr);	// CCS Routing scheme (controls rainfall if enabled; can be used without rainfall)

    // If Roe == ON used UpdateHRoe else use the normal UpdateH (JN/IV)
    // if(Statesptr->Roe==ON) UpdateHRoe(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr);
    // else UpdateH(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr);
    // change to CalcFPQxRoe and CalcFPQyRoe... now return Q's for UpdateH (JCN)
    UpdateH(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr);

    if(Statesptr->hazard==ON) UpdateV(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr);

    BoundaryFlux(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr, ChannelSegmentsVecPtr);
    if(Statesptr->acceleration==ON) UpdateQs(Parptr,Arrptr); // Old value of Q updated with current timestep values
    if(Statesptr->Roe==ON) UpdateQsRoe(Parptr,Solverptr,Arrptr); // If Roe==on UpdateQsRoe(JN/IV)
  }

  // Update t with final Tstep calculated
  if(Solverptr->t>0.0) Solverptr->MinTstep=getmin(Solverptr->MinTstep,Solverptr->Tstep);
  Solverptr->t+=Solverptr->Tstep;
  Solverptr->itCount+=1;
  // Update time of initial flood inundation (in hours) and total inundation time (in seconds)
  if((Statesptr->reset_timeinit==ON) & (Solverptr->t>Parptr->reset_timeinit_time))
  { //reset the time of initial inundation if called for in parameter file
    Statesptr->reset_timeinit=OFF;

//#pragma omp parallel for private(j,ptr) // irrelevent benefit JCN
    for(int i=0;i<Parptr->xsz;i++) for(int j=0;j<Parptr->ysz;j++)
                               {
                                 size_t ptr=i+j*Parptr->xsz;
                                 Arrptr->initHtm[ptr]=(NULLVAL);
                               }
    if(*verbose==ON) printf("\n Time of initial inundation reset \n");
  }
  // Update maxH, maxHtm, totalHtm and initHtm at the mass interval if mint_hk is specifed in the .par file OR at every time step if not
  if(Solverptr->t>=Parptr->MassTotal || Statesptr->mint_hk==OFF)
  {
#pragma omp parallel for private (j, ptr)
    for(i=0;i<Parptr->xsz;i++) for(j=0;j<Parptr->ysz;j++)
                               {
                                 ptr=i+j*Parptr->xsz;
                                 if((Arrptr->initHtm[ptr]==(NULLVAL)) && (Arrptr->H[ptr]>0.01)) Arrptr->initHtm[ptr]=Solverptr->t/3600.0;
                                 if(Arrptr->H[ptr]>0.01) Arrptr->totalHtm[ptr]+=Solverptr->Tstep/3600.0;
                                 // Update maximum water depths, and time of maximum (in hours)
                                 if(Arrptr->H[ptr]>Arrptr->maxH[ptr])
                                 {
                                   Arrptr->maxH[ptr]=Arrptr->H[ptr];
                                   Arrptr->maxHtm[ptr]=Solverptr->t/3600.0;
                                 }
                               }
  }

  // Calculate mass balance error
  if(Solverptr->t>=Parptr->MassTotal)
  {
    Solverptr->vol2=DomainVol(Statesptr,Parptr,ChannelSegments,Arrptr,ChannelSegmentsVecPtr); // CCS

    // calc losses for this mass interval
    double loss = (Parptr->InfilTotalLoss - Parptr->InfilLoss) + (Parptr->EvapTotalLoss - Parptr->EvapLoss) - (Parptr->RainTotalLoss - Parptr->RainLoss);

    //Solverptr->Qerror=BCptr->Qin-BCptr->Qout-(Solverptr->vol2+loss-Solverptr->vol1)/Parptr->MassInt;
    // New version using VolInMT and VolOutMT
    // volume error
    Solverptr->Verror=BCptr->VolInMT-BCptr->VolOutMT-(Solverptr->vol2+loss-Solverptr->vol1);
    // Q error
    Solverptr->Qerror=Solverptr->Verror/Parptr->MassInt;
    // reset to 0.0
    BCptr->VolInMT = 0.0;
    BCptr->VolOutMT = 0.0;

    // record cumulative loss for next time.
    Parptr->InfilLoss = Parptr->InfilTotalLoss;
    Parptr->EvapLoss  = Parptr->EvapTotalLoss;
    Parptr->RainLoss  = Parptr->RainTotalLoss;

    // Calculate flood area
    double FloodArea=0.0;
	double dA = Parptr->dA;
#pragma omp parallel for private(j,ptr) reduction ( + : FloodArea)
	for(i=0;i<Parptr->xsz;i++) for(j=0;j<Parptr->ysz;j++)
      {
        ptr=i+j*Parptr->xsz;
		if(Statesptr->latlong == ON) dA = Arrptr->dA[ptr]; // if latlong is on change dA to local cell area
        if(Statesptr->porosity==ON)
        {
          if(Arrptr->H[ptr]>0.01) FloodArea+=dA*Arrptr->paerial[ptr]; // If porosity used, scale flooded area by porosity (TJF)
        }
        else if (Statesptr->SGC == ON)
        {
		  if(Arrptr->H[ptr]-Arrptr->SGCbfH[ptr] > Solverptr->DepthThresh) FloodArea+=dA; // If sub-grid used remove channel depth
		}
		else
		{
		  if(Arrptr->H[ptr] > 0.01) FloodArea+=dA; // standard ara calculation
		}
	  }
    Solverptr->FArea=FloodArea;

    fprintf(Fptr->mass_fp,"%-12.3f %-10.4f %-10.4f %-10li %12.4e %12.4e  %-11.3f %-10.3f %-11.3f %12.4e %12.4e %12.4e\n",Solverptr->t,Solverptr->Tstep,Solverptr->MinTstep,Solverptr->itCount,Solverptr->FArea,Solverptr->vol2,BCptr->Qin,Solverptr->Hds,BCptr->Qout,Solverptr->Qerror,Solverptr->Verror,Parptr->RainTotalLoss-(Parptr->InfilTotalLoss+Parptr->EvapTotalLoss));
    fflush(Fptr->mass_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.

    Solverptr->vol1=Solverptr->vol2;
    Parptr->MassTotal+=Parptr->MassInt;

    //stage output
    if(Statesptr->save_stages==ON)
    {
      fprintf(Fptr->stage_fp,"%12.3f",Solverptr->t);
      for(i=0;i<Locptr->Nstages;i++)
      {
        if(Locptr->stage_check[i]==1) fprintf(Fptr->stage_fp,"%10.4f",Arrptr->H[Locptr->stage_grid_x[i]+Locptr->stage_grid_y[i]*Parptr->xsz]);
        else fprintf(Fptr->stage_fp,"-\t");
      }
      fprintf(Fptr->stage_fp,"\n");
      fflush(Fptr->stage_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
      // added to export scalar velocity
      if (Statesptr->voutput==ON)
      {
        fprintf(Fptr->vel_fp,"%12.3f",Solverptr->t);
        for(i=0;i<Locptr->Nstages;i++)
        {
          if(Locptr->stage_check[i]==1) fprintf(Fptr->vel_fp,"%10.4f", sqrt(pow(getmax(  fabs(Arrptr->Vx[Locptr->vel_grid_xy[i]]),fabs(Arrptr->Vx[Locptr->vel_grid_xy[i]+1]) ),2) +  pow(getmax(fabs(Arrptr->Vy[Locptr->vel_grid_xy[i]]),fabs(Arrptr->Vy[Locptr->vel_grid_xy[i]+Parptr->xsz+1]) ),2))   );
          else fprintf(Fptr->vel_fp,"-\t");
        }
        fprintf(Fptr->vel_fp,"\n");
        fflush(Fptr->vel_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
      }
    }
	//virtual gauge output
      if(Statesptr->gsection==ON)
      {
        fprintf(Fptr->gau_fp,"%12.2f",Solverptr->t); // print tiem to file
        for(i=0;i<Locptr->Ngauges;i++) // loop through each virtual gauge
        {
			// call discharge calculation function
			double discharge = CalcVirtualGauge(i, Parptr, Arrptr, Locptr);
			fprintf(Fptr->gau_fp," %10.3f",discharge); // Print discharge to file
        }
        fprintf(Fptr->gau_fp,"\n"); // print end of line
		fflush(Fptr->gau_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
	  }

    // Checkpointing
    if(Statesptr->checkpoint==ON)
    {
      //iteration time
      time(&Solverptr->time_check);
      Solverptr->itrn_time_now = Solverptr->itrn_time + difftime(Solverptr->time_check, Solverptr->time_start);
      if(Solverptr->itrn_time_now>=Parptr->nextcheck)
      {
        WriteCheckpoint(Fnameptr, Statesptr, Parptr, Solverptr, BCptr,ChannelSegments, Arrptr, verbose);
        Parptr->nextcheck=Solverptr->itrn_time_now+(Parptr->checkfreq*3600);
      }
    }
  }

  // Regular output

  if(Solverptr->t>=Parptr->SaveTotal)
  {

    double Comp_time, Model_Comp_Ratio, Model_time_left, Est_Time_Tot, Est_Time_Fin;

    time(&Solverptr->time_check);
    Comp_time = difftime(Solverptr->time_check, Solverptr->time_start)/60;
    if (Comp_time !=0 && Statesptr->comp_out==ON) // only of t is not zero (can't divide by zero)
    {
      Model_Comp_Ratio = ((Solverptr->t/60)/Comp_time);
      Model_time_left = (Solverptr->Sim_Time - Solverptr->t)/60;
      Est_Time_Fin = (Model_time_left/Model_Comp_Ratio);
      Est_Time_Tot = Comp_time + Est_Time_Fin;
      printf("T(mins): M: %.1lf, C: %.1lf, M/C: %.2lf, ETot: %.1lf, EFin: %.1lf\n", (Solverptr->t/60.0), Comp_time, Model_Comp_Ratio, Est_Time_Tot, Est_Time_Fin);
    }

      write_regular_output(Fnameptr, Solverptr, Statesptr, Parptr, Arrptr, SGCptr);

	  //regular profile output, if requested in param file
      if(Statesptr->profileoutput==ON)
      {
        write_profile(Fnameptr->resrootname,Parptr->SaveNo,".profile",Statesptr,ChannelSegments,Arrptr,Parptr, RiversIndexVecPtr, RiversIndexPtr); // CCS
      }

      // update interval counter
      Parptr->SaveTotal+=Parptr->SaveInt;
      Parptr->SaveNo+=1;

    }

    // Single overpass
    if(Statesptr->single_op==ON && Solverptr->t>=Parptr->op)
    {
      // write rasters
	  if (Statesptr->binary_out==ON) write_binrasterfile(Fnameptr->resrootname,-1,".opb",Arrptr->H,Arrptr->DEM,0,Statesptr,Parptr);
      else write_ascfile(Fnameptr->resrootname,-1,".op",Arrptr->H,Arrptr->DEM,0,Statesptr,Parptr);
      // write profiles
	  if (Statesptr->profileoutput==ON) write_profile(Fnameptr->resrootname,-1,".profile",Statesptr,ChannelSegments,Arrptr,Parptr, RiversIndexVecPtr, RiversIndexPtr);
      if(*verbose==ON) printf("Writing overpass at %lf seconds\n",Solverptr->t);
      Statesptr->single_op=OFF;

	  // raster elevation output
    if(Statesptr->save_elev==ON)
		{
		  if (Statesptr->SGC==ON) // SGC output
		  {
			if (Statesptr->binary_out==ON) write_binrasterfile(Fnameptr->resrootname,-1,".opelevb",Arrptr->H,Arrptr->SGCz,3,Statesptr,Parptr);
			else write_ascfile(Fnameptr->resrootname,-1,".opelev",Arrptr->H,Arrptr->SGCz,3,Statesptr,Parptr);
		  }
		  else // standard model output
		  {
			if (Statesptr->binary_out==ON) write_binrasterfile(Fnameptr->resrootname,-1,".opelevb",Arrptr->H,Arrptr->DEM,3,Statesptr,Parptr);
			else write_ascfile(Fnameptr->resrootname,-1,".opelev",Arrptr->H,Arrptr->DEM,3,Statesptr,Parptr);
		  }
		}
	}

  // Multiple overpasses
  if(Statesptr->multi_op==ON)
  {
    for(i=0;i<Parptr->op_multinum;i++)
    {
      if(Solverptr->t>=Parptr->op_multisteps[i] && Parptr->op_multiswitch[i]==0)
      {
		  // raster depth output
          if (Statesptr->binary_out==ON) write_binrasterfile(Fnameptr->resrootname,i,"-T.opb",Arrptr->H,Arrptr->DEM,0,Statesptr,Parptr);
		  else write_ascfile(Fnameptr->resrootname,i,"-T.op",Arrptr->H,Arrptr->DEM,0,Statesptr,Parptr);
		  // profiles output
		  if (Statesptr->profileoutput==ON) write_profile(Fnameptr->resrootname,i,"-T.profile",Statesptr,ChannelSegments,Arrptr,Parptr, RiversIndexVecPtr, RiversIndexPtr);

        Parptr->op_multiswitch[i]=1;
        if(*verbose==ON) printf("Writing overpass %d at %lf seconds\n",i,Parptr->op_multisteps[i]);
				  // raster elevation output
		  if(Statesptr->save_elev==ON)
	      {
			  if (Statesptr->SGC==ON) // SGC output
			  {
				if (Statesptr->binary_out==ON) write_binrasterfile(Fnameptr->resrootname,i,"-T.opelevb",Arrptr->H,Arrptr->SGCz,3,Statesptr,Parptr);
				else write_ascfile(Fnameptr->resrootname,i,"-T.opelev",Arrptr->H,Arrptr->SGCz,3,Statesptr,Parptr);
			  }
			  else // standard model output
			  {
				if (Statesptr->binary_out==ON) write_binrasterfile(Fnameptr->resrootname,i,"-T.opelevb",Arrptr->H,Arrptr->DEM,3,Statesptr,Parptr);
				else write_ascfile(Fnameptr->resrootname,i,"-T.opelev",Arrptr->H,Arrptr->DEM,3,Statesptr,Parptr);
			  }
		  }
      }
    }
  }

  // If requested on command line, check whether we should kill this simulation...
  if(Statesptr->killsim==ON)
  {
    //iteration time
    time(&Solverptr->time_check);
    Solverptr->itrn_time_now = Solverptr->itrn_time + difftime(Solverptr->time_check, Solverptr->time_start);
    // check if we have reached the kill time
    if(Solverptr->itrn_time_now>=Parptr->killsim_time) {
      if(*verbose==ON) printf("Simulation kill time reached... ");
      return;
    }
  }

  // If an auto steady-state has been requested:
  if(Statesptr->steadycheck==ON && Solverptr->t>=Parptr->steadyTotal)
  {
    Parptr->steadyQdiff=BCptr->Qin-BCptr->Qout;
    if(fabs(Parptr->steadyQdiff)<Parptr->steadyQtol) {
      steadyCount+=1; // keep going a bit to make sure...
      if(steadyCount==10) return;
    } else {
      steadyCount=0;
    }
    Parptr->steadyTotal+=Parptr->steadyInt;
  }
}
