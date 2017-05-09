/*
*****************************************************************************
CHECKPOINTING
---------------------

On the command line, -checkpoint turns on checkpointing with the default features
enabled: every 1 hour of computation time, filename based on the output file names.
Or, in the parameter file: "checkpoint" is followed by the frequency in hours of
computation time. When restarting a job, the checkpoint file (.chkpnt) is automatically
found in the output folder, or a different file to load may be specified using
-loadcheck on the command line or "loadcheck" in the parameter file, both followed
by the file name. If it is not found, the job starts from the beginning. A checkpoint
is made at the end of the simulation as well as during it - this makes it possible
to, for example, run the model in steady state for a period, then run multiple
different hydrographs from that point - the new hydrograph should include the period
of steady state in the timings. MDW


Read and write routines for checkpointing of simulations. Read function uses a
previously checkpointed simulation to start the current simulation. Write function
creates checkpoint files every user specified frequency (in hours).

*****************************************************************************
*/

#include "lisflood.h"
#include "VersionHistory.h"

//-----------------------------------------------------------------------------
// CHECKPOINTING: READ
// Added by Matt Wilson, 25 Nov 2004
// MT rewrite open file logic 28 Oct 2007
// Updated by Tim Fewtrell, February 2009
void ReadCheckpoint(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, BoundCs *BCptr,ChannelSegmentType *ChannelSegments, Arrays *Arrptr, int *verbose)
{
  int xchk,ychk;
  FILE *check_fp;
  ChannelSegmentType *csp;  //local pointer to channel segment
  int chseg;	// channel segment loop counter
  int checkv;
  int version,LFv;

  LFv = (int(VersionMajor)*100)+(int(VersionMinor)*10)+int(VersionInc);

  // #### file opening sequence before reading data

  // try and open the default checkpoint file, because if it exists it means a run was started and was halted partway.
  if((check_fp=fopen(Fnameptr->checkpointfilename,"rb"))==NULL) // open file and check it exists
  {
    // message so user knows file does not exist
    if(*verbose==ON) printf("\nUnable to find default checkpoint file: %s\n",Fnameptr->checkpointfilename);

    // check if user wants to load a different starting checkpoint file 
    if(Statesptr->checkfile == ON)
    {
      // try and open the file
      if((check_fp=fopen(Fnameptr->loadCheckpointFilename,"rb"))==NULL) // open file and check it exists
      {
        // message so user knows file does not exist
        if(*verbose==ON) printf("\nUnable to find alternative starting checkpoint file: %s\n",Fnameptr->loadCheckpointFilename);
      }			 
      else
      {
        // message so user knows alternative file opened ok
        if(*verbose==ON) printf("\nAlternative checkpoint file opened: %s",Fnameptr->loadCheckpointFilename);
      }
    }
  }
  else
  {
    if(Statesptr->checkfile == ON) // let user know program is using newer default checkpoint rather then alternative file
    {
      printf("\nAlternative checkpoint file: %s not opened as newer default file exists",Fnameptr->loadCheckpointFilename);
      printf("\n - please delete: %s if you want to start from the Alternative file\n",Fnameptr->checkpointfilename);
    }
    // message so user knows default file opened ok
    if(*verbose==ON) printf("\nDefault checkpoint file opened: %s",Fnameptr->checkpointfilename);
  }


  // #### reading data from a file if opened correctly sequence

  if(check_fp==NULL) // no file opened in initial sequence
  {
    // message so user knows no files opened and program will start from scratch
    if(*verbose==ON) printf("\nNo checkpoint file opened: Starting from scratch\n");
  }
  else		 // go ahead and read data if a file was successfully opened
  {
    // message so user knows file exists and is being read in
    if(*verbose==ON) printf("\nReading checkpoint file\n");
	
	// read checkpointing file version - for future compatibility
    fread(&checkv,sizeof(int),1,check_fp);
	if(feof(check_fp)) {
		if(*verbose==ON) printf("\nUnable to read from checkpoint file (%s) zero size: Starting from scratch\n",Fnameptr->checkpointfilename);
	}
	else if(checkv!=int(CheckVersion)) {
		if(*verbose==ON) printf("\nFile checkpoint version differs from current version: Starting from scratch\n");
	}
	else {
		// file exists and is not of zero length so go ahead and read data
		
		// check LISFLOOD-FP version
		fread(&version,sizeof(int),1,check_fp);
		if(version!=LFv) {
			if(*verbose==ON) printf("\nWARNING: LISFLOOD-FP version differs from current version\n");
		}
		fread(&Solverptr->itrn_time,sizeof(double),1,check_fp);
		fread(&Solverptr->t,sizeof(double),1,check_fp);
		fread(&Solverptr->itCount,sizeof(long),1,check_fp);

		fread(&Parptr->MassTotal,sizeof(double),1,check_fp);
		fread(&Solverptr->Tstep,sizeof(double),1,check_fp);
		fread(&Solverptr->MinTstep,sizeof(double),1,check_fp);

	    fread(&Statesptr->single_op,sizeof(int),1,check_fp);
		fread(Parptr->op_multiswitch,sizeof(int),Parptr->op_multinum,check_fp);
		fread(&Parptr->SaveNo,sizeof(int),1,check_fp);
		fread(&Parptr->SaveTotal,sizeof(double),1,check_fp);

		//Check domain dimensions are as expected to prevent memory overflows/other odd errors
		fread(&xchk,sizeof(int),1,check_fp);
		fread(&ychk,sizeof(int),1,check_fp);
		if(xchk!=Parptr->xsz || ychk!=Parptr->ysz) 
		{
			if(*verbose==ON)
			{
				printf("\nxchk %i and xsz %i\n",xchk,Parptr->xsz);
				printf("ychk %i and ysz %i\n",ychk,Parptr->ysz);
				printf("Domain dimensions do not match those in the checkpoint file!\n");
			}
			exit(0);
		}

		//Malloc dependent variables - see LoadDEM
		fread(Arrptr->H,sizeof(double),Parptr->xsz*Parptr->ysz,check_fp);
		fread(Arrptr->Qx,sizeof(double),(Parptr->xsz+1)*(Parptr->ysz+1),check_fp);
		fread(Arrptr->Qy,sizeof(double),(Parptr->xsz+1)*(Parptr->ysz+1),check_fp);
		fread(Arrptr->maxH,sizeof(double),Parptr->xsz*Parptr->ysz,check_fp);
		fread(Arrptr->maxHtm,sizeof(double),Parptr->xsz*Parptr->ysz,check_fp);
		fread(Arrptr->initHtm,sizeof(double),Parptr->xsz*Parptr->ysz,check_fp);
		if(checkv>1) fread(Arrptr->totalHtm,sizeof(double),Parptr->xsz*Parptr->ysz,check_fp);

		//fread(&ChannelSegments->N_Channel_Segments,sizeof(int),1,check_fp);
		//fread(ChannelSegments,sizeof(ChannelSegmentType),ChannelSegments->N_Channel_Segments,check_fp);

		fread(Arrptr->ChanMask,sizeof(int),Parptr->xsz*Parptr->ysz,check_fp);
		fread(Arrptr->SegMask,sizeof(int),Parptr->xsz*Parptr->ysz,check_fp);

		fread(&BCptr->Qin,sizeof(double),1,check_fp);
		fread(&BCptr->Qout,sizeof(double),1,check_fp);
		fread(&BCptr->QChanOut,sizeof(double),1,check_fp);
		fread(&Solverptr->Hds,sizeof(double),1,check_fp);
		fread(&BCptr->Qpoint,sizeof(double),1,check_fp);

		fread(&Parptr->dx,sizeof(double),1,check_fp);
		fread(&Parptr->dy,sizeof(double),1,check_fp);
		fread(&Parptr->dA,sizeof(double),1,check_fp);

		fread(&Parptr->tlx,sizeof(double),1,check_fp);
		fread(&Parptr->tly,sizeof(double),1,check_fp);
		fread(&Parptr->blx,sizeof(double),1,check_fp);
		fread(&Parptr->bly,sizeof(double),1,check_fp);

		fread(&Solverptr->FArea,sizeof(double),1,check_fp);
		fread(&Solverptr->vol1,sizeof(double),1,check_fp);
		fread(&Solverptr->vol2,sizeof(double),1,check_fp);

		fread(&Solverptr->Qerror,sizeof(double),1,check_fp);

		if(checkv>2) {
			fread(&Parptr->InfilTotalLoss,sizeof(double),1,check_fp);
			fread(&Parptr->EvapTotalLoss,sizeof(double),1,check_fp);
		}
		
		if(Statesptr->ChannelPresent==ON) {
		// main loop for channel segments
			for(chseg=0;chseg<ChannelSegments->N_Channel_Segments;chseg++)
			{
				// set up local pointer to this segment
				csp=ChannelSegments+chseg;
				// read in flows
				fread(csp->ChanQ,sizeof(double),csp->chsz,check_fp);
				// read in A
				fread(csp->A,sizeof(double),csp->chsz,check_fp);
				fread(csp->NewA,sizeof(double),csp->chsz,check_fp);
				// read in junction water depth
				fread(&csp->JunctionH,sizeof(double),1,check_fp);
        
				if(chseg>0) // trib
				{
					// trib outflow
					fread(&ChannelSegments[csp->Next_Segment].Q_Val[csp->Next_Segment_Loc],sizeof(double),1,check_fp);
				}
			}
		}

		if(*verbose==ON) printf("\n - Computation time so far: %.4lf hrs, Sim time: %.3lf of %.3lf secs\n",(Solverptr->itrn_time/3600.0),Solverptr->t,Solverptr->Sim_Time);
		fclose(check_fp);
		}
	}

  return;
}
//-----------------------------------------------------------------------------
// CHECKPOINTING: WRITE
// Added by Matt Wilson, 25 Nov 2004
// Updated by Tim Fewtrell, February 2009
void WriteCheckpoint(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, BoundCs *BCptr,ChannelSegmentType *ChannelSegments, Arrays *Arrptr, int *verbose)
{
  FILE *check_fp;
  ChannelSegmentType *csp;  //local pointer to channel segment
  int chseg;	// channel segment loop counter
  int checkv;
  int version;

  checkv = int(CheckVersion);
  version = (int(VersionMajor)*100)+(int(VersionMinor)*10)+int(VersionInc);

  //File written in binary rather than ASCII to maintain model precision
  if((check_fp=fopen(Fnameptr->checkpointfilename,"wb"))==NULL) {
    if(*verbose==ON) printf("Unable to open checkpoint file: checkpointing off");
    Statesptr->checkpoint=OFF;
    return;
  } else {
    //checkpointing file version - for future compatibility
    fwrite(&checkv,sizeof(int),1,check_fp);
	fwrite(&version,sizeof(int),1,check_fp);

    fwrite(&Solverptr->itrn_time_now,sizeof(double),1,check_fp);
    fwrite(&Solverptr->t,sizeof(double),1,check_fp);
    fwrite(&Solverptr->itCount,sizeof(long),1,check_fp);

    fwrite(&Parptr->MassTotal,sizeof(double),1,check_fp);
    fwrite(&Solverptr->Tstep,sizeof(double),1,check_fp);
    fwrite(&Solverptr->MinTstep,sizeof(double),1,check_fp);

    fwrite(&Statesptr->single_op,sizeof(int),1,check_fp);
    fwrite(Parptr->op_multiswitch,sizeof(int),Parptr->op_multinum,check_fp);
    fwrite(&Parptr->SaveNo,sizeof(int),1,check_fp);
    fwrite(&Parptr->SaveTotal,sizeof(double),1,check_fp); 

    fwrite(&Parptr->xsz,sizeof(int),1,check_fp);
    fwrite(&Parptr->ysz,sizeof(int),1,check_fp);

    fwrite(Arrptr->H,sizeof(double),Parptr->xsz*Parptr->ysz,check_fp);
    fwrite(Arrptr->Qx,sizeof(double),(Parptr->xsz+1)*(Parptr->ysz+1),check_fp);
    fwrite(Arrptr->Qy,sizeof(double),(Parptr->xsz+1)*(Parptr->ysz+1),check_fp);
    fwrite(Arrptr->maxH,sizeof(double),Parptr->xsz*Parptr->ysz,check_fp);
    fwrite(Arrptr->maxHtm,sizeof(double),Parptr->xsz*Parptr->ysz,check_fp);
    fwrite(Arrptr->initHtm,sizeof(double),Parptr->xsz*Parptr->ysz,check_fp);
    if(checkv>1) fwrite(Arrptr->totalHtm,sizeof(double),Parptr->xsz*Parptr->ysz,check_fp);

	//fwrite(&ChannelSegments->N_Channel_Segments,sizeof(int),1,check_fp);
    //fwrite(ChannelSegments,sizeof(ChannelSegmentType),ChannelSegments->N_Channel_Segments,check_fp);

    fwrite(Arrptr->ChanMask,sizeof(int),Parptr->xsz*Parptr->ysz,check_fp);
    fwrite(Arrptr->SegMask,sizeof(int),Parptr->xsz*Parptr->ysz,check_fp);

    fwrite(&BCptr->Qin,sizeof(double),1,check_fp);
    fwrite(&BCptr->Qout,sizeof(double),1,check_fp);
    fwrite(&BCptr->QChanOut,sizeof(double),1,check_fp);
    fwrite(&Solverptr->Hds,sizeof(double),1,check_fp);
    fwrite(&BCptr->Qpoint,sizeof(double),1,check_fp);

    fwrite(&Parptr->dx,sizeof(double),1,check_fp);
    fwrite(&Parptr->dy,sizeof(double),1,check_fp);
    fwrite(&Parptr->dA,sizeof(double),1,check_fp);

    fwrite(&Parptr->tlx,sizeof(double),1,check_fp);
    fwrite(&Parptr->tly,sizeof(double),1,check_fp);
    fwrite(&Parptr->blx,sizeof(double),1,check_fp);
    fwrite(&Parptr->bly,sizeof(double),1,check_fp);

	fwrite(&Solverptr->FArea,sizeof(double),1,check_fp);
	fwrite(&Solverptr->vol1,sizeof(double),1,check_fp);
    fwrite(&Solverptr->vol2,sizeof(double),1,check_fp);
    fwrite(&Solverptr->Qerror,sizeof(double),1,check_fp);

    if(checkv>2) {
      fwrite(&Parptr->InfilTotalLoss,sizeof(double),1,check_fp);
      fwrite(&Parptr->EvapTotalLoss,sizeof(double),1,check_fp);
    }

	if(Statesptr->ChannelPresent==ON) {
	// main loop for channel segments
    for(chseg=0;chseg<ChannelSegments->N_Channel_Segments;chseg++)
    {
      // set up local pointer to this segment
      csp=ChannelSegments+chseg;
      // write out flows
      fwrite(csp->ChanQ,sizeof(double),csp->chsz,check_fp);
	  // write out A
	  fwrite(csp->A,sizeof(double),csp->chsz,check_fp);
	  fwrite(csp->NewA,sizeof(double),csp->chsz,check_fp);
      // read in junction water depth
      fwrite(&csp->JunctionH,sizeof(double),1,check_fp);

      if(chseg>0) // trib
      {
        // trib outflow
        fwrite(&ChannelSegments[csp->Next_Segment].Q_Val[csp->Next_Segment_Loc],sizeof(double),1,check_fp);;
      }
	}
    }
	fclose(check_fp);
	printf("Checkpointed at %.4lf hours computation time\n",(Solverptr->itrn_time_now/3600.0));
  }

  return;
}
