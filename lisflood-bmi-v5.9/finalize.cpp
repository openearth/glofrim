#include "finalize.h"

void final_iterateq()
{
	//output regular files
	fileoutput(Fnameptr,Statesptr,Parptr,Arrptr); 
  
  // Must be last because maxH changed to Maximum elevation !!!
  // Write maximum elevation
  if (Statesptr->SGC==ON) // SGC output
  {
	for(int i=0;i<Parptr->xsz;i++) for(int j=0;j<Parptr->ysz;j++)
	{
	  size_t ptr=i+j*Parptr->xsz;
	  if(Arrptr->maxH[ptr]>1e-3) Arrptr->maxH[ptr]+=Arrptr->SGCz[ptr]; else Arrptr->maxH[ptr]=NULLVAL;
    }
	if (Statesptr->binary_out==ON) write_binrasterfile(Fnameptr->resrootname,-1,".mxeb",Arrptr->maxH,Arrptr->SGCz,0,Statesptr,Parptr);
	else write_ascfile(Fnameptr->resrootname,-1,".mxe",Arrptr->maxH,Arrptr->SGCz,0,Statesptr,Parptr);
  }
  else
  {
	for(int i=0;i<Parptr->xsz;i++) for(int j=0;j<Parptr->ysz;j++)
    {
      size_t ptr=i+j*Parptr->xsz;
      if(Arrptr->maxH[ptr]>1e-3) Arrptr->maxH[ptr]+=Arrptr->DEM[ptr]; else Arrptr->maxH[ptr]=NULLVAL;
    }
	if (Statesptr->binary_out==ON) write_binrasterfile(Fnameptr->resrootname,-1,".mxeb",Arrptr->maxH,Arrptr->DEM,0,Statesptr,Parptr);
	else write_ascfile(Fnameptr->resrootname,-1,".mxe",Arrptr->maxH,Arrptr->DEM,0,Statesptr,Parptr);
  }

  if(*verbose==ON) printf("Finished.\n\n");

}

void final()
{
  time(&Solverptr->time_finish);

  //Final checkpoint
  if(Statesptr->checkpoint==ON) WriteCheckpoint(Fnameptr,Statesptr,Parptr,Solverptr,BCptr,CSTypePtr,Arrptr,verbose);

  // get system time and echo for user
  if(*verbose==ON){
    time_t tf = time(0);
    tm timeF = *localtime(&tf);
    printf("\nFinish Date: %d/%d/%d \n",timeF.tm_mday,timeF.tm_mon + 1,timeF.tm_year + 1900);
    printf("Finish Time: %d:%d:%d \n\n",timeF.tm_hour, timeF.tm_min, timeF.tm_sec);
  }

  //iteration time
  Solverptr->itrn_time = Solverptr->itrn_time + difftime(Solverptr->time_finish, Solverptr->time_start);
  if(*verbose==ON) printf("\n  Total computation time: %.2lf mins\n\n", (Solverptr->itrn_time/60.0));

  if(Statesptr->logfile==ON)
  {
    freopen( "CON", "w", stdout );
    printf("\nLisflood run finished see log file for run details");
  }

  if(Statesptr->save_stages==ON) fclose(Fps.stage_fp);
  sprintf(t1,"%s%s",Fnameptr->resrootname,".stage");
  if(Statesptr->call_gzip==ON) {
    sprintf(tmp_sys_com,"%s%s","gzip -9 -f ",t1);
    system(tmp_sys_com);
  }

  fclose(Fps.mass_fp);
  sprintf(t1,"%s%s",Fnameptr->resrootname,".mass");
  if(Statesptr->call_gzip==ON) {
    sprintf(tmp_sys_com,"%s%s","gzip -9 -f ",t1);
    system(tmp_sys_com);
  }

  delete Arrptr;
  delete Fnameptr;
  delete Statesptr;
  delete Parptr;
  delete Solverptr;
  delete BCptr;
  delete Stageptr;
  delete SGCptr;

  delete verbose;


}
