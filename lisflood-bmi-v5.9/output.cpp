/*
*****************************************************************************
FILE OUTPUT
---------------------

A number of functions used to control the file output from LISFLOOD-FP. A short
definition of each is outlined below.

void write_ascii():	general purpose write routine for ascii raster files used generically
to pass arrays and necessary parameters.
void fileoutput():	Calls the write_ascii routine to output various model statistic
in raster format.

*****************************************************************************
*/

#include "lisflood.h"


//-----------------------------------------------------------------------------
// GENERAL ASCII RASTER FILE WRITE ROUTINE
// MT new multi purpose version - eliminates a lot of repetition
// and implements new filenaming
void write_ascfile(char *root, int SaveNumber, char *extension, double *data, double *dem, int outflag, States *Statesptr,Pars *Parptr)
/*
Purpose: writes data to a standard ascii file format
Parameters:
char *root		-	root of filename
int SaveNumber	-	save number to add to filename (if <0, will be ignored)
char *extension	-	filename extension text
double *data		-	pointer to data
double *dem		-	pointer to dem
int outflag		-	flag (if 0 = normal, if 1 or 2 indicated fluxes, if 3 indicates special option Water elev ouput DEM+H)
States *Statesptr - pointer to States structure
Pars *Parptr - pointer to Parameters structure
*/
{
  int i,j;
  double temp;
  FILE *fp;
  char fnam[800];
  char tmp_sys_com[255];

  // check if there is a savenumber to add and create filename
  if(SaveNumber>=0 && SaveNumber<= 9999)	sprintf(fnam,"%s-%.4d%s",root,SaveNumber,extension);
  else if  (SaveNumber>9999)				sprintf(fnam,"%s-%d%s",root,SaveNumber,extension);
  else 										sprintf(fnam,"%s%s",root,extension);

  // open file
  fp=fopen(fnam,"wb");

  // check file opened ok, if not give warning and exit function
  if(fp==NULL)
  {
    printf("Problems writing to file %s\n",fnam);
    return;
  }

  //Output alternative header or default header
  if(Statesptr->alt_ascheader==ON)
  {
    for(i=0;i<6;i++) 
    {
      fprintf(fp,"%s",Parptr->ascheader[i]);
    }
  } 
  else if(outflag==0 || outflag==3) 
  {
    fprintf(fp,"ncols         %i\n",Parptr->xsz);
    fprintf(fp,"nrows         %i\n",Parptr->ysz);
    fprintf(fp,"xllcorner     %lf\n",Parptr->blx);
    fprintf(fp,"yllcorner     %lf\n",Parptr->bly);
    fprintf(fp,"cellsize      %lf\n",Parptr->dx);
    fprintf(fp,"NODATA_value  %lf\n",NULLVAL);
  }

  // output data switched by type
  switch (outflag)
  {
  case 0 :	// normal cell output
    for(j=0;j<Parptr->ysz;j++)
    {
      for(i=0;i<Parptr->xsz;i++)
      {
        fprintf(fp,"%.3f\t",data[i+j*Parptr->xsz]);
      }
      // output end of line
      fprintf(fp,"\n");
    }
    break;
  case 1 :	// flux ouput - ie cell boundaries not cells
    // Edited to output the correct number of rows as Qx includes fluxes across boundaries (TJF)
    // Origin offset by dx*0.5 to output fluxes at boundaries when read into GIS (TJF)
    fprintf(fp,"ncols         %i\n",Parptr->xsz+1);
    fprintf(fp,"nrows         %i\n",Parptr->ysz);
    fprintf(fp,"xllcorner     %lf\n",Parptr->blx-(Parptr->dx/2.0));
    fprintf(fp,"yllcorner     %lf\n",Parptr->bly);
    fprintf(fp,"cellsize      %lf\n",Parptr->dx);
    fprintf(fp,"NODATA_value  %lf\n",NULLVAL);

    for(j=0;j<Parptr->ysz;j++)
    {
      for(i=0;i<Parptr->xsz+1;i++)
      {
        fprintf(fp,"%.3f\t",data[i+j*(Parptr->xsz+1)]);
      }
      // output end of line
      fprintf(fp,"\n");
    }
    break;
  case 2 :	// flux ouput - ie cell boundaries not cells
    // Edited to output the correct number of rows as Qy includes fluxes across boundaries (TJF)
    // Origin offset by dy*0.5 to output fluxes at boundaries when read into GIS (TJF)
    fprintf(fp,"ncols         %i\n",Parptr->xsz);
    fprintf(fp,"nrows         %i\n",Parptr->ysz+1);
    fprintf(fp,"xllcorner     %lf\n",Parptr->blx);
    fprintf(fp,"yllcorner     %lf\n",Parptr->bly-(Parptr->dx/2.0));
    fprintf(fp,"cellsize      %lf\n",Parptr->dx);
    fprintf(fp,"NODATA_value  %lf\n",NULLVAL);

    for(j=0;j<Parptr->ysz+1;j++)
    {
      for(i=0;i<Parptr->xsz;i++)
      {
        fprintf(fp,"%.3f\t",data[i+j*(Parptr->xsz+1)]);
      }
      // output end of line
      fprintf(fp,"\n");
    }
    break;
  case 3 :	// special output for water elevation (add DEM to data)
    for(j=0;j<Parptr->ysz;j++)
    {
      for(i=0;i<Parptr->xsz;i++)
      {

        if(data[i+j*Parptr->xsz]<=0.01) // set to null for very shallow depth
        {
          // PB: Fixed bug in generation of elev files so these are only created when H<=0.01 (not H<=0).  
          // Without this you get significant lateral curvature caused by thin films of water ponding in drying elements.
          temp=NULLVAL;
        }
        else
        {
          temp=dem[i+j*Parptr->xsz]+data[i+j*Parptr->xsz];
        }
        fprintf(fp,"%11.3f",temp);
      }
      // output end of line
      fprintf(fp,"\n");
    }
    break;
  }

  // close file
  fclose(fp);

  // check if we need to zip the file up
  if(Statesptr->call_gzip==ON) 
  {
    sprintf(tmp_sys_com,"%s%s","gzip -9 -f ",fnam);
    system(tmp_sys_com);
  }

  return;
}


//-----------------------------------------------------------------------------
// PROFILE FILE WRITE ROUTINE
void write_profile(char *root, int SaveNumber, char *extension, States *Statesptr,ChannelSegmentType *ChannelSegments, Arrays *Arrptr,Pars *Parptr, vector<int> *RiversIndexVecPtr, int *RiversIndexPtr) // CCS
/*
Purpose: writes profile data to a file
Parameters:
char *root		    -	root of filename
int SaveNumber	  -	save number to add to filename (if <0, will be ignored)
char *extension	  -	filename extension text
States *Statesptr - pointer to States structure
ChannelSegmentType *ChannelSegments - pointer to Channel structure
Arrays *Arrptr    - pointer to Array structure
Pars *Parptr      - pointer to Parameters structure
vector<int> *RiversIndexVecPtr, int *RiversIndexPtr	- CCS River index vectors for multiple river loop 
*/
{
  int chseg,j,pi,pj;
  FILE *fp;
  char fnam[800];
  char tmp_sys_com[255];
  ChannelSegmentType *csptr; // temporary pointer to one segment of channel
  int nriv, low, high; // CCS for multiple river loop.
  

   for(nriv=0; nriv<(int)RiversIndexVecPtr->size(); nriv++) // CCS 
	  {
		  high = RiversIndexPtr[nriv]-1;
		  if(nriv==0)
		  {
			  low = 0;
		  }
		  else
		  {
			  low = RiversIndexPtr[nriv-1];
		  }

		  // check if there is a savenumber to add and create filename
		  if(SaveNumber>=0 && SaveNumber<=9999)	sprintf(fnam,"%s-river%i-%.4d%s",root,nriv,SaveNumber,extension);
		  else if  (SaveNumber>9999)			sprintf(fnam,"%s-river%i-%d%s",root,nriv,SaveNumber,extension);
		  else 									sprintf(fnam,"%s-river%i-%s",root,nriv,extension);

		  // open file
		  fp=fopen(fnam,"wb");

		  // check file opened ok, if not give warning and exit function
		  if(fp==NULL)
		  {
		  printf("WARNING: Problems writing to file %s\n",fnam);
		  return;
		  }

		  // write profiles - looping through segments
		  for(chseg=low; chseg<=high; chseg++)
		  {
			  csptr=ChannelSegments+chseg;
			  fprintf(fp,"Channel_segment: %d\n",chseg);
			  fprintf(fp,"ChanX ChanY Chainage Width Mannings Slope ");
			  fprintf(fp,"BankZ BedElev WaterElev WaterDepth Flow");
			  if(Statesptr->debugmode==ON) fprintf(fp," Q_Ident Q_Val");
			  fprintf(fp,"\n");
			  for(j=0;j<csptr->chsz;j++)
				{
				  pi=csptr->ChanX[j];
				  pj=csptr->ChanY[j];
				  if(chseg!=0 && j==csptr->chsz-1)
				  {
					// last node of trib is junction with main - show dummy calc node
					fprintf(fp,"%-10.3f %10.3f", Parptr->blx + (csptr->ChanX[j] * Parptr->dx) + (Parptr->dx/2), Parptr->bly + (Parptr->ysz*Parptr->dx) - ((csptr->ChanY[j] * Parptr->dx) + (Parptr->dx/2))); 
					fprintf(fp," %10.3f %10.3f", csptr->Chainage[j], csptr->ChanWidth[j]); 
					fprintf(fp," %10.5f %10.8f", csptr->ChanN[j], csptr->Shalf[j]*fabs(csptr->Shalf[j])); 
					fprintf(fp," %10.3f %10.3f", csptr->BankZ[j], csptr->JunctionDEM);
					fprintf(fp," %10.3f %10.3f", csptr->JunctionDEM+csptr->JunctionH, csptr->JunctionH);
					fprintf(fp," %10.3f", csptr->ChanQ[j]); 
					if(Statesptr->debugmode==ON) fprintf(fp," %10i %10.3f",csptr->Q_Ident[j],csptr->Q_Val[j]);
					fprintf(fp,"\n"); 
				  }
				  else                       
				  {
					// all other points
					// note x,y position is centre of cell.
					fprintf(fp,"%-10.3f %10.3f", Parptr->blx + (csptr->ChanX[j] * Parptr->dx) + (Parptr->dx/2), Parptr->bly + (Parptr->ysz*Parptr->dx) - ((csptr->ChanY[j] * Parptr->dx) + (Parptr->dx/2))); 
					fprintf(fp," %10.3f %10.3f", csptr->Chainage[j], csptr->ChanWidth[j]); 
					fprintf(fp," %10.5f %10.8f", csptr->ChanN[j], csptr->Shalf[j]*fabs(csptr->Shalf[j])); 
					fprintf(fp," %10.3f %10.3f", csptr->BankZ[j], Arrptr->DEM[pi+pj*Parptr->xsz]);
					fprintf(fp," %10.3f %10.3f", Arrptr->DEM[pi+pj*Parptr->xsz]+Arrptr->H[pi+pj*Parptr->xsz], Arrptr->H[pi+pj*Parptr->xsz]);
					fprintf(fp," %10.3f", csptr->ChanQ[j]); 
					if(Statesptr->debugmode==ON) fprintf(fp," %10i %10.3f",csptr->Q_Ident[j],csptr->Q_Val[j]);
					fprintf(fp,"\n");
				  }
				}
		  }
   }
				  

  // close file
  fclose(fp);

  // check if we need to zip the file up
  if(Statesptr->call_gzip==ON) 
  {
    sprintf(tmp_sys_com,"%s%s","gzip -9 -f ",fnam);
    system(tmp_sys_com);
  }


  return;
}


//-----------------------------------------------------------------------------
// FILE OUTPUT
void fileoutput(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,Arrays *Arrptr)
{
  // output binary or ascii rasters
  if (Statesptr->binary_out==ON) // output binary of ascii rasters
  {
	// Write time of initial flood inundation
	write_binrasterfile(Fnameptr->resrootname,-1,".inittmb",Arrptr->initHtm,Arrptr->DEM,0,Statesptr,Parptr);
	// Write total inundation time
	write_binrasterfile(Fnameptr->resrootname,-1,".totaltmb",Arrptr->totalHtm,Arrptr->DEM,0,Statesptr,Parptr);
	// Write maximum depth
	write_binrasterfile(Fnameptr->resrootname,-1,".maxb",Arrptr->maxH,Arrptr->DEM,0,Statesptr,Parptr);
	// Write time of maximum depth
	write_binrasterfile(Fnameptr->resrootname,-1,".maxtmb",Arrptr->maxHtm,Arrptr->DEM,0,Statesptr,Parptr);
	if(Statesptr->voutput==ON)
	{
		write_binrasterfile(Fnameptr->resrootname,-1,".maxVxb",Arrptr->maxVx,Arrptr->DEM,1,Statesptr,Parptr);
		write_binrasterfile(Fnameptr->resrootname,-1,".maxVyb",Arrptr->maxVy,Arrptr->DEM,2,Statesptr,Parptr);
	}
	if(Statesptr->hazard==ON)
	{
		// Write maximum V Vd and Hazard
		write_binrasterfile(Fnameptr->resrootname,-1,".maxVcb",Arrptr->maxVc,Arrptr->DEM,0,Statesptr,Parptr);
		write_binrasterfile(Fnameptr->resrootname,-1,".maxVcdb",Arrptr->maxVcH,Arrptr->DEM,0,Statesptr,Parptr);
		write_binrasterfile(Fnameptr->resrootname,-1,".maxHazb",Arrptr->maxHaz,Arrptr->DEM,0,Statesptr,Parptr);
	}  
  }
  else
  {	  
	// Write time of initial flood inundation
	write_ascfile(Fnameptr->resrootname,-1,".inittm",Arrptr->initHtm,Arrptr->DEM,0,Statesptr,Parptr);
	// Write total inundation time
	write_ascfile(Fnameptr->resrootname,-1,".totaltm",Arrptr->totalHtm,Arrptr->DEM,0,Statesptr,Parptr);
	// Write maximum depth
	write_ascfile(Fnameptr->resrootname,-1,".max",Arrptr->maxH,Arrptr->DEM,0,Statesptr,Parptr);
	// Write time of maximum depth
	write_ascfile(Fnameptr->resrootname,-1,".maxtm",Arrptr->maxHtm,Arrptr->DEM,0,Statesptr,Parptr);
	if(Statesptr->voutput==ON)
	{
		write_ascfile(Fnameptr->resrootname,-1,".maxVx",Arrptr->maxVx,Arrptr->DEM,1,Statesptr,Parptr);
		write_ascfile(Fnameptr->resrootname,-1,".maxVy",Arrptr->maxVy,Arrptr->DEM,2,Statesptr,Parptr);
	}
	if(Statesptr->hazard==ON)
	{
		// Write maximum V Vd and Hazard
		write_ascfile(Fnameptr->resrootname,-1,".maxVc",Arrptr->maxVc,Arrptr->DEM,0,Statesptr,Parptr);
		write_ascfile(Fnameptr->resrootname,-1,".maxVcd",Arrptr->maxVcH,Arrptr->DEM,0,Statesptr,Parptr);
		write_ascfile(Fnameptr->resrootname,-1,".maxHaz",Arrptr->maxHaz,Arrptr->DEM,0,Statesptr,Parptr);
	}
  }
}

//-----------------------------------------------------------------------------
// Debug option file output (currently the modified DEM and the channel and trib masks and subgrid channels
void debugfileoutput(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,Arrays *Arrptr)
{
  double *TempChanMask;
  int i,j;

  // Write out final DEM if debug requested
  write_ascfile(Fnameptr->resrootname,-1,".dem",Arrptr->DEM,Arrptr->DEM,0,Statesptr,Parptr);

  // Write out channel mask if channel exists and debug requested
  if(Statesptr->ChannelPresent==ON) 
  {
    TempChanMask=new double[Parptr->xsz*Parptr->ysz];

    for(j=0;j<Parptr->ysz;j++)
    {
      for(i=0;i<Parptr->xsz;i++)
      {
        TempChanMask[i+j*Parptr->xsz]= Arrptr->ChanMask[i+j*Parptr->xsz];
      }
    }

    write_ascfile(Fnameptr->resrootname,-1,".chmask",TempChanMask,Arrptr->DEM,0,Statesptr,Parptr);

    for(j=0;j<Parptr->ysz;j++)
    {
      for(i=0;i<Parptr->xsz;i++)
      {
        TempChanMask[i+j*Parptr->xsz]= Arrptr->SegMask[i+j*Parptr->xsz];
      }
    }

    write_ascfile(Fnameptr->resrootname,-1,".segmask",TempChanMask,Arrptr->DEM,0,Statesptr,Parptr);

  }

  // (MT) Write out SGC arrays if used and debug requested
  if(Statesptr->SGC==ON) 
  {
	// Write out bed elevation
	write_ascfile(Fnameptr->resrootname,-1,"_SGC_bedZ.asc",Arrptr->SGCz,Arrptr->DEM,0,Statesptr,Parptr);
	// Write out width
	write_ascfile(Fnameptr->resrootname,-1,"_SGC_width.asc",Arrptr->SGCwidth,Arrptr->DEM,0,Statesptr,Parptr);
	// Write out bankfull depth
	write_ascfile(Fnameptr->resrootname,-1,"_SGC_bfdepth.asc",Arrptr->SGCbfH,Arrptr->DEM,0,Statesptr,Parptr);
  }


}

//-----------------------------------------------------------------------------
// GENERAL BINARY RASTER FILE WRITE ROUTINE
// MT new multi purpose version - eliminates a lot of repetition
// and implements new filenaming
void write_binrasterfile(char *root, int SaveNumber, char *extension, double *data, double *dem, int outflag, States *Statesptr,Pars *Parptr)
/*
Purpose: writes data that would be in an ascii raster to a binary file format
Parameters:
char *root		-	root of filename
int SaveNumber	-	save number to add to filename (if <0, will be ignored)
char *extension	-	filename extension text
double *data		-	pointer to data
double *dem		-	pointer to dem
int outflag		-	flag (if 0 = normal, if 1 or 2 indicated fluxes, if 3 indicates special option Water elev ouput DEM+H)
States *Statesptr - pointer to States structure
Pars *Parptr - pointer to Parameters structure
*/
{
  int i,j,tempi;
  double temp, no_data=-9999;
  FILE *fp;
  char fnam[800];
  //char tmp_sys_com[255];

  // check if there is a savenumber to add and create filename
  if(SaveNumber>=0 && SaveNumber<=9999)	sprintf(fnam,"%s-%.4d%s",root,SaveNumber,extension);
  else if  (SaveNumber>9999)			sprintf(fnam,"%s-%d%s",root,SaveNumber,extension);
  else 									sprintf(fnam,"%s%s",root,extension);

  // open file
  fp=fopen(fnam,"wb");

  // check file opened ok, if not give warning and exit function
  if(fp==NULL)
  {
    printf("Problems writing to file %s\n",fnam);
    return;
  }

  //Output alternative header or default header
  //if(Statesptr->alt_ascheader==ON)
  //{
    //for(i=0;i<6;i++) 
    //{
     // fprintf(fp,"%s",Parptr->ascheader[i]);
	 //   size_t fwrite ( const void * ptr, size_t size, size_t count, FILE * stream );
	//    fwrite (ptr, double , sizeof(buffer) , fp );
    //}
  //} 
  if(outflag==0 || outflag==3) 
  {
	fwrite (&Parptr->xsz, sizeof(int)   , 1, fp );
	fwrite (&Parptr->ysz, sizeof(int)   , 1, fp );
	fwrite (&Parptr->blx, sizeof(double), 1, fp );
	fwrite (&Parptr->bly, sizeof(double), 1, fp );
	fwrite (&Parptr->dx , sizeof(double), 1, fp );
	fwrite (&no_data    , sizeof(double), 1, fp );
  }

  // output data switched by type
  switch (outflag)
  {
  case 0 :	// normal cell output
	fwrite(data,sizeof(double),Parptr->xsz*Parptr->ysz,fp);
  
    break;
  case 1 :	// flux ouput - ie cell boundaries not cells
    // Edited to output the correct number of rows as Qx includes fluxes across boundaries (TJF)
    // Origin offset by dx*0.5 to output fluxes at boundaries when read into GIS (TJF)
    tempi = Parptr->xsz+1;
    fwrite (&tempi, sizeof(int) , 1 , fp );
    
    fwrite (&Parptr->ysz, sizeof(int) , 1, fp );
    temp = Parptr->blx-(Parptr->dx/2.0);
    fwrite (&temp, sizeof(double) , 1, fp );
    fwrite (&Parptr->bly, sizeof(double) , 1, fp );
    fwrite (&Parptr->dx, sizeof(double) , 1, fp );
    fwrite (&no_data, sizeof(double) , 1, fp );

    fwrite (data,sizeof(double),(Parptr->xsz+1)*(Parptr->ysz+1),fp);

    break;
  case 2 :	// flux ouput - ie cell boundaries not cells
    // Edited to output the correct number of rows as Qy includes fluxes across boundaries (TJF)
    // Origin offset by dy*0.5 to output fluxes at boundaries when read into GIS (TJF)
    fwrite (&Parptr->xsz, sizeof(int) , 1 , fp );
    tempi = Parptr->ysz+1;
    fwrite (&tempi, sizeof(int) , 1, fp );
    fwrite (&Parptr->blx, sizeof(double) , 1, fp );
    temp = Parptr->bly-(Parptr->dx/2.0);
    fwrite (&temp, sizeof(double) , 1, fp );
    fwrite (&Parptr->dx, sizeof(double) , 1, fp );
    fwrite (&no_data, sizeof(double) , 1, fp );  

    fwrite (data,sizeof(double),(Parptr->xsz+1)*(Parptr->ysz+1),fp);

    break;
  case 3 :	// special output for water elevation (add DEM to data)
    for(j=0;j<Parptr->ysz;j++)
    {
      for(i=0;i<Parptr->xsz;i++)
      {

        if(data[i+j*Parptr->xsz]<=0.01) // set to null for very shallow depth
        {
          // PB: Fixed bug in generation of elev files so these are only created when H<=0.01 (not H<=0).  
          // Without this you get significant lateral curvature caused by thin films of water ponding in drying elements.
          temp=NULLVAL;
        }
        else
        {
          temp= dem[i+j*Parptr->xsz]+ data[i+j*Parptr->xsz];
        }
        fwrite(&temp,sizeof(double),1,fp);
      }
    }
    break;
  }

  // close file
  fclose(fp);

  return;
}
// WRITE SGC FLOODPLAIN DEPTHS e.g it takes out the channel depths
void write_ascfile_SGCf(char *root, int SaveNumber, char *extension, double *data, double *SGCbfH, States *Statesptr,Pars *Parptr)
/*
Purpose: writes data to a standard ascii file format
Parameters:
char *root		-	root of filename
int SaveNumber	-	save number to add to filename (if <0, will be ignored)
char *extension	-	filename extension text
double *data		-	pointer to data
double *SGCbfH		-	pointer to SGC bankful depth
int outflag		-	flag (if 0 = normal, if 1 or 2 indicated fluxes, if 3 indicates special option Water elev ouput DEM+H)
States *Statesptr - pointer to States structure
Pars *Parptr - pointer to Parameters structure
*/
{
  int i,j;
  double temp;
  FILE *fp;
  char fnam[800];
  char tmp_sys_com[255];

  // check if there is a savenumber to add and create filename
  if(SaveNumber>=0 && SaveNumber<=9999)	sprintf(fnam,"%s-%.4d%s",root,SaveNumber,extension);
  else if (SaveNumber>9999)				sprintf(fnam,"%s-%d%s",root,SaveNumber,extension);
  else 									sprintf(fnam,"%s%s",root,extension);

  // open file
  fp=fopen(fnam,"wb");

  // check file opened ok, if not give warning and exit function
  if(fp==NULL)
  {
    printf("Problems writing to file %s\n",fnam);
    return;
  }
  //Output alternative header or default header
  if(Statesptr->alt_ascheader==ON)
  {
    for(i=0;i<6;i++) 
    {
      fprintf(fp,"%s",Parptr->ascheader[i]);
    }
  } 
  else
  {
    fprintf(fp,"ncols         %i\n",Parptr->xsz);
    fprintf(fp,"nrows         %i\n",Parptr->ysz);
    fprintf(fp,"xllcorner     %lf\n",Parptr->blx);
    fprintf(fp,"yllcorner     %lf\n",Parptr->bly);
    fprintf(fp,"cellsize      %lf\n",Parptr->dx);
    fprintf(fp,"NODATA_value  %lf\n",NULLVAL);
  }
    for(j=0;j<Parptr->ysz;j++)
    {
      for(i=0;i<Parptr->xsz;i++)
      {
		// calculate depth on floodplain
        temp = data[i+j*Parptr->xsz] - SGCbfH[i+j*Parptr->xsz];
		// if water level is below floodplain set depth to zero
		if (temp <= 0.001) temp = 0.0;
        fprintf(fp,"%.3f\t",temp);
      }
      // output end of line
      fprintf(fp,"\n");
    }
  // close file
  fclose(fp);
  // check if we need to zip the file up
  if(Statesptr->call_gzip==ON) 
  {
    sprintf(tmp_sys_com,"%s%s","gzip -9 -f ",fnam);
    system(tmp_sys_com);
  }
  return;
}
// WRITE FLOODPLAIN DEPTHS FOR SGC MODEL
void write_binrasterfile_SGCf(char *root, int SaveNumber, char *extension, double *data, double *SGCbfH, States *Statesptr,Pars *Parptr)
/*
Purpose: writes data that would be in an ascii raster to a binary file format
Parameters:
char *root		-	root of filename
int SaveNumber	-	save number to add to filename (if <0, will be ignored)
char *extension	-	filename extension text
double *data		-	pointer to data
double *dem		-	pointer to dem
int outflag		-	flag (if 0 = normal, if 1 or 2 indicated fluxes, if 3 indicates special option Water elev ouput DEM+H)
States *Statesptr - pointer to States structure
Pars *Parptr - pointer to Parameters structure
*/
{
  int i,j;
  double temp, no_data=-9999;
  FILE *fp;
  char fnam[800];
  //char tmp_sys_com[255];

  // check if there is a savenumber to add and create filename
  if(SaveNumber>=0 && SaveNumber<=9999)	sprintf(fnam,"%s-%.4d%s",root,SaveNumber,extension);
  else if (SaveNumber>9999)				sprintf(fnam,"%s-%d%s",root,SaveNumber,extension);
  else 									sprintf(fnam,"%s%s",root,extension);
  // open file
  fp=fopen(fnam,"wb");
  // check file opened ok, if not give warning and exit function
  if(fp==NULL)
  {
    printf("Problems writing to file %s\n",fnam);
    return;
  }
  fwrite (&Parptr->xsz, sizeof(int)   , 1, fp );
  fwrite (&Parptr->ysz, sizeof(int)   , 1, fp );
  fwrite (&Parptr->blx, sizeof(double), 1, fp );
  fwrite (&Parptr->bly, sizeof(double), 1, fp );
  fwrite (&Parptr->dx , sizeof(double), 1, fp );
  fwrite (&no_data    , sizeof(double), 1, fp );
  for(j=0;j<Parptr->ysz;j++)
  {
    for(i=0;i<Parptr->xsz;i++)
    {
		// calculate depth on floodplain
        temp = data[i+j*Parptr->xsz] - SGCbfH[i+j*Parptr->xsz];
		// if water level is below floodplain set depth to zero
		if (temp <= 0.0) temp = 0.0;
        fwrite(&temp,sizeof(double),1,fp);
	}
  }
  // close file
  fclose(fp);
  return;
}

void write_regular_output( Fnames *Fnameptr, Solver *Solverptr, States *Statesptr, Pars *Parptr, Arrays *Arrptr, SGCprams *SGCptr  )
{
	int ptr, p0, p1, pq0, i, j, gr0, gr1;
	double w0, w1, A0, A1, zb0, zb1, cn, bf0, bf1, h0, h1, g, dt, m, dh, Sf, R, qc, hflow, Qc;
	double *SGCVx, *SGCVy, *SGCvoutput;

	// write out water depth
    if(Statesptr->save_depth==ON)
    {
		// output binary of ascii rasters
		if (Statesptr->binary_out==ON) write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".wdb",Arrptr->H,Arrptr->DEM,0,Statesptr,Parptr);
		else write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".wd",Arrptr->H,Arrptr->DEM,0,Statesptr,Parptr);
		if (Statesptr->SGC==ON)
		{
			// export floodplain depth for SGC model only
			if (Statesptr->binary_out==ON) write_binrasterfile_SGCf(Fnameptr->resrootname,Parptr->SaveNo, ".wdfpd", Arrptr->H, Arrptr->SGCbfH, Statesptr,Parptr); 
			else write_ascfile_SGCf(Fnameptr->resrootname,Parptr->SaveNo, ".wdfp", Arrptr->H, Arrptr->SGCbfH, Statesptr,Parptr);
		}
	}

    // write out elevation data (note flag=3 so water H is added to DEM)
    if(Statesptr->save_elev==ON)
    {
	  // if using SGC the elevation neds to be from the channel bed
	  if (Statesptr->SGC==ON)
	  {
		  if (Statesptr->binary_out==ON) write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".elevb",Arrptr->H,Arrptr->SGCz,3,Statesptr,Parptr);
		  else write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".elev",Arrptr->H,Arrptr->SGCz,3,Statesptr,Parptr);
	  }
	  else
	  {
		  // output binary of ascii rasters
		  if (Statesptr->binary_out==ON) write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".elevb",Arrptr->H,Arrptr->DEM,3,Statesptr,Parptr);
		  else write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".elev",Arrptr->H,Arrptr->DEM,3,Statesptr,Parptr);
	  }
	}

    // write out x&y "flow fluxes" if required (note flag=1 for Qx and flag=2 for Qy)
    if(Statesptr->save_Qs==ON)
	{
		if (Statesptr->binary_out==ON) // output binary of ascii rasters
		{
			write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".Qxb",Arrptr->Qx,Arrptr->DEM,1,Statesptr,Parptr);
			write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".Qyb",Arrptr->Qy,Arrptr->DEM,2,Statesptr,Parptr);
		}
		else
		{
			write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".Qx",Arrptr->Qx,Arrptr->DEM,1,Statesptr,Parptr);
			write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".Qy",Arrptr->Qy,Arrptr->DEM,2,Statesptr,Parptr);
		}
		if (Statesptr->SGC==ON)
		{
			// additional sub grid output of flow width
			for(i=0;i<Parptr->xsz;i++) for(j=0;j<Parptr->ysz;j++) 
			{
				ptr=i+j*Parptr->xsz;
				w0 = Arrptr->SGCwidth[ptr];
				if (w0 > 0.0 && Arrptr->H[ptr] > 0.0)
				{
					CalcSGC_A(Arrptr->SGCgroup[ptr], Arrptr->H[ptr], Arrptr->SGCbfH[ptr], &A0,  &w0, SGCptr);
					Arrptr->SGCFlowWidth[ptr] = w0;
				}
				else Arrptr->SGCFlowWidth[ptr] = 0.0;
			}
			if (Statesptr->binary_out==ON) 
			{
				write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".Fwidthb",Arrptr->SGCFlowWidth,Arrptr->DEM,0,Statesptr,Parptr);
				// channel Q's
				write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".Qcxb",Arrptr->QxSGold,Arrptr->DEM,1,Statesptr,Parptr);
				write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".Qcyb",Arrptr->QySGold,Arrptr->DEM,2,Statesptr,Parptr);
			}
			else 
			{
				write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".Fwidth",Arrptr->SGCFlowWidth,Arrptr->DEM,0,Statesptr,Parptr);
				// channel Q's
				write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".Qcx",Arrptr->QxSGold,Arrptr->DEM,1,Statesptr,Parptr);
				write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".Qcy",Arrptr->QySGold,Arrptr->DEM,2,Statesptr,Parptr);
			}
		}
	}

    // write out x&y adaptive timesteps if required
    if(Statesptr->voutput==ON)
    {
		if (Statesptr->binary_out==ON) // output binary of ascii rasters
		{
			write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".Vxb",Arrptr->Vx,Arrptr->DEM,1,Statesptr,Parptr);
			write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".Vyb",Arrptr->Vy,Arrptr->DEM,2,Statesptr,Parptr);
		}
		else
		{
			write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".Vx",Arrptr->Vx,Arrptr->DEM,1,Statesptr,Parptr);
			write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".Vy",Arrptr->Vy,Arrptr->DEM,2,Statesptr,Parptr);
		}
	}
	// output SGC channel velocity
	if (Statesptr->SGC==ON && Statesptr->SGCvoutput == ON)
	{
		// create variables 
		SGCVx=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
        SGCVy=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
		SGCvoutput=new double[Parptr->xsz*Parptr->ysz]();
		// calculate hflow in x
		for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz-1;i++)
		{
			p0=i+j*Parptr->xsz;
			p1=i+1+j*Parptr->xsz;
			pq0=i+j*(Parptr->xsz+1)+1;
			w0=Arrptr->SGCwidth[p0]; w1=Arrptr->SGCwidth[p1];
			zb0=Arrptr->SGCz[p0]; zb1=Arrptr->SGCz[p1];
			h0=Arrptr->H[p0]; h1=Arrptr->H[p1];
			g=Solverptr->g; dt=Solverptr->Tstep;
			bf0 = Arrptr->SGCbfH[p0];    bf1 = Arrptr->SGCbfH[p1];
			if (w0>0.0 && w1>0.0)
			{
				// calculate sub grid channel Qc
				hflow=getmax(zb0+h0,zb1+h1)-getmax(zb0,zb1);
				if(hflow>Solverptr->DepthThresh) 
				{
					qc=Arrptr->QxSGold[pq0]; // Get old q in m3/s
					gr0 = Arrptr->SGCgroup[p0];	gr1 = Arrptr->SGCgroup[p1]; // get channel groups
					m = 0.5*(SGCptr->SGCm[gr0]+SGCptr->SGCm[gr1]); // get meander coefficient for sub-grid channel
					if(Arrptr->SGCManningsn!=NULL)   cn=pow(0.5*(Arrptr->SGCManningsn[p0]+Arrptr->SGCManningsn[p1]),2); // this could be pre-calculated
					else cn=SGCptr->SGCn[gr0];
					// calculate Sfch
					dh=zb0+h0-zb1-h1; 
					// Sf=-dh/Parptr->dx;  //CCS_deletion
					Sf=-dh/(Arrptr->dx[p0]*m); 
					// calculate area of cross sectional flow
					CalcSGC_A(gr0, hflow, bf0, &A0,  &w0, SGCptr);
					CalcSGC_A(gr1, hflow, bf1, &A1, &w1, SGCptr);
					// always use smallest flow area
					if (A0 < A1) // select the appropriate slope and area based on the smallest area
					{
						// calculate hydraulic radius for SGC
						R = CalcSGC_R(gr0, hflow, bf0, w0, Arrptr->SGCwidth[p0], A0, SGCptr);
						// Calculate channel flow using inertial wave equation --- in this case Qxold will be in m3s-1 not m2s-1
						Qc =  (qc-g*A0*dt*Sf) / (1+dt*g*cn*fabs(qc) / (pow(R,4.0/3.0)*A0) );
						SGCVx[pq0] = Qc/A0;
					}
					else 
					{
						// calculate hydraulic radius for SGC
						R = CalcSGC_R(gr1, hflow, bf1, w1, Arrptr->SGCwidth[p1], A1, SGCptr);
						// Calculate channel flow using inertial wave equation --- in this case Qxold will be in m3s-1 not m2s-1
						Qc =  (qc-g*A1*dt*Sf) / (1+dt*g*cn*fabs(qc) / (pow(R,4.0/3.0)*A1) );
						SGCVx[pq0] = Qc/A1;
					}
				}
				else SGCVx[pq0] = 0.0; // make sure previous flux is set to zero
			}
			
		}
		//calculate hflow in y
		for(j=0;j<Parptr->ysz-1;j++) for(i=0;i<Parptr->xsz;i++)
		{
			p0=i+j*Parptr->xsz;
			p1=i+(j+1)*Parptr->xsz;
			pq0=i+(j+1)*(Parptr->xsz+1);
			w0=Arrptr->SGCwidth[p0]; w1=Arrptr->SGCwidth[p1];
			zb0=Arrptr->SGCz[p0]; zb1=Arrptr->SGCz[p1];
			h0=Arrptr->H[p0]; h1=Arrptr->H[p1];
			g=Solverptr->g; dt=Solverptr->Tstep;
			bf0 = Arrptr->SGCbfH[p0];    bf1 = Arrptr->SGCbfH[p1];
      		if (w0>0.0 && w1>0.0)
			{
				// calculate sub grid channel Qc
				hflow=getmax(zb0+h0,zb1+h1)-getmax(zb0,zb1);
				if(hflow>Solverptr->DepthThresh) 
				{
					qc=Arrptr->QySGold[pq0]; // Get old q in m3/s
					gr0 = Arrptr->SGCgroup[p0];	gr1 = Arrptr->SGCgroup[p1]; // get channel groups
					m = 0.5*(SGCptr->SGCm[gr0]+SGCptr->SGCm[gr1]); // get meander coefficient for sub-grid channel
					// check if using distributed Mannings n
					if(Arrptr->SGCManningsn!=NULL)   cn=pow(0.5*(Arrptr->SGCManningsn[p0]+Arrptr->SGCManningsn[p1]),2); // this could be pre-calculated
					else cn=SGCptr->SGCn[gr0];
					// calculate Sfch
					dh=zb0+h0-zb1-h1; 
					//Sf=-dh/Parptr->dx; //CCS_deletion
					Sf=-dh/(Arrptr->dy[p0]*m);
					// calculate area of cross sectional flow
					CalcSGC_A(gr0, hflow, bf0, &A0,  &w0, SGCptr);
					CalcSGC_A(gr1, hflow, bf1, &A1, &w1, SGCptr);
					// always use smallest flow area
					if (A0 < A1) // select the appropriate slope and area based on the smallest area
					{
						// calculate hydraulic radius for SGC
						R = CalcSGC_R(gr0, hflow, bf0, w0, Arrptr->SGCwidth[p0], A0, SGCptr);
						// Calculate channel flow using inertial wave equation --- in this case Qxold will be in m3s-1 not m2s-1
						Qc =  (qc-g*A0*dt*Sf) / (1+dt*g*cn*fabs(qc) / (pow(R,4.0/3.0)*A0) );
						SGCVy[pq0] = Qc/A0;
					}
					else 
					{
						// calculate hydraulic radius for SGC
						R = CalcSGC_R(gr1, hflow, bf1, w1, Arrptr->SGCwidth[p1], A1, SGCptr);
						// Calculate channel flow using inertial wave equation --- in this case Qxold will be in m3s-1 not m2s-1
						Qc =  (qc-g*A1*dt*Sf) / (1+dt*g*cn*fabs(qc) / (pow(R,4.0/3.0)*A1) );
						SGCVy[pq0] = Qc/A1;
					}
				}
				else SGCVy[pq0] = 0.0; // make sure previous flux is set to zero
			}
		}
		// Get largest velocity from cell faces and place at cell centre.
		for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++) 
		{
			SGCvoutput[i+j*Parptr->xsz] = getmax(getmax(fabs(SGCVx[i+j*(Parptr->xsz+1)]), fabs(SGCVx[i+1+j*(Parptr->xsz+1)])), getmax(fabs(SGCVy[i+j*(Parptr->xsz+1)]), fabs(SGCVy[i+(j+1)*(Parptr->xsz+1)])));
		}

		if (Statesptr->binary_out==ON) // output binary of ascii rasters
		{
			write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".SGCVxb",SGCVx,Arrptr->DEM,1,Statesptr,Parptr);
			write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".SGCVyb",SGCVy,Arrptr->DEM,2,Statesptr,Parptr);
			write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".SGCVcb",SGCvoutput,Arrptr->DEM,0,Statesptr,Parptr);
		}
		else
		{
			write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".SGCVx",SGCVx,Arrptr->DEM,1,Statesptr,Parptr);
			write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".SGCVy",SGCVy,Arrptr->DEM,2,Statesptr,Parptr);
			write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".SGCVc",SGCvoutput,Arrptr->DEM,0,Statesptr,Parptr);
		}

	}
    // write out x&y Q limits if required
    if(Statesptr->save_QLs==ON)
    {
	    if (Statesptr->binary_out==ON) // output binary of ascii rasters
		{
			write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".QLxb",Arrptr->LimQx,Arrptr->DEM,1,Statesptr,Parptr);
			write_binrasterfile(Fnameptr->resrootname,Parptr->SaveNo,".QLyb",Arrptr->LimQy,Arrptr->DEM,2,Statesptr,Parptr);
		}
		else
		{
			write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".QLx",Arrptr->LimQx,Arrptr->DEM,1,Statesptr,Parptr);
			write_ascfile(Fnameptr->resrootname,Parptr->SaveNo,".QLy",Arrptr->LimQy,Arrptr->DEM,2,Statesptr,Parptr);
		}
	}

    // output regular files
    fileoutput(Fnameptr,Statesptr,Parptr,Arrptr);

	return;
}