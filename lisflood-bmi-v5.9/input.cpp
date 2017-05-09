/*
*****************************************************************************
FILE INPUT
---------------------

A number of functions used to control the file input to LISFLOOD-FP. A short
definition of each is outlined below.

*****************************************************************************
*/
#include "lisflood.h"
//-----------------------------------------------------------------------------
// LOAD STAGE DATA
void LoadStages(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,Stage *Locptr,int *verbose)
{
  //Added by Matt Wilson, 5 Apr 2004
  //Provides functionality to output regular point measurements of water stage

  int i;
  FILE *fp;

  fp=fopen(Fnameptr->stagefilename,"r");
  if(fp==NULL)
  {
    if(*verbose==ON) printf("Stages off\n");
    Statesptr->save_stages=OFF;
    return;
  }
  if(*verbose==ON) printf("\nLoading stage information:\t%s\n",Fnameptr->stagefilename);

  fscanf(fp,"%d",&Locptr->Nstages);
  fgetc(fp); // Retrieve closing EOL

  Locptr->stage_loc_x=new double[Locptr->Nstages];
  Locptr->stage_loc_y=new double[Locptr->Nstages];
  Locptr->stage_grid_x=new int[Locptr->Nstages];
  Locptr->stage_grid_y=new int[Locptr->Nstages];
  Locptr->stage_check=new int[Locptr->Nstages];

  //scan x,y locations from file
  for(i=0;i<Locptr->Nstages;i++)
  {
    fscanf(fp,"%lf",&Locptr->stage_loc_x[i]);
    fscanf(fp,"%lf",&Locptr->stage_loc_y[i]);
  }
  //convert coordinates to grid cell numbers
  for(i=0;i<Locptr->Nstages;i++)
  {
    Locptr->stage_grid_x[i]=int(floor((Locptr->stage_loc_x[i]-Parptr->blx)/Parptr->dx));
    Locptr->stage_grid_y[i]=Parptr->ysz-1-(int(floor((Locptr->stage_loc_y[i]-Parptr->bly)/Parptr->dy)));
    Locptr->stage_check[i]=1;
  }
  //check for off-image values
  for(i=0;i<Locptr->Nstages;i++)
  {
    if(Locptr->stage_grid_x[i]<0 || Locptr->stage_grid_x[i]>=Parptr->xsz || Locptr->stage_grid_y[i]<0 || Locptr->stage_grid_y[i]>=Parptr->ysz)
    {
      Locptr->stage_check[i]=0;
      if(*verbose==ON) printf("WARNING: Stage off-image: %d\n",i+1);
    }
  }
  // optional velocity output
  if (Statesptr->voutput==ON)
  {
    //convert coordinates to grid edge numbers same in x and y!
    Locptr->vel_grid_xy=new int[Locptr->Nstages];
    for(i=0;i<Locptr->Nstages;i++)
    {
      Locptr->vel_grid_xy[i]= Locptr->stage_grid_x[i] + Locptr->stage_grid_y[i]*(Parptr->xsz+1);
    }
  }
  if(*verbose==ON) printf("Done.\n\n");

  fclose(fp);
  return;
}
//-----------------------------------------------------------------------------
// LOAD WEIR DATA
void LoadWeir(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,Arrays *Arrptr,int *verbose)
{
  //One-directional flow functionality added by Matt Wilson 13 Feb 2504

  FILE *fp;
  int nw,i,j,xi,yi, p0, p1;
  double x,y, z0,z1;
  char char_tmp[10],buff[80];
  char tag_w1[] = "W";
  char tag_w2[] = "w";
  char tag_e1[] = "E";
  char tag_e2[] = "e";
  char tag_s1[] = "S";
  char tag_s2[] = "s";
  char tag_n1[] = "N";
  char tag_n2[] = "n";
  //tags for one-directional flow (culverts)
  char tag_wf1[] = "WF";
  char tag_wf2[] = "wf";
  char tag_ef1[] = "EF";
  char tag_ef2[] = "ef";
  char tag_sf1[] = "SF";
  char tag_sf2[] = "sf";
  char tag_nf1[] = "NF";
  char tag_nf2[] = "nf";
  //tags for bridge/culvert
  char tag_wb1[] = "WB";
  char tag_wb2[] = "wb";
  char tag_eb1[] = "EB";
  char tag_eb2[] = "eb";
  char tag_sb1[] = "SB";
  char tag_sb2[] = "sb";
  char tag_nb1[] = "NB";
  char tag_nb2[] = "nb";

  fp=fopen(Fnameptr->weirfilename,"r");
  if(fp==NULL)
  {
    if(*verbose==ON) printf("Weirs off\n");
    return;
  }
  if(*verbose==ON) printf("Loading weir information:\t%s\n",Fnameptr->weirfilename);

  Statesptr->weirs=ON;

  j=0;
  do{buff[j]=fgetc(fp);} while(buff[j++]!='\n');
  buff[j]='\0';
  sscanf(buff,"%i",&nw);

  Arrptr->Weir_hc=new double[nw];
  Arrptr->Weir_Cd=new double[nw];
  Arrptr->Weir_m=new double[nw];
  Arrptr->Weir_w=new double[nw];
  Arrptr->Weir_Typ=new int[nw]; // type of structure... weir = 0, bridge = 1;

  Arrptr->Weir_Fixdir=new int[nw];   // Fixed flow directions
  Arrptr->Weir_Identx=new int[(Parptr->xsz+1)*(Parptr->ysz+1)];
  Arrptr->Weir_Identy=new int[(Parptr->xsz+1)*(Parptr->ysz+1)];

  // Defalut to -1 for no weir link, and 0 for fixed flow direction
  for(i=0;i<=Parptr->xsz;i++) for(j=0;j<=Parptr->ysz;j++)
                              {
                                Arrptr->Weir_Identx[i+j*(Parptr->xsz+1)]=-1;
                                Arrptr->Weir_Identy[i+j*(Parptr->xsz+1)]=-1;
                              }

  for(i=0;i<nw;i++)
  {
    j=0;	   // load buffer until EOL
    do{buff[j]=fgetc(fp);} while(buff[j++]!='\n');
    buff[j]='\0';  // Finish off string

    if(sscanf(buff,"%lf %lf %s %lf %lf %lf %lf",
              &x,&y,char_tmp,Arrptr->Weir_Cd+i,Arrptr->Weir_hc+i,Arrptr->Weir_m+i,Arrptr->Weir_w+i)!=7)
      Arrptr->Weir_w[i]=Parptr->dx;

    xi=(int)((x-Parptr->blx)/Parptr->dx);
    yi=(int)((Parptr->tly-y)/Parptr->dy);

    // Unfixed flow direction = 0.
    if(strcmp(char_tmp,tag_w1)==0 || strcmp(char_tmp,tag_w2)==0)
    {
      Arrptr->Weir_Identx[xi+1+yi*(Parptr->xsz+1)]=i;
      Arrptr->Weir_Fixdir[i]=0;
      Arrptr->Weir_Typ[i]=0;
    }
    if(strcmp(char_tmp,tag_e1)==0 || strcmp(char_tmp,tag_e2)==0)
    {
      Arrptr->Weir_Identx[xi+yi*(Parptr->xsz+1)]=i;
      Arrptr->Weir_Fixdir[i]=0;
      Arrptr->Weir_Typ[i]=0;
    }
    if(strcmp(char_tmp,tag_s1)==0 || strcmp(char_tmp,tag_s2)==0)
    {
      Arrptr->Weir_Identy[xi+yi*(Parptr->xsz+1)]=i;
      Arrptr->Weir_Fixdir[i]=0;
      Arrptr->Weir_Typ[i]=0;
    }
    if(strcmp(char_tmp,tag_n1)==0 || strcmp(char_tmp,tag_n2)==0)
    {
      Arrptr->Weir_Identy[xi+(yi+1)*(Parptr->xsz+1)]=i;
      Arrptr->Weir_Fixdir[i]=0;
      Arrptr->Weir_Typ[i]=0;
    }
    // Control tags for one-directional flow (culverts)
    // Fixed flow directions: N = 1, E = 2, S = 3, W = 4.
    if(strcmp(char_tmp,tag_wf1)==0 || strcmp(char_tmp,tag_wf2)==0)
    {
      Arrptr->Weir_Identx[xi+1+yi*(Parptr->xsz+1)]=i;
      Arrptr->Weir_Fixdir[i]=4;
      Arrptr->Weir_Typ[i]=0;
    }
    if(strcmp(char_tmp,tag_ef1)==0 || strcmp(char_tmp,tag_ef2)==0)
    {
      Arrptr->Weir_Identx[xi+yi*(Parptr->xsz+1)]=i;
      Arrptr->Weir_Fixdir[i]=2;
      Arrptr->Weir_Typ[i]=0;
    }
    if(strcmp(char_tmp,tag_sf1)==0 || strcmp(char_tmp,tag_sf2)==0)
    {
      Arrptr->Weir_Identy[xi+yi*(Parptr->xsz+1)]=i;
      Arrptr->Weir_Fixdir[i]=3;
      Arrptr->Weir_Typ[i]=0;
    }
    if(strcmp(char_tmp,tag_nf1)==0 || strcmp(char_tmp,tag_nf2)==0)
    {
      Arrptr->Weir_Identy[xi+(yi+1)*(Parptr->xsz+1)]=i;
      Arrptr->Weir_Fixdir[i]=1;
      Arrptr->Weir_Typ[i]=0;
    }
    // control tags for bridge
    if(strcmp(char_tmp,tag_wb1)==0 || strcmp(char_tmp,tag_wb2)==0)
    {
      Arrptr->Weir_Identx[xi+1+yi*(Parptr->xsz+1)]=i;
      Arrptr->Weir_Fixdir[i]=0;
      Arrptr->Weir_Typ[i]=1;
    }
    if(strcmp(char_tmp,tag_eb1)==0 || strcmp(char_tmp,tag_eb2)==0)
    {
      Arrptr->Weir_Identx[xi+yi*(Parptr->xsz+1)]=i;
      Arrptr->Weir_Fixdir[i]=0;
      Arrptr->Weir_Typ[i]=1;
    }
    if(strcmp(char_tmp,tag_sb1)==0 || strcmp(char_tmp,tag_sb2)==0)
    {
      Arrptr->Weir_Identy[xi+yi*(Parptr->xsz+1)]=i;
      Arrptr->Weir_Fixdir[i]=0;
      Arrptr->Weir_Typ[i]=1;
    }
    if(strcmp(char_tmp,tag_nb1)==0 || strcmp(char_tmp,tag_nb2)==0)
    {
      Arrptr->Weir_Identy[xi+(yi+1)*(Parptr->xsz+1)]=i;
      Arrptr->Weir_Fixdir[i]=0;
      Arrptr->Weir_Typ[i]=1;
    }

    // now a check to make sure that Arrptr->Weir_hc is greater than the ground elevation,
    // this is especally important for the SGC model where bed elevations may change.

    // first index the cell of the weir and get z0
    p0=xi+yi*Parptr->xsz;
    if (Statesptr->SGC==ON && Arrptr->SGCwidth[p0] > 0.0) z0 = Arrptr->SGCz[p0];
    else z0 = Arrptr->DEM[p0];
    // Then index the cell it flows from and get the elevation
    if (strcmp(char_tmp,tag_w1)==0 || strcmp(char_tmp,tag_w2)==0 || strcmp(char_tmp,tag_wf1)==0 || strcmp(char_tmp,tag_wf2)==0 || strcmp(char_tmp,tag_wb1)==0 || strcmp(char_tmp,tag_wb2)==0)
    {
      p1=xi+1+yi*Parptr->xsz;
    }
    if (strcmp(char_tmp,tag_e1)==0 || strcmp(char_tmp,tag_e2)==0 || strcmp(char_tmp,tag_ef1)==0 || strcmp(char_tmp,tag_ef2)==0 || strcmp(char_tmp,tag_eb1)==0 || strcmp(char_tmp,tag_eb2)==0)
    {
      p1=xi-1+yi*Parptr->xsz;
    }
    if (strcmp(char_tmp,tag_s1)==0 || strcmp(char_tmp,tag_s2)==0 || strcmp(char_tmp,tag_sf1)==0 || strcmp(char_tmp,tag_sf2)==0 || strcmp(char_tmp,tag_sb1)==0 || strcmp(char_tmp,tag_sb2)==0)
    {
      p1=xi+(yi-1)*Parptr->xsz;
    }
    if (strcmp(char_tmp,tag_n1)==0 || strcmp(char_tmp,tag_n2)==0 || strcmp(char_tmp,tag_nf1)==0 || strcmp(char_tmp,tag_nf2)==0 || strcmp(char_tmp,tag_nb1)==0 || strcmp(char_tmp,tag_nb2)==0)
    {
      p1=xi+(yi+1)*Parptr->xsz;
    }
    if (Statesptr->SGC==ON && Arrptr->SGCwidth[p1] > 0.0) z1 = Arrptr->SGCz[p1];
    else z1 = Arrptr->DEM[p1];

    // now work out if either of the elevations (z0,z1) are above the weir crest hight.
    if (Arrptr->Weir_hc[i] < z0 || Arrptr->Weir_hc[i] < z1)
    {
      if(*verbose==ON)
      {
        if (Arrptr->Weir_Typ[i] == 0)
        {
          printf("WARNING: Weir crest height is below DEM\n");
          // for sub-grid model increase the crest height
          //if(Statesptr->SGC==ON)
          //{
          Arrptr->Weir_hc[i] = getmax(z0,z1);
          if(*verbose==ON) printf("Weir number %i crest height increased to %.3f m\n", i,Arrptr->Weir_hc[i]);
          //}
        }
        //else
        //{
        //printf("WARNING: Bridge soffit height is below DEM converted to weir!!\n");
        //Arrptr->Weir_Typ[i] = 0;
        //}
      }
    }
    // need check for bridge less than or equal to subgrid width ?.......
  }

  fclose(fp);

  if(*verbose==ON) printf("Done.\n\n");
  return;
}
//-----------------------------------------------------------------------------
/* LOAD RIVER CHANNEL NETWORK - CCS
   Added by Chris Sampson January 2011
   Loads one or more river systems into the ChannelSegments vector, using the RiversIndex vector to keep track of where each river is located
   within ChannelSegents.  This is done by recording the size of ChannelSegments after each river is loaded as an int in RiversIndex. Example:

   River Thames has 5 segments; River Severn has 4 segments.

   After first LoadRiver loop:
   ChannelSegments[0][1][2][3][4]					RiversIndex[0]
   <....Thames...>							   <5>

   After second LoadRiver loop:
   ChannelSegments[0][1][2][3][4][5][6][7][8]		RiversIndex[0][1]
   <....Thames...><..Severn..>				   <5><9>

   Note: remember that the size recorded in RiversIndex is always 1 greater than the max index of ChannelSegments because size starts at 1
   whereas the index starts at 0!

   Just as with the old system, we use a pointer called CSTypePtr throughout the rest of the model code to point at ChannelSegments[0].  We also have
   a new pointer called RiversIndexPtr to point RiversIndex[0].  Therefore CSTypePtr+5 points at ChannelSegments[5] etc.

   An extra trick with the vectors is that we can have a pointer to the vector itself (rather than a particular element of the vector as above).  These are
   called ChannelSegmentsVecPtr and RiversIndexVecPtr, and are used primarily to construct the vectors using member functions such as push.back.  They are
   also sometimes employed to determine the size loops need to be using the size() member function.  For example, if we have loaded 5 rivers, then:
   RiversIndexVecPtr->size() will equal 5. If each of these rivers contained 5 channel segments then ChannelSegmentsVecPtr->size() would equal 25.

*/

void LoadRiverNetwork(Fnames *Fnameptr,States *Statesptr,Pars *Parptr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr,Arrays *Arrptr,vector<QID7_Store> *QID7_Vec_Ptr, vector<int> *RiversIndexVecPtr, int *verbose)
{
  FILE *rfp;  // local file pointer
  int i,n;
  int tmp_size;

  if(Statesptr->multiplerivers==0)
  {
    LoadRiver(Fnameptr,Statesptr,Parptr,ChannelSegmentsVecPtr,Arrptr,QID7_Vec_Ptr, RiversIndexVecPtr, verbose); // Call LoadRiver once.
    tmp_size = ChannelSegmentsVecPtr->size();
    RiversIndexVecPtr->push_back(tmp_size); // CCS
  }
  else if(Statesptr->multiplerivers==1)
  {
    rfp=fopen(Fnameptr->multiriverfilename,"r");
    fscanf(rfp,"%i",&n);
    if(*verbose==ON) printf("Loading %i Rivers\n\n",n);

    for(i=0; i < n; i++)
    {
      fscanf(rfp,"%s",Fnameptr->rivername); //Scan the next .river filename from the .rivers file and asign it to Fnameptr->rivername.
      LoadRiver(Fnameptr,Statesptr,Parptr,ChannelSegmentsVecPtr,Arrptr,QID7_Vec_Ptr, RiversIndexVecPtr, verbose); // Call LoadRiver in loop.
      tmp_size = ChannelSegmentsVecPtr->size();
      RiversIndexVecPtr->push_back(tmp_size); /* Builds the RiversIndex vector so we know where one river stops and the next starts within
                                                 the ChannelSegments vector. */
    }

    if(*verbose==ON) printf("%i Rivers Loaded Successfully.\n\n",n);
  }

  return;
}
//-----------------------------------------------------------------------------
// LOAD RIVER DATA FROM FILE
void LoadRiver(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,vector<ChannelSegmentType> *ChannelSegmentsVecPtr,Arrays *Arrptr,vector<QID7_Store> *QID7_Vec_Ptr, vector<int> *RiversIndexVecPtr, int *verbose)
{
  FILE *fp;  // local file pointer
  int npoints,*xpi,*ypi;
  int *trib; // temp xs array to record any trib connections
  double *xp,*yp,*wp,*np,*hp,*cp,*rp,*ap, total_length=0.0,*qp;
  char buff[800],buff2[800],buff3[800];
  int i,j,pi,pj,ni,nj,oldpi,oldpj,i1,i2;
  double tmp1,tmp2,tmp3;
  double grad; // temporary gradient calculation variable eventually stored in csp->Shalf
  char *Q_Name_tmp;
  char buffer[80];
  int *Q_Ident_tmp;
  int count,chseg,SegOutNo,tmp_int;

  // MSH: csp is a utility pointer,

  fp=fopen(Fnameptr->rivername,"r");
  if(fp==NULL) return;

  Statesptr->ChannelPresent=ON;

  if(*verbose==ON) printf("Loading channel information:\t%s\n",Fnameptr->rivername);

  fscanf(fp,"%s",buffer);
  if(!strcmp(buffer,"Tribs")||!strcmp(buffer,"tribs")||!strcmp(buffer,"TRIBS"))
  {
    Statesptr->TribsPresent=ON;

    // MSH: Since we haven't allocated the memory yet, we can't assign the number of channel segments in
    // the first element of the ChannelSegments array - so read into a temp variable
    fscanf(fp,"%i",&tmp_int);
    if(*verbose==ON) printf("%i Channel Segments\n",tmp_int);
  }
  else
  {
    rewind(fp);
    tmp_int=1;
    if(*verbose==ON) printf("%i Channel Segment\n",tmp_int);
  }

  for(chseg=0;chseg<tmp_int;chseg++) // CCS Slight reorganisation of old LoadRiver function but fundamentally unchanged.
  {
    ChannelSegmentType tmp_chan; // CCS
    tmp_chan.Next_Segment=tmp_chan.Next_Segment_Loc=-1; // CCS from old code
    tmp_chan.N_Channel_Segments=tmp_int; // CCS
    ChannelSegmentType *csp=&tmp_chan; // CCS

    fscanf(fp,"%i",&npoints);
    if(*verbose==ON) printf("%i points in channel segment %i\n",npoints,chseg);

    //setup local temporary arrays, Note the () at the end ensures all elements are initialised to zero
    xp=new double[npoints] ();
    yp=new double[npoints] ();
    wp=new double[npoints] ();
    np=new double[npoints] ();
    hp=new double[npoints] ();
    cp=new double[npoints] ();
    qp=new double[npoints] ();
    rp=new double[npoints] (); // chainage ratio between entered cross sections and cell chainage
    ap=new double[npoints] (); // actual entered xs chainage
    trib=new int[npoints] ();
    Q_Name_tmp=new char[npoints*80] ();
    Q_Ident_tmp=new int[npoints] ();

    for(i=0;i<npoints;i++)
    {
      // set default qp and trib value (all the other arrays are pre zeroed with new command and end brackets () )
      qp[i]=-1;
      trib[i]=-1;

      fscanf(fp,"%lf %lf",xp+i,yp+i); // Load x,y values.

      // load buffer until EOL
      j=0;
      do{buff[j]=fgetc(fp);} while(buff[j++]!='\n');
      buff[j]='\0';									// Finish off string
      if(sscanf(buff,"%lf%lf%lf",&tmp1,&tmp2,&tmp3)==3)
      {        										// Only store values if 3 reads successful
        wp[i]=tmp1;
        np[i]=tmp2;
        hp[i]=tmp3;
        if(*verbose==ON)
          printf("Xsec %4i\tw=%8.3f n=%5.3f z=%6.3f\n",i,wp[i],np[i],hp[i]);
      }

      if(sscanf(buff,"%lf%lf%lf%s%s",&tmp1,&tmp1,&tmp1,buff2,buff3)==5  // 4+5th item found must be Q_Ident
         || sscanf(buff,"%s%s",buff2,buff3)==2) // OR No channel info - just Q in
      {
        if(!strcmp(buff2,"QFIX")||!strcmp(buff2,"qfix"))
        {
          Q_Ident_tmp[i]=4;
          sscanf(buff3,"%lf",qp+i);
          if(*verbose==ON) printf("Xsec %4i\tQFIX at %7.2f\n",i,qp[i]);
        }
        if(!strcmp(buff2,"QVAR")||!strcmp(buff2,"qvar"))
        {
          Q_Ident_tmp[i]=5;
          strcpy(Q_Name_tmp+i*80,buff3);
          if(*verbose==ON) printf("Xsec %4i\tQVAR from bdy file %s\n",i,Q_Name_tmp+i*80);
        }
        if(!strcmp(buff2,"QOUT")||!strcmp(buff2,"qout"))
        {
          Q_Ident_tmp[i]=6;
          sscanf(buff3,"%i",&SegOutNo);
          if(*verbose==ON) printf("Xsec %4i\tQOUT from segment %i discharges into segment %i\n",i,chseg,SegOutNo);
        }
        if(!strcmp(buff2,"TRIB")||!strcmp(buff2,"trib"))
        {
          Q_Ident_tmp[i]=7;
          sscanf(buff3,"%i",&SegOutNo);
          if(*verbose==ON) printf("Xsec %4i\tQin from segment %i\n",i,SegOutNo);
          trib[i]=SegOutNo;
        }
        if(!strcmp(buff2,"FREE")||!strcmp(buff2,"free"))
          // normal depth based on slope, if -1 then use last channel segment slope else use slope supplied
          // NOT fully working for Diffusive - stability issues
        {
          Q_Ident_tmp[i]=1;
          sscanf(buff3,"%lf",qp+i);
          if(qp[i] < -0.999) // ie -1 (done like this as double)
          {
            if(*verbose==ON) printf("Xsec %4i\tFREE using end slope\n",i);
          }
          else
          {
            if(*verbose==ON) printf("Xsec %4i\tFREE using slope %7.4f\n",i,qp[i]);
          }
        }
        if(!strcmp(buff2,"HFIX")||!strcmp(buff2,"hfix"))
        {
          Q_Ident_tmp[i]=2;
          sscanf(buff3,"%lf",qp+i);
          if(*verbose==ON) printf("Xsec %4i\tHFIX at %7.2f\n",i,qp[i]);
        }
        if(!strcmp(buff2,"HVAR")||!strcmp(buff2,"hvar"))
        {
          Q_Ident_tmp[i]=3;
          strcpy(Q_Name_tmp+i*80,buff3);
          if(*verbose==ON) printf("Xsec %4i\tHVAR from bdy file %s\n",i,Q_Name_tmp+i*80);
        }
        if(!strcmp(buff2,"RATE")||!strcmp(buff2,"rate"))
          // NOT fully working for Diffusive - stability issues
        {
          Q_Ident_tmp[i]=8;
          strcpy(Q_Name_tmp+i*80,buff3);
          if(*verbose==ON) printf("Xsec %4i\tRATE from bdy file %s\n",i,Q_Name_tmp+i*80);
        }
      }
    }

    if(*verbose==ON) printf("Channel data read for segment %i - interpolating values.\n",chseg);



    // Estimate number of channel pixels - to ensure we allocate enough memory for temporary arrays xpi,ypi
    // total length divided by cell size and then double this value.
    ap[0]=0;
    for(i=0;i<npoints-1;i++)
    {
      // calc straight line chainage between entered cross sections - used for cell independent chainage calcs. Add up chainage
      ap[i+1]= ap[i]+sqrt((xp[i+1]-xp[i])*(xp[i+1]-xp[i])+(yp[i+1]-yp[i])*(yp[i+1]-yp[i]));
    }
    xpi=new int[int(2.0*ap[npoints-1]/Parptr->dx)] ();
    ypi=new int[int(2.0*ap[npoints-1]/Parptr->dx)] ();



    // Insert channel into DEM grid
    oldpi=(int)((xp[0]-Parptr->tlx)/Parptr->dx);
    oldpj=(int)((Parptr->tly-yp[0])/Parptr->dy);
    total_length=0.0;

    count=0;
    for(i=1;i<npoints;i++)
    {
      pi=(int)((xp[i]-Parptr->tlx)/Parptr->dx);
      pj=(int)((Parptr->tly-yp[i])/Parptr->dy);
      for(j=0;j<=1000;j++)					// Take very small steps and insert channel
      {															// whenever x,y position changes
        ni=oldpi+((pi-oldpi)*j/1000) ;
        nj=oldpj+((pj-oldpj)*j/1000) ;
        if(ni>=0 && ni<Parptr->xsz && nj>=0 && nj<Parptr->ysz) // check it stays within DEM
        {
          if(count==0)								// Always insert first point
          {
            xpi[count]=ni;
            ypi[count]=nj;
            Arrptr->ChanMask[ni+nj*Parptr->xsz]=1; // mark mask with value of 1 - will renumber later in order
            count++;
          }
          else if(ni!=xpi[count-1] || nj!=ypi[count-1]) // if grid location changes
          {
            if (Arrptr->ChanMask[ni+nj*Parptr->xsz]==-1)   // channel mask not set
            {
              xpi[count]=ni;
              ypi[count]=nj;
              Arrptr->ChanMask[ni+nj*Parptr->xsz]=1; // mark mask with value of 1 - will renumber later in order
              total_length+=Parptr->dx*sqrt(pow(ni-xpi[count-1],(2.0))+pow(nj-ypi[count-1],(2.0)));
              count++;
            }
            else
            {
              // channel mask set, so likely that it is trib junction point.
              // NOTE, cannot have crossing channels !!!
              // DO NOT mark mask with value of 1 - as this is the junction
              // of the trib with main channel so is already marked for main channel
              xpi[count]=ni;
              ypi[count]=nj;
              total_length+=Parptr->dx*sqrt(pow( (double) (ni-xpi[count-1]),2)+pow( (double) (nj-ypi[count-1]),2));
              count++;
            }
          }
        }
      }
      oldpi=pi;
      oldpj=pj;
      cp[i]=total_length;
      rp[i]=(cp[i]-cp[i-1])/(ap[i]-ap[i-1]);
    }
    csp->chsz=count;

    if (count==0) printf("\nWARNING: no overlap with DEM cells for channel %i.\n",chseg);

    // Set up other channel rasters and fill in values
    csp->ChanX=new int[csp->chsz] ();
    csp->ChanY=new int[csp->chsz] ();
    csp->Chandx=new double[csp->chsz] ();
    csp->Shalf=new double[csp->chsz] ();
    csp->A=new double[csp->chsz] ();
    csp->NewA=new double[csp->chsz] ();
    csp->Chainage=new double[csp->chsz] ();
    csp->ChanQ=new double[csp->chsz] (); // only used to record Q values for output in profile - not used in calc
    csp->ChanWidth=new double[csp->chsz] ();
    csp->ChanN=new double[csp->chsz] ();
    csp->Q_Val=new double[csp->chsz] ();
    csp->BankZ=new double[csp->chsz] ();
    csp->Q_Name=new char[csp->chsz*80] ();
    csp->Q_Ident=new int[csp->chsz] ();

    for(i=0;i<csp->chsz;i++)
    {
      csp->ChanX[i]=xpi[i];
      csp->ChanY[i]=ypi[i];
    }


    // Find chainage and dx along channel
    for(i=0;i<csp->chsz-1;i++)
      csp->Chandx[i]=sqrt(pow(Parptr->dx*(csp->ChanX[i]-csp->ChanX[i+1]),2)+
                          pow(Parptr->dy*(csp->ChanY[i]-csp->ChanY[i+1]),2));

    // assume dx for last cell is same as last segment
    csp->Chandx[csp->chsz-1]=sqrt(pow(Parptr->dx*(csp->ChanX[csp->chsz-1]-csp->ChanX[csp->chsz-2]),2)+
                                  pow(Parptr->dy*(csp->ChanY[csp->chsz-1]-csp->ChanY[csp->chsz-2]),2));
    csp->Chainage[0]=0;

    // add up dx to get chainage
    for(i=1;i<csp->chsz;i++)csp->Chainage[i]=csp->Chainage[i-1]+csp->Chandx[i-1];

    // Fill in channel mask
    for(i=0;i<csp->chsz;i++)
    {
      pi=csp->ChanX[i];
      pj=csp->ChanY[i];

      // renumber channel mask
      Arrptr->ChanMask[pi+pj*Parptr->xsz]=i;
      // set bank level
      csp->BankZ[i]=Arrptr->DEM[pi+pj*Parptr->xsz];
      // mark segment mask
      Arrptr->SegMask[pi+pj*Parptr->xsz]=chseg;
    }

    // Adjust chainage calcs so that it is independent of cell size. ####
    if(Statesptr->chainagecalc == ON)
    {
      if(*verbose==ON) printf("Cell size independent channel chainage calculations are ON.\n");
      // adjust each dx by ratio calculated previously
      for(i=0;i<npoints;i++)
      {
        for(j=0;j<csp->chsz;j++)
        {
          if(csp->Chainage[j]>=cp[i] && csp->Chainage[j]<(cp[i+1]))
          {
            csp->Chandx[j]=csp->Chandx[j]/rp[i+1]; // adjust by ratio calculated previously
          }
        }
      }
      // calc last seg dx as same as penultimate
      csp->Chandx[csp->chsz-1]=csp->Chandx[csp->chsz-2];

      // add up chainage again
      for(i=1;i<csp->chsz;i++)csp->Chainage[i]=csp->Chainage[i-1]+csp->Chandx[i-1];

      // also adjust original xs chainage as this was based on centre cell distance
      for(i=1;i<npoints;i++) cp[i]=ap[i];
    }
    else if(*verbose==ON) printf("Cell size independent channel chainage calculations are OFF.\n");




    // Interpolate width, Mannings n
    i1=0;i2=0;
    while(i2<npoints-1)
    {
      for(i2=i1+1;i2<npoints;i2++)
      {
        if(wp[i2]>0.0) break; // loop until next non zero value found
      }
      for(i=0;i<csp->chsz;i++)
      {
        if(csp->Chainage[i]>=cp[i1] && csp->Chainage[i]<=(cp[i2]+0.1)) // add 0.1m to cp[] to get last point
        {
          pi=csp->ChanX[i];
          pj=csp->ChanY[i];
          csp->ChanWidth[i]=wp[i1]+(wp[i2]-wp[i1])*(csp->Chainage[i]-cp[i1])/(cp[i2]-cp[i1]);
          csp->ChanN[i]=np[i1]+(np[i2]-np[i1])*(csp->Chainage[i]-cp[i1])/(cp[i2]-cp[i1]);
          if (i==csp->chsz-1 && chseg != 0)
          {
            // special case where trib junction with main channel. Do not overwrite main channel bed elevation but
            // otherwise record channel info for trib end point in dummy node (allows channel BC link for diffusive
            // and correct slope calc for kinematic).
            csp->JunctionDEM=hp[i1]+(hp[i2]-hp[i1])*(csp->Chainage[i]-cp[i1])/(cp[i2]-cp[i1]);
          }
          else
          {
            Arrptr->DEM[pi+pj*Parptr->xsz]=hp[i1]+(hp[i2]-hp[i1])*(csp->Chainage[i]-cp[i1])/(cp[i2]-cp[i1]);
          }
        }
      }
      i1=i2;
    }


    // Interpolate slope
    for(i=0;i<csp->chsz;i++)
    {
      // find gradient of segment
      if(chseg==0) // main channel
      {
        if(i==csp->chsz-1) // last point
        {
          // special case for last point as we can only know the slope of the segment behind it
          grad=(Arrptr->DEM[csp->ChanX[i]+csp->ChanY[i]*Parptr->xsz]-Arrptr->DEM[csp->ChanX[i-1]+csp->ChanY[i-1]*Parptr->xsz])
            /csp->Chandx[i-1];
        }
        else // all other points
        {
          grad=(Arrptr->DEM[csp->ChanX[i+1]+csp->ChanY[i+1]*Parptr->xsz]-Arrptr->DEM[csp->ChanX[i]+csp->ChanY[i]*Parptr->xsz])
            /csp->Chandx[i];
        }
      }
      else // trib special case at end due to junction
      {
        if(i==csp->chsz-2) // last but one point
        {
          // tribs use dummy node for last point - ie junction ##
          grad=(csp->JunctionDEM - Arrptr->DEM[csp->ChanX[i]+csp->ChanY[i]*Parptr->xsz])
            /csp->Chandx[i];
        }
        else if(i==csp->chsz-1) // last point
        {
          // tribs use dummy node for last point - ie junction ##
          grad=( csp->JunctionDEM - Arrptr->DEM[csp->ChanX[i-1]+csp->ChanY[i-1]*Parptr->xsz])
            /csp->Chandx[i-1];
        }
        else // all other points as main channel
        {
          grad=(Arrptr->DEM[csp->ChanX[i+1]+csp->ChanY[i+1]*Parptr->xsz]-Arrptr->DEM[csp->ChanX[i]+csp->ChanY[i]*Parptr->xsz])
            /csp->Chandx[i];
        }
      }


      // check if slope is positive or negative
      if(grad>=0)
      {
        if (Statesptr->diffusive==ON)
        {
          // only keep uphill slopes for diffusive
          csp->Shalf[i]=-1*sqrt(fabs(grad));
        }
        else
        {
          // for kinematic we just pretend it is a downhill
          csp->Shalf[i]=sqrt(fabs(grad));
          // also warn user, in case they don't know
          if(*verbose==ON) printf("\nWARNING: Kinematic solver BUT uphill slope at point %i for channel %i.\n",i-1,chseg);
        }
      }
      else
      {
        csp->Shalf[i]=sqrt(fabs(grad));
      }
    }

    // Fill in Q boundary conditions from _tmp arrays
    for(i=0;i<npoints;i++) // loop through the input cross-sections
    {
      if(Q_Ident_tmp[i]!=0)  // if there is a boundary condition
      {
        for(j=0;j<csp->chsz;j++) // loop through the cross-sections mapped onto the dem space
        {
          if((cp[i]+0.1) >= csp->Chainage[j] && (cp[i]+0.1) < (csp->Chainage[j]+csp->Chandx[j]))
            // make sure we only apply the BC to one point
          {
            csp->Q_Ident[j]=Q_Ident_tmp[i]; // copy type across
            csp->Q_Val[j]=qp[i];            // copy value across
            if(Q_Ident_tmp[i]==3 || Q_Ident_tmp[i]==5 || Q_Ident_tmp[i]==8)	// check if BC type has name
            {
              strcpy(csp->Q_Name+j*80,Q_Name_tmp+i*80); // copy name across
            }
            if(Q_Ident_tmp[i]==7) // only for tribs
            {
              /*
                record in the trib data the location and segment it links to
                ChannelSegments[trib[i]].Next_Segment_Loc=j; // CCS
                ChannelSegments[trib[i]].Next_Segment=chseg; // CCS

                ^^ The above code is left commented out to show why the following QID7_Store vector is needed.
                As these terms need to be written to an instance of ChannelSegmentType not pointed to by csp
                (ChannelSegments[trib[i]]), we store them in the QID7 vector and use the UpdateChannelVector function
                to move the contents to the correct place after LoadRiver has finished. // CCS
              */

              QID7_Store QID7_tmp; // CCS create temp instance of QID7_Store and populate struct:
              QID7_tmp.chseg=chseg;
              QID7_tmp.Next_Segment_Loc=j;
              QID7_tmp.trib=trib[i];
              if(Statesptr->multiplerivers==ON)
              {
                QID7_tmp.RiverID=RiversIndexVecPtr->size(); /*#CCS# this allows us to keep track of which river we are in when later using QID7.
                                                              1st river ID will be 0, 2nd will be 1, etc.*/
              }
              QID7_Vec_Ptr->push_back(QID7_tmp); // CCS push_back temp instance into external vector
            }
          }
        }
      }
    }

    // release memory used for temporary variables
    delete[] xp;
    delete[] yp;
    delete[] wp;
    delete[] np;
    delete[] hp;
    delete[] xpi;
    delete[] ypi;
    delete[] cp;
    delete[] Q_Name_tmp;

    if (Statesptr->diffusive==1 && csp->Q_Ident[csp->chsz-1]==0)
    {
      if(*verbose==ON) printf("\nWARNING: Channel %i has no d/s BC using FREE with channel slope\n",chseg); // warn user that for diffusive no BC is set so using free
      csp->Q_Ident[csp->chsz-1]=1;  // copy type across
      csp->Q_Val[csp->chsz-1]=-1;   // copy value across
    }

    ChannelSegmentsVecPtr->push_back(tmp_chan);

  } // end of channel segment loop


  fclose(fp);
  if(*verbose==ON) printf("Done.\n\n");

  return;
}
//-----------------------------------------------------------------------------
/*UPDATE CHANNEL VECTOR FUNCTION // CCS
  This function is needed to update the ChannelSegments vector after LoadRiver has run.  It simply copies the  terms held
  in QID7_Store to the correct place in ChannelSegments.  It solves the problem of the LoadRiver function attempting to write
  to an index of ChannelSegemnts that hasn't yet been created when tagging tributary junctions.*/

void UpdateChannelsVector(States *Statesptr,ChannelSegmentType *CSTypePtr, vector<QID7_Store> *QID7_Vec_Ptr,QID7_Store *QID7Ptr, int *RiversIndexPtr)
{
  int vecsize, i, n;
  if(Statesptr->ChannelPresent==OFF) return;
  vecsize = QID7_Vec_Ptr->size();
  for(i = 0; i < vecsize; i++)
  {
    if(Statesptr->multiplerivers==OFF)
    {
      n = QID7Ptr[i].trib;
      CSTypePtr[n].Next_Segment=QID7Ptr[i].chseg;
      CSTypePtr[n].Next_Segment_Loc=QID7Ptr[i].Next_Segment_Loc;
    }
    else if(Statesptr->multiplerivers==ON)
    {
      if(QID7Ptr[i].RiverID == 0)
      {
        n = QID7Ptr[i].trib;
        CSTypePtr[n].Next_Segment=QID7Ptr[i].chseg;
        CSTypePtr[n].Next_Segment_Loc=QID7Ptr[i].Next_Segment_Loc;
      }
      else
      {
        n = QID7Ptr[i].trib + RiversIndexPtr[QID7Ptr[i].RiverID-1]; //Make sure we write to the correct CSTypePtr index when not in first river.
        CSTypePtr[n].Next_Segment=QID7Ptr[i].chseg+RiversIndexPtr[QID7Ptr[i].RiverID-1]; //Again, when not in first river we need to ensure correct ID of next segment.
        CSTypePtr[n].Next_Segment_Loc=QID7Ptr[i].Next_Segment_Loc;
      }
    }
  }

  return;

}
//-----------------------------------------------------------------------------
// LOAD FLOODPLAIN FRICTION FROM FILE
void LoadManningsn(Fnames *Fnameptr,Pars *Parptr,Arrays *Arrptr,int *verbose)
{
  FILE *fp;
  int i,j;
  char dum[800];
  double no_data_value=-9999;

  fp=fopen(Fnameptr->nfilename,"r");
  if(fp==NULL) return;

  if(*verbose==ON) printf("Loading floodplain Manning's n data:\t%s\t",Fnameptr->nfilename);

  for(i=0;i<5;i++) fscanf(fp,"%s %s",dum,dum);
  fscanf(fp,"%s %lf",dum,&no_data_value);

  Arrptr->Manningsn=new double[Parptr->xsz*Parptr->ysz];
  for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                             {
                               fscanf(fp,"%lf",Arrptr->Manningsn+i+j*Parptr->xsz);
                               if((int)Arrptr->Manningsn[i+j*Parptr->xsz]==no_data_value) Arrptr->Manningsn[i+j*Parptr->xsz]=Parptr->FPn;
                             }
  fclose(fp);

  if(*verbose==ON) printf("Done.\n\n");
  return;
}
void LoadSGCManningsn(Fnames *Fnameptr,Pars *Parptr,Arrays *Arrptr,int *verbose)
{
  FILE *fp;
  int i,j;
  char dum[800];
  double no_data_value=-9999;

  fp=fopen(Fnameptr->SGCnfilename,"r");
  if(fp==NULL) return;

  if(*verbose==ON) printf("Loading SGC Manning's n data:\t%s\t",Fnameptr->SGCnfilename);

  for(i=0;i<5;i++) fscanf(fp,"%s %s",dum,dum);
  fscanf(fp,"%s %lf",dum,&no_data_value);

  Arrptr->SGCManningsn=new double[Parptr->xsz*Parptr->ysz];
  for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                             {
                               fscanf(fp,"%lf",Arrptr->SGCManningsn+i+j*Parptr->xsz);
                               if((int)Arrptr->SGCManningsn[i+j*Parptr->xsz]==no_data_value) Arrptr->SGCManningsn[i+j*Parptr->xsz]=Parptr->SGC_n;
                             }
  fclose(fp);

  if(*verbose==ON) printf("Done.\n\n");
  return;
}
//-----------------------------------------------------------------------------
// LOAD POROSITY FROM FILE, TJF
void LoadPor(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,Arrays *Arrptr, int *verbose)
{
  FILE *fp;
  int i,j,k,m;
  double incr;
  char dum[800], buff[800];
  double no_data_value=-9999;

  fp=fopen(Fnameptr->porfilename,"r");
  if(fp==NULL)
  {
    if(*verbose==ON)
    {
      Statesptr->porosity=OFF;
      printf("Porosity off\n");
    }
    return;
  }
  if(*verbose==ON) printf("Loading porosity data:\t%s\n",Fnameptr->porfilename);

  Statesptr->porosity=ON;

  fscanf(fp,"%s",buff);

  // Determining the porosity method
  if((!strcmp(buff,"PFIX")||!strcmp(buff,"pfix")))
  {
    Parptr->Por_Ident=1;
    if(*verbose==ON) printf("Fixed Aerial Porosity Method\n");
  }

  if((!strcmp(buff,"PVAR")||!strcmp(buff,"pvar")))
  {
    Parptr->Por_Ident=2;
    if(*verbose==ON) printf("Water Height Dependent Aerial Porosity Method\n");
  }

  if((!strcmp(buff,"PBOUND")||!strcmp(buff,"pbound")))
  {
    Parptr->Por_Ident=3;
    if(*verbose==ON) printf("Fixed Boundary Porosity Method\n");
  }

  if((!strcmp(buff,"PBVAR")||!strcmp(buff,"pbvar")))
  {
    Parptr->Por_Ident=4;
    if(*verbose==ON) printf("Water Height Dependent Boundary Porosity Method\n");
  }

  // Loading the porosity values
  // Por_Ident=1
  if(Parptr->Por_Ident==1)
  {
    Arrptr->paerial=new double[Parptr->xsz*Parptr->ysz];
    for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                               {
                                 fscanf(fp,"%lf",Arrptr->paerial+i+j*Parptr->xsz);
                                 if((int)Arrptr->paerial[i+j*Parptr->xsz]==no_data_value) Arrptr->paerial[i+j*Parptr->xsz]=1;
                               }
    fclose(fp);
  }

  // Por_Ident=2
  else if(Parptr->Por_Ident==2)
  {
    fscanf(fp,"%s %i",dum,&Parptr->zsz);
    fscanf(fp,"%s %lf",dum,&Parptr->maxelev);
    fscanf(fp,"%s %lf",dum,&Parptr->zlev);

    Arrptr->paerial=new double[Parptr->xsz*Parptr->ysz*Parptr->zsz];

    for(k=0;k<Parptr->zsz;k++)
    {
      fscanf(fp,"%lf",&incr);
      for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                                 {
                                   fscanf(fp,"%lf",Arrptr->paerial+i+j*Parptr->xsz+k*Parptr->xsz*Parptr->ysz);
                                   if((int)Arrptr->paerial[i+j*Parptr->xsz+k*Parptr->xsz*Parptr->ysz]==no_data_value) Arrptr->paerial[i+j*Parptr->xsz+k*Parptr->xsz*Parptr->ysz]=1;
                                 }
    }
    fclose(fp);
  }


  // Por_Ident=3
  else if(Parptr->Por_Ident==3)
  {
    fscanf(fp,"%s",dum);

    Arrptr->paerial=new double[Parptr->xsz*Parptr->ysz];
    Arrptr->pbound=new double[Parptr->xsz*Parptr->ysz*4];

    for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                               {
                                 fscanf(fp,"%lf",Arrptr->paerial+i+j*Parptr->xsz);
                                 if((int)Arrptr->paerial[i+j*Parptr->xsz]==no_data_value) Arrptr->paerial[i+j*Parptr->xsz]=1;
                               }

    fscanf(fp,"%s",dum);

    for(k=0;k<4;k++)
    {
      for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                                 {
                                   fscanf(fp,"%lf",Arrptr->pbound+i+j*Parptr->xsz+k*Parptr->xsz*Parptr->ysz);
                                   if((int)Arrptr->pbound[i+j*Parptr->xsz+k*Parptr->xsz*Parptr->ysz]==no_data_value) Arrptr->pbound[i+j*Parptr->xsz+k*Parptr->xsz*Parptr->ysz]=1;
                                 }
    }
    fclose(fp);
  }

  // Por_Ident=4
  else if(Parptr->Por_Ident==4)
  {
    fscanf(fp,"%s %i",dum,&Parptr->zsz);
    fscanf(fp,"%s %lf",dum,&Parptr->maxelev);
    fscanf(fp,"%s %lf",dum,&Parptr->zlev);
    fscanf(fp,"%s",dum);

    Arrptr->paerial=new double[Parptr->xsz*Parptr->ysz*Parptr->zsz];
    Arrptr->pbound=new double[Parptr->xsz*Parptr->ysz*Parptr->zsz*4];

    for(k=0;k<Parptr->zsz;k++)
    {
      fscanf(fp,"%lf",&incr);
      for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                                 {
                                   fscanf(fp,"%lf",Arrptr->paerial+i+j*Parptr->xsz+k*Parptr->xsz*Parptr->ysz);
                                   if((int)Arrptr->paerial[i+j*Parptr->xsz+k*Parptr->xsz*Parptr->ysz]==no_data_value) Arrptr->paerial[i+j*Parptr->xsz+k*Parptr->xsz*Parptr->ysz]=1;
                                 }
    }

    fscanf(fp,"%s",dum);

    for(m=0;m<Parptr->zsz;m++)
    {
      fscanf(fp,"%lf",&incr);
      for(k=0;k<4;k++) for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                                                  {
                                                    fscanf(fp,"%lf",Arrptr->pbound+i+j*Parptr->xsz+m*Parptr->xsz*Parptr->ysz+k*Parptr->xsz*Parptr->ysz);
                                                    if((int)Arrptr->pbound[i+j*Parptr->xsz+m*Parptr->xsz*Parptr->ysz+k*Parptr->xsz*Parptr->ysz]==no_data_value)
                                                    {
                                                      Arrptr->pbound[i+j*Parptr->xsz+m*Parptr->xsz*Parptr->ysz+k*Parptr->xsz*Parptr->ysz]=1;
                                                    }
                                                  }
    }
    fclose(fp);
  }

  if(*verbose==ON) printf("Done.\n\n");

  return;
}
//-----------------------------------------------------------------------------
// LOAD INITIAL DEPTHS FROM FILE
void LoadStart(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,Arrays *Arrptr,SGCprams *SGCptr, int *verbose)
{
  FILE *fp;
  int i,j,gr;
  char dum[800];
  double no_data_value=-9999;

  fp=fopen(Fnameptr->startfilename,"r");
  if(fp==NULL)
  {
    if(*verbose==ON) printf("\nWARNING: Unable to load startfile:\t%s\t\n",Fnameptr->startfilename);
    return;
  }

  if(*verbose==ON) printf("Loading initial depths:\t%s\t",Fnameptr->startfilename);


  for(i=0;i<5;i++) fscanf(fp,"%s %s",dum,dum);
  fscanf(fp,"%s %lf",dum,&no_data_value);

  for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                             {
                               fscanf(fp,"%lf",Arrptr->H+i+j*Parptr->xsz);
                               // if no_data set depth to zero
                               if((int)Arrptr->H[i+j*Parptr->xsz]==no_data_value) Arrptr->H[i+j*Parptr->xsz]=0.0;
                               else if (Statesptr->startelev == ON) // convert water surface elevation to depth is this is being used
                               {
                                 // check to see if SGC is on
                                 if (Statesptr->SGC == ON)
                                 {
                                   gr = Arrptr->SGCgroup[i+j*Parptr->xsz]; // channel group number
                                   // is SGC is on calculate both a depth from the channel bed and the domain volume
                                   Arrptr->H[i+j*Parptr->xsz] = getmax(Arrptr->H[i+j*Parptr->xsz] - Arrptr->SGCz[i+j*Parptr->xsz], 0.0 );
                                   if (Arrptr->H[i+j*Parptr->xsz] <= Arrptr->SGCbfH[i+j*Parptr->xsz] || Arrptr->SGCwidth[i+j*Parptr->xsz] > Parptr->dx)	Arrptr->SGCVol[i+j*Parptr->xsz] = CalcSGC_UpV(SGCptr->SGCchantype[gr], Arrptr->H[i+j*Parptr->xsz], SGCptr->SGCs[gr], Arrptr->SGCc[i+j*Parptr->xsz]);
                                   else																													Arrptr->SGCVol[i+j*Parptr->xsz] = Arrptr->SGCbfV[i+j*Parptr->xsz] + (Arrptr->H[i+j*Parptr->xsz]-Arrptr->SGCbfH[i+j*Parptr->xsz])*Parptr->dA; // out of bank level
                                 }
                                 else Arrptr->H[i+j*Parptr->xsz] = getmax(Arrptr->H[i+j*Parptr->xsz] - Arrptr->DEM[i+j*Parptr->xsz], 0.0 );

                               }
                               else if (Statesptr->SGC == ON)
                               {
                                 gr = Arrptr->SGCgroup[i+j*Parptr->xsz]; // channel group number
                                 if (Arrptr->H[i+j*Parptr->xsz] <= Arrptr->SGCbfH[i+j*Parptr->xsz] || Arrptr->SGCwidth[i+j*Parptr->xsz] > Parptr->dx)	Arrptr->SGCVol[i+j*Parptr->xsz] = CalcSGC_UpV(SGCptr->SGCchantype[gr], Arrptr->H[i+j*Parptr->xsz], SGCptr->SGCs[gr], Arrptr->SGCc[i+j*Parptr->xsz]);
                                 else  Arrptr->SGCVol[i+j*Parptr->xsz] = Arrptr->SGCbfV[i+j*Parptr->xsz] + (Arrptr->H[i+j*Parptr->xsz]-Arrptr->SGCbfH[i+j*Parptr->xsz])*Parptr->dA; // out of bank level
                               }
                             }
  fclose(fp);

  if(*verbose==ON) printf("Done.\n\n");

  return;
}
//----------------------------------------------------------------------------
// LOAD DEM FROM FILE
// Also loads cell size, lower left corner coordinates, and new all rasters
void LoadDEM(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,Arrays *Arrptr,int *verbose)
{
  FILE *fp;
  char dum[800];
  double no_data_value=-9999;
  int i,j;

  fp=fopen(Fnameptr->demfilename,"rb");
  if(fp==NULL) {fprintf(stderr,"No DEM file found. Aborting.\n");exit(0);}

  if(*verbose==ON) printf("\nLoading DEM:\t%s\n",Fnameptr->demfilename);

  fscanf(fp,"%s %i",dum,&Parptr->xsz);
  fscanf(fp,"%s %i",dum,&Parptr->ysz);
  fscanf(fp,"%s %lf",dum, &Parptr->blx);
  fscanf(fp,"%s %lf",dum, &Parptr->bly);
  fscanf(fp,"%s %lf",dum, &Parptr->dx);

  Parptr->dx_sqrt = sqrt((double) Parptr->dx); // sqrt now for later use in flooplain calcs - small speed increase
  Parptr->dy=Parptr->dx;Parptr->dA=Parptr->dx*Parptr->dy;
  Parptr->tlx=Parptr->blx;Parptr->tly=Parptr->bly+Parptr->ysz*Parptr->dy;
  fscanf(fp,"%s %lf",dum,&no_data_value);

  if(*verbose==ON)
  {
    printf("%ix%i\nBL corner\t(%lf,%lf)\nNODATA_value\t%lf\n",
           Parptr->xsz,Parptr->ysz,Parptr->blx,Parptr->bly,no_data_value);
  }

  // allocate memory for arrays, Note the () at the end ensures all elements are initialised to zero
  Arrptr->H=new double[Parptr->xsz*Parptr->ysz]();
  // If Roe==on allocate memory (JN/IV)
  if(Statesptr->Roe==ON)
  {
    Arrptr->HU=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->HV=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->LSHU=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->RSHU=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->BSHU=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->TSHU=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->LSHV=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->RSHV=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->BSHV=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->TSHV=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->FHx=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
    Arrptr->FHUx=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
    Arrptr->FHVx=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
    Arrptr->FHy=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
    Arrptr->FHUy=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
    Arrptr->FHVy=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
  }
  Arrptr->maxH=new double[Parptr->xsz*Parptr->ysz]();
  Arrptr->totalHtm=new double[Parptr->xsz*Parptr->ysz]();
  Arrptr->Qx=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
  Arrptr->Qy=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
  // added to record Hflow (needed for velocity and hazard calculations
  //Arrptr->Hflowx=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
  //Arrptr->Hflowy=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();

  if(Statesptr->voutput==ON)
  {
    Arrptr->Vx=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
    Arrptr->Vy=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
    Arrptr->maxVx=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
    Arrptr->maxVy=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
  }
  if (Statesptr->hazard==ON)
  {
    Arrptr->maxVc=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->maxVcH=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->maxHaz=new double[Parptr->xsz*Parptr->ysz]();
  }
  if (Statesptr->SGC==ON)
  {
    // Geometric variables
    Arrptr->SGCwidth=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->SGCz=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->SGCc=new double[Parptr->xsz*Parptr->ysz]();
    // Flow variables
    Arrptr->QxSGold=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
    Arrptr->QySGold=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
    // Bank full variables
    Arrptr->SGCbfH=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->SGCbfV=new double[Parptr->xsz*Parptr->ysz]();
    // Volume variables
    Arrptr->SGCVol=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->SGCdVol=new double[Parptr->xsz*Parptr->ysz]();
    Arrptr->SGCQin=new double[Parptr->xsz*Parptr->ysz](); //JMH
    // Model parameters
    Arrptr->SGCgroup=new int[Parptr->xsz*Parptr->ysz]();
    if(Statesptr->save_Qs==ON)
    {
      Arrptr->SGCFlowWidth=new double[Parptr->xsz*Parptr->ysz]();
    }
  }

  Arrptr->Qxold=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
  Arrptr->Qyold=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();

  // allocate memory for velocity arrays U and V
  // currently only used in acceleration version - initialised under all conditions as may want them
  // for other lisflood versions (TJF)
  Arrptr->U=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();
  Arrptr->V=new double[(Parptr->xsz+1)*(Parptr->ysz+1)]();

  // allocate memory for none zero arrays
  Arrptr->maxHtm=new double[Parptr->xsz*Parptr->ysz];
  Arrptr->initHtm=new double[Parptr->xsz*Parptr->ysz];
  Arrptr->TRecx=new double[(Parptr->xsz+1)*(Parptr->ysz+1)];
  Arrptr->TRecy=new double[(Parptr->xsz+1)*(Parptr->ysz+1)];
  Arrptr->LimQx=new double[(Parptr->xsz+1)*(Parptr->ysz+1)];
  Arrptr->LimQy=new double[(Parptr->xsz+1)*(Parptr->ysz+1)];

  Arrptr->ChanMask=new int[Parptr->xsz*Parptr->ysz];
  Arrptr->SegMask=new int[Parptr->xsz*Parptr->ysz];

  Arrptr->DEM=new double[Parptr->xsz*Parptr->ysz];

  // allocate memory for flow direction array for routing very shallow flows from rainfall component CCS 13/03/2012
  if (Statesptr->routing==ON)
  {
    Arrptr->FlowDir=new int[Parptr->xsz*Parptr->ysz]();
    for(i=0;i<Parptr->xsz*Parptr->ysz;i++) Arrptr->FlowDir[i]=(int)NULLVAL;
  }

  // allocate memory for lat long arrays
  Arrptr->dx=new double[Parptr->xsz*Parptr->ysz]();
  Arrptr->dy=new double[Parptr->xsz*Parptr->ysz]();
  Arrptr->dA=new double[Parptr->xsz*Parptr->ysz]();



  // set initial values of elements for some arrays to NULLVAL
  for(i=0;i<Parptr->xsz*Parptr->ysz;i++) Arrptr->maxHtm[i]=NULLVAL;
  for(i=0;i<Parptr->xsz*Parptr->ysz;i++) Arrptr->initHtm[i]=NULLVAL;
  for(i=0;i<(Parptr->xsz+1)*(Parptr->ysz+1);i++) Arrptr->TRecx[i]=NULLVAL;
  for(i=0;i<(Parptr->xsz+1)*(Parptr->ysz+1);i++) Arrptr->TRecy[i]=NULLVAL;
  for(i=0;i<(Parptr->xsz+1)*(Parptr->ysz+1);i++) Arrptr->LimQx[i]=NULLVAL;
  for(i=0;i<(Parptr->xsz+1)*(Parptr->ysz+1);i++) Arrptr->LimQy[i]=NULLVAL;

  // set initial values of elements for mask arrays to NULLVAL
  for(i=0;i<Parptr->xsz*Parptr->ysz;i++) Arrptr->ChanMask[i]=-1;
  for(i=0;i<Parptr->xsz*Parptr->ysz;i++) Arrptr->SegMask[i]=-1;

  for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                             {
                               fscanf(fp,"%lf",Arrptr->DEM+i+j*Parptr->xsz);
                               if((int)Arrptr->DEM[i+j*Parptr->xsz]==no_data_value) Arrptr->DEM[i+j*Parptr->xsz]=1e10;
                             }
  fclose(fp);

  if(*verbose==ON) printf("Done.\n\n");

  return;
}
//-----------------------------------------------------------------------------------
// LOADS FILE GIVING IDENTIFIERS FOR EACH BOUNDARY CELL FROM .bci FILE
// (e.g. HFIX, QVAR etc)
void LoadBCs(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,BoundCs *BCptr,Arrays *Arrptr,int *verbose)
{
  int numBCs,i,j,BCi1,BCi2,tmpi;
  double start,finish;
  FILE *fp;
  char buff[800], buff2[800], buff3[800], side;
  double BC_tmp;
  int pi=-1,maxpi=MAXPI;  // increase max from 10 to 20K (MT)
  double px,py;

  int *new_xpi,*new_ypi,*new_PS_Ident;
  double *new_PS_Val, *new_PS_qold, *new_PS_qSGold;
  char *new_PS_Name;

  // POINT SOURCE STUFF
  BCptr->xpi=new int[maxpi];
  BCptr->ypi=new int[maxpi];
  for(i=0;i<maxpi;i++) BCptr->xpi[i]=BCptr->ypi[i]=-1;
  BCptr->PS_Ident=new int[maxpi];
  BCptr->PS_Val=new double[maxpi];
  BCptr->PS_Name=new char[maxpi*80];
  for(i=0;i<maxpi;i++) {BCptr->PS_Ident[i]=0;BCptr->PS_Val[i]=-1.0;BCptr->PS_Name[i]='\0';}
  BCptr->numPS=-1;

  // BOUNDARY CONDITION STUFF
  numBCs=2*Parptr->xsz+2*Parptr->ysz;
  BCptr->BC_Ident=new int[numBCs];
  BCptr->BC_Val=new double[numBCs];
  BCptr->BC_Name=new char[numBCs*80];
  for(i=0;i<numBCs;i++) {BCptr->BC_Ident[i]=0;BCptr->BC_Val[i]=-1.0;BCptr->BC_Name[i]='\0';}

  fp=fopen(Fnameptr->bcifilename,"rb");
  if(fp==NULL) return;
  if(*verbose==ON) printf("Loading boundary condition IDs:\t%s\n",Fnameptr->bcifilename);

  while(!feof(fp))
  {
    BCi1=BCi2=-1;side='\0';
    // Read NSEW and location, and determine start/finish of BC-->(BCi1,BCi2)
    fscanf(fp,"%s",buff);
    if(feof(fp)) break;
    if(buff[0]=='N')
    {
      fscanf(fp,"%lf%lf",&start,&finish);

      if(start<Parptr->blx) start=Parptr->blx;
      if(start>Parptr->blx+Parptr->xsz*Parptr->dx) start=Parptr->blx+Parptr->xsz*Parptr->dx;
      if(finish<Parptr->blx) finish=Parptr->blx;
      if(finish>Parptr->blx+Parptr->xsz*Parptr->dx) finish=Parptr->blx+Parptr->xsz*Parptr->dx;

      BCi1=(int)((start-Parptr->blx)/Parptr->dy);
      BCi2=(int)((finish-Parptr->blx)/Parptr->dy);
      if(BCi1>BCi2) {tmpi=BCi1;BCi1=BCi2;BCi2=tmpi;}
      BCi2--;side='N';
    }

    if(buff[0]=='W')
    {
      fscanf(fp,"%lf%lf",&start,&finish);

      if(start<Parptr->bly) start=Parptr->bly;
      if(start>Parptr->tly) start=Parptr->tly;
      if(finish<Parptr->bly) finish=Parptr->bly;
      if(finish>Parptr->tly) finish=Parptr->tly;

      BCi1=(int)(2*Parptr->xsz+2*Parptr->ysz-(Parptr->tly-start)/Parptr->dy);
      BCi2=(int)(2*Parptr->xsz+2*Parptr->ysz-(Parptr->tly-finish)/Parptr->dy);
      if(BCi1>BCi2) {tmpi=BCi1;BCi1=BCi2;BCi2=tmpi;}
      BCi2--;side='W';
    }

    if(buff[0]=='S')
    {
      fscanf(fp,"%lf%lf",&start,&finish);

      if(start<Parptr->blx) start=Parptr->blx;
      if(start>Parptr->blx+Parptr->xsz*Parptr->dx) start=Parptr->blx+Parptr->xsz*Parptr->dx;
      if(finish<Parptr->blx) finish=Parptr->blx;
      if(finish>Parptr->blx+Parptr->xsz*Parptr->dx) finish=Parptr->blx+Parptr->xsz*Parptr->dx;

      BCi1=(int)(2*Parptr->xsz+Parptr->ysz-(start-Parptr->blx)/Parptr->dy);
      BCi2=(int)(2*Parptr->xsz+Parptr->ysz-(finish-Parptr->blx)/Parptr->dy);

      if(BCi1>BCi2) {tmpi=BCi1;BCi1=BCi2;BCi2=tmpi;}
      BCi2--;side='S';
    }

    if(buff[0]=='E')
    {
      fscanf(fp,"%lf%lf",&start,&finish);

      if(start<Parptr->bly) start=Parptr->bly;
      if(start>Parptr->tly) start=Parptr->tly;
      if(finish<Parptr->bly) finish=Parptr->bly;
      if(finish>Parptr->tly) finish=Parptr->tly;

      BCi1=(int)(Parptr->xsz+(Parptr->tly-start)/Parptr->dy);
      BCi2=(int)(Parptr->xsz+(Parptr->tly-finish)/Parptr->dy);
      if(BCi1>BCi2) {tmpi=BCi1;BCi1=BCi2;BCi2=tmpi;}
      BCi2--;side='E';
    }

    // Read locations of point sources or point free
    if(buff[0]=='P' || buff[0]=='F')
    {
      pi++;
      fscanf(fp,"%lf%lf",&px,&py);
      BCptr->xpi[pi]=(int)((px-Parptr->blx)/Parptr->dx);
      BCptr->ypi[pi]=(int)((Parptr->tly-py)/Parptr->dy);
    }

    // Read free boundary condition locations
    // load buffer until EOL
    j=0;
    do{buff2[j]=fgetc(fp);} while(buff2[j++]!='\n' && !feof(fp));
    buff2[j-1]='\0';               // Finish off string
    // get buff so you know boundary type
    BC_tmp = -1;
    sscanf(buff2,"%s%lf",buff,&BC_tmp);

    // If a FREE surface boundary condition
    if((!strcmp(buff,"FREE")||!strcmp(buff,"free"))&&BCi1>-1)
    {
      // if BC_tmp is -1 there is no slope specifed... use local slope from elevation model (origional mehod)
      if(BC_tmp < -0.999) // ie -1 (done like this as double)
      {
        if(*verbose==ON) printf("FREE on %c side start %lf end %lf\n",side,start,finish);
      }
      else
      {
        if(*verbose==ON) printf("FREE on %c side start %lf end %lf using slope %.5f\n",side,start,finish,BC_tmp);
      }
      for(i=BCi1;i<=BCi2;i++)
      {
        BCptr->BC_Ident[i]=1;
        // store floodplain slope in BCptr->BC_Val[i]
        if(Statesptr->adaptive_ts==ON || Statesptr->qlim==ON)
        {
          if (BC_tmp < -0.999) BCptr->BC_Val[i] = BC_tmp; // make BC_Val equal to -1 to use local water surface slope (origional lisflood)
          else BCptr->BC_Val[i] = sqrt(BC_tmp); // sqrt of user specified slope for diffusive version (jcn)
        }
        else
        {
          BCptr->BC_Val[i] = BC_tmp; // user specified slope or -1(for local slope) for accelleration or any other version (jcn)
        }
      }
    }
    // Read in fixed values of H for boundary
    else if((!strcmp(buff,"HFIX")||!strcmp(buff,"hfix")) &&BCi1>-1)
    {
      sscanf(buff2,"%s%lf",buff,&BC_tmp);
      if(*verbose==ON) printf("HFIX at %lf on %c side start %lf end %lf\n",BC_tmp,side,start,finish);

      for(i=BCi1;i<=BCi2;i++)
      {
        BCptr->BC_Val[i]=BC_tmp;
        BCptr->BC_Ident[i]=2;
      }
    }
    // Read in fixed values if Q for boundary
    else if((!strcmp(buff,"QFIX")||!strcmp(buff,"qfix")) &&BCi1>-1)
    {
      sscanf(buff2,"%s%lf",buff,&BC_tmp);
      if(*verbose==ON) printf("QFIX at %lf on %c side start %lf end %lf\n",BC_tmp,side,start,finish);

      for(i=BCi1;i<=BCi2;i++)
      {
        BCptr->BC_Val[i]=BC_tmp;
        BCptr->BC_Ident[i]=4;
      }
    }

    //	Read boundary names for varying values
    else if((!strcmp(buff,"QVAR")||!strcmp(buff,"qvar"))&&BCi1>-1)
    {
      //fscanf(fp,"%s",buff);
      sscanf(buff2,"%s%s",buff3,buff);
      if(*verbose==ON)
        printf("QVAR from bdy file %s on %c side start %lf end %lf\n",
               buff,side,start,finish);
      for(i=BCi1;i<=BCi2;i++)
      {
        strcpy(BCptr->BC_Name+i*80,buff);
        BCptr->BC_Ident[i]=5;
      }
    }
    else if((!strcmp(buff,"HVAR")||!strcmp(buff,"hvar"))&&BCi1>-1)
    {
      //fscanf(fp,"%s",buff);
      sscanf(buff2,"%s%s",buff3,buff);
      if(*verbose==ON)
        printf("HVAR from bdy file %s on %c side start %lf end %lf\n",buff,side,start,finish);
      for(i=BCi1;i<=BCi2;i++)
      {
        strcpy(BCptr->BC_Name+i*80,buff);
        BCptr->BC_Ident[i]=3;
      }
    }
    // Fixed/Varying values/names for point sources
    // Note these need to come after the boundary conditions in the code else both will get implemented!
    // I'm not convinced &&BCptr->xpi[pi]>-1&&pi>-1 is strict enough (JCN)
    else if((!strcmp(buff,"HFIX")||!strcmp(buff,"hfix"))&&BCptr->xpi[pi]>-1&&pi>-1)
    {
      //fscanf(fp,"%lf",&BC_tmp);
      sscanf(buff2,"%s%s",buff3,buff);
      if(*verbose==ON) printf("HFIX at point [%lf,%lf] %lf\n",px,py,BC_tmp);
      BCptr->PS_Val[pi]=BC_tmp;
      BCptr->PS_Ident[pi]=2;
    }
    else if((!strcmp(buff,"QFIX")||!strcmp(buff,"qfix"))&&BCptr->xpi[pi]>-1&&pi>-1)
    {
      //fscanf(fp,"%lf",&BC_tmp);
      if(*verbose==ON) printf("QFIX at point [%lf,%lf] %lf\n",px,py,BC_tmp);
      BCptr->PS_Val[pi]=BC_tmp;
      BCptr->PS_Ident[pi]=4;
    }
    else if((!strcmp(buff,"QVAR")||!strcmp(buff,"qvar"))&&BCptr->xpi[pi]>-1&&pi>-1)
    {
      //fscanf(fp,"%s",buff);
      sscanf(buff2,"%s%s",buff3,buff);
      if(*verbose==ON)
        printf("QVAR at point [%lf,%lf] %s\n",px,py,buff);
      strcpy(BCptr->PS_Name+pi*80,buff);
      BCptr->PS_Ident[pi]=5;
    }
    else if((!strcmp(buff,"HVAR")||!strcmp(buff,"hvar"))&&BCptr->xpi[pi]>-1&&pi>-1)
    {
      //fscanf(fp,"%s",buff);
      sscanf(buff2,"%s%s",buff3,buff);
      if(*verbose==ON) printf("HVAR at point [%lf,%lf] %s\n",px,py,buff);
      strcpy(BCptr->PS_Name+pi*80,buff);
      BCptr->PS_Ident[pi]=3;
    }
    else if((!strcmp(buff,"FREE")||!strcmp(buff,"free"))&&BCptr->xpi[pi]>-1&&pi>-1&&Statesptr->SGC==ON)
    {
      // point FREE boundary for internal SGC boundaries
      if(*verbose==ON) printf("FREE at point [%lf,%lf] with slope %lf\n",px,py,BC_tmp);
      BCptr->PS_Val[pi]=BC_tmp;
      BCptr->PS_Ident[pi]=6;
    }
    else
    {
      if(*verbose==ON) printf("WARNING: Incorrect boundary condition in .bci file\n");
    }

  }

  if(pi>-1)
  {
    pi++;
    new_xpi=new int[pi];
    new_ypi=new int[pi];
    new_PS_Ident=new int[pi];
    new_PS_Val=new double[pi];
    new_PS_qold=new double[pi]();
    new_PS_qSGold=new double[pi]();
    new_PS_Name=new char[pi*80];

    for(i=0;i<pi;i++)
    {
      new_xpi[i]=BCptr->xpi[i];
      new_ypi[i]=BCptr->ypi[i];
      new_PS_Ident[i]=BCptr->PS_Ident[i];
      new_PS_Val[i]=BCptr->PS_Val[i];
      for(j=0;j<80;j++) new_PS_Name[i*80+j]=BCptr->PS_Name[i*80+j];
    }
    BCptr->xpi=new_xpi;
    BCptr->ypi=new_ypi;
    BCptr->PS_Ident=new_PS_Ident;
    BCptr->PS_Val=new_PS_Val;
    BCptr->PS_qold=new_PS_qold;
    BCptr->PS_qSGold=new_PS_qSGold;
    BCptr->PS_Name=new_PS_Name;

    BCptr->numPS=pi;
  }

  if(*verbose==ON) printf("Done.\n\n");
  //  LoadBCVar();

  fclose(fp);

  return;
}
//-----------------------------------------------------------------------------
// LOAD TIME VARYING BOUNDARY CONDITIONS FROM .bdy FILE
void LoadBCVar(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,BoundCs *BCptr,ChannelSegmentType *ChannelSegments,Arrays *Arrptr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr, int *verbose)
{
  FILE *fp;
  int i,j,nbdy=0,ndata,numBCs,chseg;
  char buff[255],units[80];
  double **peterpointer;

  numBCs=2*Parptr->xsz+2*Parptr->ysz;
  BCptr->BCVarlist=new double*[numBCs];
  for(i=0;i<numBCs;i++) BCptr->BCVarlist[i]=NULL;

  if(Statesptr->ChannelPresent==ON) {
    for(chseg=0;chseg<(int) ChannelSegmentsVecPtr->size();chseg++) // CCS
    {
      ChannelSegments[chseg].QVarlist=new double*[ChannelSegments[chseg].chsz];
      for(i=0;i<ChannelSegments[chseg].chsz;i++) ChannelSegments[chseg].QVarlist[i]=NULL;
    }
  }

  fp=fopen(Fnameptr->bdyfilename,"r");
  if(fp==NULL) return;
  if(*verbose==ON) printf("Loading time varying boundary conditions:\t%s\n",Fnameptr->bdyfilename);

  while(!feof(fp))
  {
    j=0;  // skip 1st comment line

    if(nbdy==0){
      do{
        buff[j]=fgetc(fp);
        if(feof(fp)) break;
      } while(buff[j++]!='\n');

    }

    fscanf(fp,"%s",buff);
    if(feof(fp) || buff[0]=='\n') break;

    // Check through list (2d domain) of boundary names, if a match found, assign
    // pixel to this boundary and set peterpointer to point to
    // relevant structure. If none found set to NULL, boundary condition data is
    // unassigned but file is still read through to get to next set of data.
    peterpointer=NULL;
    for(i=0;i<numBCs;i++)
    {
      if(!strcmp(buff,(BCptr->BC_Name+i*80)))
      {
        BCptr->BC_Val[i]=nbdy;
        peterpointer=BCptr->BCVarlist;
      }
    }

    // Check through list (river channel) of boundary names, if a match found, assign
    // channel node to this boundary and set peterpointer to point to
    // relevant structure.
    if(Statesptr->ChannelPresent==ON)
    {
      for(chseg=0;chseg<(int)ChannelSegmentsVecPtr->size();chseg++) // CCS
      {
        for(i=0;i<ChannelSegments[chseg].chsz;i++)
        {
          if(!strcmp(buff,(ChannelSegments[chseg].Q_Name+i*80)))
          {
            ChannelSegments[chseg].Q_Val[i]=nbdy;
            peterpointer=ChannelSegments[chseg].QVarlist;
          }
        }
      }
    }

    // same for point sources
    for(i=0;i<BCptr->numPS;i++)
    {
      if(!strcmp(buff,(BCptr->PS_Name+i*80)))
      {
        BCptr->PS_Val[i]=nbdy;
        peterpointer=BCptr->BCVarlist;
      }
    }

    // Set up arrays, 1 element per time step for BCVarlist, 1 per data
    // point for Vargiven and Tgiven
    fscanf(fp,"%i%s",&ndata,units);
    if(peterpointer!=NULL)
    {
      peterpointer[nbdy]=new double[ndata*2+2];
      peterpointer[nbdy][ndata*2+1]=-1;

      // Go through the motions even if peterpointer==NULL
      for(i=0;i<ndata;i++) fscanf(fp,"%lf%lf",peterpointer[nbdy]+i*2,peterpointer[nbdy]+i*2+1);
      if(!strcmp(units,"hours")) for(i=0;i<ndata;i++) *(peterpointer[nbdy]+i*2+1)*=3600;
      else if(!strcmp(units,"days")) for(i=0;i<ndata;i++) *(peterpointer[nbdy]+i*2+1)*=(3600*24);

      fgetc(fp);		// Get end of line character
    }

    //    for(i=0;i<ndata*2+2;i++){
    //      printf("\nLoadBCVar: i=%d, peterpointer[0]+i*2=%lf, peterpointer[nbdy]+i*2+1=%lf",i,*(peterpointer[0]+i*2),*(peterpointer[0]+i*2+1));
    //    }

    if(peterpointer==NULL)
    {
      printf("WARNING: bdy %s is unreferenced - data ignored.\n",buff);
      continue;
    }

    nbdy++;

    if(*verbose==ON) printf("bdy %s read.\n",buff);
  }

  // Check for bdy names not found in bdy file
  if(Statesptr->ChannelPresent==ON) {
    for(chseg=0;chseg<(int)ChannelSegmentsVecPtr->size();chseg++) for(i=0;i<ChannelSegments[chseg].chsz;i++) if(ChannelSegments[chseg].Q_Ident[i]==5 && ChannelSegments[chseg].Q_Val[i]<0) // CCS
                                                                                                             {
                                                                                                               printf("WARNING: bdy %s in river file not found in bdy file - ignored.\n",ChannelSegments[chseg].Q_Name+i*80);
                                                                                                               ChannelSegments[chseg].Q_Ident[i]=0;
                                                                                                             }
  }
  for(i=0;i<numBCs;i++) if((BCptr->BC_Ident[i]==5 || BCptr->BC_Ident[i]==3) && BCptr->BC_Val[i]<0)
                        {
                          printf("WARNING: bdy %s in bci file not found in bdy file - ignored.\n",BCptr->BC_Name+i*80);
                          for(j=0;j<numBCs;j++) if(!strcmp(BCptr->BC_Name+i*80,BCptr->BC_Name+j*80)) BCptr->BC_Ident[j]=0;
                        }

  for(i=0;i<BCptr->numPS;i++) if((BCptr->PS_Ident[i]==5 || BCptr->PS_Ident[i]==3) && BCptr->PS_Val[i]<0)
                              {
                                printf("WARNING: bdy %s in bci file not found in bdy file - ignored.\n",BCptr->BC_Name+i*80);
                                for(j=0;j<BCptr->numPS;j++) if(!strcmp(BCptr->PS_Name+i*80,BCptr->PS_Name+j*80)) BCptr->PS_Ident[j]=0;
                              }

  fclose(fp);

  if(*verbose==ON) printf("Done.\n\n");
  return;
}

//-----------------------------------------------------------------------------
// LOAD TIME EVAPORATION FROM .evap FILE
void LoadEvap(Fnames *Fnameptr,Arrays *Arrptr, int *verbose)
{
  FILE *fp;
  int i,j,ndata;
  char buff[255],units[80];

  fp=fopen(Fnameptr->evapfilename,"r");
  if(fp==NULL) return;
  if(*verbose==ON) printf("\nLoading time varying evaporation:\t%s\n",Fnameptr->evapfilename);

  //while(!feof(fp)) // removed this because it used to write over the evaportaion data if you had a carage return at the end of the file (JCN)
  // this will need to be change if you ever what to import multiple evaporation time series (can't think why you would do spatially varying this way though)
  //{
  j=0; // skip 1st comment line
  do{
    buff[j]=fgetc(fp);
    if(feof(fp)) break;
  } while(buff[j++]!='\n');
  //if(feof(fp) || buff[0]=='\n') break; // removed because of removal of while loop see comment above

  fscanf(fp,"%i%s",&ndata,units);

  Arrptr->evap=new double[ndata*2+2]();
  Arrptr->evap[ndata*2+1]=-1;

  for(i=0;i<ndata;i++) fscanf(fp,"%lf%lf",Arrptr->evap+i*2,Arrptr->evap+i*2+1);

  // convert time to seconds
  if(!strcmp(units,"hours")) for(i=0;i<ndata;i++) *(Arrptr->evap+i*2+1)*=3600;
  else if(!strcmp(units,"days")) for(i=0;i<ndata;i++) *(Arrptr->evap+i*2+1)*=(3600*24);

  // convert evaporation rate from mm/day to m/second
  for(i=0;i<ndata;i++) *(Arrptr->evap+i*2)/=(1000*24*3600);
  //}
  fclose(fp);

  if(*verbose==ON) printf("Done.\n\n");
  return;
}
//-----------------------------------------------------------------------------
// LOAD TIME VARYING RAINFALL FROM .rain FILE
void LoadRain(Fnames *Fnameptr,Arrays *Arrptr, int *verbose)
{
  FILE *fp;
  int i,j,ndata;
  char buff[255],units[80];

  fp=fopen(Fnameptr->rainfilename,"r");
  if(fp==NULL) return;
  if(*verbose==ON) printf("\nLoading time varying rainfall:\t%s\n",Fnameptr->rainfilename);

  //while(!feof(fp))// removed this because it used to write over the evaportaion data if you had a carage return at the end of the file (JCN)
  // this will need to be change if you ever what to import multiple evaporation time series (can't think why you would do spatially varying this way though)
  //{
  j=0; // skip 1st comment line
  do{
    buff[j]=fgetc(fp);
    if(feof(fp)) break;
  } while(buff[j++]!='\n');

  //if(feof(fp) || buff[0]=='\n') break;

  fscanf(fp,"%i%s",&ndata,units);

  Arrptr->rain=new double[ndata*2+2];
  Arrptr->rain[ndata*2+1]=-1;

  for(i=0;i<ndata;i++) fscanf(fp,"%lf%lf",Arrptr->rain+i*2,Arrptr->rain+i*2+1);

  // convert time to seconds
  if(!strcmp(units,"hours")) for(i=0;i<ndata;i++) *(Arrptr->rain+i*2+1)*=3600;
  else if(!strcmp(units,"days")) for(i=0;i<ndata;i++) *(Arrptr->rain+i*2+1)*=(3600*24);

  // convert rainfall rate from mm/hr to m/second
  for(i=0;i<ndata;i++) *(Arrptr->rain+i*2)/=(1000*3600);

  //}
  fclose(fp);

  if(*verbose==ON) printf("Done.\n\n");
  return;
}
//-----------------------------------------------------------------------------
// LOAD INITIAL DEPTHS FROM FILE
void LoadSGC(Fnames *Fnameptr,Pars *Parptr,Arrays *Arrptr,States *Statesptr, SGCprams *SGCptr, int *verbose)
{
  FILE *fp;
  int i,j;
  char dum[80];
  double no_data_value=-9999, tmp;

  fp=fopen(Fnameptr->SGCbankfilename,"r");
  if(fp==NULL) return;
  if(*verbose==ON) printf("Loading SGCbank:\t%s\t",Fnameptr->SGCbankfilename);
  for(i=0;i<5;i++) fscanf(fp,"%s %s",dum,dum);
  fscanf(fp,"%s %lf",dum,&no_data_value);
  for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                             {
                               fscanf(fp,"%lf",Arrptr->SGCz+i+j*Parptr->xsz);
                               if (Arrptr->SGCz[i+j*Parptr->xsz] == no_data_value) Arrptr->SGCz[i+j*Parptr->xsz] = Arrptr->DEM[i+j*Parptr->xsz]; // In the case of no_data_value set the bank height to the DEM
                             }
  fclose(fp);
  if(*verbose==ON) printf("Done.\n");

  fp=fopen(Fnameptr->SGCwidthfilename,"r");
  if(fp==NULL) return;
  if(*verbose==ON) printf("Loading SGCwidth:\t%s\t",Fnameptr->SGCwidthfilename);
  for(i=0;i<5;i++) fscanf(fp,"%s %s",dum,dum);
  fscanf(fp,"%s %lf",dum,&no_data_value);
  for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                             {
                               fscanf(fp,"%lf",Arrptr->SGCwidth+i+j*Parptr->xsz);
                               if (Arrptr->SGCwidth[i+j*Parptr->xsz] == no_data_value) Arrptr->SGCwidth[i+j*Parptr->xsz] = 0.0; // In the case of no_data_value set the width to zero
                             }
  fclose(fp);
  if(*verbose==ON) printf("Done.\n\n");

  // This loads the distributed SGC group information
  if (Statesptr->SGCchangroup == ON)
  {
    fp=fopen(Fnameptr->SGCchangroupfilename,"r");
    if(fp==NULL) return;
    if(*verbose==ON) printf("Loading SGC channel group:\t%s\t",Fnameptr->SGCchangroupfilename);
    for(i=0;i<5;i++) fscanf(fp,"%s %s",dum,dum);
    fscanf(fp,"%s %lf",dum,&no_data_value);
    for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
                               {
                                 fscanf(fp,"%lf",&tmp); // read a double in case someone doesn't put intergers in the ascii file.
                                 Arrptr->SGCgroup[i+j*Parptr->xsz] = (int)tmp;
                                 if (Arrptr->SGCgroup[i+j*Parptr->xsz] == (int)no_data_value || Arrptr->SGCgroup[i+j*Parptr->xsz] < 0) Arrptr->SGCgroup[i+j*Parptr->xsz] = 0; // In the case of no_data_value set the chan group to zero
                               }
    fclose(fp);
    if(*verbose==ON) printf("Done.\n\n");
  }
  return;
}
//-----------------------------------------------------------------------------
// LOAD GAUGE DATA
void LoadGauges(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,Stage *Locptr,int *verbose)
{
  //Added by Jeff Neal, 22 Jul 2011
  //Provides functionality to output regular section measurements of discharge

  int i;
  char dum[10];
  FILE *fp;

  fp=fopen(Fnameptr->gaugefilename,"r");
  Statesptr->gsection=ON;
  if(fp==NULL)
  {
    if(*verbose==ON) printf("Gauges off\n");
    Statesptr->gsection=OFF;
    return;
  }
  if(*verbose==ON) printf("\nGauge section information:\t%s\n",Fnameptr->gaugefilename);

  fscanf(fp,"%d",&Locptr->Ngauges);
  fgetc(fp); // Retrieve closing EOL

  Locptr->gauge_loc_x=new double[Locptr->Ngauges]();
  Locptr->gauge_loc_y=new double[Locptr->Ngauges]();
  Locptr->gauge_grid_x=new int[Locptr->Ngauges]();
  Locptr->gauge_grid_y=new int[Locptr->Ngauges]();
  Locptr->gauge_grid_xy=new int[Locptr->Ngauges]();
  Locptr->gauge_dir=new int[Locptr->Ngauges]();
  Locptr->gauge_dist=new double[Locptr->Ngauges]();
  Locptr->gauge_cells=new int[Locptr->Ngauges]();

  //scan x,y locations from file
  for(i=0;i<Locptr->Ngauges;i++)
  {
    fscanf(fp,"%lf",&Locptr->gauge_loc_x[i]);
    fscanf(fp,"%lf",&Locptr->gauge_loc_y[i]);
    fscanf(fp,"%s",dum);
    if(!strcmp(dum,"N")) Locptr->gauge_dir[i] = 1;
    if(!strcmp(dum,"E")) Locptr->gauge_dir[i] = 2;
    if(!strcmp(dum,"S")) Locptr->gauge_dir[i] = 3;
    if(!strcmp(dum,"W")) Locptr->gauge_dir[i] = 4;
    if(!strcmp(dum,"n")) Locptr->gauge_dir[i] = 1;
    if(!strcmp(dum,"e")) Locptr->gauge_dir[i] = 2;
    if(!strcmp(dum,"s")) Locptr->gauge_dir[i] = 3;
    if(!strcmp(dum,"w")) Locptr->gauge_dir[i] = 4;
    fscanf(fp,"%lf",&Locptr->gauge_dist[i]);
    Locptr->gauge_cells[i] = int(ceil(Locptr->gauge_dist[i]/Parptr->dx)); // work out number of cells
  }
  for(i=0;i<Locptr->Ngauges;i++)
  {
    // convert coordinates to cells
    Locptr->gauge_grid_x[i]=int(floor((Locptr->gauge_loc_x[i]-Parptr->blx)/Parptr->dx));
    Locptr->gauge_grid_y[i]=Parptr->ysz-1-(int(floor((Locptr->gauge_loc_y[i]-Parptr->bly)/Parptr->dy)));

    // check for off-image values and set to domain edge
    if (Locptr->gauge_grid_x[i]<0) Locptr->gauge_grid_x[i] = 0;
    if (Locptr->gauge_grid_x[i]>=Parptr->xsz) Locptr->gauge_grid_x[i] = Parptr->xsz-1;
    if (Locptr->gauge_grid_y[i]<0) Locptr->gauge_grid_y[i] = 0;
    if (Locptr->gauge_grid_y[i]>=Parptr->ysz) Locptr->gauge_grid_y[i] = Parptr->ysz-1;

    // work out the location in the xy vector
    if (Locptr->gauge_dir[i] == 1) Locptr->gauge_grid_xy[i]= Locptr->gauge_grid_x[i]   + Locptr->gauge_grid_y[i]     *(Parptr->xsz+1);
    if (Locptr->gauge_dir[i] == 2) Locptr->gauge_grid_xy[i]= Locptr->gauge_grid_x[i]+1 + Locptr->gauge_grid_y[i]     *(Parptr->xsz+1);
    if (Locptr->gauge_dir[i] == 3) Locptr->gauge_grid_xy[i]= Locptr->gauge_grid_x[i]   + (Locptr->gauge_grid_y[i]+1) *(Parptr->xsz+1);
    if (Locptr->gauge_dir[i] == 4) Locptr->gauge_grid_xy[i]= Locptr->gauge_grid_x[i]   + Locptr->gauge_grid_y[i]     *(Parptr->xsz+1);

    // adjust distances if these will go off the domain - check these
    if (Locptr->gauge_dir[i] == 1 && Locptr->gauge_grid_y[i]-Locptr->gauge_cells[i] < 0) Locptr->gauge_cells[i] = Locptr->gauge_grid_y[i];
    else if (Locptr->gauge_dir[i] == 2 && Locptr->gauge_grid_x[i]+Locptr->gauge_cells[i] > Parptr->xsz-1) Locptr->gauge_cells[i] = Parptr->xsz-1-Locptr->gauge_grid_x[i];
    else if (Locptr->gauge_dir[i] == 3 && Locptr->gauge_grid_y[i]+Locptr->gauge_cells[i] > Parptr->ysz-1) Locptr->gauge_cells[i] = Parptr->ysz-1-Locptr->gauge_grid_y[i];
    else if (Locptr->gauge_dir[i] == 4 && Locptr->gauge_grid_x[i]-Locptr->gauge_cells[i] < 0) Locptr->gauge_cells[i] = Locptr->gauge_grid_x[i];
  }

  if(*verbose==ON) printf("Done.\n\n");

  fclose(fp);
  return;
}
//-----------------------------------------------------------------------------
// LOAD INITIAL DEPTHS FROM BINARY FILE
void LoadBinaryStart(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,Arrays *Arrptr, int *verbose)
{
  FILE *fp;
  int i,j, tmpi;
  double tmp, no_data_value=-9999;

  fp=fopen(Fnameptr->startfilename,"rb");
  if(fp==NULL)
  {
    if(*verbose==ON) printf("\nWARNING: Unable to load initial binary startfile:\t%s\t\n",Fnameptr->startfilename);
    return;
  }

  if(*verbose==ON) printf("Loading initial depths:\t%s\t",Fnameptr->startfilename);

  // read and dump header information, but check data are compatable
  fread(&tmpi, sizeof(int), 1, fp);
  if (tmpi != Parptr->xsz) printf("\nWARNING: incorrect number of cells in .start file\n");
  fread(&tmpi, sizeof(int), 1, fp);
  if (tmpi != Parptr->ysz) printf("\nWARNING: incorrect number of cells in .start file\n");

  fread(&tmp, sizeof(double), 1, fp);
  fread(&tmp, sizeof(double), 1, fp);
  fread(&tmp, sizeof(double), 1, fp);
  fread(&no_data_value, sizeof(double), 1, fp);

  for(j=0;j<Parptr->ysz;j++) //for(i=0;i<Parptr->xsz;i++) no need to loop x with fread.
  {
    fread(Arrptr->H+j*Parptr->xsz, sizeof(double), Parptr->xsz, fp);
    // loop through x
    for(i=0;i<Parptr->xsz;i++)
    {
      // Set depth to zero if no_data
      if((int)Arrptr->H[i+j*Parptr->xsz]==no_data_value) Arrptr->H[i+j*Parptr->xsz]=0.0;
      else if (Statesptr->startelev == ON) // if elevation file used for initial depths this needs to be converted to depths
	    {
        // check to see if SGC is on
        if (Statesptr->SGC == ON) Arrptr->H[i+j*Parptr->xsz] = getmax(Arrptr->H[i+j*Parptr->xsz] - Arrptr->SGCz[i+j*Parptr->xsz], 0.0 );
        else Arrptr->H[i+j*Parptr->xsz] = getmax(Arrptr->H[i+j*Parptr->xsz] - Arrptr->DEM[i+j*Parptr->xsz], 0.0 );
      }
    }
  }
  fclose(fp);

  if(*verbose==ON) printf("Done.\n\n");

  return;
}
// LOAD SGC PARAMETER DATA
void LoadSGCChanPrams(Fnames *Fnameptr,States *Statesptr,Pars *Parptr,SGCprams *SGCptr, int *verbose)
{
  //Added by Jeff Neal, 06 Aug 2012
  //Provides functionality to import difstrbuted sub-grid channel parameters

  int i, j, tmp, buff_size=800;
  char buff[800];
  FILE *fp;

  fp=fopen(Fnameptr->SGCchanpramsfilename,"r");
  if(fp==NULL)
  {
    if(*verbose==ON) printf("Failed to read SGC channel parameters\n");
    Statesptr->SGCchanprams=OFF;
    Statesptr->SGCchangroup=OFF;
    return;
  }
  if(*verbose==ON) printf("\nSGC channel parameter information:\t%s\n",Fnameptr->SGCchanpramsfilename);

  for(j=0;j<buff_size;j++)
  {
    buff[j]=fgetc(fp);
    if(buff[j] =='\r\n' || buff[j]=='\n' || buff[j]==EOF) break;
  }
  buff[j]='\0';									// Finish off string
  sscanf(buff,"%i",&SGCptr->NSGCprams);

  // create new variables
  SGCptr->SGCchantype=new int   [SGCptr->NSGCprams]();
  SGCptr->SGCp       =new double[SGCptr->NSGCprams]();
  SGCptr->SGCr       =new double[SGCptr->NSGCprams]();
  SGCptr->SGCs       =new double[SGCptr->NSGCprams]();
  SGCptr->SGCn       =new double[SGCptr->NSGCprams]();
  SGCptr->SGCm       =new double[SGCptr->NSGCprams]();
  SGCptr->SGCa       =new double[SGCptr->NSGCprams]();

  if(*verbose==ON) printf("Num   Type  p     r     sl    n     m     a    \n");
  //scan x,y locations from file
  for(i=0;i<SGCptr->NSGCprams;i++)
  {
    // initalise with defaults
    SGCptr->SGCchantype[i]	=	1;
    SGCptr->SGCp[i]			=	Parptr->SGC_p;
    SGCptr->SGCr[i]			=	Parptr->SGC_r;
    SGCptr->SGCs[i]			=	Parptr->SGC_s;
    SGCptr->SGCn[i]			=	Parptr->SGC_n;
    SGCptr->SGCm[i]			=	Parptr->SGC_m;
    SGCptr->SGCa[i]			=	Parptr->SGC_a;
    // load buffer until EOL
    for(j=0;j<buff_size;j++)
    {
      buff[j]=fgetc(fp);
      if(buff[j] =='\r\n' || buff[j]=='\n' || buff[j]==EOF) break;
    }
    buff[j]='\0';									// Finish off string
    sscanf(buff,"%i%i%lf%lf%lf%lf%lf%lf",&tmp,&SGCptr->SGCchantype[i],&SGCptr->SGCp[i],&SGCptr->SGCr[i],&SGCptr->SGCs[i],&SGCptr->SGCn[i],&SGCptr->SGCm[i],&SGCptr->SGCa[i]);

    if (*verbose==ON && SGCptr->SGCchantype[i] == 2 && SGCptr->SGCs[i] > 20)  printf("Warning channel shape exponent above recomended value");
    if (*verbose==ON && SGCptr->SGCchantype[i] == 2 && SGCptr->SGCs[i] < 1.3) printf("Warning channel shape exponent below recomended value");
    if (*verbose==ON && SGCptr->SGCs[i] < 0.0) printf("ERROR SGC meander coefficient is too low!");
    if (*verbose==ON) printf("%i     %i     %.3f %.3f %.3f %.3f %.3f %.3f\n", i, SGCptr->SGCchantype[i],SGCptr->SGCp[i],SGCptr->SGCr[i],SGCptr->SGCs[i],SGCptr->SGCn[i],SGCptr->SGCm[i],SGCptr->SGCa[i]);
  }
  if(*verbose==ON) printf("Done.\n\n");

  fclose(fp);
  return;
}
