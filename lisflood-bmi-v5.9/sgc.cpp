/*
*****************************************************************************
SUB GRID CHANNELS
---------------------

Simple sub grid channels method that calculates channel bed from bank elevations
and a width depth function. JCN

Keyword in parameter file is "SGC" and filename "SGCwidth" "SGCbank".

*****************************************************************************
*/

#include "lisflood.h"
#include "VersionHistory.h"

//---------------------------------------------------------------------------
// CALCULATES SGC FLOODPLAIN
void SGC_FloodplainQ(States *Statesptr,Pars *Parptr,Solver *Solverptr,Arrays *Arrptr,SGCprams *SGCptr)
{
  int i,j,*wiptr;
  double h0,h1;
  double *hptr0,*hptr1,*qptr;

  // SGC uses time step calculated in update H and stored in Solverptr->SGCtmpTstep
  Solverptr->Tstep = Solverptr->SGCtmpTstep;

  // Calculate Qx
#pragma omp parallel for private( i, h0, h1, hptr0,qptr,wiptr)
  for(j=0;j<Parptr->ysz;j++)
  {
    hptr0=Arrptr->H+j*Parptr->xsz;
    qptr=Arrptr->Qx+j*(Parptr->xsz+1)+1;
    wiptr=Arrptr->Weir_Identx+j*(Parptr->xsz+1)+1;
    for(i=0;i<Parptr->xsz-1;i++)
    {
      h0=*hptr0;
      h1=*(hptr0+1);
      *qptr=0.0;
      if(h0>Solverptr->DepthThresh || h1>Solverptr->DepthThresh)
      {
        if(Statesptr->weirs==ON && *wiptr!=-1) *qptr=CalcWeirQx(i,j,Parptr, Arrptr, Solverptr, Statesptr, SGCptr);
        else *qptr=CalcFPQxSGC(i,j,Statesptr, Parptr, Solverptr, Arrptr, SGCptr);
      }
      qptr++;
      hptr0++;
      wiptr++;
    }
  }
  // Calculate Qy
#pragma omp parallel for private( i, h0, h1, hptr0,hptr1,qptr,wiptr)
  for(j=0;j<Parptr->ysz-1;j++)
  {
    hptr0=Arrptr->H+j*Parptr->xsz;
    hptr1=Arrptr->H+(j+1)*Parptr->xsz;
    qptr=Arrptr->Qy+(j+1)*(Parptr->xsz+1);
    wiptr=Arrptr->Weir_Identy+(j+1)*(Parptr->xsz+1);
    for(i=0;i<Parptr->xsz;i++)
    {
      h0=*hptr0;
      h1=*hptr1;
      *qptr=0.0;
      if(h0>Solverptr->DepthThresh || h1>Solverptr->DepthThresh)
      {
        if(Statesptr->weirs==ON && *wiptr!=-1) *qptr=CalcWeirQy(i,j, Parptr, Arrptr, Solverptr, Statesptr, SGCptr);
        else *qptr=CalcFPQySGC(i,j,Statesptr, Parptr, Solverptr, Arrptr, SGCptr);
      }
      hptr0++;
      hptr1++;
      qptr++;
      wiptr++;
    }
  }
  return;
}
// CALCULATE FLOODPLAIN AND SGC FLOWS
double CalcFPQxSGC(int i,int j,States *Statesptr,Pars *Parptr,Solver *Solverptr,Arrays *Arrptr,SGCprams *SGCptr)
{
  double z0,z1,h0,h1,Sf=0.0,hflow,fn,dh=0.0,Q=0.0,Qc=0.0,q0,g,dt,we=0.0,w0,w1,zb0,zb1,qc, A0,A1, R,bf0,bf1, cn, m;
  int p0,p1,pq0,gr0,gr1;

  p0=i+j*Parptr->xsz;
  p1=i+1+j*Parptr->xsz;
  pq0=i+j*(Parptr->xsz+1)+1;
  // Geometry variables
  z0=Arrptr->DEM[p0]; z1=Arrptr->DEM[p1];
  w0=Arrptr->SGCwidth[p0]; w1=Arrptr->SGCwidth[p1];
  zb0=Arrptr->SGCz[p0]; zb1=Arrptr->SGCz[p1];
  // State variables
  h0=Arrptr->H[p0]; h1=Arrptr->H[p1];
  // Parameters
  g=Solverptr->g; dt=Solverptr->Tstep;
  // Bankfull depths
  bf0 = Arrptr->SGCbfH[p0];    bf1 = Arrptr->SGCbfH[p1];
  if(Arrptr->Manningsn!=NULL)   fn=0.5*(Arrptr->Manningsn[p0]+Arrptr->Manningsn[p1]); // this could be pre-calculated
  else fn=Parptr->FPn;
  //////////////////SG flows ////////////////////////////////
  // if both cells are a SGC
  if (w0>0.0 && w1>0.0)
  {
    // calculate sub grid channel Qc
    hflow=getmax(zb0+h0,zb1+h1)-getmax(zb0,zb1);
    if(hflow>Solverptr->DepthThresh)
    {
      qc=Arrptr->QxSGold[pq0]; // Get old q in m3/s
      gr0 = Arrptr->SGCgroup[p0]; gr1 = Arrptr->SGCgroup[p1]; // get channel groups
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
        Arrptr->QxSGold[pq0] =  (qc-g*A0*dt*Sf) / (1+dt*g*cn*fabs(qc) / (pow(R,4.0/3.0)*A0) );
        we = w0;
      }
      else
      {
        // calculate hydraulic radius for SGC
        R = CalcSGC_R(gr1, hflow, bf1, w1, Arrptr->SGCwidth[p1], A1, SGCptr);
        // Calculate channel flow using inertial wave equation --- in this case Qxold will be in m3s-1 not m2s-1
        Arrptr->QxSGold[pq0] =  (qc-g*A1*dt*Sf) / (1+dt*g*cn*fabs(qc) / (pow(R,4.0/3.0)*A1) );
        we = w1; // set effective width
      }
      Qc = Arrptr->QxSGold[pq0]; // Update channel flows
    }
    else Arrptr->QxSGold[pq0] = 0.0; // make sure previous flux is set to zero
  }
  //////////////////FP flows ////////////////////////////////
  // don't calculate floodplain flow if we (channel width) is greater than cell width, note we is initalised as zero
  //if (we < Parptr->dx) //CCS_deletion
  if (we < Arrptr->dy[p0]) // CCS dy[p0] for cell width relative to flow in x
  {
    // Correct h for any SGC depth
    if (w0>0.0) h0 = getmax(h0-bf0, 0.0);
    if (w1>0.0) h1 = getmax(h1-bf1, 0.0);
    // Calculating hflow based on floodplain levels
    hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
    if(hflow>Solverptr->DepthThresh)
    {
      q0=Arrptr->Qxold[pq0]; // in m2/s
      // calculate Qfp
      dh=z0+h0-z1-h1;
      //Sf=-dh/Parptr->dx; //CCS_deletion
      Sf=-dh/Arrptr->dx[p0];
      /* inertial model.  If routing scheme is off we always calc Q.
         If routing scheme is on, we do not calc Q where water surface gradient is steep*/
      if (Statesptr->routing==OFF || fabs(Sf)<Parptr->RouteSfThresh)
      {
        Arrptr->Qxold[pq0]=(q0-(g*dt*hflow*Sf))/(1+g*dt*fn*fn*fabs(q0)/pow(hflow,(7.0/3.0)));//*Parptr->dx;
        //Q= Arrptr->Qxold[pq0]*(Parptr->dx-we); //CCS_deletion
        Q= Arrptr->Qxold[pq0]*(Arrptr->dx[p0]-we);
      }
      else Arrptr->Qxold[pq0]=0.0; // Set Qxold to zero
    }
    else Arrptr->Qxold[pq0]=0.0; // Set Qxold to zero
  }
  // sub grid channel combuined Q
  Q = Qc + Q;
  // option to save V's
  if (Statesptr->voutput==ON)
  {
    //if ((Q-Qc) == 0.0 || we >= Parptr->dx) Arrptr->Vx[pq0]=0.0; //CCS_deletion
    if ((Q-Qc) == 0.0 || we >= Arrptr->dy[p0]) Arrptr->Vx[pq0]=0.0; // CCS dy[p0] for cell width relative to channel flow in x
    else
    {
      //Arrptr->Vx[pq0]=(Q-Qc)/(Parptr->dx-we)/hflow; //CCS_deletion
      Arrptr->Vx[pq0]=(Q-Qc)/(Arrptr->dx[p0]-we)/hflow;
      Arrptr->maxVx[pq0]=getmax(Arrptr->maxVx[pq0],fabs(Arrptr->Vx[pq0]));
    }
  }
  return(Q);
}
//-----------------------------------------------------------------------------------
// CALCULATE FLOODPLAIN AND SGC FLOWS
double CalcFPQySGC(int i,int j,States *Statesptr,Pars *Parptr,Solver *Solverptr,Arrays *Arrptr,SGCprams *SGCptr)
{
  double z0,z1,h0,h1,Sf=0.0,hflow,fn,dh=0.0,Q=0.0,Qc=0.0,q0,g,dt,we=0.0,w0,w1,zb0,zb1,qc, A0,A1, R,bf0,bf1, cn, m;
  int p0,p1,pq0,gr0,gr1;

  p0=i+j*Parptr->xsz;
  p1=i+(j+1)*Parptr->xsz;
  pq0=i+(j+1)*(Parptr->xsz+1);
  // Geometry variables
  z0=Arrptr->DEM[p0]; z1=Arrptr->DEM[p1];
  w0=Arrptr->SGCwidth[p0]; w1=Arrptr->SGCwidth[p1];
  zb0=Arrptr->SGCz[p0]; zb1=Arrptr->SGCz[p1];
  // State variables
  h0=Arrptr->H[p0]; h1=Arrptr->H[p1];
  // Parameters
  g=Solverptr->g; dt=Solverptr->Tstep;
  // Bankfull depths
  bf0 = Arrptr->SGCbfH[p0];    bf1 = Arrptr->SGCbfH[p1];
  if(Arrptr->Manningsn!=NULL)   fn=0.5*(Arrptr->Manningsn[p0]+Arrptr->Manningsn[p1]); // this could be pre-calculated
  else fn=Parptr->FPn;
  //////////////////SG flows ////////////////////////////////
  // if both cells are a SGC
  if (w0>0.0 && w1>0.0)
  {
    // calculate sub grid channel Qc
    hflow=getmax(zb0+h0,zb1+h1)-getmax(zb0,zb1);
    if(hflow>Solverptr->DepthThresh)
    {
      qc=Arrptr->QySGold[pq0]; // Get old q in m3/s
      gr0 = Arrptr->SGCgroup[p0]; gr1 = Arrptr->SGCgroup[p1]; // get channel groups
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
        Arrptr->QySGold[pq0] =  (qc-g*A0*dt*Sf) / (1+dt*g*cn*fabs(qc) / (pow(R,4.0/3.0)*A0) );
        we = w0;
      }
      else
      {
        // calculate hydraulic radius for SGC
        R = CalcSGC_R(gr1, hflow, bf1, w1, Arrptr->SGCwidth[p1], A1, SGCptr);
        // Calculate channel flow using inertial wave equation --- in this case Qxold will be in m3s-1 not m2s-1
        Arrptr->QySGold[pq0] =  (qc-g*A1*dt*Sf) / (1+dt*g*cn*fabs(qc) / (pow(R,4.0/3.0)*A1) );
        we = w1; // set effective width
      }
      Qc = Arrptr->QySGold[pq0]; // Update channel flows
    }
    else Arrptr->QySGold[pq0] = 0.0; // make sure previous flux is set to zero
  }
  /////////////////FP flows ////////////////////////////////
  // don't calculate floodplain flow if we (channel width) is greater than cell width, note we is initalised as zero
  //if (we < Parptr->dx) //CCS_deletion
  if (we < Arrptr->dx[p0]) // CCS use dx[p0] for cell width relative to flow in y
  {
    // Correct h for any SGC depth
    if (w0>0) h0 = getmax(h0-bf0, 0.0);
    if (w1>0) h1 = getmax(h1-bf1, 0.0);
    // calculate Qfp
    hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
    if(hflow>Solverptr->DepthThresh)
    {
      q0=Arrptr->Qyold[pq0]; // in m2/s
      dh=z0+h0-z1-h1;
      //Sf=-dh/Parptr->dx; //CCS_deletion
      Sf=-dh/Arrptr->dy[p0];
      /* inertial model.  If routing scheme is off we always calc Q.
         If routing scheme is on, we do not calc Q where water surface gradient is steep*/
      if (Statesptr->routing==OFF || fabs(Sf)<Parptr->RouteSfThresh)
      {
        Arrptr->Qyold[pq0]=(q0-(g*dt*hflow*Sf))/(1+g*dt*fn*fn*fabs(q0)/pow(hflow,(7.0/3.0)));
        //Q= Arrptr->Qyold[pq0]*(Parptr->dx-we); //CCS_deletion
        Q=Arrptr->Qyold[pq0]*(Arrptr->dy[p0]-we);
      }
      else Arrptr->Qyold[pq0]=0.0; // Set Qxold to zero
    }
    else Arrptr->Qyold[pq0]=0.0; // make sure previous flux is set to zero
  }
  // sub grid channel combined Q
  Q = Qc + Q;
  // option to save V's
  if (Statesptr->voutput==ON)
  {
    //if ((Q-Qc) == 0.0 || we >= Parptr->dx) Arrptr->Vy[pq0]=0.0;// check Qfp is above CCS_deletion
    if ((Q-Qc) == 0.0 || we >= Arrptr->dx[p0]) Arrptr->Vy[pq0]=0.0; // check Qfp is above  CCS dx[p0] for cell width relative to channel flow in y
    else
    {
      //Arrptr->Vy[pq0]=(Q-Qc)/(Parptr->dx-we)/hflow; //CCS_deletion
      Arrptr->Vy[pq0]=(Q-Qc)/(Arrptr->dy[p0]-we)/hflow;
      Arrptr->maxVy[pq0]=getmax(Arrptr->maxVy[pq0],fabs(Arrptr->Vy[pq0]));
    }
  }
  return(Q);
}
//-----------------------------------------------------------------------------------
// CALCULATE SGC bed elevation
void CalcSGCz(Fnames *Fnameptr,States *Statesptr,Pars *Parptr, Arrays *Arrptr,SGCprams *SGCptr, int *verbose)
{
  FILE *fpb, *fpa;
  int i,j,p0, gr, ct;
  double w0, z0, r, p, sl, m, a, zf, dx, dy, chan_length;
  char dum[80];
  double no_data=-9999;

  if (Statesptr->SGCchangroup == OFF) // spatially uniform channel parameters
  {
    SGCptr->NSGCprams               = 1; // there is only one global channel type
    // create variables
    SGCptr->SGCchantype             = new int   [SGCptr->NSGCprams]();
    SGCptr->SGCr                    = new double[SGCptr->NSGCprams](); SGCptr->SGCp                 = new double[SGCptr->NSGCprams]();
    SGCptr->SGCs                    = new double[SGCptr->NSGCprams](); SGCptr->SGCn                 = new double[SGCptr->NSGCprams]();
    SGCptr->SGCm                    = new double[SGCptr->NSGCprams](); SGCptr->SGCa                 = new double[SGCptr->NSGCprams]();
    // set uniform channel parameters
    SGCptr->SGCchantype[0]  = ct = Parptr->SGCchan_type;
    SGCptr->SGCr[0]                 = r  = Parptr->SGC_r;
    SGCptr->SGCp[0]                 = p  = Parptr->SGC_p;
    SGCptr->SGCs[0]                 = sl = Parptr->SGC_s;
    SGCptr->SGCn[0]                 =      Parptr->SGC_n; // Note converted to n^2 below
    SGCptr->SGCm[0]                 = m  = Parptr->SGC_m; // meander coefficient for SGC model, default 1, else set in .par file
    SGCptr->SGCa[0]                 = a  = Parptr->SGC_a; // Calchment area (or other external explanitory variable) default is -1 and not used.
  }
  // Open upstream catchment area for channel if it exsists
  if (Statesptr->SGCcat_area==ON)
  {
    fpa=fopen(Fnameptr->SGCbedfilename,"r");
    if(fpa==NULL)
    {
      if(*verbose==ON) printf("Failed to Load SGCcat_area:\t%s\t using channel width\n",Fnameptr->SGCbedfilename);
      fclose(fpa); Statesptr->SGCcat_area=OFF;
    }
    else
    {
      if(*verbose==ON) printf("Loading SGCcat_area:\t%s\t",Fnameptr->SGCbedfilename);
      for(i=0;i<5;i++) fscanf(fpa,"%s %s",dum,dum);
      fscanf(fpa,"%s %lf",dum,&no_data);
    }
  }
  // Open bed and channel area files if they exsist
  if (Statesptr->SGCbed==ON)
  {
    fpb=fopen(Fnameptr->SGCbedfilename,"r");
    if(fpb==NULL)
    {
      if(*verbose==ON) printf("Failed to Load SGCbed:\t%s\t using hydraulic geometry\n",Fnameptr->SGCbedfilename);
      fclose(fpb); Statesptr->SGCbed=OFF;
    }
    else
    {
      if(*verbose==ON) printf("Loading SGCbed:\t%s\t",Fnameptr->SGCbedfilename);
      for(i=0;i<5;i++) fscanf(fpb,"%s %s",dum,dum);
      fscanf(fpb,"%s %lf",dum,&no_data);
    }
  }
  // Populate SGCptr with further parameters needed by certain channel models (such as the power shaped channel)
  SGC_wp_prams (SGCptr);
  // changes n to  n^2 to save on computation
  for(i=0;i<SGCptr->NSGCprams;i++) SGCptr->SGCn[i] = SGCptr->SGCn[i]*SGCptr->SGCn[i];
  // set Parptr->SGC_m to a maximum of 1 to prevent time step increasing for meander coefficients above 1;
  Parptr->SGC_m = getmin(1.0,Parptr->SGC_m);

  for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++) // loop through model domain
                             {
                               p0=i+j*Parptr->xsz;      // cell index
                               w0=Arrptr->SGCwidth[p0]; // channel width
                               z0=Arrptr->SGCz[p0];     // bank height

                               // Work out the channel type and parameters
                               if (Statesptr->SGCchangroup == ON) // spatially varying parameters
                               {
                                 // Which group does the channel belong to
                                 gr = Arrptr->SGCgroup[p0];
                                 ct = SGCptr->SGCchantype[gr];
                                 r  = SGCptr->SGCr[gr];
                                 p  = SGCptr->SGCp[gr];
                                 a  = SGCptr->SGCa[gr];
                                 sl = SGCptr->SGCs[gr];                    // e.g. trapazoidal cahnnel slope or power law coefficient
                                 m  = SGCptr->SGCm[gr];                     // sub-grid channel meander coefficient.
                                 Parptr->SGC_m = getmin(m, Parptr->SGC_m); // keep track of minimum meander coefficient in Parptr->SGC_m for time step calculation
                               }

                               // Get bed elevation and catchement area from file if they exsist... must read every cell so needs to be before if(w0>0.0)
                               if (Statesptr->SGCcat_area==ON)  fscanf(fpa,"%lf ",&a); // read catchment upstream area from file, if not present will be -1
                               if (Statesptr->SGCbed==ON)  fscanf(fpb,"%lf ",&zf);     // set zf to bed elevation from file if it exsists
                               else zf = no_data;                                      // else set bed elevation from file to no_data

                               if(w0 > 0.0) // if there is a sub-grid cahnnel
                               {
                                 dx = Arrptr->dx[p0];
                                 dy = Arrptr->dy[p0];

                                 // Calculate the length of sub-grid channel in a cell
                                 chan_length = 0.0; // set length of channel to zer
                                 if (i == 0               || Arrptr->SGCwidth[p0-1] > 0.0)           chan_length += 0.5*dx;    // check for domain edge or west side trib
                                 if (i == (Parptr->xsz-1) || Arrptr->SGCwidth[p0+1] > 0.0)           chan_length += 0.5*dx;    // check for domain edge or east side trib
                                 if (j == 0               || Arrptr->SGCwidth[p0-Parptr->xsz] > 0.0) chan_length += 0.5*dy;    // check for domain edge or north side trib
                                 if (j == (Parptr->ysz-1) || Arrptr->SGCwidth[p0+Parptr->xsz] > 0.0) chan_length += 0.5*dy;    // check for domain edge or south side trib
                                 Arrptr->SGCc[p0] = m*getmax(chan_length,getmin(dx,dy)); // make sure channel is at least as long as the shortest cell width then apply meander coefficient m

                                 // check that we*SGCtribs is not greater than dA to prevent overdoing the cell volume... unless w0 is greater than the cell width already in which case carry on with the larger cell surface area
                                 // note from above that if w0 is greater than dx then SGCc will be dx
                                 if (w0 < getmax(dx,dy) && w0*Arrptr->SGCc[p0] > Arrptr->dA[p0]) Arrptr->SGCc[p0] = Arrptr->dA[p0]/w0;

                                 if (Statesptr->SGCbed==ON && zf != no_data)      Arrptr->SGCz[p0] = getmin (z0 - 0.01       , zf);                                               // If bed elevation is specified in the bedfile set bed elevation bed to the bed file elevation and if above the DEM set to 0.01 below the DEM
                                 else if (Statesptr->SGCbfh_mode == ON)           Arrptr->SGCz[p0] = getmin (z0 - p , Arrptr->DEM[p0]-0.01);                              // If the model is told to use p as depth use p as depth e.g. SGCbfh_mode
                                 else if (Statesptr->SGCA_mode == ON)                     Arrptr->SGCz[p0] = getmin (z0 - p/w0 , Arrptr->DEM[p0]-0.01);                   // If the model is told to use p as channel area use p as area and divide by width e.g. SGCA_mode
                                 else if (a > 0)                                                          Arrptr->SGCz[p0] = getmin (z0 - r*pow(a,p) , Arrptr->DEM[p0]-0.01);             // Calculate z using the hydraulic geometry and catchement area, and if above the DEM set to 0.01 below the DEM
                                 else                                                                                     Arrptr->SGCz[p0] = getmin (z0 - r*pow(w0,p), Arrptr->DEM[p0]-0.01);             // Calculate z using the hydraulic geometry and channel width, and if above the DEM set to 0.01 below the DEM

                                 // Build the channel information into Arrptr from its parameters
                                 switch (ct)
                                 {
                                 case 1: // Rectangular channel (default) - This model has a top width and bed elevation and top
                                   Arrptr->SGCbfH[p0] = Arrptr->DEM[p0]-Arrptr->SGCz[p0]; // Calculate SGC bank full depth
                                   Arrptr->SGCc[p0]   = Arrptr->SGCc[p0]*Arrptr->SGCwidth[p0]; // calculate SGCc
                                   Arrptr->SGCbfV[p0] = Arrptr->SGCc[p0]*Arrptr->SGCbfH[p0]; // calculate bank full volume
                                   break;       // Break terminates the switch statement

                                 case 2: // h = x^c channel
                                   Arrptr->SGCbfH[p0] = Arrptr->DEM[p0]-Arrptr->SGCz[p0]; // Calculate SGC bank full depth
                                   Arrptr->SGCc[p0]   = Arrptr->SGCc[p0]*Arrptr->SGCwidth[p0]*sl*pow(Arrptr->SGCbfH[p0],-1.0/sl)/(sl+1.0); // calculate SGCc
                                   Arrptr->SGCbfV[p0] = Arrptr->SGCc[p0]*pow(Arrptr->SGCbfH[p0],1.0/sl +1.0); // calculate bank full volume
                                   break;  // Break terminates the switch statement

                                 case 3: // linear slope - This model has a bed elevation, top width and top elevation.
                                   Arrptr->SGCbfH[p0] = Arrptr->DEM[p0]-Arrptr->SGCz[p0]; // Calculate SGC bank full depth
                                   sl                 = w0/(getmin(z0,Arrptr->DEM[p0])-Arrptr->SGCz[p0]);  // calculate section slope using the minimum of DEM and z0 to preserve channel width at bank full, this channel will be rectangular with a linear slop bottom.
                                   if (Arrptr->DEM[p0] > z0) Arrptr->SGCwidth[p0] = (Arrptr->DEM[p0] - Arrptr->SGCz[p0])*sl; // Calculate a new we if DEM is above bank height (e.g. we needs to be bigger)
                                   Arrptr->SGCc[p0]   = 0.5*Arrptr->SGCc[p0]*Arrptr->SGCwidth[p0]; // calculate SGCc
                                   Arrptr->SGCbfV[p0] = Arrptr->SGCc[p0]*Arrptr->SGCbfH[p0]; // calculate bank full volume
                                   break;       // Break terminates the switch statement

                                 case 4: // triangular channel - This model has the bed elevation, top width and top elevation
                                   Arrptr->SGCbfH[p0] = Arrptr->DEM[p0]-Arrptr->SGCz[p0]; // Calculate SGC bank full depth
                                   sl = w0/(2*(getmin(z0,Arrptr->DEM[p0])-Arrptr->SGCz[p0]));  // calculate section slope using the minimum of DEM and z0 to preserve channel width at bank full, this channel will be rectangular with a linear slop bottom.
                                   // Calculate a new we if DEM is above bank height (e.g. we needs to be bigger)
                                   if (Arrptr->DEM[p0] > z0) Arrptr->SGCwidth[p0] = 2.0*(Arrptr->DEM[p0] - Arrptr->SGCz[p0])*sl;
                                   Arrptr->SGCc[p0]   = 0.5*Arrptr->SGCc[p0]*Arrptr->SGCwidth[p0]; // calculate SGCc
                                   Arrptr->SGCbfV[p0] = Arrptr->SGCc[p0]*Arrptr->SGCbfH[p0]; // calculate bank full volume
                                   break;       // Break terminates the switch statement

                                 case 5: // parabolic channel
                                   Arrptr->SGCbfH[p0] = Arrptr->DEM[p0] -Arrptr->SGCz[p0]; // Calculate SGC bank full depth
                                   Arrptr->SGCc[p0]   = Arrptr->SGCc[p0]*Arrptr->SGCwidth[p0]*2.0*pow(Arrptr->SGCbfH[p0],-0.5)/(2.0+1.0); // calculate SGCc
                                   Arrptr->SGCbfV[p0] = Arrptr->SGCc[p0]*pow(Arrptr->SGCbfH[p0],3.0/2.0); // calculate bank full volume
                                   break;       // Break terminates the switch statement

                                 case 6: // No banks
                                   Arrptr->SGCbfH[p0] = Arrptr->DEM[p0]-Arrptr->SGCz[p0]; // Calculate SGC bank full depth
                                   Arrptr->SGCc[p0]   = Arrptr->SGCc[p0]*Arrptr->SGCwidth[p0]; // calculate SGCc
                                   Arrptr->SGCbfV[p0] = Arrptr->SGCc[p0]*Arrptr->SGCbfH[p0]; // calculate bank full volume
                                   break;       // Break terminates the switch statement

                                   /*
                                     case 7: // trapazoidal channel
                                             // sub-grid channels pre-processor for trapazoidal channel
                                             // NOTE that in the case of the trapazoidal channel width at bankfull is changed to width at bed by this pre-processor
                                             // Calculate z using the hydraulic geometry
                                             Arrptr->SGCz[p0] = getmin (z0 - r*pow(w0,p), Arrptr->DEM[p0] - 0.01);
                                             // check if DEM is channel bed and if so make the channel bed 0.01 m below the DEM such that the channel is essentially rectangular at DEM height
                                             //if (Arrptr->SGCz[p0] == Arrptr->DEM[p0]) Arrptr->SGCz[p0]-=0.01;
                                             Arrptr->SGCsl[p0] = sl;  // for the trapazoidal channel SGCsl is an input parameter
                                             // Calculate a new for trapazoidal channel bed
                                             Arrptr->SGCwidth[p0] = Arrptr->SGCwidth[p0] - 2*(Arrptr->DEM[p0]-Arrptr->SGCz[p0])*Arrptr->SGCsl[p0];
                                             break;     // Break terminates the switch statement
                                   */
                                 default: // its all gone wrong
                                   printf("Should not be here! Something is wrong with the SGC channel model in CalcSGCz");
                                   break;
                                 }
                               }
                               else  Arrptr->SGCz[p0] = Arrptr->DEM[p0]; // if no SGC set SGCz to floodplain elevation to allow easy water surface elevation outputs
                               // finally calculate SGC bank full depth from DEM and bed elevation
                               Arrptr->SGCbfH[p0] = Arrptr->DEM[p0]-Arrptr->SGCz[p0];
                             }

  // Close SGCbed and SCGcat_area files
  if (Statesptr->SGCbed==ON) fclose(fpb);
  if (Statesptr->SGCcat_area==ON) fclose(fpa);
  if(*verbose==ON) printf("SGC calc z done.\n");

  return;
}
//-----------------------------------------------------------------------
// SGC - SUM Qs INTO A CELL AND UPDATE DEPTH ACCORDINGLY
void SGC_UpdateH(States *Statesptr, Pars *Parptr, Solver *Solverptr,BoundCs *BCptr,ChannelSegmentType *ChannelSegments,Arrays *Arrptr,SGCprams *SGCptr)
{
  int i,j,pi,pTH,p0,gr;
  double *qxptr0,*qyptr0,*qyptr1;
  double dV=0.0, himp,ThreadMax=0.0,Hmax=0.0, Qpoint=0.0, Qfree, cn, fn;
  double Q_multiplier=1.0;

  // CCS Multiplier for Q inputs. If using regular grid, Qs are specified as m^2 and need to be multiplied by dx;
  // if using lat-long Qs are specified in m^3 and therefore multiplier is 1. Note intialised above as 1.0.
  if(Statesptr->latlong==OFF) Q_multiplier=Parptr->dx;

  // SGC Insert point sources ((MT) moved before dV check, in case flow out of cell results in negative H reset and the inflow would prevent this)
  // (JCN) moved HFIX and HVAR point sources after update h to prevent potential incorrect water surface elevation due dh in update H.
#pragma omp parallel for private(pTH, dV, himp, Qfree) reduction (+:Qpoint) // Note: should not be more than one point source per cell else further reduction on SGCdVol needed.
  for(pi=0;pi<BCptr->numPS;pi++)
  {
    // Set initial dV and himp as zero
    dV=0.0;
    himp=0.0;
    // location in vector
    pTH = BCptr->xpi[pi]+BCptr->ypi[pi]*Parptr->xsz;
    // different boundary conditions
    if(BCptr->PS_Ident[pi]==4) dV = BCptr->PS_Val[pi]*Q_multiplier*Solverptr->Tstep; // QFIX // Calculate change in volume
    else if(BCptr->PS_Ident[pi]==5) dV = InterpBC(BCptr->BCVarlist[(int)BCptr->PS_Val[pi]],Solverptr->t)*Q_multiplier*Solverptr->Tstep; // QVAR // Calculate change in volume
    else if (BCptr->PS_Ident[pi]==6) // point free boundary
    {
      if(Arrptr->H[pTH]>Solverptr->DepthThresh)
      {
        if(Arrptr->SGCManningsn!=NULL)   cn=Arrptr->SGCManningsn[pTH]*Arrptr->SGCManningsn[pTH]; // this could be pre-calculated
        else cn=SGCptr->SGCn[Arrptr->SGCgroup[pTH]];
        if(Arrptr->Manningsn!=NULL) fn=Arrptr->Manningsn[pTH];
        else fn=Parptr->FPn;
        CalcSGC_pointFREE(Arrptr->H[pTH], Arrptr->SGCwidth[pTH], BCptr->PS_Val[pi], Solverptr->DepthThresh, Solverptr->Tstep, getmin(Arrptr->dx[pTH], Arrptr->dy[pTH]), Solverptr->g, cn, fn, Arrptr->SGCbfH[pTH], Arrptr->SGCgroup[pTH], -1, &BCptr->PS_qold[pi], &Qfree, &BCptr->PS_qSGold[pi], SGCptr);
        dV+=Qfree*Solverptr->Tstep;
      }

    }
    // update the cell volume change and in point source Q
    Arrptr->SGCdVol[pTH]+=dV; // Add volume to SGCdVol for use later by update H e.g wait for main update H before calculating H
    Qpoint+=dV/Solverptr->Tstep; // Update Qpoint
  }

  // Calculate dV (-ve => outflow) and update NewH
#pragma omp parallel for private( i,qxptr0,qyptr0,qyptr1,dV,ThreadMax, p0, gr)
  for(j=0;j<Parptr->ysz;j++)
  {
    qxptr0=Arrptr->Qx+j*(Parptr->xsz+1);
    qyptr0=Arrptr->Qy+j*(Parptr->xsz+1);
    qyptr1=Arrptr->Qy+(j+1)*(Parptr->xsz+1);
    ThreadMax=0.0;
    for(i=0;i<Parptr->xsz;i++)
    {
      p0 = i+j*Parptr->xsz; // location in Arrptr
      dV = Solverptr->Tstep*(*qxptr0-*(qxptr0+1)+*qyptr0-*qyptr1); // compute volume change in cell
      dV+= Arrptr->SGCdVol[p0]; // add on all other volume changes in cell
      dV+=Arrptr->SGCQin[p0]*Solverptr->Tstep; // JMH
      Arrptr->SGCdVol[p0] = 0.0; // rest SGCdVol to zero
      if (dV != 0.0) // if volume has changed
      {
        Arrptr->SGCVol[p0] += dV; // Add dV to Vol. This is now the total volume.
        Arrptr->SGCVol[p0] =  getmax(Arrptr->SGCVol[p0],0.0); // if volume has gone negative set to zero, should not happen so will result in mass balance error
        if (Arrptr->SGCVol[p0] > 0.0)
        {
          gr = Arrptr->SGCgroup[p0]; // channel group number
          /* //CCS_deletion
             if (Arrptr->SGCwidth[p0] == 0.0)                                                                                                         Arrptr->H[p0] = Arrptr->SGCVol[p0]/Parptr->dA; // there is no sub-grid channel so normal updateH
             else if (Arrptr->SGCVol[p0] >= Arrptr->SGCbfV[p0] && Arrptr->SGCwidth[p0] < Parptr->dx)  Arrptr->H[p0] = Arrptr->SGCbfH[p0]+(Arrptr->SGCVol[p0]-Arrptr->SGCbfV[p0])/Parptr->dA; // there is a sub-grid channel but its over bank // Note need to include check for SGCwidth < cell width to make sure out of bank volume is not accounted for in cells with channels wider than a cell
             else                                                                                                                                                                     Arrptr->H[p0] = CalcSGC_UpH(SGCptr->SGCchantype[gr],  Arrptr->SGCVol[p0], SGCptr->SGCs[gr], Arrptr->SGCc[p0]); // there is a sub-grid channel and its within bank
          */
          if (Arrptr->SGCwidth[p0] == 0.0)                                                                                                                                                                    Arrptr->H[p0] = Arrptr->SGCVol[p0]/Arrptr->dA[p0]; // there is no sub-grid channel so normal updateH
          else if (Arrptr->SGCVol[p0] >= Arrptr->SGCbfV[p0] && Arrptr->SGCwidth[p0] < 0.5*(Arrptr->dx[p0]+Arrptr->dy[p0])     )       Arrptr->H[p0] = Arrptr->SGCbfH[p0]+(Arrptr->SGCVol[p0]-Arrptr->SGCbfV[p0])/Arrptr->dA[p0]; // there is a sub-grid channel but its over bank // Note need to include check for SGCwidth < cell width to make sure out of bank volume is not accounted for in cells with channels wider than a cell
          else                                                                                                                                                                                                                                Arrptr->H[p0] = CalcSGC_UpH(SGCptr->SGCchantype[gr],  Arrptr->SGCVol[p0], SGCptr->SGCs[gr], Arrptr->SGCc[p0]); // there is a sub-grid channel and its within bank
        }
        else Arrptr->H[p0] = 0.0; // make sure H is zero if Vol is zero
      }
      // record maximum h in the thread after UpdateH for Tstep calc
      ThreadMax=getmax(Arrptr->H[p0],ThreadMax);
      qxptr0++;
      qyptr0++;
      qyptr1++;
    }
    // find the maximum depth in all the threads
#pragma omp critical
    {
      Hmax=getmax(Hmax,ThreadMax);
    }
  }
#pragma omp parallel for private(pTH, dV, himp, gr) reduction (+:Qpoint) // Note: should not be more than one point source per cell else further reduction on SGCdVol
  for(pi=0;pi<BCptr->numPS;pi++)
  {
    // Set initial dV and himp as zero
    dV=0.0;
    himp=0.0;
    // location in vector
    pTH = BCptr->xpi[pi]+BCptr->ypi[pi]*Parptr->xsz;
    // different boundary conditions
    if(BCptr->PS_Ident[pi]==2 || BCptr->PS_Ident[pi]==3) // HFIX or HVAR
    {
      if(BCptr->PS_Ident[pi]==2) himp = BCptr->PS_Val[pi]; // HFIX
      if(BCptr->PS_Ident[pi]==3) himp = InterpBC(BCptr->BCVarlist[(int)BCptr->PS_Val[pi]],Solverptr->t); // HVAR
      if (Arrptr->SGCwidth[pTH] > 0.0)
      {
        // sub-grid channel
        himp=himp-Arrptr->SGCz[pTH];
        if(himp<0.0) himp=0.0;
        // Calculate volume after update and subtract from before update
        gr = Arrptr->SGCgroup[pTH]; // channel group number
        if (himp <= Arrptr->SGCbfH[pTH])                                dV = CalcSGC_UpV(SGCptr->SGCchantype[gr], himp, SGCptr->SGCs[gr], Arrptr->SGCc[pTH]); //Calculate channel volume
        else                                                                                    dV = Arrptr->SGCbfV[pTH] + (himp-Arrptr->SGCbfH[pTH])*Arrptr->dA[pTH]; // out of bank level
        // calculate volume before update
        if (Arrptr->H[pTH] <= Arrptr->SGCbfH[pTH])          dV -= CalcSGC_UpV(SGCptr->SGCchantype[gr], Arrptr->H[pTH], SGCptr->SGCs[gr], Arrptr->SGCc[pTH]); //Calculate channel volume
        else                                                                                    dV -= Arrptr->SGCbfV[pTH] + (Arrptr->H[pTH]-Arrptr->SGCbfH[pTH])*Arrptr->dA[pTH]; // out of bank level
        Arrptr->H[pTH] = himp;
      }
      else
      {
        // floodplain only cell
        himp=himp-Arrptr->DEM[pTH]; // Get depth on floodplain
        if(himp<0.0) himp=0.0;
        dV = (himp-Arrptr->H[pTH])*Arrptr->dA[pTH];
        Arrptr->H[pTH] = himp;

      }
      Arrptr->SGCVol[pTH] += dV;
      Qpoint+=dV/Solverptr->Tstep; // Update Qpoint
#pragma omp critical
      {
        Hmax=getmax(Hmax,himp);
      }
    }
  }

#pragma omp parallel for private(pTH, dV, himp, gr) reduction (+:Qpoint) // Note: should not be more than one point source per cell else further reduction on SGCdVol
  for(pi=0;pi<BCptr->numPS;pi++)
  {
    // Set initial dV and himp as zero
    dV=0.0;
    himp=0.0;
    // location in vector
    pTH = BCptr->xpi[pi]+BCptr->ypi[pi]*Parptr->xsz;
    // different boundary conditions
    if(BCptr->PS_Ident[pi]==2 || BCptr->PS_Ident[pi]==3) // HFIX or HVAR
    {
      if(BCptr->PS_Ident[pi]==2) himp = BCptr->PS_Val[pi]; // HFIX
      if(BCptr->PS_Ident[pi]==3) himp = InterpBC(BCptr->BCVarlist[(int)BCptr->PS_Val[pi]],Solverptr->t); // HVAR
      if (Arrptr->SGCwidth[pTH] > 0.0)
      {
        // sub-grid channel
        himp=himp-Arrptr->SGCz[pTH];
        if(himp<0.0) himp=0.0;
        // Calculate volume after update and subtract from before update
        gr = Arrptr->SGCgroup[pTH]; // channel group number
        if (himp <= Arrptr->SGCbfH[pTH])                                dV = CalcSGC_UpV(SGCptr->SGCchantype[gr], himp, SGCptr->SGCs[gr], Arrptr->SGCc[pTH]); //Calculate channel volume
        else                                                                                    dV = Arrptr->SGCbfV[pTH] + (himp-Arrptr->SGCbfH[pTH])*Arrptr->dA[pTH]; // out of bank level
        // calculate volume before update
        if (Arrptr->H[pTH] <= Arrptr->SGCbfH[pTH])          dV -= CalcSGC_UpV(SGCptr->SGCchantype[gr], Arrptr->H[pTH], SGCptr->SGCs[gr], Arrptr->SGCc[pTH]); //Calculate channel volume
        else                                                                                    dV -= Arrptr->SGCbfV[pTH] + (Arrptr->H[pTH]-Arrptr->SGCbfH[pTH])*Arrptr->dA[pTH]; // out of bank level
        Arrptr->H[pTH] = himp;
      }
      else
      {
        // floodplain only cell
        himp=himp-Arrptr->DEM[pTH]; // Get depth on floodplain
        if(himp<0.0) himp=0.0;
        dV = (himp-Arrptr->H[pTH])*Arrptr->dA[pTH];
        Arrptr->H[pTH] = himp;

      }
      Arrptr->SGCVol[pTH] += dV;
      Qpoint+=dV/Solverptr->Tstep; // Update Qpoint
#pragma omp critical
      {
        Hmax=getmax(Hmax,himp);
      }
    }
  }
  // save Qpoint
  BCptr->Qpoint = Qpoint;

  // now calulate the next time step to be used if it is samller than inital time step
  //Solverptr->SGCtmpTstep = getmin(Solverptr->cfl*Parptr->dx/(sqrt(Solverptr->g*Hmax)),Solverptr->InitTstep); //CCS_deletion
  Solverptr->SGCtmpTstep = getmin(Solverptr->cfl*Parptr->min_dx_dy*Parptr->SGC_m/(sqrt(Solverptr->g*Hmax)),Solverptr->InitTstep);

  return;
}
//-----------------------------------------------------------------------------------
// BOUNDARY CONDITIONS
// Calculate Qx and Qy at edges of the domain in response to boundary
// conditions
void SGC_BCs(States *Statesptr,Pars *Parptr,Solver *Solverptr,BoundCs *BCptr,ChannelSegmentType *ChannelSegments,Arrays *Arrptr,SGCprams *SGCptr)
{
  int BCi,numBCs,p0,p1,sign,dir,edge,pTQ,gr;
  double h0,h1,z1,hflow,fn,dh,Sf,g,w0,zb0,Qfp = 0.0, A, R, cn, cell_length, cell_width;
  double *qptr,*qoldptr,*qSGoldptr;
  double Q_multiplier=1.0;

  // CCS Multiplier for Q boundaries. If using regular grid, Qs are specified as m^2 and need to be multiplied by dx;
  // if using lat-long Qs are specified in m^3 and therefore multiplier is 1. Note intialised above as 1.0.
  if(Statesptr->latlong==OFF) Q_multiplier=Parptr->dx;

  g=Solverptr->g;
  numBCs=2*Parptr->xsz+2*Parptr->ysz;
  for(BCi=0;BCi<numBCs;BCi++)
  {
    if (BCptr->BC_Ident[BCi]!=0) // if BCi = 0 do nothing
    {
      // First for each edge number work out where it is on the boundary,
      // the associated edge pixels and whether it's facing in the x or y
      // direction
      if(BCi<=Parptr->xsz-1)
      {   // N(j=0) edge
        p0=BCi;
        p1=BCi+Parptr->xsz;
        pTQ=BCi;
        qptr=Arrptr->Qy+BCi;
        qoldptr=Arrptr->Qyold+BCi;
        qSGoldptr=Arrptr->QySGold+BCi;
        sign=-1; dir=1; edge=1;
      }
      else if(BCi>=Parptr->xsz && BCi<=Parptr->xsz+Parptr->ysz-1)
      {  // E edge
        p0=Parptr->xsz-1+(BCi-Parptr->xsz)*Parptr->xsz;
        p1=Parptr->xsz-2+(BCi-Parptr->xsz)*Parptr->xsz;
        pTQ=Parptr->xsz+(BCi-Parptr->xsz)*(Parptr->xsz+1);
        qptr=Arrptr->Qx+Parptr->xsz+(BCi-Parptr->xsz)*(Parptr->xsz+1);
        qoldptr=Arrptr->Qxold+Parptr->xsz+(BCi-Parptr->xsz)*(Parptr->xsz+1);
        qSGoldptr=Arrptr->QxSGold+Parptr->xsz+(BCi-Parptr->xsz)*(Parptr->xsz+1);
        sign=1; dir=2; edge=2;
      }
      else if(BCi>=Parptr->xsz+Parptr->ysz && BCi<=2*Parptr->xsz+Parptr->ysz-1)
      {  // S(j=ysz-1) edge
        p0=2*Parptr->xsz+Parptr->ysz-1-BCi+(Parptr->ysz-1)*Parptr->xsz;
        p1=2*Parptr->xsz+Parptr->ysz-1-BCi+(Parptr->ysz-2)*Parptr->xsz;
        pTQ=2*Parptr->xsz+Parptr->ysz-1-BCi+(Parptr->ysz)*(Parptr->xsz+1);
        qptr=Arrptr->Qy+2*Parptr->xsz+Parptr->ysz-1-BCi+(Parptr->ysz)*(Parptr->xsz+1);
        qoldptr=Arrptr->Qyold+2*Parptr->xsz+Parptr->ysz-1-BCi+(Parptr->ysz)*(Parptr->xsz+1);
        qSGoldptr=Arrptr->QySGold+2*Parptr->xsz+Parptr->ysz-1-BCi+(Parptr->ysz)*(Parptr->xsz+1);
        sign=1; dir=1; edge=3;
      }
      else
      {   // W edge
        p0=0+(numBCs-1-BCi)*Parptr->xsz;
        p1=1+(numBCs-1-BCi)*Parptr->xsz;
        pTQ=(numBCs-1-BCi)*(Parptr->xsz+1);
        qptr=Arrptr->Qx+(numBCs-1-BCi)*(Parptr->xsz+1);
        qoldptr=Arrptr->Qxold+(numBCs-1-BCi)*(Parptr->xsz+1);
        qSGoldptr=Arrptr->QxSGold+(numBCs-1-BCi)*(Parptr->xsz+1);
        sign=-1; dir=2; edge=4;
      }

      // CCS Record cell length and cell width relative to direction of channel (for lat-long grids)
      if(edge==1 || edge==3) {
        cell_length=Arrptr->dy[p0];
        cell_width=Arrptr->dx[p0]; // N or S boundary; assume flow is N-S or S-N
      }
      if(edge==2 || edge==4) {
        cell_length=Arrptr->dx[p0];
        cell_width=Arrptr->dy[p0]; // E or W boundary; assume flow is E-W or W-E
      }

      // Now calculate flows
      if(BCptr->BC_Ident[BCi]==1) // FREE boundary
      {
        if(Arrptr->H[p0]>Solverptr->DepthThresh)
        {
          if(Arrptr->SGCManningsn!=NULL)   cn=Arrptr->SGCManningsn[p0]*Arrptr->SGCManningsn[p0]; // this could be pre-calculated
          else cn=SGCptr->SGCn[Arrptr->SGCgroup[p0]];
          if(Arrptr->Manningsn!=NULL) fn=Arrptr->Manningsn[p0];
          else fn=Parptr->FPn;
          // calcQ
          CalcSGC_pointFREE(Arrptr->H[p0], Arrptr->SGCwidth[p0], BCptr->BC_Val[BCi], Solverptr->DepthThresh, Solverptr->Tstep, cell_width, g, cn, fn, Arrptr->SGCbfH[p0], Arrptr->SGCgroup[p0], sign, qoldptr, qptr, qSGoldptr, SGCptr);
          /*// use sub grid implmentation if there is a channel presents
            hflow = Arrptr->H[p0];
            w0 = Arrptr->SGCwidth[p0];
            Sf = BCptr->BC_Val[BCi];

            if(w0>0.0) // channel present
            {
            gr = Arrptr->SGCgroup[p0];
            // channel friction
            if(Arrptr->SGCManningsn!=NULL)   cn=Arrptr->SGCManningsn[p0]*Arrptr->SGCManningsn[p0]; // this could be pre-calculated
            else cn=SGCptr->SGCn[gr];
            CalcSGC_A(gr, hflow, Arrptr->SGCbfH[p0], &A, &w0, SGCptr); // calculate channel area for SGC
            R = CalcSGC_R(gr, hflow, Arrptr->SGCbfH[p0], w0, Arrptr->SGCwidth[p0], A, SGCptr); // calculate hydraulic radius for SGC
            // Calculate channel flow using inertial wave equation --- in this case Qxold will be in m3s-1 not m2s-1
            *qSGoldptr =  sign* (fabs(*qSGoldptr)+fabs(g*Solverptr->Tstep*A*Sf)) / (1+Solverptr->Tstep*g*cn*fabs(*qSGoldptr) / (pow(R,4.0/3.0)*A) );
            hflow = getmax(hflow-Arrptr->SGCbfH[p0], 0.0); // reduce hflow to account for channel
            }
            else *qSGoldptr = 0.0;
            // multiply flux by -sign and use absolute value of q0 to get flux directions correctly assigned at boundaries
            // fabs on Sf and q0 always results in positive or no flow... sign then sorts out the direction(jcn)
            //if(hflow>Solverptr->DepthThresh && w0 < Parptr->dx) // only calculate floodpplain flow if the depth is above bank and the channel is narrower than a cell //CCS_deletion
            if(hflow>Solverptr->DepthThresh && w0 < cell_width) // only calculate floodpplain flow if the depth is above bank and the channel is narrower than a cell
            {
            // floodplain friction
            if(Arrptr->Manningsn!=NULL) fn=Arrptr->Manningsn[p0];
            else fn=Parptr->FPn;
            // calculate FP flow
            *qoldptr=sign*(fabs(*qoldptr)+fabs(g*Solverptr->Tstep*hflow*Sf))/ (1+g*Solverptr->Tstep*fn*fn*fabs(*qoldptr)/(pow(hflow,(7./3.))));
            }
            else *qoldptr=0.0;
            // *qptr= (*qSGoldptr) + (*qoldptr) * (Parptr->dx-w0); //Combine SGC and floodplain flows //CCS_deletion
            *qptr= (*qSGoldptr) + (*qoldptr) * (cell_width-w0);*/
        }
        else
        {
          *qptr=0.0;
          *qoldptr=0.0;
          *qSGoldptr = 0.0;
        }
      }
      // HFIX & HVAR boundary
      if(BCptr->BC_Ident[BCi]==2 || BCptr->BC_Ident[BCi]==3)
      {
        if (Arrptr->H[p0]>Solverptr->DepthThresh)
        {
          if (BCptr->BC_Ident[BCi]==2) h0=BCptr->BC_Val[BCi];   // boundary depth for HFIX
          else h0=InterpBC(BCptr->BCVarlist[(int)BCptr->BC_Val[BCi]],Solverptr->t); // boundary depth for HVAR
          h1=Arrptr->H[p0];        // cell depth
          z1=Arrptr->DEM[p0];      // FP elevation
          w0=Arrptr->SGCwidth[p0]; // SGC width

          if(w0>0.0) //  check for sub-grid channel
          {
            // channel present
            zb0=Arrptr->SGCz[p0]; // channel bed
            hflow=getmax(h1,h0-zb0); // use max of cell depth and boundary depth
            dh=h0-(zb0+h1);
            gr = Arrptr->SGCgroup[p0];
            //if (edge == 1 || edge == 4) Sf=-dh/Parptr->dx; //CCS_deletion
            //else Sf=dh/Parptr->dx; //CCS_deletion
            Sf=dh/(cell_length*SGCptr->SGCm[gr]);
            if (edge == 1 || edge == 4) Sf=-Sf;

            // channel friction
            if(Arrptr->SGCManningsn!=NULL)   cn=Arrptr->SGCManningsn[p0]*Arrptr->SGCManningsn[p0]; // this could be pre-calculated
            else cn=SGCptr->SGCn[gr]; // FP channel... note this one is already squared
            CalcSGC_A(gr, hflow, Arrptr->SGCbfH[p0], &A, &w0, SGCptr); // calculate channel area for SGC
            R = CalcSGC_R(gr, hflow, Arrptr->SGCbfH[p0], w0, Arrptr->SGCwidth[p0], A, SGCptr); // calculate hydraulic radius for SGC
            // Calculate channel flow using inertial wave equation --- in this case Qxold will be in m3s-1 not m2s-1
            *qSGoldptr           =  ((*qSGoldptr)-g*Solverptr->Tstep*A*Sf) / (1+Solverptr->Tstep*g*cn*fabs(*qSGoldptr) / (pow(R,4.0/3.0)*A) );
            hflow = getmax(hflow-Arrptr->SGCbfH[p0], 0.0); // reduce hflow to account for channel
          }
          else hflow=getmax(h1,h0-z1);
          // multiply flux by -sign and use absolute value of q0 to get flux directions correctly assigned at boundaries
          // fabs on Sf and q0 always results in positive or no flow... sign then sorts out the direction(jcn)
          //if(hflow>Solverptr->DepthThresh && w0 < Parptr->dx) //CCS_deletion
          if(hflow>Solverptr->DepthThresh && w0 < cell_width)
          {
            // floodplain friction
            if(Arrptr->Manningsn!=NULL) fn=Arrptr->Manningsn[p0];
            else fn=Parptr->FPn; // FP Manning's
            dh=h0-z1-h1;
            //if (edge == 1 || edge == 4) Sf=-dh/Parptr->dx; //CCS_deletion
            //else Sf=dh/Parptr->dx; //CCS_deletion
            Sf=dh/cell_length;
            if (edge == 1 || edge == 4) Sf=-Sf;
            // calculate FP flow
            *qoldptr=((*qoldptr)-g*Solverptr->Tstep*hflow*Sf)/(1+g*Solverptr->Tstep*fn*fn*fabs(*qoldptr)/(pow(hflow,(7./3.))));
            //Qfp = (*qoldptr) * (Parptr->dx-w0);// CCS_deletion
            Qfp = (*qoldptr) * (cell_width-w0);
          }
          else *qoldptr=0.0;
          *qptr= (*qSGoldptr) +  Qfp;//Combine SGC and floodplain flows
        }
        else
        {
          *qptr=0.0;
          *qoldptr=0.0;
          *qSGoldptr = 0.0;
        }
      }
      // QFIX boundary
      if(BCptr->BC_Ident[BCi]==4)
      {
        //*qptr=-sign*BCptr->BC_Val[BCi]*Parptr->dx; //CCS_deletion
        *qptr=-sign*BCptr->BC_Val[BCi]*Q_multiplier;
      }
      // QVAR boundary
      if(BCptr->BC_Ident[BCi]==5)
      {
        //*qptr=-sign*InterpBC(BCptr->BCVarlist[(int)BCptr->BC_Val[BCi]],Solverptr->t)*Parptr->dx; //CCS_deletion
        *qptr=-sign*InterpBC(BCptr->BCVarlist[(int)BCptr->BC_Val[BCi]],Solverptr->t)*Q_multiplier;

      }
    }
  }
  return;
}
//-----------------------------------------------------------------------------------
// Calculates q for starting SGC model from water surfaces
void SGC_hotstart(States *Statesptr,Pars *Parptr,Solver *Solverptr,Arrays *Arrptr)
{
  int p0,p1,pq0,i,j;
  double z0,z1,zb0,zb1,h0,h1,fn,w0,w1,dh,Sf,hflow;

  // Calc Q in x:
#pragma omp parallel for private(i,p0,p1,pq0,z0,z1,zb0,zb1,h0,h1,fn,w0,w1,dh,Sf,hflow)
  for(j=0;j<Parptr->ysz;j++)
  {
    for(i=0;i<Parptr->xsz-1;i++)
    {
      p0=i+j*Parptr->xsz;
      p1=i+1+j*Parptr->xsz;
      pq0=i+j*(Parptr->xsz+1)+1;

      z0=Arrptr->DEM[p0];
      z1=Arrptr->DEM[p1];
      zb0=Arrptr->SGCz[p0];
      zb1=Arrptr->SGCz[p1];
      h0=Arrptr->H[p0];
      h1=Arrptr->H[p1];
      w0=Arrptr->SGCwidth[p0];
      w1=Arrptr->SGCwidth[p1];
      if(Arrptr->Manningsn!=NULL) fn=0.5*(Arrptr->Manningsn[p0]+Arrptr->Manningsn[p1]);
      else fn=Parptr->FPn;
      //////////////////SG flows ////////////////////////////////
      // if both cells are a SGC
      if (w0>0 && w1>0)
      {
        if(zb0+h0>zb1+h1 && h0>Solverptr->DepthThresh) // Flow from 0->1
        {
          dh=zb0+h0-zb1-h1;
          //Sf=sqrt(dh/Parptr->dx); //CCS_deletion
          Sf=sqrt(dh/Arrptr->dx[p0]);
          hflow=getmax(zb0+h0,zb1+h1)-getmax(zb0,zb1);

          if(hflow>Solverptr->DepthThresh)
          {
            Arrptr->QxSGold[pq0]=(pow(hflow,(5.0/3.0))*Sf/fn)*getmin(w0,w1);
          }
        }
        else if(zb0+h0<zb1+h1 && h1>Solverptr->DepthThresh)  // Flow from 1->0
        {
          dh=zb1+h1-zb0-h0;
          //Sf=sqrt(dh/Parptr->dx); //CCS_deletion
          Sf=sqrt(dh/Arrptr->dx[p0]);
          hflow=getmax(zb0+h0,zb1+h1)-getmax(zb0,zb1);
          if(hflow>Solverptr->DepthThresh)
          {
            Arrptr->QxSGold[pq0]=(-pow(hflow,(5.0/3.0))*Sf/fn);
          }
        }
      }
      //////////////////FP flows ////////////////////////////////
      // Correct h for any SGC depth
      if (w0>0) h0 = getmax(h0-z0+zb0, 0.0);
      if (w1>0) h1 = getmax(h1-z1+zb1, 0.0);
      // calculate qfp
      if(z0+h0>z1+h1 && h0>Solverptr->DepthThresh) // Flow from 0->1
      {
        dh=z0+h0-z1-h1;
        //Sf=sqrt(dh/Parptr->dx); //CCS_deletion
        Sf=sqrt(dh/Arrptr->dx[p0]);
        hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
        if(hflow>Solverptr->DepthThresh)
        {
          Arrptr->Qxold[pq0]=(pow(hflow,(5.0/3.0))*Sf/fn);
        }
      }
      else if(z0+h0<z1+h1 && h1>Solverptr->DepthThresh)  // Flow from 1->0
      {
        dh=z1+h1-z0-h0;
        //Sf=sqrt(dh/Parptr->dx); //CCS_deletion
        Sf=sqrt(dh/Arrptr->dx[p0]);
        hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
        if(hflow>Solverptr->DepthThresh)
        {
          Arrptr->Qxold[pq0]=(-pow(hflow,(5.0/3.0))*Sf/fn);
        }
      }
    }
  }

  // Calc Q in y:
#pragma omp parallel for private(i,p0,p1,pq0, z0,z1,zb0,zb1,h0,h1,fn,w0,w1,dh,Sf,hflow)
  for(j=0;j<Parptr->ysz-1;j++)
  {
    for(i=0;i<Parptr->xsz;i++)
    {
      p0=i+j*Parptr->xsz;
      p1=i+(j+1)*Parptr->xsz;
      pq0=i+(j+1)*(Parptr->xsz+1);

      z0=Arrptr->DEM[p0];
      z1=Arrptr->DEM[p1];
      w0=Arrptr->SGCwidth[p0];
      w1=Arrptr->SGCwidth[p1];
      zb0=Arrptr->SGCz[p0];
      zb1=Arrptr->SGCz[p1];
      h0=Arrptr->H[p0];
      h1=Arrptr->H[p1];
      if(Arrptr->Manningsn!=NULL) fn=0.5*(Arrptr->Manningsn[p0]+Arrptr->Manningsn[p1]);
      else fn=Parptr->FPn;
      ////////////////SG flows ////////////////////////////////
      // if both cells are a SGC
      // if both cells are a SGC
      if (w0>0 && w1>0)
      {
        if(zb0+h0>zb1+h1 && h0>Solverptr->DepthThresh) // Flow from 0->1
        {
          dh=zb0+h0-zb1-h1;
          //Sf=sqrt(dh/Parptr->dx); //CCS_deletion
          Sf=sqrt(dh/Arrptr->dy[p0]);
          hflow=getmax(zb0+h0,zb1+h1)-getmax(zb0,zb1);
          if(hflow>Solverptr->DepthThresh)
          {
            Arrptr->QySGold[pq0]=(pow(hflow,(5.0/3.0))*Sf/fn)*getmin(w0,w1);
          }
        }
        else if(zb0+h0<zb1+h1 && h1>Solverptr->DepthThresh)  // Flow from 1->0
        {
          dh=zb1+h1-zb0-h0;
          //Sf=sqrt(dh/Parptr->dx); //CCS_deletion
          Sf=sqrt(dh/Arrptr->dy[p0]);
          hflow=getmax(zb0+h0,zb1+h1)-getmax(zb0,zb1);
          if(hflow>Solverptr->DepthThresh)
          {
            Arrptr->QySGold[pq0]=(-pow(hflow,(5.0/3.0))*Sf/fn);
          }
        }
      }
      //////////////////FP flows ////////////////////////////////
      // Correct h for any SGC depth
      if (w0>0) h0 = getmax(h0-z0+zb0, 0.0);
      if (w1>0) h1 = getmax(h1-z1+zb1, 0.0);
      // calculate qfp
      if(z0+h0>z1+h1 && h0>Solverptr->DepthThresh) // Flow from 0->1
      {
        dh=z0+h0-z1-h1;
        //Sf=sqrt(dh/Parptr->dx); //CCS_deletion
        Sf=sqrt(dh/Arrptr->dy[p0]);
        hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
        if(hflow>Solverptr->DepthThresh)
        {
          Arrptr->Qyold[pq0]=(pow(hflow,(5.0/3.0))*Sf/fn);
        }
      }
      else if(z0+h0<z1+h1 && h1>Solverptr->DepthThresh)  // Flow from 1->0
      {
        dh=z1+h1-z0-h0;
        //Sf=sqrt(dh/Parptr->dx); //CCS_deletion
        Sf=sqrt(dh/Arrptr->dy[p0]);
        hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
        if(hflow>Solverptr->DepthThresh)
        {
          Arrptr->Qyold[pq0]=(-pow(hflow,(5.0/3.0))*Sf/fn);
        }
      }
    }
  }
  return;
}
//-----------------------------------------------------------------------------
// FLOODPLAIN EVAPORATION
// with correction for sub grid channels
void SGC_Evaporation(Pars *Parptr, Solver *Solverptr, Arrays *Arrptr,SGCprams *SGCptr)
{
  int i,j, p0, gr;
  double bfh, V, cell_evap,h0,evap_rate,h_old, loc_evap_loss=0.0;

  evap_rate=InterpBC(Arrptr->evap,Solverptr->t); //constant rate across whole floodplain
  // calculate evaporation depth
  cell_evap = evap_rate * Solverptr->Tstep; //rate for depth, not area

  // Calculate Evaporation
#pragma omp parallel for private( i, gr, h0, p0, h_old, bfh, V) reduction (+:loc_evap_loss)
  for(j=0;j<Parptr->ysz;j++)
  {
    for(i=0;i<Parptr->xsz;i++)
    {
      // location of cell
      p0  = i+j*Parptr->xsz;
      // depth in cell
      h0  = Arrptr->H[p0];
      if(h0>Solverptr->DepthThresh) // There is water to evaporate
      {
        // retain old depth for mass balance calc
        h_old = h0;
        // update depth by subtracting evap depth
        h0-=cell_evap;
        //check for -ve depths
        if(h0<0.0)
        {
          // reduce evap loss to account for dry bed
          cell_evap+=h0;
          // make h zero
          h0=0;
        }
        // update h (*CCS - is this redundant?)
        Arrptr->H[p0]=h0;
        // get bankfull depth
        bfh = Arrptr->SGCbfH[p0];
        // calculate mass loss in rectangular channel or floodplain including treatment for bank transitions
        // if there is a sub-grid channel and h0 is below bankfull
        if (Arrptr->SGCwidth[p0] > 0.0 && h0 < bfh)
        {
          // sub-grid channel evap or transition evap
          //if (h_old < Arrptr->SGCbfH[p0] || Arrptr->SGCwidth[p0] > Parptr->dx) //CCS_deletion
          if (h_old < Arrptr->SGCbfH[p0] || Arrptr->SGCwidth[p0] > 0.5*(Arrptr->dx[p0]+Arrptr->dy[p0]))
          {
            // calculate loss in vol
            gr = Arrptr->SGCgroup[p0]; // channel group number
            V =  CalcSGC_UpV(SGCptr->SGCchantype[gr], h_old, SGCptr->SGCs[gr], Arrptr->SGCc[p0]); //Calculate channel volume
            V -= CalcSGC_UpV(SGCptr->SGCchantype[gr], h0,    SGCptr->SGCs[gr], Arrptr->SGCc[p0]); //Calculate channel volume
            Arrptr->SGCdVol[p0]-=V;
            loc_evap_loss+=V; // accounts for SGC channel loss
          }
          // old water level must be above bank height and the channel is smaller than a cell width
          // but the new water level is below bank height, evap mass loss for bank transition
          else
          {
            // mass lost from floodplain
            cell_evap = h_old - bfh;
            // mass lost from channel
            gr = Arrptr->SGCgroup[p0]; // channel group number
            V =     CalcSGC_UpV(SGCptr->SGCchantype[gr], bfh, SGCptr->SGCs[gr], Arrptr->SGCc[p0]); //Calculate bankfull area
            V -= CalcSGC_UpV(SGCptr->SGCchantype[gr], h0,  SGCptr->SGCs[gr], Arrptr->SGCc[p0]); //Calculate channel area
            // mass loss from floodplain
            //V+= (cell_evap*Parptr->dA); // add mass loss from floodplain to that of channel //CCS_deletion
            V+= (cell_evap*Arrptr->dA[p0]);
            Arrptr->SGCdVol[p0]-=V; // remove volume from cell
            loc_evap_loss+=V; //mass-balance for a standard cell
          }
        }
        // standard evap case of non sub-grid channel or overbank level
        else
        {
          //V = (cell_evap*Parptr->dA); //CCS_deletion
          V = (cell_evap*Arrptr->dA[p0]);
          Arrptr->SGCdVol[p0]-=V; // remove volume from cell
          loc_evap_loss+=V; //mass-balance for a standard cell
        }
      }
    }
  }
  // add time-step evap loss to the total evap loss
  Parptr->EvapTotalLoss+=loc_evap_loss;
  return;
}
//-----------------------------------------------------------------------------------
// SGC RAINFALL, CCS May 2013
// Rainfall routine WITHOUT routing scheme.  If routing is enabled, rainfall is added by SGC_Routing function.
void SGC_Rainfall(Pars *Parptr, Solver *Solverptr, Arrays *Arrptr)
{
  int i,j, p0;
  double rain_rate,V,loc_rainfall_total=0.0;

  rain_rate=InterpBC(Arrptr->rain,Solverptr->t);//constant rate across whole floodplain
  //V=(rain_rate*Solverptr->Tstep*Parptr->dA); // rainfall volume = (rainfall rate * timestep * surface area) //CCS_deletion

  // Add rainfall volume to cell volumes:
#pragma omp parallel for private(i, p0, V) reduction (+:loc_rainfall_total)
  for(j=0;j<Parptr->ysz;j++)
  {
    for(i=0;i<Parptr->xsz;i++)
    {
      p0=i+j*Parptr->xsz; // location of cell
      V=(rain_rate*Solverptr->Tstep*Arrptr->dA[p0]); // now calc inside pragma loop to handle variable dA lat-long grids

      if(Arrptr->DEM[p0]!=1e10) // if not nodata cell
      {
        Arrptr->SGCdVol[p0]+=V; // add rainfall volume to cell
        loc_rainfall_total+=V; // mass balance for local cell (cumulative)
        //Arrptr->H[p0]+=rain_rate*Solverptr->Tstep;
      }
    }
  }
  Parptr->RainTotalLoss+=loc_rainfall_total; // Update domain mass balance
  return;
}
//-----------------------------------------------------------------------------------
// SGC Routing Scheme, CCS May 2013
// Routing Scheme.  Also adds rainfall when rainfall is enabled, but can be used without rainfall to force stability on steep domains.
void SGC_Routing(States *Statesptr, Pars *Parptr, Solver *Solverptr, Arrays *Arrptr)
{
  int i,j,p0,p1;
  double cell_rain,rain_rate,h0,h1,z0,z1,dh,Sfy,Sfx,V,flow,flowV,flow_fraction;

  if(Statesptr->rainfall==ON) // Calc rainfall H and dVol for this timestep:
  {
    rain_rate=InterpBC(Arrptr->rain,Solverptr->t);//constant rate across whole floodplain
    cell_rain = rain_rate * Solverptr->Tstep; //rate for depth, not area
    //V=(rain_rate*Solverptr->Tstep*Parptr->dA); // rainfall volume = (rainfall rate * timestep * surface area) //CCS_deletion
  }

  // Main routing loop:
  for(j=0;j<Parptr->ysz;j++)
  {
    for(i=0;i<Parptr->xsz;i++)
    {
      p0=i+j*Parptr->xsz;

      if(Arrptr->DEM[p0]!=1e10)
      {
        if(Statesptr->rainfall==ON) // If rainfall is enabled, add rainfall:
        {
          V=(rain_rate*Solverptr->Tstep*Arrptr->dA[p0]);
          Arrptr->H[p0]+=cell_rain; // Update H array for cell rainfall
          Arrptr->SGCdVol[p0]+=V; // add rainfall volume to cell
          Parptr->RainTotalLoss+=V; //for mass-balance
        }

        if(Arrptr->H[p0]>0) // Only proceed with rest of loop if cell is wet:
        {
          h0=Arrptr->H[p0];
          z0=Arrptr->DEM[p0];

          //Calculate friction slopes to assess whether routing scheme is needed to force stability:
          if(i==Parptr->xsz-1) Sfx=0.0; // boundary cell, don't calc Sfx:
          else
          {
            p1=i+1+j*Parptr->xsz;
            z1=Arrptr->DEM[p1];
            h1=Arrptr->H[p1];
            dh=z0+h0-z1-h1;
            //Sfx=fabs(dh/Parptr->dx); // friction slope in x //CCS_deletion
            Sfx=fabs(dh/Arrptr->dx[p0]); // friction slope in x
          }

          if(j==Parptr->ysz-1) Sfy=0.0; // boundary cell, don't calc Sfy:
          else
          {
            p1=i+(j+1)*Parptr->xsz;
            z1=Arrptr->DEM[p1];
            h1=Arrptr->H[p1];
            dh=z0+h0-z1-h1;
            //Sfy=fabs(dh/Parptr->dx); // friction slope in y //CCS_deletion
            Sfy=fabs(dh/Arrptr->dy[p0]); // friction slope in y
          }

          /*if(Statesptr->rainfall==ON) // If rainfall is enabled, add rainfall:
            {
            Arrptr->H[p0]+=cell_rain; // Update H array for cell rainfall
            Arrptr->SGCdVol[p0]+=V; // add rainfall volume to cell
            Parptr->RainTotalLoss+=V; //for mass-balance
            }*/

          if(Arrptr->SGCwidth[p0]==0.0) // if not subgrid channel cell
          {
            h0=Arrptr->H[p0];
            if(h0<Solverptr->DepthThresh || Sfy>=Parptr->RouteSfThresh || Sfx>=Parptr->RouteSfThresh)
            {
              h1=Arrptr->H[Arrptr->FlowDir[p0]]-Arrptr->SGCbfH[Arrptr->FlowDir[p0]];
              if(h1<0) {h1=0.0;}  // If h1 negative due to below bankful SG channel cell, set h1 to zero
              z0=Arrptr->DEM[p0]; //cell DEM height
              z1=Arrptr->DEM[Arrptr->FlowDir[p0]]; //lowest neighbour cell DEM height

              flow=(z0+h0)-(z1+h1);/*calculate the maximum possible flow into lowest neigbour cell by
                                     comparing water surface elevations:*/

              if(flow>h0) {flow=h0;} /*where water surface elevation of neighbour cell is below DEM
                                       level of current cell, set maxflow to h0*/
              if(flow<0) {flow=0;} /*where water surface elevation of neighbour cell is above water
                                     surface elevation of current cell, set maxflow to 0*/

              flow_fraction=Solverptr->Tstep/Arrptr->RouteInt[p0]; // fraction of cell volume to route in this time step
              if(flow_fraction>1) flow_fraction=1; // prevent flow fraction>1 (should never happen!)

              //flowV=flow*Parptr->dA*flow_fraction; // Calc flow volume (flow depth * area * flow fraction) //CCS_deletion
              flowV=flow*Arrptr->dA[p0]*flow_fraction; // Calc flow volume (flow depth * area * flow fraction)

              Arrptr->SGCdVol[p0]-=flowV; // Remove flow volume from source cell
              Arrptr->SGCdVol[Arrptr->FlowDir[p0]]+=flowV; // Add flow volume to recipient cell

            }
          }
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------------
// This function calculates the area for a single side of a sub-grid channel interface given the channel type
void CalcSGC_A(int gr, double hflow, double bf, double *A, double *we, SGCprams *SGCptr)
{
  // This function calculates the area of flow (A) for a given flow depth (hflow), in some cases it
  // also returnes the widths of flow (We).
  double sl;
  // switch depending on the channel type
  switch (SGCptr->SGCchantype[gr])
  {
  case 1: // Rectangular channel (default) - This model has a top width and bed elevation and top
    (*A) = (*we)*hflow;
    break;  // Break terminates the switch statement

  case 2: // Power shaped channel. h = x^sl channel
    sl = SGCptr->SGCs[gr];
    if (hflow < bf)
    {
      (*we) = (*we)*pow(hflow/bf,1.0/sl);
      (*A)  = hflow*(*we)*(1.0-1.0/(sl+1.0));
    }
    else (*A)  = bf*(*we)*(1-1/(sl+1)) + (*we)*(hflow-bf); // out of bank flow area
    break;  // Break terminates the switch statement

  case 3: // linear slope.
    if (hflow < bf)
    {
      (*we) = ((*we)/bf)*hflow;
      (*A)  = (*we)*hflow*0.5;    // within bank flow area
    }
    else  (*A)  = (*we)*bf*0.5 + (*we)*(hflow-bf); // out of bank flow area
    break;  // Break terminates the switch statement

  case 4: // triangular channel
    if (hflow < bf)
    {
      (*we) = ((*we)/bf)*hflow;
      (*A)  = (*we)*hflow*0.5;    // within bank flow area
    }
    else  (*A)  = (*we)*bf*0.5   + (*we)*(hflow-bf); // out of bank flow area
    break;  // Break terminates the switch statement

  case 5: // parabolic channel
    if (hflow < bf)
    {
      (*we) = (*we)*sqrt(hflow/bf);
      (*A)  = hflow*(*we)*(2.0/3.0);
    }
    else    (*A)  = bf*   (*we)*(2.0/3.0) + (*we)*(hflow-bf); // out of bank flow area
    break;  // Break terminates the switch statement

  case 6: // Rectangular channel (default) - This model has a top width and bed elevation and top
    (*A) = (*we)*hflow;
    break;  // Break terminates the switch statement

  case 7: // trapazoidal channel
    sl = SGCptr->SGCs[gr];
    if (hflow < bf) (*A)  = ((*we)+sl*hflow)*hflow;    // within bank flow area
    else            (*A)  = ((*we)+sl*bf   )*bf   + ((*we)+sl*bf)*(hflow-bf);
    break;  // Break terminates the switch statement

  default: // its all gone wrong
    printf("should not be here! Something is wrong with the SGC channel model A calculation");
    break;
  }
  return;
}
//-----------------------------------------------------------------------------------
// This function calculates the hydraulic radius of a sub-grid channel given the channel type
double CalcSGC_R(int gr, double h, double hbf, double w, double wbf, double A, SGCprams *SGCptr)
{
  // This function calculates the hydraulic redius of the channel given the flow area (A) and depth (h)
  double R=0.0, sl, hp, b1, b2;

  switch (SGCptr->SGCchantype[gr])
  {
  case 1: // Rectangular channel (default) - This model has a top width and bed elevation and top
    R = getmin(h, hbf); // don't exceed bankfull
    R = A / (w+2*R); // calculate hydraulic radius
    break; // Break terminates the switch statement

  case 2: // Exponent channel
    h = getmin(h, hbf); // don't exceed bankfull
    hp = 2*h/wbf; // use half width
    // get beta parameters for the channel group number
    if (hp <= SGCptr->SGCbetahmin)
    {
      // beta parameters for flow depths below SGCbetahmin (0.05) depth/bankfull width
      b1 = SGCptr->SGCbeta1[gr]; b2 = SGCptr->SGCbeta2[gr];
      R  = b1*hp + b2*hp*hp; // calculate wetted perimiter component
    }
    else
    {
      // beta parameters for flow depths above or equal to SGCbetahmin (0.05) depth/bankfull width
      b1 = SGCptr->SGCbeta3[gr]; b2 = SGCptr->SGCbeta4[gr];
      hp -= SGCptr->SGCbetahmin; // decrease hp to account for for wetted perimiter fraction below 0.05 depth/bankful depth
      R  = SGCptr->SGCbeta5[gr] + b1*hp + b2*hp*hp; // calculate wetted perimiter component
    }
    R = A / (w + R*wbf);  // calculate hydraulic radius
    break;  // Break terminates the switch statement

  case 3: // linear slope - This model has a bed elevation, slope, top width and top elevation.
    if (h < hbf)    R = A / (h   + sqrt(h  *h  +w*w));      // within bank flow hydraulic radius
    else                    R = A / (hbf + sqrt(hbf*hbf+w*w));   // out of bank hydraulic radius (wetted perimeter is actually constant !!)
    break;  // Break terminates the switch statement

  case 4: // triangular channel - This model has the bed elevation, slope, top width and top elevation
    w = 0.5*w;
    if (h < hbf)    R = A / (2*sqrt(h  *h  +w*w));   // within bank flow hydraulic radius
    else                    R = A / (2*sqrt(hbf*hbf+w*w));   // out of bank hydraulic radius (wetted perimeter is actually constant !!)
    break;  // Break terminates the switch statement

  case 5: // parabolic channel
    h = getmin(h, hbf); // don't exceed bankfull
    w = w/2.0; // half width
    R = sqrt(w*w + 16.0*h*h);
    R = 0.5*R + (w*w)/(8.0*h) * log((4.0*h+R)/w);
    R = 2.0*R;
    R = A / R;  // calculate hydraulic radius  */
    break;  // Break terminates the switch statement

  case 6: // Rectangular channel (no banks) - This model has a top width and bed elevation and top
    R = A / w; // calculate hydraulic radius
    break; // Break terminates the switch statement

  case 7: // trapazoidal channel
    sl = SGCptr->SGCs[gr];
    if (h < hbf)    R = A / (w + 2* h   * sqrt(1+sl*sl));   // within bank flow hydraulic radius
    else                    R = A / (w + 2* hbf * sqrt(1+sl*sl));   // out of bank hydraulic radius (wetted perimeter is actually constant !!)
    break;  // Break terminates the switch statement

  default: // its all gone wrong
    printf("should not be here! Something is wrong with the SGC channel model R calculation");
    break;
  }
  return(R);
}
/*

  case 7: // trapazoidal channel
  // calculate hydraulic radius
  if (hflow < bf)       R = A / (we + 2* hflow * sqrt(1+sl*sl));   // within bank flow hydraulic radius
  else                  R = A / (we + 2* bf    * sqrt(1+sl*sl));   // out of bank hydraulic radius (wetted perimeter is actually constant !!)
  break;        // Break terminates the switch statement

*/

//-----------------------------------------------------------------------------------
// This function updates H for a sub-grid channel given the channel type
double CalcSGC_UpH(int SGCchan_type, double V, double sl, double c)
{
  double h=0.0;

  // switch to the correct sub-grid channel, default 1 is the rectangular
  switch (SGCchan_type)
  {
  case 1: // Rectangular channel (default) - This model has a top width and bed elevation and top
    h = V/c;
    break;

  case 2: // y = x^sl channel
    h = pow(V/c,sl/(sl+1.0));
    break;  // Break terminates the switch statement

  case 3: // Rectangular channel (default) - This model has a top width and bed elevation and top
    h = V/c;
    break;

  case 4: // Rectangular channel (default) - This model has a top width and bed elevation and top
    h = V/c;
    break;

  case 5: // parabiolic channel
    h = pow(V/c,2.0/3.0);
    break;  // Break terminates the switch statement

  case 6: // Rectangular channel (no banks) - This model has a top width and bed elevation and top
    h = V/c;
    break;

  default: // its all gone wrong
    printf("should not be here! Something is wrong with the SGC channel model in SGC_UpdateH");
    break;
  } // end of switch statement
  return(h);
}
/*
  case 7: // trapazoidal channel

  // channel flow is within bank and the channel is not as wide as the cell
  if (*hptr < zbf ||  we > dx )
  {
  //trapazoidal channel calculate new h
  A = (dV/chanx) + (we+sl0*(*hptr))*(*hptr); //Calculate channel area
  (*hptr) = (sl0*sqrt((we*we+4*sl0*A)/(sl0*sl0))-we)/(2*sl0); // convert area to h
  ///// bank transitions //////
  // the channel must have started within bank if it has now gone overbank and width is less than the cell width correct h for overbank transition
  if (*hptr > zbf && we < dx)
  {
  // calculate bankfull area
  a = zbf*we +zbf*zbf*sl0;
  // Calculate bank full A then take this away from A, devide by dx then add to (z0-zb) to get h
  (*hptr) = zbf + (A - a)*chanx/dA;
  }
  }
  else // there is a flooded sub grid channel
  {
  (*hptr)+= dV/dA; // water level above sub grid channel use normal updateH
  ///// bank transitions //////
  // the channel must have started overbank if the water level has decreases below bank height correct for return within bank
  if (*hptr < zbf)
  {
  // bank height has been crossed... decreasing
  if ((*hptr)<0)  *hptr=0.0;
  else
  {
  // loss of area below bank full based Parptr->dA
  A = ((zbf-(*hptr))*dA)/chanx;
  // bank full A (a) minus loss of A (A)
  a = zbf*we + zbf*zbf*sl0;
  A = a - A;
  // calculate h given A and bank slope
  (*hptr) = (sl0*sqrt((we*we+4*sl0*A)/(sl0*sl0))-we)/(2*sl0);
  }
  }
  }
  break;        // Break terminates the switch statement
  }*/


double CalcSGC_UpV(int SGCchan_type, double h, double sl, double c)
{
  // This function calculates the volume of a sub-grid channel given a depth
  double v=0.0;
  // switch to the correct sub-grid channel, default 1 is the rectangular
  switch (SGCchan_type)
  {
  case 1:
    v = h*c;
    break;

  case 2:
    v = c*pow(h,1.0/sl+1.0);
    break;

  case 3:
    v = h*c;
    break;

  case 4:
    v = h*c;
    break;

  case 5:
    v = c*pow(h,3.0/2.0);
    break;

  case 6: // Rectangular channel (no banks)
    v = h*c;
    break;

  default: // its all gone wrong
    printf("should not be here! Something is wrong with the SGC channel model in CalcSGC_UpV");
    break;
  } // end of switch statement
  return(v);
}
void CalcSGC_pointFREE(double hflow, double w0, double Sf, double DT, double Tstep, double cell_width, double g, double cn, double fn, double SGCbfH, int gr, int sign, double *qoldptr, double *qptr, double *qSGoldptr, SGCprams *SGCptr)
{
  double A, R, wbfH;

  wbfH = w0; // save bank full width
  if(w0>0.0) // channel present
  {
    CalcSGC_A(gr, hflow, SGCbfH, &A, &w0, SGCptr); // calculate channel area for SGC
    R = CalcSGC_R(gr, hflow, SGCbfH, w0, wbfH, A, SGCptr); // calculate hydraulic radius for SGC
    // Calculate channel flow using inertial wave equation --- in this case Qxold will be in m3s-1 not m2s-1
    *qSGoldptr =  sign* (fabs(*qSGoldptr)+fabs(g*Tstep*A*Sf)) / (1+Tstep*g*cn*fabs(*qSGoldptr) / (pow(R,4.0/3.0)*A) );
    hflow = getmax(hflow-SGCbfH, 0.0); // reduce hflow to account for channel
  }
  else *qSGoldptr = 0.0;
  // multiply flux by -sign and use absolute value of q0 to get flux directions correctly assigned at boundaries
  // fabs on Sf and q0 always results in positive or no flow... sign then sorts out the direction(jcn)
  //if(hflow>Solverptr->DepthThresh && w0 < Parptr->dx) // only calculate floodpplain flow if the depth is above bank and the channel is narrower than a cell //CCS_deletion
  if(hflow>DT && w0 < cell_width) // only calculate floodpplain flow if the depth is above bank and the channel is narrower than a cell
  {
    // calculate FP flow
    *qoldptr=sign*(fabs(*qoldptr)+fabs(g*Tstep*hflow*Sf))/ (1+g*Tstep*fn*fn*fabs(*qoldptr)/(pow(hflow,(7./3.))));
  }
  else *qoldptr=0.0;
  //*qptr= (*qSGoldptr) + (*qoldptr) * (Parptr->dx-w0); //Combine SGC and floodplain flows //CCS_deletion
  *qptr= (*qSGoldptr) + (*qoldptr) * (cell_width-w0);

}
void SGC_wp_prams (SGCprams *SGCptr)
{
  int i;
  double s;
  /* Populate SGC gamma parmaters for power shaped channel (chantype 2)
     The equation to be used is beta = gamma + gamma*1/sl + gamma*(1/sl)^2 + gamma*(1/sl)^3 + gamma*sl + gamma*sl^2 + gamma*sl^3 + gamma*sl^0.5;
     or a simplification of
     This is repeated for beta's 1 to 4
     where beta 1 and 2 are for values less the
     Where b is beta that will be calulated and g is gamme the structure of SGCgamma is as follows
     SGCgamma = [b1g0, b1g1, b1g2, b1g3, b1g4, b2g0, b2g1...
  */
  SGCptr->SGCgamma       =new double[32]();
  // beta1
  SGCptr->SGCgamma[0] = 0.328523411998739;
  SGCptr->SGCgamma[1] = -0.136629598093381;
  SGCptr->SGCgamma[2] = 0.226782890677302;
  SGCptr->SGCgamma[3] = 0.407317315366405;
  SGCptr->SGCgamma[4] = 0.0573936111784235;
  SGCptr->SGCgamma[5] = -0.000827027398309836;
  SGCptr->SGCgamma[6] = 7.73363755705211e-06;
  SGCptr->SGCgamma[7] = -0.218593761377831;


  //SGCptr->SGCgamma[0] = -0.00382640781129463;
  //SGCptr->SGCgamma[1] = 0.192705773900624;
  //SGCptr->SGCgamma[2] = -0.0332165808200565;
  //SGCptr->SGCgamma[3] = 0.00253291121412471;
  //SGCptr->SGCgamma[4] = 0.0125545811978110;
  //SGCptr->SGCgamma[5] = -5.64900319095837e-05;
  // beta2
  SGCptr->SGCgamma[8] = -1.73499845293949;
  SGCptr->SGCgamma[9] = 1.86961902684339;
  SGCptr->SGCgamma[10] = -0.0787095731609086;
  SGCptr->SGCgamma[11] = -1.70280956920105;
  SGCptr->SGCgamma[12] = -0.236031863429058;
  SGCptr->SGCgamma[13] = 0.00114799994407101;
  SGCptr->SGCgamma[14] = 3.43885058357733e-06;
  SGCptr->SGCgamma[15] = 1.73575973862887;

  //SGCptr->SGCgamma[6] = 0.693473047865080;
  //SGCptr->SGCgamma[7] = 0.0998814928897974;
  //SGCptr->SGCgamma[8] = -0.112446326301890;
  //SGCptr->SGCgamma[9] = 0.0179336144901384;
  //SGCptr->SGCgamma[10] = 0.106947061485761;
  //SGCptr->SGCgamma[11] = -0.00634898682856126;
  //SGCptr->SGCgamma[12] = 0.000113534050231379;
  //SGCptr->SGCgamma[13] = 0.154047342257332;
  // beta3
  SGCptr->SGCgamma[16] = -0.0914107976436178;
  SGCptr->SGCgamma[17] = 0.0400444500138610;
  SGCptr->SGCgamma[18] = 0.387741337558434;
  SGCptr->SGCgamma[19] = -0.152454059221963;
  SGCptr->SGCgamma[20] = -0.0859683978597510;
  SGCptr->SGCgamma[21] = 0.00121037617690324;
  SGCptr->SGCgamma[22] = -1.24573143624190e-05;
  SGCptr->SGCgamma[23] = 0.515523126751616;

  //SGCptr->SGCgamma[14] = 0.0786876032500024;
  //SGCptr->SGCgamma[15] = 0.210468008584605;
  //SGCptr->SGCgamma[16] = -0.0562968512037699;
  //SGCptr->SGCgamma[17] = 0.00619136021801920;
  //SGCptr->SGCgamma[18] = -0.0552653878505992;
  //SGCptr->SGCgamma[19] = 0.000514162069550961;
  //SGCptr->SGCgamma[20] = -2.31911554664107e-06;
  //SGCptr->SGCgamma[21] = 0.383834975591764;
  // beta4
  SGCptr->SGCgamma[24] = 0.440010360908978;
  SGCptr->SGCgamma[25] = -0.291924093276058;
  SGCptr->SGCgamma[26] = 0.00781944059584046;
  SGCptr->SGCgamma[27] = 0.0164084553274571;
  SGCptr->SGCgamma[28] = 0.0337017710333279;
  SGCptr->SGCgamma[29] = -0.000464008615742078;
  SGCptr->SGCgamma[30] = 4.69370047950512e-06;
  SGCptr->SGCgamma[31] = -0.204804996998889;

  //SGCptr->SGCgamma[22] = 0.141643432781455;
  //SGCptr->SGCgamma[23] = -0.0722734103471909;
  //SGCptr->SGCgamma[24] = 0.0195083451137352;
  //SGCptr->SGCgamma[25] = -0.00206587357521111;
  //SGCptr->SGCgamma[26] = -0.0131507972423559;
  //SGCptr->SGCgamma[27] = 0.000641712101745405;
  //SGCptr->SGCgamma[28] = -1.23372890042230e-05;

  // now for each velaue os SGC1 calculate the beta parameters
  // create arrays for the number of model parameters
  SGCptr->SGCbeta1    =new double[SGCptr->NSGCprams]();
  SGCptr->SGCbeta2    =new double[SGCptr->NSGCprams]();
  SGCptr->SGCbeta3    =new double[SGCptr->NSGCprams]();
  SGCptr->SGCbeta4    =new double[SGCptr->NSGCprams]();
  SGCptr->SGCbeta5    =new double[SGCptr->NSGCprams]();

  // Calculate beta's for each value of SGC1 ... these will now be used by the wetted perimiter function of the radius function
  for(i=0;i<SGCptr->NSGCprams;i++)
  {
    s = SGCptr->SGCs[i];
    SGCptr->SGCbeta1[i] =   SGCptr->SGCgamma[0]  + SGCptr->SGCgamma[1] *(1/s) + SGCptr->SGCgamma[2]  *((1/s)*(1/s)) + SGCptr->SGCgamma[3]  *((1/s)*(1/s)*(1/s)) + SGCptr->SGCgamma[4]*s   + SGCptr->SGCgamma[5]*s*s  + SGCptr->SGCgamma[6]*s*s*s + SGCptr->SGCgamma[7] *sqrt(s) ;
    SGCptr->SGCbeta2[i] =   SGCptr->SGCgamma[8]  + SGCptr->SGCgamma[9] *(1/s) + SGCptr->SGCgamma[10]  *((1/s)*(1/s)) + SGCptr->SGCgamma[11]  *((1/s)*(1/s)*(1/s)) + SGCptr->SGCgamma[12]*s  + SGCptr->SGCgamma[13]*s*s + SGCptr->SGCgamma[14]*s*s*s + SGCptr->SGCgamma[15] *sqrt(s) ;
    SGCptr->SGCbeta3[i] =   SGCptr->SGCgamma[16] + SGCptr->SGCgamma[17]*(1/s) + SGCptr->SGCgamma[18] *((1/s)*(1/s)) + SGCptr->SGCgamma[19] *((1/s)*(1/s)*(1/s)) + SGCptr->SGCgamma[20]*s  + SGCptr->SGCgamma[21]*s*s + SGCptr->SGCgamma[22]*s*s*s + SGCptr->SGCgamma[23] *sqrt(s) ;
    SGCptr->SGCbeta4[i] =   SGCptr->SGCgamma[24] + SGCptr->SGCgamma[25]*(1/s) + SGCptr->SGCgamma[26] *((1/s)*(1/s)) + SGCptr->SGCgamma[27] *((1/s)*(1/s)*(1/s)) + SGCptr->SGCgamma[28]*s  + SGCptr->SGCgamma[29]*s*s + SGCptr->SGCgamma[30]*s*s*s + SGCptr->SGCgamma[31] *sqrt(s) ;
    // this is the wetted perimiter fraction at SGC beta min threshold
    SGCptr->SGCbeta5[i] =   SGCptr->SGCbeta1[i]*SGCptr->SGCbetahmin + SGCptr->SGCbeta2[i]*SGCptr->SGCbetahmin*SGCptr->SGCbetahmin;
  }
}
