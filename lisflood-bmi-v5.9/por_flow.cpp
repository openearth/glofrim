/*
*****************************************************************************
POROSITY
---------------------

Simple porosity scaling algorithm used to represent floodplain structures and
topographic complexity. A number of different porosity approaches are incorporated
and each are outlined in the full documentation. Porosity maps currently derived
from porosity.cpp (but may be implemented in LISFLOOD-FP at a later date.) TJF

Keyword in parameter file is "porfile" followed by the filename.

*****************************************************************************
*/

#include "lisflood.h"

//-----------------------------------------------------------------------------------
// CALCULATE FLOODPLAIN (ie NON WIER) FLOW BETWEEN A POINT AND ITS W NEIGHBOUR
// SCALING THE FLOW BASED ON THE POROSITY, TJF
double CalcFPQxPor(int i,int j,States *Statesptr,Pars *Parptr, Solver *Solverptr, Arrays *Arrptr)
{
  double z0,z1,h0,h1,Sf,hflow,fn,dh,Qlim,Q,alpha;
  double por0,por1,por;
  int p0,p1;
  int pH0,pH1;
  int zsz,maxelev,zlev;

  zsz=Parptr->zsz;
  maxelev=(int)Parptr->maxelev;
  zlev=(int)Parptr->zlev;

  p0=i+j*Parptr->xsz;
  p1=i+1+j*Parptr->xsz;

  z0=Arrptr->DEM[p0];
  z1=Arrptr->DEM[p1];
  h0=Arrptr->H[p0];
  h1=Arrptr->H[p1];


  if(Arrptr->Manningsn!=NULL) fn=0.5*(Arrptr->Manningsn[p0]+Arrptr->Manningsn[p1]);
  else fn=Parptr->FPn;

  if(z0+h0>z1+h1 && h0>Solverptr->DepthThresh) { // Flow from 0->1
    dh=z0+h0-z1-h1;
    Sf=sqrt(dh/Parptr->dx);
    hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
    hflow=getmax(hflow,0);
    hflow=getmin(hflow,Solverptr->MaxHflow);

    if(MaskTest(Arrptr->ChanMask[p0],Arrptr->ChanMask[p1]) && hflow>Solverptr->DepthThresh) 
    {
      if(dh<Solverptr->dhlin) 
      {
        Sf=sqrt(Parptr->dx/Solverptr->dhlin)*(dh/Parptr->dx);
        alpha=(pow(hflow,(5./3.))*Parptr->dx_sqrt)/(fn*sqrt(Solverptr->dhlin));
      } else alpha=pow(hflow,(5./3.))/(2.*fn*Sf);

      if(Parptr->Por_Ident==1)
      {
        por0=Arrptr->paerial[p0];
        por1=Arrptr->paerial[p1];
        por=getmin(por0, por1);
        Q=(pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }
      else if(Parptr->Por_Ident==2)
      {
        pH0=(int)(h0/zlev);
        pH1=(int)(h1/zlev);
        if(pH0 > (maxelev/zlev)) pH0=(maxelev/zlev);
        if(pH1 > (maxelev/zlev)) pH1=(maxelev/zlev);
        por0=Arrptr->paerial[i+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz];
        por1=Arrptr->paerial[i+j*Parptr->xsz+pH1*Parptr->xsz*Parptr->ysz];
        por=getmin(por0, por1);
        Q=(pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }
      else if(Parptr->Por_Ident==3)
      {
        por0=Arrptr->pbound[i+j*Parptr->xsz+1*Parptr->xsz*Parptr->ysz];
        por1=Arrptr->pbound[i+1+j*Parptr->xsz+3*Parptr->xsz*Parptr->ysz];
        por=getmin(por0, por1);
        Q=(pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }
      else if(Parptr->Por_Ident==4)
      {
        pH0=(int)(h0/zlev);
        pH1=(int)(h1/zlev);
        if(pH0 > (maxelev/zlev)) pH0=(maxelev/zlev);
        if(pH1 > (maxelev/zlev)) pH1=(maxelev/zlev);
        por0=Arrptr->pbound[i+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz+1*Parptr->xsz*Parptr->ysz];
        por1=Arrptr->pbound[i+1+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz+3*Parptr->xsz*Parptr->ysz];
        por=getmin(por0, por1);
        Q=(pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }

      if(Statesptr->adaptive_ts==ON) 
      {
        Solverptr->Tstep=getmin(Solverptr->Tstep,(0.25*Parptr->dy*Parptr->dy/alpha));
        // MT: added to record Tstep
        Arrptr->TRecx[p0]=Solverptr->Tstep;
      } 
      else 
      {
        // flow limiter
        Qlim=Parptr->dA*fabs(dh)/(8*Solverptr->Tstep);
        if(fabs(Q)>Qlim) 
        {
          if(Q>0) Q=Qlim;
          if(Q<0) Q=-Qlim;
          // MT added to record Qlim
          Arrptr->LimQx[p0]=Q;
        }
      }
    }
    else Q=0.0;

  } 
  else if(z0+h0<z1+h1 && h1>Solverptr->DepthThresh) 
  { // Flow from 1->0
    dh=z1+h1-z0-h0;
    Sf=sqrt(dh/Parptr->dx);
    hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
    hflow=getmax(hflow,0);
    hflow=getmin(hflow,Solverptr->MaxHflow);

    if(MaskTest(Arrptr->ChanMask[p0],Arrptr->ChanMask[p1]) && hflow>Solverptr->DepthThresh) 
    {
      if(dh<Solverptr->dhlin) 
      {
        Sf=sqrt(Parptr->dx/Solverptr->dhlin)*(dh/Parptr->dx);
        alpha=(pow(hflow,(5./3.))*Parptr->dx_sqrt)/(fn*sqrt(Solverptr->dhlin));
      } else alpha=pow(hflow,(5./3.))/(2.*fn*Sf);

      if(Parptr->Por_Ident==1)
      {
        por0=Arrptr->paerial[p0];
        por1=Arrptr->paerial[p1];
        por=getmin(por0, por1);
        Q=(-pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }
      else if(Parptr->Por_Ident==2)
      {
        pH0=(int)(h0/zlev);
        pH1=(int)(h1/zlev);
        if(pH0 > (maxelev/zlev)) pH0=(maxelev/zlev);
        if(pH1 > (maxelev/zlev)) pH1=(maxelev/zlev);
        por0=Arrptr->paerial[i+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz];
        por1=Arrptr->paerial[i+j*Parptr->xsz+pH1*Parptr->xsz*Parptr->ysz];
        por=getmin(por0, por1);
        Q=(-pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }
      else if(Parptr->Por_Ident==3)
      {
        por0=Arrptr->pbound[i+j*Parptr->xsz+1*Parptr->xsz*Parptr->ysz];
        por1=Arrptr->pbound[i+1+j*Parptr->xsz+3*Parptr->xsz*Parptr->ysz];
        por=getmin(por0, por1);
        Q=(-pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }
      else if(Parptr->Por_Ident==4)
      {
        pH0=(int)(h0/zlev);
        pH1=(int)(h1/zlev);
        if(pH0 > (maxelev/zlev)) pH0=(maxelev/zlev);
        if(pH1 > (maxelev/zlev)) pH1=(maxelev/zlev);
        por0=Arrptr->pbound[i+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz+1*Parptr->xsz*Parptr->ysz];
        por1=Arrptr->pbound[i+1+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz+3*Parptr->xsz*Parptr->ysz];
        por=getmin(por0, por1);
        Q=(-pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }

      if(Statesptr->adaptive_ts==ON) 
      {
        Solverptr->Tstep=getmin(Solverptr->Tstep,(0.25*Parptr->dy*Parptr->dy/alpha));
        // MT: added to record Tstep
        Arrptr->TRecx[p0]=Solverptr->Tstep;
      } 
      else 
      {
        // flow limiter
        Qlim=Parptr->dA*fabs(dh)/(8*Solverptr->Tstep);
        if(fabs(Q)>Qlim) 
        {
          if(Q>0) Q=Qlim;
          if(Q<0) Q=-Qlim;
          // MT added to record Qlim
          Arrptr->LimQx[p0]=Q;
        }
      }
    }
    else Q=0.0;

  }
  else Q=0.0;

  return(Q);
}
//-----------------------------------------------------------------------------------
// CALCULATE FLOODPLAIN (ie NON WIER) FLOW BETWEEN A POINT AND ITS S NEIGHBOUR
// SCALING THE FLOW BASED ON THE POROSITY, TJF
double CalcFPQyPor(int i,int j,States *Statesptr,Pars *Parptr, Solver *Solverptr, Arrays *Arrptr)
{
  double z0,z1,h0,h1,Sf,hflow,fn,dh,Qlim,Q,alpha;
  double por0,por1,por;
  int p0,p1;
  int pH0,pH1;
  int zsz,maxelev,zlev;
  //ChannelSegmentType *csp;

  zsz=Parptr->zsz;
  maxelev=(int)Parptr->maxelev;
  zlev=(int)Parptr->zlev;

  p0=i+j*Parptr->xsz;
  p1=i+(j+1)*Parptr->xsz;

  z0=Arrptr->DEM[p0];
  z1=Arrptr->DEM[p1];
  h0=Arrptr->H[p0];
  h1=Arrptr->H[p1];

  if(Arrptr->Manningsn!=NULL) fn=0.5*(Arrptr->Manningsn[p0]+Arrptr->Manningsn[p1]);
  else fn=Parptr->FPn;


  if(z0+h0>z1+h1 && h0>Solverptr->DepthThresh) 
  {
    dh=z0+h0-z1-h1;
    Sf=sqrt(dh/Parptr->dx);
    hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
    hflow=getmax(hflow,0);
    hflow=getmin(hflow,Solverptr->MaxHflow);

    if(MaskTest(Arrptr->ChanMask[p0],Arrptr->ChanMask[p1]) && hflow>Solverptr->DepthThresh) 
    {
      if(dh<Solverptr->dhlin) 
      {
        Sf=sqrt(Parptr->dx/Solverptr->dhlin)*(dh/Parptr->dx);
        alpha=(pow(hflow,(5./3.))*Parptr->dx_sqrt)/(fn*sqrt(Solverptr->dhlin));
      } else alpha=pow(hflow,(5./3.))/(2.*fn*Sf);

      if(Parptr->Por_Ident==1)
      {
        por0=Arrptr->paerial[p0];
        por1=Arrptr->paerial[p1];
        por=getmin(por0, por1);
        Q=(pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }
      else if(Parptr->Por_Ident==2)
      {
        pH0=(int)(h0/zlev);
        pH1=(int)(h1/zlev);
        if(pH0 > (maxelev/zlev)) pH0=(maxelev/zlev);
        if(pH1 > (maxelev/zlev)) pH1=(maxelev/zlev);
        por0=Arrptr->paerial[i+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz];
        por1=Arrptr->paerial[i+j*Parptr->xsz+pH1*Parptr->xsz*Parptr->ysz];
        por=getmin(por0, por1);
        Q=(pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }
      else if(Parptr->Por_Ident==3)
      {
        por0=Arrptr->pbound[i+j*Parptr->xsz+2*Parptr->xsz*Parptr->ysz];
        por1=Arrptr->pbound[i+1+j*Parptr->xsz+0*Parptr->xsz*Parptr->ysz];
        por=getmin(por0, por1);
        Q=(pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }
      else if(Parptr->Por_Ident==4)
      {
        pH0=(int)(h0/zlev);
        pH1=(int)(h1/zlev);
        if(pH0 > (maxelev/zlev)) pH0=(maxelev/zlev);
        if(pH1 > (maxelev/zlev)) pH1=(maxelev/zlev);
        por0=Arrptr->pbound[i+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz+2*Parptr->xsz*Parptr->ysz];
        por1=Arrptr->pbound[i+1+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz+0*Parptr->xsz*Parptr->ysz];
        por=getmin(por0, por1);
        Q=(pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }

      if(Statesptr->adaptive_ts==ON) 
      {
        Solverptr->Tstep=getmin(Solverptr->Tstep,(0.25*Parptr->dy*Parptr->dy/alpha));
        // MT: added to record Tstep
        Arrptr->TRecy[p0]=Solverptr->Tstep;
      } 
      else 
      {
        // flow limiter
        Qlim=Parptr->dA*fabs(dh)/(8*Solverptr->Tstep);
        if(fabs(Q)>Qlim) 
        {
          //printf("Q limited %lf->%lf at t=%i\n",fabs(Q),Qlim,ts);
          if(Q>0) Q=Qlim;
          if(Q<0) Q=-Qlim;
          // MT added to record Qlim
          Arrptr->LimQy[p0]=Q;
        }
      }
    }
    else Q=0.0;

  } 
  else if(z0+h0<z1+h1 && h1>Solverptr->DepthThresh) 
  {
    dh=z1+h1-z0-h0;
    Sf=sqrt(dh/Parptr->dx);
    hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
    hflow=getmax(hflow,0);
    hflow=getmin(hflow,Solverptr->MaxHflow);

    if(MaskTest(Arrptr->ChanMask[p0],Arrptr->ChanMask[p1]) && hflow>Solverptr->DepthThresh) 
    {
      if(dh<Solverptr->dhlin) 
      {
        Sf=sqrt(Parptr->dx/Solverptr->dhlin)*(dh/Parptr->dx);
        alpha=(pow(hflow,(5./3.))*Parptr->dx_sqrt)/(fn*sqrt(Solverptr->dhlin));
      } else alpha=pow(hflow,(5./3.))/(2.*fn*Sf);

      if(Parptr->Por_Ident==1)
      {
        por0=Arrptr->paerial[p0];
        por1=Arrptr->paerial[p1];
        por=getmin(por0, por1);
        Q=(-pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }
      else if(Parptr->Por_Ident==2)
      {
        pH0=(int)(h0/zlev);
        pH1=(int)(h1/zlev);
        if(pH0 > (maxelev/zlev)) pH0=(maxelev/zlev);
        if(pH1 > (maxelev/zlev)) pH1=(maxelev/zlev);
        por0=Arrptr->paerial[i+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz];
        por1=Arrptr->paerial[i+j*Parptr->xsz+pH1*Parptr->xsz*Parptr->ysz];
        por=getmin(por0, por1);
        Q=(-pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }
      else if(Parptr->Por_Ident==3)
      {
        por0=Arrptr->pbound[i+j*Parptr->xsz+0*Parptr->xsz*Parptr->ysz];
        por1=Arrptr->pbound[i+1+j*Parptr->xsz+2*Parptr->xsz*Parptr->ysz];
        por=getmin(por0, por1);
        Q=(-pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }
      else if(Parptr->Por_Ident==4)
      {
        pH0=(int)(h0/zlev);
        pH1=(int)(h1/zlev);
        if(pH0 > (maxelev/zlev)) pH0=(maxelev/zlev);
        if(pH1 > (maxelev/zlev)) pH1=(maxelev/zlev);
        por0=Arrptr->pbound[i+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz+0*Parptr->xsz*Parptr->ysz];
        por1=Arrptr->pbound[i+1+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz+2*Parptr->xsz*Parptr->ysz];
        por=getmin(por0, por1);
        Q=(-pow(hflow,(5./3.))*Sf*Parptr->dy/fn)*por;
      }

      if(Statesptr->adaptive_ts==ON) 
      {
        Solverptr->Tstep=getmin(Solverptr->Tstep,(0.25*Parptr->dy*Parptr->dy/alpha));
        // MT: added to record Tstep
        Arrptr->TRecy[p0]=Solverptr->Tstep;
      } 
      else 
      {
        // flow limiter
        Qlim=Parptr->dA*fabs(dh)/(8*Solverptr->Tstep);
        if(fabs(Q)>Qlim) 
        {
          if(Q>0) Q=Qlim;
          if(Q<0) Q=-Qlim;
          // MT added to record Qlim
          Arrptr->LimQy[p0]=Q;
        }
      }
    }
    else Q=0.0;

  }
  else Q=0.0;

  return(Q);
}

//-----------------------------------------------------------------------------------
// CALCULATE VOLUME IN CELL WHEN POROSITY SCALES THE AREA AVAILABLE FOR STORAGE
double PorArea(int i,int j,Pars *Parptr, Arrays *Arrptr)
{
  double h0,por0;
  int p0,pH0;
  double dAPor;
  int zsz,maxelev,zlev;

  zsz=Parptr->zsz;
  maxelev=(int)Parptr->maxelev;
  zlev=(int)Parptr->zlev;

  p0=i+j*Parptr->xsz;

  h0=Arrptr->H[p0];


  // Calculate area based on the porosity value
  if(Parptr->Por_Ident==1 || Parptr->Por_Ident==3)
  {
    por0=Arrptr->paerial[p0];
    dAPor=Parptr->dA*por0;
  }
  else if(Parptr->Por_Ident==2 || Parptr->Por_Ident==4)
  {
    pH0=(int)(h0/zlev);
    if(pH0 > (maxelev/zlev)) pH0=(maxelev/zlev);
    por0=Arrptr->paerial[i+j*Parptr->xsz+pH0*Parptr->xsz*Parptr->ysz];
    dAPor=Parptr->dA*por0;

  }


  return(dAPor);
}
//-----------------------------------------------------------------------------------
