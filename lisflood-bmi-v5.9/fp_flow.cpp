/*
*****************************************************************************
FLOODPLAIN FLOW
---------------------

A series of function used to calculate the flow between flodoplain cells based
on one or a combination of analytical approaches (eg. Manning's equation, weir
equations).

*****************************************************************************
*/

#include "lisflood.h"

//-----------------------------------------------------------------------------
// CALCULATES FLOODPLAIN, OVERBANK AND WEIR FLOWS
void FloodplainQ(States *Statesptr,Pars *Parptr,Solver *Solverptr, Arrays *Arrptr,SGCprams *SGCptr)
{
  int i,j;
  double h0,h1,ThreadTS, TmpTstep;
  double *hptr0,*hptr1,*qptr, *TSptr;
  int *wiptr;

    //Calculate maximum water depth and time step for acceleration version
  if(Statesptr->acceleration==ON) CalcT(Parptr, Solverptr, Arrptr);
  if(Statesptr->Roe==ON) CalcTRoe(Parptr, Solverptr, Arrptr);
  TmpTstep = Solverptr->Tstep;

	// Calculate Qx
#pragma omp parallel for private( i, h0, h1, hptr0,qptr,wiptr,TSptr,ThreadTS) 
  for(j=0;j<Parptr->ysz;j++)
  {
    hptr0=Arrptr->H+j*Parptr->xsz;
    qptr=Arrptr->Qx+j*(Parptr->xsz+1)+1;
    wiptr=Arrptr->Weir_Identx+j*(Parptr->xsz+1)+1;
	// initialise thread time step for openMP
	ThreadTS = TmpTstep;
	TSptr = &ThreadTS;
    for(i=0;i<Parptr->xsz-1;i++)
    {
      h0=*hptr0;
      h1=*(hptr0+1);
      *qptr=0.0;

      if(h0>Solverptr->DepthThresh || h1>Solverptr->DepthThresh)
      {
        if(Statesptr->weirs==ON && *wiptr!=-1) *qptr=CalcWeirQx(i,j,Parptr, Arrptr, Solverptr, Statesptr, SGCptr); // timestep update needed here TSptr
        else if(Statesptr->porosity==ON) *qptr=CalcFPQxPor(i,j,Statesptr, Parptr, Solverptr, Arrptr); // timestep update needed here TSptr
		else if(Statesptr->acceleration==ON) *qptr=CalcFPQxAcc(i,j,Statesptr, Parptr, Solverptr, Arrptr);
        else if(Statesptr->Roe==ON) *qptr=CalcFPQxRoe(i,j,Statesptr, Parptr, Solverptr, Arrptr);
		else if(Statesptr->adaptive_ts==ON || Statesptr->qlim==ON) *qptr=CalcFPQx(i,j,Statesptr, Parptr, Solverptr, Arrptr, TSptr);
      }
      qptr++;
      hptr0++;
      wiptr++;
    }
	if(Statesptr->acceleration==ON || Statesptr->Roe==ON)
	{
		// do nothing
	}
	else
	{
		#pragma omp critical
		{
			Solverptr->Tstep=getmin(Solverptr->Tstep,ThreadTS);
		}
	}
  }
	// Calculate Qy
//#pragma omp section
#pragma omp parallel for private( i, h0, h1, hptr0,hptr1,qptr,wiptr,TSptr,ThreadTS)
  for(j=0;j<Parptr->ysz-1;j++)
  {
    hptr0=Arrptr->H+j*Parptr->xsz;
    hptr1=Arrptr->H+(j+1)*Parptr->xsz;
    qptr=Arrptr->Qy+(j+1)*(Parptr->xsz+1);
    wiptr=Arrptr->Weir_Identy+(j+1)*(Parptr->xsz+1);
	// initialise thread time step for openMP
	ThreadTS = TmpTstep;
	TSptr = &ThreadTS;
    for(i=0;i<Parptr->xsz;i++)
    {
      h0=*hptr0;
      h1=*hptr1;
      *qptr=0.0;

      if(h0>Solverptr->DepthThresh || h1>Solverptr->DepthThresh)
      {
        if(Statesptr->weirs==ON && *wiptr!=-1) *qptr=CalcWeirQy(i,j, Parptr, Arrptr, Solverptr, Statesptr, SGCptr); // timestep update needed here TSptr
        else if(Statesptr->porosity==ON) *qptr=CalcFPQyPor(i,j,Statesptr, Parptr, Solverptr, Arrptr); // timestep update needed here TSptr
		else if(Statesptr->acceleration==ON) *qptr=CalcFPQyAcc(i,j,Statesptr, Parptr, Solverptr, Arrptr);
        else if(Statesptr->Roe==ON) *qptr=CalcFPQyRoe(i,j,Statesptr, Parptr, Solverptr, Arrptr);
		else if(Statesptr->adaptive_ts==ON || Statesptr->qlim==ON) *qptr=CalcFPQy(i,j,Statesptr, Parptr, Solverptr, Arrptr, TSptr);
      }
      hptr0++;
      hptr1++;
      qptr++;
      wiptr++;
    }
	if(Statesptr->acceleration==ON || Statesptr->Roe==ON)
	{
		// do nothing
	}
	else
	{
		#pragma omp critical
		{
			Solverptr->Tstep=getmin(Solverptr->Tstep,ThreadTS);
		}
	}
  }
  return;
}
//-----------------------------------------------------------------------------------
// CALCULATE FLOODPLAIN (ie NON WIER) FLOW BETWEEN A POINT AND ITS W NEIGHBOUR
double CalcFPQx(int i,int j,States *Statesptr,Pars *Parptr,Solver *Solverptr,Arrays *Arrptr, double * TSptr)
{
  double z0,z1,h0,h1,Sf,hflow,fn,dh,Qlim,Q,alpha;
  int p0,p1,pTQ;

  p0=i+j*Parptr->xsz;
  p1=i+1+j*Parptr->xsz;
  pTQ=i+j*(Parptr->xsz+1)+1;

  z0=Arrptr->DEM[p0];
  z1=Arrptr->DEM[p1];
  h0=Arrptr->H[p0];
  h1=Arrptr->H[p1];

  if(Arrptr->Manningsn!=NULL) fn=0.5*(Arrptr->Manningsn[p0]+Arrptr->Manningsn[p1]);
  else fn=Parptr->FPn;

  if(z0+h0>z1+h1 && h0>Solverptr->DepthThresh) // Flow from 0->1
  { 
    dh=z0+h0-z1-h1;
    Sf=sqrt(dh/Parptr->dx);
    hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
    hflow=getmax(hflow,0);
    hflow=getmin(hflow,Solverptr->MaxHflow);
	// added to record Hflow
	//Arrptr->Hflowx[pTQ] = hflow;

    if(MaskTest(Arrptr->ChanMask[p0],Arrptr->ChanMask[p1]) && hflow>Solverptr->DepthThresh) {
      if(dh<Solverptr->dhlin) 
      {
        Sf=sqrt(Parptr->dx/Solverptr->dhlin)*(dh/Parptr->dx);
        alpha=(pow(hflow,(5.0/3.0))*Parptr->dx_sqrt)/(fn*sqrt(Solverptr->dhlin));
      } 
      else alpha=pow(hflow,(5.0/3.0))/(2.*fn*Sf);

      Q=(pow(hflow,(5.0/3.0))*Sf*Parptr->dy/fn);
      if(Statesptr->adaptive_ts==ON) 
      {
        *TSptr=getmin(*TSptr,(0.25*Parptr->dy*Parptr->dy/alpha));
        // MT: added to record Tstep
        Arrptr->TRecx[pTQ]=*TSptr;
      } 
      else 
      {
        // flow limiter
        Qlim=Solverptr->Qlimfact*Parptr->dA*fabs(dh)/(8*Solverptr->Tstep);
        if(fabs(Q)>Qlim) 
        {
          if(Q>0) Q=Qlim;
          if(Q<0) Q=-Qlim;
          // MT added to record Qlim
          Arrptr->LimQx[pTQ]=Q;
        }
      }
    }
    else Q=0.0;
  } 
  else if(z0+h0<z1+h1 && h1>Solverptr->DepthThresh)  // Flow from 1->0
  {
    dh=z1+h1-z0-h0;
    Sf=sqrt(dh/Parptr->dx);
    hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
    hflow=getmax(hflow,0);
    hflow=getmin(hflow,Solverptr->MaxHflow);
	// added to record Hflow
	//Arrptr->Hflowx[pTQ] = hflow;

    if(MaskTest(Arrptr->ChanMask[p0],Arrptr->ChanMask[p1]) && hflow>Solverptr->DepthThresh) 
    {
      if(dh<Solverptr->dhlin) 
      {
        Sf=sqrt(Parptr->dx/Solverptr->dhlin)*(dh/Parptr->dx);
        alpha=(pow(hflow,(5.0/3.0))*Parptr->dx_sqrt)/(fn*sqrt(Solverptr->dhlin));
      } 
      else alpha=pow(hflow,(5.0/3.0))/(2.*fn*Sf);

      Q=(-pow(hflow,(5.0/3.0))*Sf*Parptr->dy/fn);
      if(Statesptr->adaptive_ts==ON) 
      {
        *TSptr=getmin(*TSptr,(0.25*Parptr->dy*Parptr->dy/alpha));
        // MT: added to record Tstep
        Arrptr->TRecx[pTQ]=*TSptr;
      } 
      else 
      {
        // flow limiter
        Qlim=Solverptr->Qlimfact*Parptr->dA*fabs(dh)/(8*Solverptr->Tstep);
        if(fabs(Q)>Qlim) 
        {
          if(Q>0) Q=Qlim;
          if(Q<0) Q=-Qlim;
          // MT added to record Qlim
          Arrptr->LimQx[pTQ]=Q;
        }
      }
    }
    else Q=0.0;

  }
  else Q=0.0;
  // option to save V's
  if (Statesptr->voutput==ON)
  {
	  if (Q!=0)
	  {
	    Arrptr->Vx[pTQ]=Q/Parptr->dx/hflow;
	    Arrptr->maxVx[pTQ]=getmax(Arrptr->maxVx[pTQ],fabs(Arrptr->Vx[pTQ]));
	  }
	  else Arrptr->Vx[pTQ]=0.0;
  }
  return(Q);
}
//-----------------------------------------------------------------------------------
// CALCULATE FLOODPLAIN (ie NON WIER) FLOW BETWEEN A POINT AND ITS S NEIGHBOUR
double CalcFPQy(int i,int j,States *Statesptr,Pars *Parptr,Solver *Solverptr,Arrays *Arrptr, double * TSptr)
{
  double z0,z1,h0,h1,Sf,hflow,fn,dh,Qlim,Q,alpha;
  int p0,p1,pTQ;

  p0=i+j*Parptr->xsz;
  p1=i+(j+1)*Parptr->xsz;
  pTQ=i+(j+1)*(Parptr->xsz+1);

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
	// added to record Hflow
	//Arrptr->Hflowy[pTQ] = hflow;

    if(MaskTest(Arrptr->ChanMask[p0],Arrptr->ChanMask[p1]) && hflow>Solverptr->DepthThresh) 
    {
      if(dh<Solverptr->dhlin) 
      {
        Sf=sqrt(Parptr->dx/Solverptr->dhlin)*(dh/Parptr->dx);
        alpha=(pow(hflow,(5.0/3.0))*Parptr->dx_sqrt)/(fn*sqrt(Solverptr->dhlin));
      } 
      else alpha=pow(hflow,(5.0/3.0))/(2.*fn*Sf);

      Q=(pow(hflow,(5.0/3.0))*Sf*Parptr->dy/fn);
      if(Statesptr->adaptive_ts==ON) 
      {
        *TSptr=getmin(*TSptr,(0.25*Parptr->dy*Parptr->dy/alpha));
        // MT: added to record Tstep
        Arrptr->TRecy[pTQ]=*TSptr;
      } 
      else 
      {
        // flow limiter
        Qlim=Solverptr->Qlimfact*Parptr->dA*fabs(dh)/(8*Solverptr->Tstep);
        if(fabs(Q)>Qlim) 
        {
          if(Q>0) Q=Qlim;
          if(Q<0) Q=-Qlim;
          // MT added to record Qlim
          Arrptr->LimQy[pTQ]=Q;
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
	// added to record Hflow
	// Arrptr->Hflowy[pTQ] = hflow;

    if(MaskTest(Arrptr->ChanMask[p0],Arrptr->ChanMask[p1]) && hflow>Solverptr->DepthThresh) 
    {
      if(dh<Solverptr->dhlin) 
      {
        Sf=sqrt(Parptr->dx/Solverptr->dhlin)*(dh/Parptr->dx);
        alpha=(pow(hflow,(5.0/3.0))*Parptr->dx_sqrt)/(fn*sqrt(Solverptr->dhlin));
      } 
      else alpha=pow(hflow,(5.0/3.0))/(2.*fn*Sf);

      Q=(-pow(hflow,(5.0/3.0))*Sf*Parptr->dy/fn);
      if(Statesptr->adaptive_ts==ON) 
      {
        *TSptr=getmin(*TSptr,(0.25*Parptr->dy*Parptr->dy/alpha));
        // MT: added to record Tstep
        Arrptr->TRecy[pTQ]=*TSptr;
      } 
      else 
      {
        // flow limiter
        Qlim=Solverptr->Qlimfact*Parptr->dA*fabs(dh)/(8*Solverptr->Tstep);
        if(fabs(Q)>Qlim) 
        {
          if(Q>0) Q=Qlim;
          if(Q<0) Q=-Qlim;
          // MT added to record Qlim
          Arrptr->LimQy[pTQ]=Q;
        }
      }
    }
    else Q=0.0;
  }
  else Q=0.0;
  // option to save V's
  if (Statesptr->voutput==ON)
  {
	  if (Q!=0)
	  {
	    Arrptr->Vy[pTQ]=Q/Parptr->dx/hflow;
	    Arrptr->maxVy[pTQ]=getmax(Arrptr->maxVy[pTQ],fabs(Arrptr->Vy[pTQ]));
	  }
	  else Arrptr->Vy[pTQ]=0.0;
  }
  return(Q);
}

//----------------------------------------------------------------------------
int MaskTest(int m1,int m2)
{
  if(m1==-1 && m2==-1) return(1);
  if(m1==-1 && m2>0) return(1);
  if(m1>0 && m2==-1) return(1);
  return(0);
}
