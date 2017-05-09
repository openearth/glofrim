/*
*****************************************************************************
FLOODPLAIN FLOW WITH ACCELERATION
---------------------------------

Calculate flow between floodplain cells using acceleration formulation.

*****************************************************************************
*/

#include "lisflood.h"

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------------

// CALCULATE FLOODPLAIN (ie NON WIER) FLOW BETWEEN A POINT AND ITS W NEIGHBOUR USING ACCELERATION FORMULATION
double CalcFPQxAcc(int i,int j,States *Statesptr,Pars *Parptr,Solver *Solverptr,Arrays *Arrptr)
{
  double z0,z1,h0,h1,Sf=0.0,hflow,fn,dh=0.0,Q,q0,g,dt;
  int p0,p1,pq0;
  int pqy1,pqy2, pqy3,pqy4;
  double qy_avg,qvect;
  double qup,qdown, y0, y1;

  p0=i+j*Parptr->xsz;
  p1=i+1+j*Parptr->xsz;



  pq0=i+j*(Parptr->xsz+1)+1;

  z0=Arrptr->DEM[p0];
  z1=Arrptr->DEM[p1];
  h0=Arrptr->H[p0];
  h1=Arrptr->H[p1];
  q0=Arrptr->Qxold[pq0]; // in m2/s
  g=Solverptr->g;
  dt=Solverptr->Tstep;

  if(Solverptr->fricSolver2D==ON){
	    pqy1=i+j*(Parptr->xsz+1);
	    pqy2=(i+1)+j*(Parptr->xsz+1);
	    pqy3=i+(j+1)*(Parptr->xsz+1);
	    pqy4=(i+1)+(j+1)*(Parptr->xsz+1);
	    qy_avg=(Arrptr->Qyold[pqy1]+Arrptr->Qyold[pqy2]+Arrptr->Qyold[pqy3]+Arrptr->Qyold[pqy4])/4;
	    qvect=sqrt(q0*q0+qy_avg*qy_avg);
  }
  else{
	  qvect=q0;
  }


  if(Arrptr->Manningsn!=NULL) fn=0.5*(Arrptr->Manningsn[p0]+Arrptr->Manningsn[p1]);
  else fn=Parptr->FPn;
	y0=z0+h0;
	y1=z1+h1;

	dh=y0-y1;
  
    Sf=-dh/Parptr->dx;
    hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
    hflow=getmax(hflow,0);
    hflow=getmin(hflow,Solverptr->MaxHflow);
	// added to record Hflow
	//Arrptr->Hflowx[pq0] = hflow;

    qup=Arrptr->Qxold[pq0-1];
    qdown=Arrptr->Qxold[pq0+1];

    if(MaskTest(Arrptr->ChanMask[p0],Arrptr->ChanMask[p1]) && hflow>Solverptr->DepthThresh)
	{
		if (Statesptr->routing==OFF || fabs(Sf)<Parptr->RouteSfThresh) // CCS If routing scheme is enabled, only calc Q where Sf is below threshold value
		{
			//Q-centred scheme
    		Q=((Solverptr->theta*q0+0.5*(1-Solverptr->theta)*(qup+qdown))-(g*dt*hflow*Sf))/(1+g*dt*hflow*fn*fn*fabs(qvect)/(pow(hflow,(10.0/3.0))))*Parptr->dx;

			//Correction for situations where centred scheme introduces mass errors
			//if(copysign(1,Q)!= copysign(1,dh)){
			if (Q*dh < 0.0){ // version of line above that will compile on windows machine
				//Semi-implicit scheme of Bates et al (2010)
				Q=(q0-(g*dt*hflow*Sf))/(1+g*dt*hflow*fn*fn*fabs(qvect)/(pow(hflow,(10.0/3.0))))*Parptr->dx;
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
	    Arrptr->Vx[pq0]=Q/Parptr->dx/hflow;
	    Arrptr->maxVx[pq0]=getmax(Arrptr->maxVx[pq0],fabs(Arrptr->Vx[pq0]));
	  }
	  else Arrptr->Vx[pq0]=0.0;
  }
  return(Q);
}

//-----------------------------------------------------------------------------------
// CALCULATE FLOODPLAIN (ie NON WIER) FLOW BETWEEN A POINT AND ITS S NEIGHBOUR USING ACCELERATION FORMULATION
double CalcFPQyAcc(int i,int j,States *Statesptr,Pars *Parptr,Solver *Solverptr,Arrays *Arrptr)
{
  double z0,z1,h0,h1,Sf=0.0,hflow,fn,dh=0.0,Q=0.0,q0,g,dt;
  int p0,p1,pq0;
  int pqx1,pqx2, pqx3,pqx4;
  double qx_avg,qvect;
  double qup,qdown, y0, y1;

  p0=i+j*Parptr->xsz;
  p1=i+(j+1)*Parptr->xsz;
  
  pq0=i+(j+1)*(Parptr->xsz+1);

  z0=Arrptr->DEM[p0];
  z1=Arrptr->DEM[p1];
  h0=Arrptr->H[p0];
  h1=Arrptr->H[p1];
  q0=Arrptr->Qyold[pq0]; // in m2/s
  g=Solverptr->g;
  dt=Solverptr->Tstep;

  if(Solverptr->fricSolver2D==ON){
	  pqx1=i+j*(Parptr->xsz+1);
	  pqx2=(i+1)+j*(Parptr->xsz+1);
	  pqx3=i+(j+1)*(Parptr->xsz+1);
	  pqx4=(i+1)+(j+1)*(Parptr->xsz+1);
	  qx_avg=(Arrptr->Qxold[pqx1]+Arrptr->Qxold[pqx2]+Arrptr->Qxold[pqx3]+Arrptr->Qxold[pqx4])/4;
	  qvect=sqrt(q0*q0+qx_avg*qx_avg);
  }
  else{
	  qvect=q0;
  }


  if(Arrptr->Manningsn!=NULL) fn=0.5*(Arrptr->Manningsn[p0]+Arrptr->Manningsn[p1]);
  else fn=Parptr->FPn;
	y0=z0+h0;
	y1=z1+h1;
	dh=y0-y1;

    Sf=-dh/Parptr->dx;
    hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
    hflow=getmax(hflow,0);
    hflow=getmin(hflow,Solverptr->MaxHflow);
	// added to record Hflow
	//Arrptr->Hflowy[pq0] = hflow;

    qup=Arrptr->Qyold[i+(j)*(Parptr->xsz+1)];
    qdown=Arrptr->Qyold[i+(j+2)*(Parptr->xsz+1)];

    if(MaskTest(Arrptr->ChanMask[p0],Arrptr->ChanMask[p1]) && hflow>Solverptr->DepthThresh) 
    {
		if (Statesptr->routing==OFF || fabs(Sf)<Parptr->RouteSfThresh) // CCS If routing scheme is enabled, only calc Q where Sf is below threshold value
		{
    		// q-centred scheme
			Q=((Solverptr->theta*q0+0.5*(1-Solverptr->theta)*(qup+qdown))-(g*dt*hflow*Sf))/(1+g*dt*hflow*fn*fn*fabs(qvect)/(pow(hflow,(10.0/3.0))))*Parptr->dx;

			//Correction for situations where centred scheme introduces mass errors
			//if(copysign(1,Q)!= copysign(1,dh)){
			if (Q*dh < 0.0){ // version of line above that will compile on windows machine
				//Semi-implicit scheme of Bates et al (2010)
				Q=(q0-(g*dt*hflow*Sf))/(1+g*dt*hflow*fn*fn*fabs(qvect)/(pow(hflow,(10.0/3.0))))*Parptr->dx;
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
	    Arrptr->Vy[pq0]=Q/Parptr->dx/hflow;
	    Arrptr->maxVy[pq0]=getmax(Arrptr->maxVy[pq0],fabs(Arrptr->Vy[pq0]));
	  }
	  else Arrptr->Vy[pq0]=0.0;
  }
  return(Q);
}
////-----------------------------------------------------------------------------------
//// Calculate max(H) for each timestep
double CalcMaxH(Pars *Parptr, Arrays *Arrptr)
{
	int i,j,p0;
	double h0, Hmax=0.0,ThreadMax;


// Loop through to calculate maximum water depth
#pragma omp parallel for private(i,p0,h0,ThreadMax)
	for(j=0;j<Parptr->ysz;j++) {
		ThreadMax=0.0;
		for(i=0;i<Parptr->xsz;i++)
		{
			p0=i+j*Parptr->xsz;
			h0=Arrptr->H[p0];
			if(MaskTestAcc(Arrptr->ChanMask[p0])) 
			{
				ThreadMax=getmax(h0,ThreadMax);
			}
		}
		#pragma omp critical
		{
			Hmax=getmax(Hmax,ThreadMax);
		}
	}

return(Hmax);
}
//-----------------------------------------------------------------------------------
// Calculate timestep for acceleration version based on 2D shallow water CFL condition
void CalcT(Pars *Parptr, Solver *Solverptr, Arrays *Arrptr)
{
	double cfl=0.7, g;
	double MH=0.0;
	double locT;

	g=Solverptr->g;
	cfl=Solverptr->cfl;
	
	// Calculate maximum water depth every timestep
	MH=CalcMaxH(Parptr, Arrptr);
	
	// Apply local 2D CFL condition and use global minimum timestep to ensure global stability
	if(MH >	Solverptr->DepthThresh) //Don't apply where h=0 as h appears in equation denominator
	{
		// Old time step control based on Brett's code (obsolete/wrong)
		/*
		for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
		{
			p0=i+j*Parptr->xsz;
			pqx=i+j*(Parptr->xsz+1)+1;
			pqy=i+(j+1)*(Parptr->xsz+1);
			locH=Arrptr->H[p0];
			locV=Arrptr->Qyold[pqx]/Parptr->dx;
			locU=Arrptr->Qxold[pqy]/Parptr->dx;
			if(locH > Solverptr->DepthThresh)
			{
				//locT=cfl*Parptr->dx/(fabs(locU)+fabs(locV)+(2*sqrt(g*locH)));
				locT=cfl*Parptr->dx/(sqrt(g*locH));
				Solverptr->Tstep=getmin(Solverptr->Tstep,locT);
			}
		}
		*/
		// Time step control for stability from actual equations (MSH, implemented by TJF)
		locT=cfl*Parptr->dx/(sqrt(g*MH));
		Solverptr->Tstep=getmin(Solverptr->Tstep,locT);
	}
	else // Set to initial timestep if h=0
	{
		Solverptr->Tstep=Solverptr->InitTstep;
	}
	
return;
}
//----------------------------------------------------------------------------
// MaskTest for channel for MaxH calculation
int MaskTestAcc(int m1)
{
  if(m1==-1) return(1);
  
  return(0);
}

//----------------------------------------------------------------------------
// Update all the old qs with new qs for acceleration version
void UpdateQs(Pars *Parptr,Arrays *Arrptr)
{
	int pqptr,i,j;
	double dxinv;
	
	dxinv=1/Parptr->dx;

// Loop over all the old qx cells to update with the new qxs (in m2/s)
#pragma omp parallel for private(i,pqptr)
	for(j=0;j<Parptr->ysz;j++)
    {
      for(i=0;i<Parptr->xsz+1;i++)
      {
		  pqptr=i+j*(Parptr->xsz+1);
		  // Divide Q by dx to return m2/s
		  Arrptr->Qxold[pqptr]=Arrptr->Qx[pqptr]*dxinv;
      }
    }

// Loop over all the old qy cells to update with the new qys (in m2/s)
#pragma omp parallel for private(i,pqptr)
	for(j=0;j<Parptr->ysz+1;j++)
    {
      for(i=0;i<Parptr->xsz;i++)
      {
		  pqptr=i+j*(Parptr->xsz+1);
		  // Divide Q by dx to return m2/s
		  Arrptr->Qyold[pqptr]=Arrptr->Qy[pqptr]*dxinv;
      }
    }
	
	return;
}
//----------------------------------------------------------------------------

