/*
*****************************************************************************
CHANNEL FLOW
---------------------

Newton-Raphson non-linear solver for channel flow discretised using a 1D kinematic
wave approximation.

void   ChannelQ():	Explicit non-linear solver to calculate channel flows one step
forward in time using Newton-Rapson (double Newton_Rapshon()).
double CalcQ,CalcA:	Calcualtes discharge in a rectangular channel and calculates
the cross sectional area.
double BankQ:		Calculate the flow out of the channel by summing appropriate
floodplain flows.

*****************************************************************************
*/

#include "lisflood.h"

//-----------------------------------------------------------------------------
// Set start water depths for channels
void SetChannelStartH(States *Statesptr, Pars *Parptr, Arrays *Arrptr,ChannelSegmentType *ChannelSegments, vector<int> *RiversIndexVecPtr, int *RiversIndexPtr)
{

  int i, chseg, pi, pj;
  int nriv, high, low; // CCS For mulitple river loop
  ChannelSegmentType *csp;

  for(nriv=0; nriv<(int)RiversIndexVecPtr->size(); nriv++) // CCS Multiple river loop
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

	  for(chseg=low; chseg<=high; chseg++) // CCS
	  {
		csp=ChannelSegments+chseg;

		// Setup starting water depth
		for(i=0;i<csp->chsz;i++)
		{
		  pi=csp->ChanX[i];
		  pj=csp->ChanY[i];

		  // set all water depths to start value (either default or defined in par file)
		  if(chseg!=low && i==csp->chsz-1) // trib and last point - ie dummy junction node // CCS
		  {
			csp->JunctionH=Parptr->ch_start_h; 
		  }
		  else // all points on main channel and all but last one on tribs
		  {
			Arrptr->H[pi+pj*Parptr->xsz]=Parptr->ch_start_h; 
		  }
		}
	  }
  }

}

//-----------------------------------------------------------------------------
// calculate initial channel flows assuming steady state and no fp flow
void CalcChannelStartQ(States *Statesptr, Pars *Parptr, Arrays *Arrptr,ChannelSegmentType *ChannelSegments, vector<int> *RiversIndexVecPtr, int *RiversIndexPtr)
{

  int i, chseg, pi, pj;
  int nriv, high, low; // CCS For multiple river loop
  double Qstart;
  ChannelSegmentType *csp;

  // reverse loop to get trib flows to pass on to d/s channel
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

	  for(chseg=high; chseg>=low ; chseg--) // CCS
	  {
		csp=ChannelSegments+chseg;

		//set Q to zero at start
		Qstart = 0;

		// loop through and pickup any inflows adding them up as we go downstream and calculating the flow depth as we go
		for(i=0;i<csp->chsz;i++)
		{
		  pi=csp->ChanX[i];
		  pj=csp->ChanY[i];

		  if(csp->Q_Ident[i]==4 ) 
		  {
			// fixed flow 
			Qstart+=csp->Q_Val[i];	
		  }
		  if(csp->Q_Ident[i]==5) 
		  {
			// interpolate from hydrograph
			Qstart+=InterpBC(csp->QVarlist[(int)csp->Q_Val[i]],0); 
		  }
		  if(csp->Q_Ident[i]==7) 
		  {
			// tributary connection
			Qstart+=csp->Q_Val[i];	
		  }

		  csp->ChanQ[i] = Qstart;

		}

		if(chseg!=low) // only for tribs record flow in d/s channel // CCS
		{
		  ChannelSegments[csp->Next_Segment].Q_Val[csp->Next_Segment_Loc]=Qstart;
		}
	  }
  }

  return;
}
//-----------------------------------------------------------------------------
// Set start water depths for channels based on initial channel flows
void SetChannelStartHfromQ(States *Statesptr, Pars *Parptr, Arrays *Arrptr,ChannelSegmentType *ChannelSegments, Solver *Solverptr, vector<int> *RiversIndexVecPtr, int *RiversIndexPtr)
{

  int i, chseg, pi, pj,iter,pid,pjd;
  double divg,Elevus,Elevds,dh,err_d,Jus=0.0,Jds=0.0,Eus=0.0,Eds=0.0,dx,eps;
  double WSbc;  // value used to hold last channel element flow h (d/s boundary condition)
  int HoutFREE=OFF; // flag to indicate calc h BC from slope on  the fly later.
  int nodes;  // temporary variable to hold number of cells/nodes in segment
  int nriv, low, high; // CCS For multiple channel loop
  ChannelSegmentType *csp;
  divg=Solverptr->divg;

  // Diffusive channel solver start H from Q
  if(Statesptr->diffusive==ON) 
  {
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
		  for(chseg=low; chseg<=high; chseg++) // CCS
		  {
		  csp=ChannelSegments+chseg;
	    // loop through and calc H from previously worked out Q (reverse loop)
		  for(i=csp->chsz-1;i>=0;i--)
		  {
			pi=csp->ChanX[i];
			pj=csp->ChanY[i];
			// setup temporary variable - just because it is easier to read code
			nodes = csp->chsz;
	  
			if(chseg!=low && i==csp->chsz-1) // trib and last point - ie dummy junction node // CCS
			{
				//csp->JunctionH=CalcA(csp->ChanN[i],csp->Shalf[i],csp->ChanWidth[i],csp->ChanQ[i])/csp->ChanWidth[i]; // calc new area from this flow and calc h
				csp->JunctionH=Arrptr->H[pi+pj*Parptr->xsz];
			}
			else if (chseg==low && i==csp->chsz-1) // last point in main channel // CCS
			{
				if(csp->Q_Ident[nodes-1]==1) //free boundary, calc h from slope ## This BC needs some stability work
				{
					// set flag so we can calc h from the flow and slope "on the fly" in CalcF() function
					WSbc=0; // set dummy value as will be calculated later
					HoutFREE=ON;
				}
				else if(csp->Q_Ident[nodes-1]==2) // fixed H out 
				{
					// note - we will need to subtract DEM Elev to get h from water elevation entered
					WSbc=csp->Q_Val[nodes-1];  
				}
				else if(csp->Q_Ident[nodes-1]==3) //interpolate H from stage hydrograph
				{
					// note - we will need to subtract DEM Elev to get h from water elevation entered
					WSbc=InterpBC(csp->QVarlist[(int)csp->Q_Val[nodes-1]],Solverptr->t); 
				}
				else if(csp->Q_Ident[nodes-1]==8) //interpolate H from rating curve ## This BC needs some stability work
				{
					// set flag to 2 so we can calc h from the stage discharge curve "on the fly" in CalcF() function
					HoutFREE=2;
					WSbc = InterpBC(csp->QVarlist[(int)csp->Q_Val[nodes-1]],csp->Q_Val[nodes-1]);
				}
				else 
				{
					printf("WARNING: Main Channel has no d/s BC\n"); 
					// already checked in input function, and set to FREE by default, but leave this here just in case of any odd bugs
				}
				if (HoutFREE==ON) // only if HoutFREE flagged ON
				{
					// calc h "on the fly" from slope for FREE BC
					if(csp->Q_Val[csp->chsz-1] < -0.999) // if -1 then use channel slope - use bed slope of last segment
					{
						Arrptr->H[csp->ChanX[nodes-1]+csp->ChanY[nodes-1]*Parptr->xsz]=CalcA(csp->ChanN[nodes-1],csp->Shalf[nodes-1],csp->ChanWidth[nodes-1],csp->ChanQ[nodes-1])
							/csp->ChanWidth[nodes-1]; 				  
					}
					else // else use user supplied slope
					{
						Arrptr->H[csp->ChanX[nodes-1]+csp->ChanY[nodes-1]*Parptr->xsz]=CalcA(csp->ChanN[nodes-1],sqrt(csp->Q_Val[nodes-1]),csp->ChanWidth[nodes-1],csp->ChanQ[nodes-1])
							/ csp->ChanWidth[nodes-1];
					}
				}
				else if (HoutFREE==2) // HoutFREE flagged for stage discharge curve
				{
					// subtract dem to get water depth
					Arrptr->H[csp->ChanX[nodes-1]+csp->ChanY[nodes-1]*Parptr->xsz] = WSbc - Arrptr->DEM[csp->ChanX[nodes-1]+csp->ChanY[nodes-1]*Parptr->xsz];
					// This BC needs some stability work
				}
				else // HoutFREE not flagged
				{
					// subtract dem to get water depth
					if(chseg==low) Arrptr->H[csp->ChanX[nodes-1]+csp->ChanY[nodes-1]*Parptr->xsz] = WSbc - Arrptr->DEM[csp->ChanX[nodes-1]+csp->ChanY[nodes-1]*Parptr->xsz]; // CCS
					else         Arrptr->H[csp->ChanX[nodes-1]+csp->ChanY[nodes-1]*Parptr->xsz] = WSbc - csp->JunctionDEM; //trib uses dummy junction node
				}
			}
	  
			/*
			Full dynamic steady state solution to incorporate dys/dx and u/g du/dx terms to improve stability of
			initial water surface profile
			*/
			else // all points on main channel and all but last one on trib
			{
				iter=0;
				pid=csp->ChanX[i+1];
				pjd=csp->ChanY[i+1];
				Arrptr->H[pi+pj*Parptr->xsz]=Arrptr->H[pid+pjd*Parptr->xsz];
				Elevds=Arrptr->DEM[pid+pjd*Parptr->xsz]+Arrptr->H[pid+pjd*Parptr->xsz];
				Jds=CalcEnergySlope(csp->ChanN[i+1],csp->ChanWidth[i+1],Arrptr->H[pid+pjd*Parptr->xsz],csp->ChanQ[i+1]);
				Eds=Elevds+pow(csp->ChanQ[i+1]/(csp->ChanWidth[i+1]*Arrptr->H[pid+pjd*Parptr->xsz]),2)*divg*Solverptr->dynsw;
				dx=(csp->Chainage[i+1]-csp->Chainage[i]);
				dh=0.00001;
				eps=1.5*dh;
				err_d=-2.0*dh;
				do {
					iter=iter+1;
					Arrptr->H[pi+pj*Parptr->xsz]=Arrptr->H[pi+pj*Parptr->xsz]+dh*signR(err_d);
					Elevus=Arrptr->DEM[pi+pj*Parptr->xsz]+Arrptr->H[pi+pj*Parptr->xsz];
					Jus=CalcEnergySlope(csp->ChanN[i],csp->ChanWidth[i],Arrptr->H[pi+pj*Parptr->xsz],csp->ChanQ[i]);
					Eus=Elevus+pow(csp->ChanQ[i]/(csp->ChanWidth[i]*Arrptr->H[pi+pj*Parptr->xsz]),2)*divg*Solverptr->dynsw;
					err_d=0.5*(Jus+Jds)*dx+Eds-Eus;
				} while (fabs(err_d)>eps);
			}
		}
	}
	}
	}

  // Kinematic channel solver start H from Q
  else 
  {
	  for(nriv=0; nriv<(int)RiversIndexVecPtr->size(); nriv++) // CCS Multiple channel loop
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
    
		// reverse loop 
		for(chseg=high; chseg>=low; chseg--) // CCS
		{
			csp=ChannelSegments+chseg;

			// loop through and calc H from previously worked out Q
			for(i=0;i<csp->chsz;i++)
			{
			pi=csp->ChanX[i];
			pj=csp->ChanY[i];

			if(chseg!=low && i==csp->chsz-1) // trib and last point - ie dummy junction node // CCS
			{
				csp->JunctionH=CalcA(csp->ChanN[i],csp->Shalf[i],csp->ChanWidth[i],csp->ChanQ[i])/csp->ChanWidth[i]; // calc new area from this flow and calc h
			}
			else // all points on main channel and all but last one on tribs
			{
				Arrptr->H[pi+pj*Parptr->xsz]=CalcA(csp->ChanN[i],csp->Shalf[i],csp->ChanWidth[i],csp->ChanQ[i])/csp->ChanWidth[i]; // calc new area from this flow and calc h
			}
		}
	}
	}
	 
  }
	 

  
  //write_profile("FREE",-1,".profile",Statesptr,ChannelSegments,Arrptr,Parptr);
  
  return;
}
//----------------------------------------------------------------------------
// Calculate energy slope using mannings, slope^1/2, channel width, water depth
// and initial Q
double CalcEnergySlope(double n,double w,double h,double Q)
{
  double s;
  s=pow((Q*n)/(w*pow(h,(5./3.))),2);
  return(s);
}
//---------------------------------------------------------------------------
void ChannelQ(double deltaT, States *Statesptr, Pars *Parptr,Solver *Solverptr,BoundCs *BCptr,ChannelSegmentType *ChannelSegments,Arrays *Arrptr, vector<int> *RiversIndexVecPtr, int *RiversIndexPtr)
{
  // #######  kinematic solver #######

  int i;			// loop counter
  int pi;			// x position of the next element
  int pj;			// y position of the next element
  int chseg;	// channel segment loop counter
  double alpha0; // alpha for current element
  double alpha1; // alpha for next element (alpha is the constant bit of the mannings equation without the Area term)
  double dxdt; // channel segment length / timestep
  double C;		// flow at current element
  double ct;		// mannings for current element
  double phi;  // ratio of cell size to channel width (dx/w)
  double Qtmp;  //temporary flow value used to hold first channel element flow area at start
  double qc;    //temporary flow value used to hold channel segment outflow
  ChannelSegmentType *csp;  //local pointer to channel segment
  int nriv, low, high; // CCS For multiple river loop
  double qc_total; // CCS for storing total domain channel Qout

  double *qc_store = new double[RiversIndexVecPtr->size()]; //CCS preallocate array for storing parallel channel solver output 

  #pragma omp parallel for private(high, low, chseg, csp, i, Qtmp, alpha0, alpha1, ct, pi, pj, phi, dxdt, C, qc) // parallelised by CCS
  for(nriv=0; nriv<(int)RiversIndexVecPtr->size(); nriv++) // CCS 
	  {
		  // Main River Loop CCS
		  high = RiversIndexPtr[nriv]-1;
		  if(nriv==0)
		  {
			  low = 0;
		  }
		  else
		  {
			  low = RiversIndexPtr[nriv-1];
		  }

		  // Main channel Loop CCS
		  for(chseg=low; chseg<=high; chseg++)
		  {  
			// set up local pointer to this segment
			csp=ChannelSegments+chseg;

			// calculate Areas from water height and channel width
			#pragma omp parallel for
			for(i=0;i<csp->chsz;i++) 
			{
				if(chseg!=low && i==csp->chsz-1) csp->A[i]=csp->JunctionH*csp->ChanWidth[i];  // trib and last point - ie dummy junction node // CCS
				else                           csp->A[i]=Arrptr->H[csp->ChanX[i]+csp->ChanY[i]*Parptr->xsz]*csp->ChanWidth[i]; // all other points
			}

			// Fix first channel element with inflow area
			if(csp->Q_Ident[0]==4) Qtmp=csp->Q_Val[0];  // fixed Q inflow
			if(csp->Q_Ident[0]==5) Qtmp=InterpBC(csp->QVarlist[(int)csp->Q_Val[0]],Solverptr->t); //interpolate Q from hydrograph
			csp->NewA[0]=CalcA(csp->ChanN[0],csp->Shalf[0],csp->ChanWidth[0],Qtmp); // calc new area from this flow

			// #### Solver #### - Use explicit non-linear scheme to fill in A values

			// calculate alpha1 for first element
			alpha1=csp->Shalf[0]/(csp->ChanN[0]*pow(csp->ChanWidth[0],(2./3.)));

			// loop through channel elements and fill in A values
			for(i=0;i<csp->chsz-1;i++)
			{
				// get mannings 
				ct=csp->ChanN[i+1];
				// get the x position 
				pi=csp->ChanX[i+1];
				// get the y position 
				pj=csp->ChanY[i+1];

				// if all the following are true together (AND)
				if(Statesptr->NCFS // what is NCFS? default set to 1 but not reset anywhere - also checked in channel vol calc - see DomainVol()
				&& Arrptr->DEM[pi+pj*Parptr->xsz]+Arrptr->H[pi+pj*Parptr->xsz]>csp->BankZ[i+1] //water elevation is above bank elevation
				&& csp->ChanWidth[i+1]<Parptr->dx)																							//channel width is less than cell size
				{
				// set phi to ratio of cell size to channel width (dx/w)
				phi=Parptr->dx/csp->ChanWidth[i+1];
				}
				else phi=1.0; // set ratio to 1 

				// set alpha0 to previous iteration's alpha1
				alpha0=alpha1;

				// calculate new alpha1
				// alpha = ((1/n)*S^(1/2))/w^(2/3) i.e. the constant bit of the mannings equation without the Area term Q=VA=alpha*A^(5/3)
				alpha1=csp->Shalf[i+1]/(ct*pow(csp->ChanWidth[i+1],(2./3.)));

				// calculate dx/dt i.e. channel segment length / timestep
				dxdt=csp->Chandx[i]/deltaT;

				// calculate flow i.e. Q = alpha0 * NewA^(5/3) + dx/dt term
				C=alpha0*pow(csp->NewA[i],(5./3.))+phi*dxdt*csp->A[i+1];

				// remove flow that goes out of bank
				if(chseg!=low && i+1==csp->chsz-1) // CCS
				{
				// if a tributary and node is last one - ie dummy node - we do not calculate bank flows 
				// - else we double count because it is already done for main channel node.
				}
				else
				{ // all others we do
				C-=BankQ(i+1,csp,Parptr,Arrptr);
				}

				// check to see if there are any inflows and add them
				if(csp->Q_Ident[i+1]==4) C+=csp->Q_Val[i+1];														    // fixed flow
				if(csp->Q_Ident[i+1]==5) C+=InterpBC(csp->QVarlist[(int)csp->Q_Val[i+1]],Solverptr->t); // interpolate from hydrograph
				if(csp->Q_Ident[i+1]==7) C+=csp->Q_Val[i+1];															 // tributary connection

				// use Newton Raphson solver to calculate new A for next element
				csp->NewA[i+1]=Newton_Raphson(csp->A[i+1],dxdt*phi,alpha0,alpha1,C,Solverptr);
			}

			// update H based on new areas calculated, divided by channel width
			for(i=0;i<csp->chsz;i++) 
			{
				if(chseg!=low && i==csp->chsz-1) csp->JunctionH=csp->NewA[i]/csp->ChanWidth[i];  // trib and last point - ie dummy junction node // CCS
				else                       Arrptr->H[csp->ChanX[i]+csp->ChanY[i]*Parptr->xsz]=csp->NewA[i]/csp->ChanWidth[i]; // all other points

				// record Q for profile output
				csp->ChanQ[i]=CalcQ(csp->ChanN[i], csp->Shalf[i],
				csp->ChanWidth[i],Arrptr->H[csp->ChanX[i]+csp->ChanY[i]*Parptr->xsz]); 
			}

			// calculate outflow from last channel segment
			if(chseg==low) // main channel // CCS
			{
				qc=CalcQ(csp->ChanN[csp->chsz-1], csp->Shalf[csp->chsz-1], csp->ChanWidth[csp->chsz-1],
				Arrptr->H[csp->ChanX[csp->chsz-1]+csp->ChanY[csp->chsz-1]*Parptr->xsz]);
			}
			else // trib dummy node
			{
				qc=CalcQ(csp->ChanN[csp->chsz-1], csp->Shalf[csp->chsz-1], csp->ChanWidth[csp->chsz-1],
				csp->JunctionH);
			}

			if(chseg!=low) // trib // CCS
			{
				// if not main channel, set the QVal of the next segment to the outflow from the end of this segment
				ChannelSegments[csp->Next_Segment].Q_Val[csp->Next_Segment_Loc]=qc;
			}
			else // main channel
			{
				// store the main channel outflow CCS
				qc_store[nriv]=qc;

			    // update downstream water depth
				Solverptr->Hds=Arrptr->H[csp->ChanX[csp->chsz-1]+csp->ChanY[csp->chsz-1]*Parptr->xsz];
			}
			
		}

	} // end of river channel segment loop
	//sum all the main channel qouts CCS
	qc_total=0;
	for(int n=0; n<(int)RiversIndexVecPtr->size(); n++)
	{
		qc_total=qc_total+qc_store[n];
	}
	delete [] qc_store; // delete qc store array

	BCptr->QChanOut=qc_total;

  return;
}
//---------------------------------------------------------------------------
// NON-LINEAR SOLVER FOR CHANNEL FLOW
double Newton_Raphson(double Ai,double dx,double a0,double a1,double c,Solver *Solverptr)
{
  double res,A,f;

  do
  {
    f=dx*Ai+a1*pow(Ai,(5./3.));
    A=Ai-(f-c)/(dx+(5./3)*a1*pow(Ai,(2./3.)));
    res=fabs((f-c)/c);
    Ai=A;
  } while(res>Solverptr->SolverAccuracy);

  return A;
}

//---------------------------------------------------------------------------
// CALCULATE BANK FLOWS FOR ChannelQ()
double BankQ(int chani, ChannelSegmentType *ChannelSegments, Pars *Parptr, Arrays *Arrptr)
{
  double qbank0,q0,q1,q2,q3;
  int pi0,pj0;

  pi0=ChannelSegments->ChanX[chani];
  pj0=ChannelSegments->ChanY[chani];

  // Find Qbank0 from 4-neighbours, masking out channel flows
  q0=Arrptr->Qx[pi0+1+pj0*(Parptr->xsz+1)];

  q1=-Arrptr->Qx[pi0+pj0*(Parptr->xsz+1)];

  q2=Arrptr->Qy[pi0+(pj0+1)*(Parptr->xsz+1)];

  q3=-Arrptr->Qy[pi0+pj0*(Parptr->xsz+1)];

  qbank0=q0+q1+q2+q3;

  return(qbank0);
}

//----------------------------------------------------------------------------
// calculate flow using mannings, slope^1/2, channel width and water depth
double CalcQ(double n,double s,double w,double h)
{
  double Q;
  Q=fabs(s)*pow(w*h,(5./3.))/(n*pow(w,(2./3.)));
  return(Q);
}
//---------------------------------------------------------------------------
// calculate flow area using mannings, slope^1/2, channel width and flow
double CalcA(double n,double s,double w,double Q)
{
  double A;
  Q=fabs(Q);
  A=pow(n*pow(w,(2./3.))*Q/fabs(s),(3./5.));

  return(A);
}

//---------------------------------------------------------------------------
// calculate channel volume
double ChannelVol(Pars *Parptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr)
{
  double v=0;
  int i,j;
#pragma omp parallel for private (j) reduction (+ : v)
  for(i=0;i<Parptr->xsz;i++) 
  {
    for(j=0;j<Parptr->ysz;j++)
    {
      if(Arrptr->ChanMask[i+j*Parptr->xsz]!=-1) v+=Arrptr->H[i+j*Parptr->xsz]*Parptr->dA;
    }
  }

  return(v);
}

//---------------------------------------------------------------------------
//void ChannelQ_Diff(double deltaT, States *Statesptr, Pars *Parptr,Solver *Solverptr,BoundCs *BCptr,ChannelSegmentType *ChannelSegments,Arrays *Arrptr
void ChannelQ_Diff(double deltaT, States *Statesptr, Pars *Parptr,Solver *Solverptr,BoundCs *BCptr,ChannelSegmentType *ChannelSegments, Arrays *Arrptr, vector<int> *RiversIndexVecPtr, int *RiversIndexPtr)
// CHANNEL DIFFUSIVE SOLVER #################
{

  // main arrays for diffusive
  double *x,*xn,*f,**J,**JinvL;
  int *indx;
  double d,res;

  int i;		// loop counter
  int chseg;	// channel segment loop counter
  int nodes;	// temporary variable to hold number of cells/nodes in segment
  double Qbc;	// flow value used to hold first channel element flow area at start (u/s boundary condition)
  double WSbc;  // value used to hold last channel element flow h (d/s boundary condition)
  int HoutFREE=OFF; // flag to indicate calc h BC from slope on  the fly later.
  double qc;    //temporary flow value used to hold channel segment outflow
  ChannelSegmentType *csp;  //local pointer to channel segment
  int exitcrit=1;	// variable for determing the exit criterion for the NR solver; see manual for details (TJF)

  // new diffusive solver declarations
  int maxiterations=200;		// maximum number of iterations for solver
  int itcount=0;				  // iteration counter

  int nriv, low, high; // CCS For multiple river loop
  double qc_total; // CCS

  double *qc_store = new double[RiversIndexVecPtr->size()]; //CCS preallocate array for storing parallel channel solver output 

  #pragma omp parallel for private(high, low, chseg, itcount, HoutFREE, csp, nodes, x, xn, f, indx, J, i, JinvL, Qbc, WSbc, res, qc) //parrallelised by CCS
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
		  // main loop for channel segments
		  for(chseg=low; chseg<=high; chseg++) // CCS
		  {
			// set iteration counter to zero for each channel
			itcount=0;				
			// set default flags
			HoutFREE=OFF;
			// set up local pointer to this segment
			csp=ChannelSegments+chseg;
			// setup temporary variable - just because it is easier to read code
			nodes = csp->chsz;

			// allocate memory for local solver arrays
			x=new double[nodes*2] ();
			xn=new double[nodes*2] ();
			f=new double[nodes*2] ();
			indx=new int[nodes*2] ();
			J=new double*[nodes*2] ();
			for(i=0;i<2*nodes;i++) J[i]=new double[5] ();
			JinvL=new double*[2*nodes] ();
			for(i=0;i<2*nodes;i++) JinvL[i]=new double[5] ();

			// calculate Areas from water height and channel width
			for(i=0;i<nodes;i++) 
			{
			  if(chseg!=low && i==nodes-1) x[2*i]=xn[2*i]=csp->JunctionH*csp->ChanWidth[i]; // trib and last point - ie dummy junction node // CCS
			  else                         x[2*i]=xn[2*i]=Arrptr->H[csp->ChanX[i]+csp->ChanY[i]*Parptr->xsz]*csp->ChanWidth[i]; // all other points
			}

			// fill in Q values
			for(i=0;i<nodes;i++)
			{
			  // use flows from last iteration as start point for this one.
			  // for first timestep, these are already filled in by CalcChannelStartQ()
			  x[2*i+1]=xn[2*i+1]=csp->ChanQ[i];
			}

			// Fix first channel element with inflow area, ie u/s BC
			if(csp->Q_Ident[0]==4) Qbc=csp->Q_Val[0];  // fixed Q inflow
			if(csp->Q_Ident[0]==5) Qbc=InterpBC(csp->QVarlist[(int)csp->Q_Val[0]],Solverptr->t); //interpolate Q from hydrograph

			// work out d/s boundary condition ie d/s H 
			if(chseg>low) //trib // CCS
			{
			  // get H from d/s channel - add DEM value to convert to elevation - to transfer as d/s BC to trib (other bc types are elevations). 
			  // Also covers the case when d/s trib bed elevation is different from main channel.
			  WSbc =   Arrptr->H[ChannelSegments[csp->Next_Segment].ChanX[csp->Next_Segment_Loc]+ChannelSegments[csp->Next_Segment].ChanY[csp->Next_Segment_Loc]*Parptr->xsz]
			  + Arrptr->DEM[ChannelSegments[csp->Next_Segment].ChanX[csp->Next_Segment_Loc]+ChannelSegments[csp->Next_Segment].ChanY[csp->Next_Segment_Loc]*Parptr->xsz];
			}
			else  // main channel d/s BC
			{
			  if(csp->Q_Ident[nodes-1]==1) //free boundary, calc h from slope ## This BC needs some stability work
			  {
				// set flag so we can calc h from the flow and slope "on the fly" in CalcF() function
				WSbc=0; // set dummy value as will be calculated later
				HoutFREE=ON;
			  }
			  else if(csp->Q_Ident[nodes-1]==2) // fixed H out 
			  {
				// note - we will need to subtract DEM Elev to get h from water elevation entered
				WSbc=csp->Q_Val[nodes-1];  
			  }
			  else if(csp->Q_Ident[nodes-1]==3) //interpolate H from stage hydrograph
			  {
				// note - we will need to subtract DEM Elev to get h from water elevation entered
				WSbc=InterpBC(csp->QVarlist[(int)csp->Q_Val[nodes-1]],Solverptr->t); 
			  }
			  else if(csp->Q_Ident[nodes-1]==8) //interpolate H from rating curve ## This BC needs some stability work
			  {
				// set flag to 2 so we can calc h from the stage discharge curve "on the fly" in CalcF() function
				HoutFREE=2;
				WSbc = InterpBC(csp->QVarlist[(int)csp->Q_Val[nodes-1]],x[2*nodes-3]);
			  }
			  else 
			  {
				printf("WARNING: Main Channel has no d/s BC\n"); 
				// already checked in input function, and set to FREE by default, but leave this here just in case of any odd bugs
			  }
			}

			// Use implicit newton raphson scheme to solve 
			// Utilises LU decomposition - Crout's method to solve the banded diagonal matrix of linear equations
			do
			{
			  res=0.0; // reset the res incase we use the max rather than the norm
			  // setup function matrix
			  calcF(x,xn,f,deltaT,csp,Parptr,Arrptr,Qbc,chseg,WSbc,HoutFREE,Solverptr,low);
			  // setup Jacobian matrix - Q determinant of function
			  calcJ(x,xn,J,deltaT,csp,chseg,HoutFREE);
	
			  bandec(J,2*nodes,2,2,JinvL,indx,d);
			  banbks(J,2*nodes,2,2,JinvL,indx,f); // On input, f is function vector from CalcF. On output, f is solution.
			  //mprove(); // Routine for iterative improvement - not currently implemented (TJF).

			  for(i=0;i<2*nodes;i++) xn[i]-=f[i];
			  for(i=0;i<nodes;i++) if(xn[2*i]<1.0) xn[2*i]=1.0;
	  
			  switch(exitcrit) // Switch for different exit criteria to improve channel model solution. See manual for details.
			  {
			  case 1:
				  res=norm(f,2*nodes); // norm of x (the solution)
				  break;
			  case 2:
				  for(i=0;i<2*nodes;i++) res=getmax(fabs(f[i]),res); // max of x (the solution)
				  break;
			  case 3:
				  calcF(x,xn,f,deltaT,csp,Parptr,Arrptr,Qbc,chseg,WSbc,HoutFREE,Solverptr,low); // recalculate f(x) for exit criteria
				  res=norm(f,2*nodes);  // norm of f(x) (the function vector)
				  break;
			  case 4:
				  calcF(x,xn,f,deltaT,csp,Parptr,Arrptr,Qbc,chseg,WSbc,HoutFREE,Solverptr,low); // recalculate f(x) for exit criteria
				  for(i=0;i<2*nodes;i++) res=getmax(fabs(f[i]),res); // max of f(x) (the function vector)
				  break;
			  }
			} 
			while( res > Solverptr->SolverAccuracy && itcount++ < maxiterations);
			if(itcount>=maxiterations) printf("WARNING: Max iterations exceeded diffusive channel at t=%.3f in channel %i.\n",Solverptr->t,chseg); //Warning for iterations

			// update H based on new areas calculated, divided by channel width
			for(i=0;i<nodes;i++) 
			{
				if(chseg!=low && i==nodes-1)	csp->JunctionH=xn[2*i]/csp->ChanWidth[i];  // trib and last point - ie dummy junction node // CCS
				else						Arrptr->H[csp->ChanX[i]+csp->ChanY[i]*Parptr->xsz]=xn[2*i]/ csp->ChanWidth[i]; // all other points

				csp->ChanQ[i]=xn[2*i+1]; // record Q for profile output
			}

			// get outflow from last channel segment
			qc=xn[2*(nodes-1)+1];

			if(chseg>low) // trib #CCS#
			{
			  // set the QVal of the next segment to the outflow from the end of this segment
			  ChannelSegments[csp->Next_Segment].Q_Val[csp->Next_Segment_Loc]=qc;
			}
			else // main channel
			{
			  // store the main channel outflow CCS
			  qc_store[nriv]=qc;
			  // update downstream water depth
			  Solverptr->Hds=Arrptr->H[csp->ChanX[nodes-1]+csp->ChanY[nodes-1]*Parptr->xsz];
			}

			// clean up memory
			for(i=0;i<2*nodes;i++) delete[] J[i];
			for(i=0;i<2*nodes;i++) delete[] JinvL[i];
			delete[] J;
			delete[] JinvL;
			delete[] x;
			delete[] xn;
			delete[] f;
			delete[] indx;

		  } // end of main river channel segment loop
		}
	//Calculate QChanOut from all rivers CCS
	qc_total=0;
	for(int n=0; n<(int)RiversIndexVecPtr->size(); n++)
	{
		qc_total=qc_total+qc_store[n];
	}
	delete [] qc_store; // delete qc store array

	BCptr->QChanOut=qc_total;

  return;
}


//---------------------------------------------------------------------------
void calcF(double *x,double *xn,double *f, double dt, ChannelSegmentType *csp, Pars *Parptr, 
           Arrays *Arrptr, double Qin, int chseg, double WSout, int HoutFREE, Solver *Solverptr, int low)
{
  double q,a,dx;
  int i;

  int nodes=csp->chsz;
  double th=1.0; 
  double w,w0,w1;
  double w43,mann;
  double nSq;
  double Hout; //h calculated at downstream boundary condition

  for(i=0;i<nodes-1;i++)
  {

    dx=csp->Chandx[i];

    // need to use both channel widths
    w0=csp->ChanWidth[i];
    w1=csp->ChanWidth[i+1];
    w=(w0+w1)/2.0;
    w43=pow(w0,(4./3.));

    nSq=1./((1./csp->ChanN[i]+1./csp->ChanN[i+1])/2.); // calc inverse average for 2 cross-sections
    nSq=nSq*nSq; // square for later use

    f[2*i+1]=(xn[2*i]+xn[2*i+2]-x[2*i]-x[2*i+2])/(2.*dt);
    f[2*i+1]+=th*(xn[2*i+3]-xn[2*i+1])/dx;
    f[2*i+1]+=(1-th)*(x[2*i+3]-x[2*i+1])/dx;  // if th=1, use fully implicit method

    // Add overbank flows - only if not in startup mode
    f[2*i+1]+=BankQ(i+1,csp,Parptr,Arrptr)/dx;

    // check to see if there are any inflows at the node and add them
    //
    if (i!=0) 
      // make sure it is not start of channel inflow BC
    {
      if(csp->Q_Ident[i+1]==4 ) 
      {
        // fixed flow 
        f[2*i+1]-=csp->Q_Val[i+1]/dx;	
      }
      if(csp->Q_Ident[i+1]==5) 
      {
        // interpolate from hydrograph
        f[2*i+1]-=InterpBC(csp->QVarlist[(int)csp->Q_Val[i+1]],Solverptr->t)/dx; 
      }
      if(csp->Q_Ident[i+1]==7) 
      {
        // tributary connection
        f[2*i+1]-=csp->Q_Val[i+1]/dx;	
      }
    }

    f[2*i+2]=pow(csp->Shalf[i],2.0);
    if(csp->Shalf[i]<0) f[2*i+2]=-f[2*i+2];

    // new variable width code
    f[2*i+2]-=(((th/w1)*xn[2*i+2])-((th/w0)*xn[2*i]))/dx;
    f[2*i+2]-=((((1-th)/w1)*x[2*i+2])-(((1-th)/w0)*x[2*i]))/dx; // if th=1, use fully implicit method

    a=th*(xn[2*i]+xn[2*i+2])/2;
    a+=(1-th)*(x[2*i]+x[2*i+2])/2;
    q=th*(xn[2*i+1]+xn[2*i+3])/2;
    q+=(1-th)*(x[2*i+1]+x[2*i+3])/2;
    f[2*i+2]-=nSq*w43*q*fabs(q)*pow( a,(-10./3.0));
  }

  // boundary conditions
  // Qdiff of first element - Qin
  f[0]=xn[1]-Qin;


  if (HoutFREE==ON) // only if HoutFREE flagged ON
  {
    // calc h "on the fly" from slope for FREE BC
    if(csp->Q_Val[csp->chsz-1] < -0.999) // if -1 then use channel slope - use bed slope of last segment
    {
		// use a combination of previous iteration (xn) and previous timestep (x) for estimate of Q - set by th
		// if th=1, use fully implicit method  
		Hout=CalcA(csp->ChanN[nodes-1],csp->Shalf[nodes-1],csp->ChanWidth[nodes-1],((1-th)*x[2*nodes-1])+th*xn[2*nodes-1])
        /csp->ChanWidth[nodes-1];
		
		// Recalculate f(x) at last section by substituting Manning's into diffusive form of momentum equation 
		f[2*nodes-2]=pow(csp->Shalf[nodes-1],2.0);
		f[2*nodes-2]-=(((th/csp->ChanWidth[nodes-1])*xn[2*nodes-2])-((th/csp->ChanWidth[nodes-2])*xn[2*nodes-4]))/csp->Chandx[nodes-2];
		a=th*(xn[2*nodes-4]+xn[2*nodes-2])/2;
		mann=(pow(xn[2*nodes-2],(5./3.))*csp->Shalf[nodes-1])/(csp->ChanN[nodes-1]*pow(csp->ChanWidth[nodes-1],(2./3.)));
		q=th*(xn[2*nodes-3]+mann)/2;
		nSq=1./((1./csp->ChanN[nodes-2]+1./csp->ChanN[nodes-1])/2.); // calc inverse average for 2 cross-sections
		nSq=nSq*nSq;
		w43=pow(csp->ChanWidth[nodes-2],(4./3.));
		f[2*nodes-2]-=nSq*w43*q*fabs(q)*pow(a,(-10./3.0));
		//f[2*nodes-2]=((1-alpha)*xn[2*nodes-3]+alpha*xn[2*nodes-1]*csp->ChanN[nodes-1]*pow((1-alpha)*csp->ChanWidth[nodes-2]+alpha*csp->ChanWidth[nodes-1],(2./3.)))/csp->Shalf[nodes-1]-pow((1-alpha)*xn[2*nodes-4]+alpha*xn[2*nodes-2],(5./3.));
    }
    else // else use user supplied slope
    {
		// use a combination of previous iteration (xn) and previous timestep (x) for estimate of Q - set by th
		// if th=1, use fully implicit method 
		Hout=CalcA(csp->ChanN[nodes-1],sqrt(csp->Q_Val[nodes-1]),csp->ChanWidth[nodes-1],((1-th)*x[2*nodes-1])+th*xn[2*nodes-1])
        / csp->ChanWidth[nodes-1];
		
		// Recalculate f(x) at last section by substituting Manning's into diffusive form of momentum equation 
		f[2*nodes-2]=csp->Q_Val[nodes-1];
		f[2*nodes-2]-=(((th/csp->ChanWidth[nodes-1])*xn[2*nodes-2])-((th/csp->ChanWidth[nodes-2])*xn[2*nodes-4]))/csp->Chandx[nodes-2];
		a=th*(xn[2*nodes-4]+xn[2*nodes-2])/2;
		mann=(pow(xn[2*nodes-2],(5./3.))*sqrt(csp->Q_Val[nodes-1]))/(csp->ChanN[nodes-1]*pow(csp->ChanWidth[nodes-1],(2./3.)));
		q=th*(xn[2*nodes-3]+mann)/2;
		nSq=1./((1./csp->ChanN[nodes-2]+1./csp->ChanN[nodes-1])/2.); // calc inverse average for 2 cross-sections
		nSq=nSq*nSq;
		w43=pow(csp->ChanWidth[nodes-2],(4./3.));
		f[2*nodes-2]-=nSq*w43*q*fabs(q)*pow(a,(-10./3.0));
		//f[2*nodes-2]=((1-alpha)*xn[2*nodes-3]+alpha*xn[2*nodes-1]*csp->ChanN[nodes-1]*pow((1-alpha)*csp->ChanWidth[nodes-2]+alpha*csp->ChanWidth[nodes-1],(2./3.)))/sqrt(csp->Q_Val[nodes-1])-pow((1-alpha)*xn[2*nodes-4]+alpha*xn[2*nodes-2],(5./3.));
    }
  }
  else if (HoutFREE==2) // HoutFREE flagged for stage discharge curve
  {
    // subtract dem to get water depth
    Hout = WSout - Arrptr->DEM[csp->ChanX[nodes-1]+csp->ChanY[nodes-1]*Parptr->xsz];
    // This BC needs some stability work
  }
  else // HoutFREE not flagged
  {
    // subtract dem to get water depth
	if(chseg==low) Hout = WSout - Arrptr->DEM[csp->ChanX[nodes-1]+csp->ChanY[nodes-1]*Parptr->xsz]; // CCS

    else         Hout = WSout - csp->JunctionDEM; //trib uses dummy junction node
  }
  
  // Set global variables to the value of Hout calculated here
  if(chseg!=0) csp->JunctionH=Hout;
  else Arrptr->H[csp->ChanX[nodes-1]+csp->ChanY[nodes-1]*Parptr->xsz]=Hout;

  // Area diff of last element - bc area
  f[2*nodes-1]=xn[2*nodes-2]-(Hout)*csp->ChanWidth[nodes-1];

  return;
}
//---------------------------------------------------------------------------
void calcJ(double *x,double *xn,double **J, double dt, ChannelSegmentType *csp, int chseg, int HoutFREE)
{
  double a,q,dx;
  int i;
  int n=csp->chsz;
  double th=1; 
  double w,mann;
  double n2w43th,fu3;
  
  fu3=5./3.; // five divided by three for later use

  for(i=0;i<5;i++) J[0][i]=0;
  J[0][3]=1;
  
  for(i=0;i<2*n-3;i+=2)
  {
    dx=csp->Chandx[i/2];

    w=csp->ChanWidth[i/2];	
    n2w43th=csp->ChanN[i/2]*csp->ChanN[i/2]*pow(w,(4./3.))*th; 
    
	// Setting the odd rows of the compact form of the Jacobian [Continuity Equation]
    J[i+1][0]=0;
    J[i+1][1]=0.5/dt;
    J[i+1][2]=-th/dx;
    J[i+1][3]=J[i+1][1];
    J[i+1][4]=-J[i+1][2];

	// Setting the even rows of the compact form of the Jacobian [Momentum Equation]
	a=th*(xn[i]+xn[i+2])/2+(1-th)*(x[i]+x[i+2])/2;
    q=th*(xn[i+1]+xn[i+3])/2+(1-th)*(x[i+1]+x[i+3])/2;

    J[i+2][0]=th/(w*dx)+10*n2w43th*fabs(q)*q/(6*pow(a, (13./3.)));
    J[i+2][1]=-n2w43th*fabs(q)/pow(a,(10./3.));
    J[i+2][2]=-J[i+2][0];
    J[i+2][3]=J[i+2][1];
    J[i+2][4]=0;
  }
  
  // if FREE boundary, reset the final lines of the Jacobian as the derivatives of Manning's equation rather than
  // the diffusive form of the de St. Venant.
  if (HoutFREE==ON) // only if HoutFREE flagged ON
  {
	  for(i=0;i<5;i++) J[2*n-2][i]=0; // reset to 0
	  if(csp->Q_Val[csp->chsz-1] < -0.999) // if -1 then use channel slope - use bed slope of last segment
	  {
		  // Combining Manning's and diffusive continuity and momentum for final section
		  // derivative of combination wrt Q
		  mann=(pow(xn[2*n-2],(5./3.))*csp->Shalf[n-1])/(csp->ChanN[n-1]*pow(csp->ChanWidth[n-1],(2./3.)));
		  J[2*n-2][2]=-pow(csp->ChanN[n-1],2.0)*pow(csp->ChanWidth[n-2],(4./3.))*pow((xn[2*n-2]+xn[2*n-4])/2,(-10./3.))*(mann+xn[2*n-3])/2;

		  // Combining Manning's and diffusive continuity and momentum for final section
		  // derivative of combination wrt A
		  J[2*n-2][3]=-(1/(csp->Chandx[n-1]*csp->ChanWidth[n-1]));
		  J[2*n-2][3]-=pow(csp->ChanN[n-1],2.0)*pow(csp->ChanWidth[n-2],(4./3.))*pow((mann+xn[2*n-3])/2,2.0)*(10./6.)*pow((xn[2*n-2]+xn[2*n-4])/2,(-13./3.));
		  J[2*n-2][3]+=pow((xn[2*n-2]+xn[2*n-4])/2,(-10./3.))*((mann+xn[2*n-3])/2)*-fu3*pow(xn[2*n-2],(2./3.))*(csp->Shalf[n-1]/(csp->ChanN[n-1]*pow(csp->ChanWidth[n-1],(2./3.))));

		  // Old version
		  //J[2*n-2][2]+=(csp->ChanN[n-1]*pow((1-alpha)*csp->ChanWidth[n-2]+alpha*csp->ChanWidth[n-1],(2./3.)))/csp->Shalf[n-1]; // derivative of Mannings wrt Q
		  //J[2*n-2][3]=-fu3*(pow((1-alpha)*xn[2*n-4]+alpha*xn[2*n-2],(2./3.))); // derivative of Mannings wrt A
	  }
	  else // user supplied slope
	  {
		  // Combining Manning's and diffusive continuity and momentum for final section
		  // derivative of combination wrt Q
		  mann=(pow(xn[2*n-2],(5./3.))*sqrt(csp->Q_Val[n-1]))/(csp->ChanN[n-1]*pow(csp->ChanWidth[n-1],(2./3.)));
		  J[2*n-2][2]=-pow(csp->ChanN[n-1],2.0)*pow(csp->ChanWidth[n-2],(4./3.))*pow((xn[2*n-2]+xn[2*n-4])/2,(-10./3.))*(mann+xn[2*n-3])/2;

		  // Combining Manning's and diffusive continuity and momentum for final section
		  // derivative of combination wrt A
		  J[2*n-2][3]=-(1/(csp->Chandx[n-1]*csp->ChanWidth[n-1]));
		  J[2*n-2][3]-=pow(csp->ChanN[n-1],2.0)*pow(csp->ChanWidth[n-2],(4./3.))*pow((mann+xn[2*n-3])/2,2.0)*(10./6.)*pow((xn[2*n-2]+xn[2*n-4])/2,(-13./3.));
		  J[2*n-2][3]+=pow((xn[2*n-2]+xn[2*n-4])/2,(-10./3.))*((mann+xn[2*n-3])/2)*-fu3*pow(xn[2*n-2],(2./3.))*(sqrt(csp->Q_Val[n-1])/(csp->ChanN[n-1]*pow(csp->ChanWidth[n-1],(2./3.))));
		  
		  //J[2*n-2][2]+=(csp->ChanN[n-1]*pow((1-alpha)*csp->ChanWidth[n-2]+alpha*csp->ChanWidth[n-1],(2./3.)))/sqrt(csp->Q_Val[n-1]); // derivative of Mannings wrt Q
		  //J[2*n-2][3]+=-fu3*(pow((1-alpha)*xn[2*n-4]+alpha*xn[2*n-2],(2./3.))); // derivative of Mannings wrt A
	  }
  }
  
  for(i=0;i<5;i++) J[2*n-1][i]=0;
  J[2*n-1][1]=1;

  return;
}
//---------------------------------------------------------------------------
double norm(double *x,int n)
{
  double res=0.;
  for(int i=0;i<n;i++) res+=x[i]*x[i];
  return sqrt(res/n);
}
//---------------------------------------------------------------------------

//====================================================================================================
// bandec(), banbks(), SWAP() and mprove() below are all from Chapter 2.4 in the book 
// "Numerical Recipes in C" p51-54. These allow solution of band diagonal 
// linear systems by LU decomposition (Crout's method). Rewritten for matrix 
// indexes starting at zero rather than 1.
//---------------------------------------------------------------------------
#define TINY 1.e-20

void bandec(double **a, int n, int m1, int m2, double **al, int indx[], double &d)
// Given an n x n band diagonal matrix A with m1 subdiagonal rows and m2 superdiagonal rows,
// compactly stored in the array a[1..n][1..m1+m2+1]. The diagonal elements are in a[1..n][m1+1]. 
// Subdiagonal elements are in a[j..n][1..m1] (with j > 1 appropriate to the number of elements 
// on each subdiagonal). Superdiagonal elements are in a[1..j][m1+2..m1+m2+1] with j < n appropriate 
// to the number of elements on each superdiagonal. This routine constructs an LU decomposition of 
// a rowwise permutation of A. The upper triangular matrix replaces a, while the lower triangular 
// matrix is returned in al[1..n][1..m1]. indx[1..n] is an output vector which records the row 
// permutation effected by the partial pivoting; d is output as +-1 depending on whether the 
// number of row interchanges was even or odd, respectively. This routine is used in combination 
// with banbks to solve band-diagonal sets of equations.
{
  int i,j,k,l;
  int mm;
  double dum;

  mm=m1+m2+1;
  l=m1;
  for (i=1;i<=m1;i++) // Rearrange the storage a bit.
  {
    for (j=m1+2-i;j<=mm;j++) a[i-1][j-l-1]=a[i-1][j-1];
    l--;
    for (j=mm-l;j<=mm;j++) a[i-1][j-1]=0.0;
  }

  d=1.0;
  l=m1;
  for (k=1;k<=n;k++) // For each row...
  {
    dum=a[k-1][0];
    i=k;
    if (l < n) l++;
    for (j=k+1;j<=l;j++) // Find the pivot element.
    {
      if (fabs(a[j-1][1-1]) > fabs(dum))
      {
        dum=a[j-1][0];
        i=j;
      }
    }
    indx[k-1]=i;
    if (dum == 0.0) a[k-1][0]=TINY;

    //  Matrix is algorithmically singular, but proceed anyway with
    //  TINY pivot (desirable in some applications).

    if (i != k) // Interchange rows.
    {
      d = -(d);
      for (j=1;j<=mm;j++) SWAP(a[k-1][j-1],a[i-1][j-1]);
    }
    for (i=k+1;i<=l;i++) // Do the elimination.
    {
      dum=a[i-1][0]/a[k-1][0];
      al[k-1][i-k-1]=dum;
      for (j=2;j<=mm;j++) a[i-1][j-2]=a[i-1][j-1]-dum*a[k-1][j-1];
      a[i-1][mm-1]=0.0;
    }
  }
  return;
}

//---------------------------------------------------------------------------
// swap a and b
void SWAP(double &a,double &b)
{
  double dum=a;
  a=b;
  b=dum;
  return;
}

//---------------------------------------------------------------------------
void banbks(double **a, int n, int m1, int m2, double **al, int indx[], double b[])
// Given the arrays a, al, and indx as returned from bandec, and given a right-hand side vector
// b[1..n], solves the band diagonal linear equations A . x = b. The solution vector x overwrites
// b[1..n]. The other input arrays are not modified, and can be left in place for successive calls
// with different right-hand sides.
{
  int i,k,l;
  int mm;
  double dum;

  mm=m1+m2+1;
  l=m1;
  for (k=1;k<=n;k++) // Forward substitution, unscrambling the permuted rows
  {                  // as we go.
    i=indx[k-1];
    if (i != k) SWAP(b[k-1],b[i-1]);
    if (l < n) l++;
    for (i=k+1;i<=l;i++)
    {
      b[i-1] -= al[k-1][i-k-1]*b[k-1];
    }
  }

  l=1;
  for (i=n;i>=1;i--) // Backsubstitution.
  {
    dum=b[i-1];
    for (k=2;k<=l;k++) dum -= a[i-1][k-1]*b[k+i-2];
    b[i-1]=dum/a[i-1][0];
    if (l < mm) l++;
  }
}
//----------------------------------------------------------------------------
// Return the sign of a number
int signR(double a)
{
	int b;

	if(a>=0) b=1;
	else b=-1;

	return(b);
}
/*
//====================================================================================================
// Iterative improvement routine to improve the solution vector x[1..n]. Not currently implemented but
// and code not fully working but may be written correctly and used in the future if channel solver
// has trouble converging. Will get round to this at some point. (TJF)
//----------------------------------------------------------------------------
void mprove(float **a, float **alud, int n, int indx[], float b[], float x[])
// Improves a solution vector x[1..n] of the linear set of equations A.x=b. 
// The matrix a[1..n][1..n] and the vectors b[1..n] and x[1..n] are input as 
// is the dimension n, Also input is alud[1..n][1..n], the LU decomposoition of
// a as returned by bandec, and the vecotr indx[1..n] also returned by that routine.
// On output, only x[1..n] is modified to an improved set of values.
{
	int j, i;
	double sdp;
	float *r;

	r=vector(1,n); // needs rewriting
	for(i=0;i<2*n;i++)
	{
		sdp=-b[i];
		for(j=0;j<2*n;j++) sdp+=a[i][j]*x[j];
		r[i]=sdp;
	}
	banbks(alud,n,indx,r);
	for(i=0;i<2*n;i++) x[i]-=r[i];
	delete r;
}
*/