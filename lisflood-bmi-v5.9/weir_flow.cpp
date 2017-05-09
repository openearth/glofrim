/*
*****************************************************************************
WEIR and BRIDGE FLOW
---------------------

A series of function used to calculate the flow between floodplain cells based
on on free flow or drowned flow weir equations (implemented by Rich Dawson,
one-directional flow added by Matt Wilson, sub-grid model flow added by Jeff Neal).
Bridge orifice flow calculation added by Jeff Neal and Mark Trigg.

*****************************************************************************
*/

#include "lisflood.h"

//-----------------------------------------------------------------------------------
// CALCULATE WIER FLOW BETWEEN A POINT AND ITS W NEIGHBOUR
double CalcWeirQx(int i,int j,Pars *Parptr,Arrays *Arrptr, Solver *Solverptr, States * Statesptr, SGCprams * SGCptr)
{
  double z0,z1,h0,h1,Q,hu,hd;
  double usVel; // MT upstream velocity for energy gradient height calc.
  double heg; // MT upstream energy gradient height
  int p0,p1,pq0,weir_id;

	double Qoc; // open channel flow
	double Qp; // pressure(orifice) flow
	double Tz; // transit zone
	double Cd; // bridge Cd
	double Soffit; // bridge soffit elevation
	double Area; // bridge open area - precaclulated in input
	double g;  // gravity accel
	double Z; //bridge opening
	double Zratio; // flow depth to bridge opening ratio
	double Width; // bridge width
	double dt; // timestep (for open channel flow calc)
	double cn; // mannings (for open channel flow calc)
	double dh; // flow head change (for open channel flow calc)
	double Sf; // friction slope (for open channel flow calc)
	double hflow; // flow depth
	double A; // open channel flow area
	double R; // open channel hydraulic radius
	int pQus; // pointer to upstream q flux boundary (between cells)
	int pQds; // pointer to downstream q flux boundary (between cells)

	g = Solverptr->g;

  p0=i+j*Parptr->xsz;
  p1=i+1+j*Parptr->xsz;
  pq0=i+j*(Parptr->xsz+1)+1;

  pQus=i+j*(Parptr->xsz+1);
  pQds=i+j*(Parptr->xsz+1)+2;

  z0=Arrptr->DEM[p0];
  z1=Arrptr->DEM[p1];
  h0=Arrptr->H[p0];
  h1=Arrptr->H[p1];
  weir_id=Arrptr->Weir_Identx[i+1+j*(Parptr->xsz+1)];

  if (Statesptr->SGC==ON) // reset z to subgrid channels if used
  {
	  if (Arrptr->SGCwidth[p0] > 0.0) z0=Arrptr->SGCz[p0];
	  if (Arrptr->SGCwidth[p1] > 0.0) z1=Arrptr->SGCz[p1];
  }

  //Weir equation implementation altered by Rich Dawson 12 Jan 2004/3 Feb 2004.
  //One-directional flow functionality added by Matt Wilson 13 Feb 2004.
  //Bridge equation added by Jeff Neal 1 Sep 2012.

  Q=0.0;

  if (Arrptr->Weir_Typ[weir_id] == 0) // simulate a weir
  {
	if((z0+h0)>(z1+h1))		// Flow in +x direction
	{
		if ((h0+z0)>Arrptr->Weir_hc[weir_id]  && h0>0) // check depth is above weir and that the cell is wet
		{
		    if(Arrptr->Weir_Fixdir[weir_id]==0 || Arrptr->Weir_Fixdir[weir_id]==2) // check for one-directional flow (culvert)
			{
				hu=h0+z0-Arrptr->Weir_hc[weir_id]; // upstream head
				hd=h1+z1-Arrptr->Weir_hc[weir_id]; // downstream head
				if((hd/hu)<Arrptr->Weir_m[weir_id]) Q=Arrptr->Weir_Cd[weir_id]*Arrptr->Weir_w[weir_id]*pow(hu,(1.5)); // Free flow
				else								Q=Arrptr->Weir_Cd[weir_id]*Arrptr->Weir_w[weir_id]*hu*(sqrt(hu-hd))/sqrt(Arrptr->Weir_m[weir_id]); // Drowned flow
			}
		}
	}

	if((z0+h0)<(z1+h1))		// Flow in -x direction
	{
		if ((h1+z1)>Arrptr->Weir_hc[weir_id] && h1>0) // check depth is above weir and that the cell is wet
		{
			if(Arrptr->Weir_Fixdir[weir_id]==0 || Arrptr->Weir_Fixdir[weir_id]==4) // check for one-directional flow (culvert)
			{
				hu=h1+z1-Arrptr->Weir_hc[weir_id]; // upstream head
				hd=h0+z0-Arrptr->Weir_hc[weir_id]; // downstream head
				if((hd/hu)<Arrptr->Weir_m[weir_id]) Q=-Arrptr->Weir_Cd[weir_id]*Arrptr->Weir_w[weir_id]*pow(hu,(1.5)); // Free flow
				else								Q=-Arrptr->Weir_Cd[weir_id]*Arrptr->Weir_w[weir_id]*hu*(sqrt(hu-hd))/sqrt(Arrptr->Weir_m[weir_id]); // Drowned flow
			}
		}
	}
  }
  else if (Arrptr->Weir_Typ[weir_id] == 1) // simulate a bridge
  {
	  hu=h0+z0;
	  hd=h1+z1;

    // get basic bridge params from global arrays;
		Tz = Arrptr->Weir_m[weir_id];
		Soffit = Arrptr->Weir_hc[weir_id];
		Cd = Arrptr->Weir_Cd[weir_id];
		Width = Arrptr->Weir_w[weir_id];

		// calculate some more bridge parameters from basic ones
		Z = min(Soffit-z1,Soffit-z0); // bridge opening (smallest opening)
		Area = Width*Z; // bridge flow area (again smallest opening)

		// get some basic paramters for the open channel flow calc
		dt=Solverptr->Tstep;
		cn = 0.5* (SGCptr->SGCn[Arrptr->SGCgroup[p0]] + SGCptr->SGCn[Arrptr->SGCgroup[p1]]); // mean mannings (note n2)
		dh=z0+h0-z1-h1; // difference in water level
		//Sf=-dh/Parptr->dx; //CCS_deletion
		Sf=-dh/Arrptr->dx[p0]; // CCS added for subgrid lat long compatibility.
		hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1); // calculate the max flow depth
		A = Width*hflow; // calc open channel flow area
		R = A / (Width+2*hflow); // calc hydraulic radius from open channel flow


		// calculate open channel flow using SGC flow
		Qoc=(Arrptr->QxSGold[pq0]-g*A*dt*Sf) / (1+dt*g*cn*fabs(Arrptr->QxSGold[pq0]) / (pow(R,4.0/3.0)*A)); // note mannings is squared by preprocessor
		
		// orifice bridge flow
		if (hu > hd) { // Positive flow
			Zratio = h0/Z; // calc current flow depth to bridge opening ratio
			usVel = Arrptr->QxSGold[pQus] / (h0*Arrptr->SGCwidth[p0]); // MT calculate upstream velocity
			heg = (usVel*usVel)/(2*g); // MT calculate upstream energy gradient height
			Qp = Cd*Area*sqrt(2*g*(hu-hd+heg)); // calc bridge orifice flow
		}
		else { // Negative flow
			Zratio = h1/Z; // calc current flow depth to bridge opening ratio
			usVel = Arrptr->QxSGold[pQds] / (h1*Arrptr->SGCwidth[p1]); // MT calculate upstream velocity
			heg = (usVel*usVel)/(2*g); // MT calculate upstream energy gradient height
			Qp = -(Cd*Area*sqrt(2*g*(hd-hu+heg))); // calc bridge orifice flow
		}


	  if (hu < Soffit && hd < Soffit)
	  {
		  // flow is below the soffit so use SGC open channel flow 
		  Q = Qoc;
	  }
	  else if(Zratio >= 1.0 && Zratio <= Tz) // transition flow between open and orifice/pressure
	  {
		  Q = (Qoc*(Tz-Zratio)/(Tz-1.0))+(Qp*(Zratio-1.0)/(Tz-1.0));
	  }
	  else if(Zratio > Tz) // pressure flow
	  {
		  Q = Qp;
	  }
	  else  // other flow - should not happen, but put here to catch other cases?
	  {
		  printf("WARNING: Unexpected Bridge flow calc fail at t=%.3f , Soffit %f m.\n",Solverptr->t,Soffit); //Warning for fail
		  Q = Qoc;
	  }

	}
  
	// update flow array
  if (Statesptr->SGC==ON) Arrptr->QxSGold[pq0] = Q;

  return(Q);
}
//-----------------------------------------------------------------------------------
// CALCULATE WIER FLOW BETWEEN A POINT AND ITS S NEIGHBOUR
double CalcWeirQy(int i,int j,Pars *Parptr,Arrays *Arrptr, Solver *Solverptr, States * Statesptr, SGCprams * SGCptr)
{
  double z0,z1,h0,h1,Q,hu,hd;
  double usVel; // MT upstream velocity for energy gradient height calc.
  double heg; // MT upstream energy gradient height
  int p0,p1,pq0,weir_id;

	double Qoc; // open channel flow
	double Qp; // pressure(orifice) flow
	double Tz; // transit zone
	double Cd; // bridge Cd
	double Soffit; // bridge soffit elevation
	double Area; // bridge open area - precaclulated in input
	double g;  // gravity accel
	double Z; //bridge opening
	double Zratio; // flow depth to bridge opening ratio
	double Width; // bridge width
	double dt; // timestep (for open channel flow calc)
	double cn; // mannings (for open channel flow calc)
	double dh; // flow head change (for open channel flow calc)
	double Sf; // friction slope (for open channel flow calc)
	double hflow; // flow depth
	double A; // open channel flow area
	double R; // open channel hydraulic radius
	int pQus; // pointer to upstream q flux boundary (between cells)
	int pQds; // pointer to downstream q flux boundary (between cells)

	g = Solverptr->g;

  pq0=i+(j+1)*(Parptr->xsz+1);
	pQus=i+j*(Parptr->xsz+1);
	pQds=i+(j+2)*(Parptr->xsz+1);

  p0=i+j*Parptr->xsz;
  p1=i+(j+1)*Parptr->xsz;
  z0=Arrptr->DEM[p0];
  z1=Arrptr->DEM[p1];
  h0=Arrptr->H[p0];
  h1=Arrptr->H[p1];
  weir_id=Arrptr->Weir_Identy[i+(j+1)*(Parptr->xsz+1)];

  if (Statesptr->SGC==ON)
  {
	  if (Arrptr->SGCwidth[p0] > 0.0) z0=Arrptr->SGCz[p0];
	  if (Arrptr->SGCwidth[p1] > 0.0) z1=Arrptr->SGCz[p1];
  }
  //Weir equation implementation altered by Rich Dawson 12 Jan 2004/3 Feb 2004.
  //One-directional flow functionality added by Matt Wilson 13 Feb 2004.

  Q=0.0;
  if (Arrptr->Weir_Typ[weir_id] == 0)
  {
    if((z0+h0)>(z1+h1))		// Flow in +y direction
    {
      if ((h0+z0)>Arrptr->Weir_hc[weir_id] && h0>0) // check depth is above weir and that the cell is wet
      {
		  if(Arrptr->Weir_Fixdir[weir_id]==0 || Arrptr->Weir_Fixdir[weir_id]==3)  // check for one-directional flow (culvert)
		  {
			hu=h0+z0-Arrptr->Weir_hc[weir_id]; // upstream head
			hd=h1+z1-Arrptr->Weir_hc[weir_id]; // downstream head
			if((hd/hu)<Arrptr->Weir_m[weir_id]) Q=Arrptr->Weir_Cd[weir_id]*Arrptr->Weir_w[weir_id]*pow(hu,(1.5)); // Free flow
			else								Q=Arrptr->Weir_Cd[weir_id]*Arrptr->Weir_w[weir_id]*(hu)*(sqrt(hu-hd))/sqrt(Arrptr->Weir_m[weir_id]); // Drowned flow
		  }
      }
    }
    if((z0+h0)<(z1+h1))			// Flow in -y direction
    {
      if ((h1+z1)>Arrptr->Weir_hc[weir_id] && h1>0) // check depth is above weir and that the cell is wet
      {
		  if(Arrptr->Weir_Fixdir[weir_id]==0 || Arrptr->Weir_Fixdir[weir_id]==1)  // check for one-directional flow (culvert)
		  {
			hu=h1+z1-Arrptr->Weir_hc[weir_id]; // upstram head
			hd=h0+z0-Arrptr->Weir_hc[weir_id]; // downstream head
			
			if((hd/hu)<Arrptr->Weir_m[weir_id]) Q=-Arrptr->Weir_Cd[weir_id]*Arrptr->Weir_w[weir_id]*pow(hu,(1.5)); // Free flow
			else								Q=-Arrptr->Weir_Cd[weir_id]*Arrptr->Weir_w[weir_id]*(hu)*(sqrt(hu-hd))/sqrt(Arrptr->Weir_m[weir_id]); // Drowned flow
		  }
      }
    }
  }
  else if (Arrptr->Weir_Typ[weir_id] == 1 && Statesptr->SGC==ON) // simulate a bridge
  {

	  hu=h0+z0;
	  hd=h1+z1;

    // get basic bridge params from global arrays;
		Tz = Arrptr->Weir_m[weir_id];
		Soffit = Arrptr->Weir_hc[weir_id];
		Cd = Arrptr->Weir_Cd[weir_id];
		Width = Arrptr->Weir_w[weir_id];

		// calculate some more bridge parameters from basic ones
		Z = min(Soffit-z1,Soffit-z0); // bridge opening (smallest opening)
		Area = Width*Z; // bridge flow area (again smallest opening)

		// get some basic paramters for the open channel flow calc
	    dt=Solverptr->Tstep;
	    cn = 0.5* (SGCptr->SGCn[Arrptr->SGCgroup[p0]] + SGCptr->SGCn[Arrptr->SGCgroup[p1]]); // mean mannings (note n2)
		dh=z0+h0-z1-h1; // difference in water level
		//Sf=-dh/Parptr->dx; //CCS_deletion
		Sf=-dh/Arrptr->dy[p0]; // CCS added for sub grid lat-long compatibility.
		hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1); // calculate the max flow depth
		A = Width*hflow; // calc open channel flow area
		R = A / (Width+2*hflow); // calc hydraulic radius from open channel flow


		// calculate open channel flow using SGC flow
		Qoc=(Arrptr->QySGold[pq0]-g*A*dt*Sf) / (1+dt*g*cn*fabs(Arrptr->QySGold[pq0]) / (pow(R,4.0/3.0)*A)); // note mannings is squared
		
		// orifice bridge flow
		if (hu > hd) { // Positive flow
			Zratio = h0/Z; // calc current flow depth to bridge opening ratio
			usVel = Arrptr->QySGold[pQus] / (h0*Arrptr->SGCwidth[p0]); // MT calculate upstream velocity
			heg = (usVel*usVel)/(2*g); // MT calculate upstream energy gradient height
			Qp = Cd*Area*sqrt(2*g*(hu-hd+heg)); // calc bridge orifice flow
		}
		else { // Negative flow
			Zratio = h1/Z; // calc current flow depth to bridge opening ratio
			usVel = Arrptr->QySGold[pQds] / (h1*Arrptr->SGCwidth[p1]); // MT calculate upstream velocity
			heg = (usVel*usVel)/(2*g); // MT calculate upstream energy gradient height
			Qp = -(Cd*Area*sqrt(2*g*(hd-hu+heg))); // calc bridge orifice flow
		}


	  if (hu < Soffit && hd < Soffit)
	  {
		  // flow is below the soffit so use SGC open channel flow 
		  Q = Qoc;
	  }
	  else if(Zratio >= 1.0 && Zratio <= Tz) // transition flow between open and orifice/pressure
	  {
		  Q = (Qoc*(Tz-Zratio)/(Tz-1.0))+(Qp*(Zratio-1.0)/(Tz-1.0));
	  }
	  else if(Zratio > Tz) // pressure flow
	  {
		  Q = Qp;
	  }
	  else  // other flow - should not happen, but put here to catch other cases?
	  {
		  printf("WARNING: Unexpected Bridge flow calc fail at t=%.3f , Soffit %f m.\n",Solverptr->t,Soffit); //Warning for fail
		  Q = Qoc;
	  }

	}
  
	// update flow array
  if (Statesptr->SGC==ON) Arrptr->QxSGold[pq0] = Q;

  return(Q);
}
