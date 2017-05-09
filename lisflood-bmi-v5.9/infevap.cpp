/*
*****************************************************************************
INFILTRATION AND EVAPORATION
---------------------


*****************************************************************************
*/

#include "lisflood.h"

//-----------------------------------------------------------------------------
// FLOODPLAIN INFILTRATION, MDW
void FPInfiltration(Pars *Parptr, Solver *Solverptr, Arrays *Arrptr)
{
  int i,j;
  double cell_inf,h0;

  // Calculate Infiltration
  for(j=0;j<Parptr->ysz;j++) 
  {
    for(i=0;i<Parptr->xsz;i++) 
    {
      if(Arrptr->ChanMask[i+j*Parptr->xsz]==-1) { //only on non-channel cells
        h0=Arrptr->H[i+j*Parptr->xsz];
        if(h0>Solverptr->DepthThresh) 
        {
          cell_inf = Parptr->InfilRate * Solverptr->Tstep; //rate for depth, not area
          h0-=cell_inf;

          //check for -ve depths
          if(h0<0) 
          {
            cell_inf+=h0;
            h0=0;
          }
          Arrptr->H[i+j*Parptr->xsz]=h0;
          //for mass-balance
          Parptr->InfilTotalLoss+=cell_inf*Parptr->dA;
        }
      }
    }
  }

  return;
}


//-----------------------------------------------------------------------------
// FLOODPLAIN EVAPORATION, MDW
void Evaporation(Pars *Parptr, Solver *Solverptr, Arrays *Arrptr)
{
  int i,j;
  double cell_evap,h0,evap_rate;

  evap_rate=InterpBC(Arrptr->evap,Solverptr->t);//constant rate across whole floodplain
  cell_evap = evap_rate * Solverptr->Tstep; //rate for depth, not area

  // Calculate Evaporation
  for(j=0;j<Parptr->ysz;j++) 
  {
    for(i=0;i<Parptr->xsz;i++) 
    {
      if(Arrptr->ChanMask[i+j*Parptr->xsz]==-1) 
      { //only on non-channel cells
        h0=Arrptr->H[i+j*Parptr->xsz];
        if(h0>Solverptr->DepthThresh) 
        {    
          h0-=cell_evap;

          //check for -ve depths 
          if(h0<0) 
          {
            cell_evap+=h0;
            h0=0;
          }
          Arrptr->H[i+j*Parptr->xsz]=h0;
          //for mass-balance
          Parptr->EvapTotalLoss+=cell_evap*Parptr->dA;
        }
      }
    }
  }

  return;
}
//-----------------------------------------------------------------------------
// PRECIPITATION INPUT, TJF
// Adds rainfall (only used when routing scheme is disabled; if routing is enabled, rainfall
// is added in Routing function)
void Rainfall(Pars *Parptr, Solver *Solverptr, Arrays *Arrptr)
{
  int i,j, p0;
  double cell_rain,h0,rain_rate,loc_rainfall_total=0.0;

  rain_rate=InterpBC(Arrptr->rain,Solverptr->t);//constant rate across whole floodplain

  // Add rainfall depth to water heights
#pragma omp parallel for private(i, p0, h0, cell_rain) reduction (+:loc_rainfall_total)  // Parallelised CCS May 2013
  for(j=0;j<Parptr->ysz;j++) 
  {
    for(i=0;i<Parptr->xsz;i++) 
    {
		p0=i+j*Parptr->xsz; // location of cell
		cell_rain=rain_rate * Solverptr->Tstep; //rate for depth, not area.  Has to be inside loop due to negative rain legacy support.

		if(Arrptr->DEM[p0]!=1e10)
		{
			h0=Arrptr->H[p0];
			h0+=cell_rain;
			
			//check for -ve depths (legacy support for pre-evap code when rainfall could be negative)
			if(h0<0) 
			{
				cell_rain+=h0;
				h0=0;
			}
			Arrptr->H[p0]=h0;			
			loc_rainfall_total+=cell_rain*Parptr->dA; // mass balance for local cell (cumulative)		
		}
	}
  }
  Parptr->RainTotalLoss+=loc_rainfall_total; // Update domain mass balance
  return;
}
//-----------------------------------------------------------------------------
// SIMPLE FLOW DIRECTION MAP BASED ON DEM, CCS March 2012; Revised May 2013
/* The function also calculates how often the routing function is run to 
obtain the default or user specified flow routing speed*/
void FlowDirDEM(Pars *Parptr, Arrays *Arrptr, States *Statesptr, BoundCs *BCptr)
{
	int i,j,p0,numBCs,BCi,pi;
	double north,south,east,west,dx,dy,Routing_Speed;

	//if (Statesptr->latlong==ON)
	//{	
		// create RouteInt array for LatLong grids
		Arrptr->RouteInt=new double[Parptr->xsz*Parptr->ysz](); // always create distributed RouteInt JCN
	//}
	//else
	//{
		// else fixed grid size so we can calculate a single value:
	//	Parptr->RouteInt=Parptr->dx/Parptr->Routing_Speed; JCN removed as option because of both latlong and dist_routing need distributed RoutInt
	//}
		
	/* Determine lowest neighbour and assign DEM index to flowdir array.*/
	for(j=0;j<Parptr->ysz;j++) 
	{
		for(i=0;i<Parptr->xsz;i++) 
		{
			p0=i+j*Parptr->xsz; // location of cell

			if(Arrptr->DEM[p0]!=1e10)
			{
				north=1e10;east=1e10;south=1e10;west=1e10;
				if((j-1)>=0)
				{
					north=Arrptr->DEM[i+(j-1)*Parptr->xsz];
				}
				if(i+1<=Parptr->xsz-1)
				{
					east=Arrptr->DEM[(i+1)+j*Parptr->xsz];
				}
				if(j+1<=Parptr->ysz-1)
				{
					south=Arrptr->DEM[i+(j+1)*Parptr->xsz];
				}
				if((i-1)>=0)
				{
					west=Arrptr->DEM[(i-1)+j*Parptr->xsz];
				}

				// calculate dx and dy for rectangular grid or latlong
				if (Statesptr->latlong==ON) {dx = Arrptr->dx[p0]; dy = Arrptr->dy[p0];}
				else {dx = Parptr->dx; dy = Parptr->dx;}

				if(north<=east && north<=south && north<=west && north<Arrptr->DEM[p0]) // flow to the north
				{
					Arrptr->FlowDir[p0]=i+(j-1)*Parptr->xsz; 
					if (Statesptr->dist_routing==OFF) Routing_Speed = Parptr->Routing_Speed; // If distributed routing is OFF use routing speed provide 
					else Routing_Speed = Parptr->Routing_Speed*sqrt((Arrptr->DEM[p0]-north)/dy); // If distributed routing is ON use slope dependent routing speed
					Arrptr->RouteInt[p0]=dy/Routing_Speed; // calculate RouteInt
				}
				else if(east<=north && east<=south && east<=west && east<Arrptr->DEM[p0]) // flow to east
				{
					Arrptr->FlowDir[p0]=(i+1)+j*Parptr->xsz;
					if (Statesptr->dist_routing==OFF) Routing_Speed = Parptr->Routing_Speed; // If distributed routing is OFF use routing speed provide 
					else Routing_Speed = Parptr->Routing_Speed*sqrt((Arrptr->DEM[p0]-east)/dx); // If distributed routing is ON use slope dependent routing speed
					Arrptr->RouteInt[p0]=dx/Routing_Speed; // calculate RouteInt
				}
				else if(south<=north && south<=east && south<=west && south<Arrptr->DEM[p0]) // flow to south
				{
					Arrptr->FlowDir[p0]=i+(j+1)*Parptr->xsz;
					if (Statesptr->dist_routing==OFF) Routing_Speed = Parptr->Routing_Speed; // If distributed routing is OFF use routing speed provide 
					else Routing_Speed = Parptr->Routing_Speed*sqrt((Arrptr->DEM[p0]-south)/dy); // If distributed routing is ON use slope dependent routing speed
					Arrptr->RouteInt[p0]=dy/Routing_Speed; // calculate RouteInt
				}
				else if(west<=north && west<=east && west<=south && west<Arrptr->DEM[p0]) // flow to west
				{
					Arrptr->FlowDir[p0]=(i-1)+j*Parptr->xsz;
					if (Statesptr->dist_routing==OFF) Routing_Speed = Parptr->Routing_Speed; // If distributed routing is OFF use routing speed provide 
					else Routing_Speed = Parptr->Routing_Speed*sqrt((Arrptr->DEM[p0]-west)/dx); // If distributed routing is ON use slope dependent routing speed
					Arrptr->RouteInt[p0]=dx/Routing_Speed; // calculate RouteInt
				}
				
				/* If the above four statements are all false then cell elevation is lower than all neighbours and 
				no routing should occur.  This is achieved by identifying the lowest neighbour as itself in FlowDir.
				By doing this, the routing scheme will compare the cell with itself and detect no difference in water elevation,
				therefore calculating the volume of water to be routed as zero.*/
				else Arrptr->FlowDir[p0]=p0; 

				/* If the subgrid model is ON then we want to disable routing in cells that contain subgrid
				channels.  We do this as above by setting such cells as their own flow direction targets:*/
				if (Statesptr->SGC==ON)
				{
					if (Arrptr->SGCwidth[p0]>0.0)
					{
						Arrptr->FlowDir[p0]=p0; 
					}
				}
			}
		}
	}

	/* Routing needs to be disabled on all non-closed boundary cells.  This avoids routing a volume
	of water that is also operated on by a boundary condition within a single time step, as this
	can result in a mass balance error*/	
	numBCs=2*Parptr->xsz+2*Parptr->ysz;
	for(BCi=0;BCi<numBCs;BCi++) 
	{
		if (BCptr->BC_Ident[BCi]!=0) // if BCi = 0 do nothing
		{
			if(BCi<=Parptr->xsz-1) // N (j=0) edge
			{	
				p0=BCi;
			}
			else if(BCi>=Parptr->xsz && BCi<=Parptr->xsz+Parptr->ysz-1) // E edge
			{  
				p0=Parptr->xsz-1+(BCi-Parptr->xsz)*Parptr->xsz;
			}
			else if(BCi>=Parptr->xsz+Parptr->ysz && BCi<=2*Parptr->xsz+Parptr->ysz-1) // S(j=ysz-1) edge
			{  
				p0=2*Parptr->xsz+Parptr->ysz-1-BCi+(Parptr->ysz-1)*Parptr->xsz;
			} 
			else // W edge
			{   
				p0=0+(numBCs-1-BCi)*Parptr->xsz;
			}
			Arrptr->FlowDir[p0]=p0; // Set recipient cell to be source cell where boundary is present
		}
	}

	// Same principle as above, but now for point source boundary conditions:
	for(pi=0;pi<BCptr->numPS;pi++)
	{
	  p0 = BCptr->xpi[pi]+BCptr->ypi[pi]*Parptr->xsz;
	  Arrptr->FlowDir[p0]=p0; // Set recipient cell to be source cell where boundary is present
	}

	/* Allocate memory for the Route_dH array.  This array is used to temporarily store the changes to
	H for each cell while the routing loop is being calculated.  H is then updated at end of scheme.  This
	is needed to prevent flow from being routed over more than one cell per call to routing scheme.  Not needed
	for subgrid model as SGCdVol is changed rather than H.*/
	if (Statesptr->routing==ON && Statesptr->SGC==OFF)
	{
		Arrptr->Route_dH=new double[Parptr->xsz*Parptr->ysz]();
		return;
	}
}
//-----------------------------------------------------------------------------
//ROUTE SHALLOW FLOWS FROM RAINFALL ACCORDING TO FLOW DIRECTION MAP, CCS March 2012; Revised May 2013	
void Routing(States *Statesptr, Pars *Parptr, Solver *Solverptr, Arrays *Arrptr)
{
	int i,j,p0,p1;
	double cell_rain,rain_rate,h0,h1,z0,z1,flow, flow_fraction;
	double dh,Sfx,Sfy;
	
	if(Statesptr->rainfall==ON) // calc rainfall dH for this timestep
	{
		rain_rate=InterpBC(Arrptr->rain,Solverptr->t);//constant rate across whole floodplain
		cell_rain = rain_rate * Solverptr->Tstep; //rate for depth, not area
	}
		
	//flow_fraction=Solverptr->Tstep/Parptr->RouteInt; // fraction of cell volume to route in this time step // JCN Not used now that RouteInt is always dirstibuted
	//if(flow_fraction>1) flow_fraction=1; // prevent flow fraction>1 (should never happen!) 

	//  Main routing loop
	for(j=0;j<Parptr->ysz;j++) 
	{
		for(i=0;i<Parptr->xsz;i++) 
		{		
			p0=i+j*Parptr->xsz;

			if(Arrptr->DEM[p0]!=1e10)
			{
				if(Statesptr->rainfall==ON) // Add rainfall if on
				{
					Arrptr->H[p0]+=cell_rain; // Add cell rainfall
					Parptr->RainTotalLoss+=cell_rain*Parptr->dA; // Mass balance
				}

				if(Arrptr->H[p0]>0) // Only proceed with rest of loop if cell is wet:
				{
					z0=Arrptr->DEM[p0];
					h0=Arrptr->H[p0];

					//Calculate friction slopes to assess whether routing scheme is needed to force stability:
					if(i==Parptr->xsz-1) Sfx=0.0; // boundary cell, don't calc Sfx:
					else
					{
						p1=i+1+j*Parptr->xsz;
						z1=Arrptr->DEM[p1]; 
						h1=Arrptr->H[p1];
						dh=z0+h0-z1-h1; 
						Sfx=fabs(dh/Parptr->dx); // friction slope in x
					}

					if(j==Parptr->ysz-1) Sfy=0.0; // boundary cell, don't calc Sfy:
					else
					{
						p1=i+(j+1)*Parptr->xsz; 
						z1=Arrptr->DEM[p1]; 
						h1=Arrptr->H[p1];
						dh=z0+h0-z1-h1; 
						Sfy=fabs(dh/Parptr->dx); // friction slope in y
					}

				
				
					h0=Arrptr->H[p0];//cell water height

					// Perform routing
					if(h0<Solverptr->DepthThresh || Sfy>=Parptr->RouteSfThresh || Sfx>=Parptr->RouteSfThresh)
					{
						h1=Arrptr->H[Arrptr->FlowDir[p0]]; //lowest neighbour cell water height 
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

						flow=flow*flow_fraction; // calc correct H to route this timestep

						Arrptr->Route_dH[p0]-=flow; // Record change in H in source cell 
						Arrptr->Route_dH[Arrptr->FlowDir[p0]]+=flow; // Record change in H in recipient cell
					}
							
				}
			}
		}
	}

	// Update H and reset Route_dH:
	for(j=0;j<Parptr->ysz;j++) 
	{
		for(i=0;i<Parptr->xsz;i++) 
		{
			p0=i+j*Parptr->xsz;
			Arrptr->H[p0]+=Arrptr->Route_dH[p0];
			Arrptr->Route_dH[p0]=0.0;
		}
	}
}