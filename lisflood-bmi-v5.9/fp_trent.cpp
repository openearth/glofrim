/*
*****************************************************************************
FLOODPLAIN FLOW WITH ROE APPROXIMATE RIEMANN SOLVER 1st ORDER SCHEME
---------------------------------

Calculate flow between floodplain cells

*****************************************************************************
*/

#include "lisflood.h"

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------------

// CALCULATE FLOODPLAIN (ie NON WIER) FLOW BETWEEN A POINT AND ITS W NEIGHBOUR USING ACCELERATION FORMULATION
double CalcFPQxRoe(int i,int j,States *Statesptr,Pars *Parptr,Solver *Solverptr,Arrays *Arrptr)
{
  double z0,z1,h0,h1,hflow,dh=0.0,Q,g,dt;
  int p0,p1,pq0;  

  double dl,hl,hr,s0,hul,hur,hvl,hvr;
  double dtol;
  double shl,shr,aux1,ul,ur,vl,vr,ubarra,vbarra,cbarra,a1,a2,a3,a1m,a2m,a3m,a1l,a1r,a3l,a3r,epsilon,dhu,dhv; 	
  double alfa1,alfa2,alfa3,e11,e12,e13,e21,e22,e23,e31,e32,e33,f1lp,f2lp,f3lp,f1rp,f2rp,f3rp,fctm,wsurface;//beta1,beta2,beta3; CCS threw unreferenced local variable warning so I guess no longer needed?
  double maximum(double a,double b,double c);

    // Indexing	
    p0=i+j*Parptr->xsz;
    p1=i+1+j*Parptr->xsz;
    pq0=i+1+j*(Parptr->xsz+1); // staggered grid for discharge
	
    // Predefined-variables	
    g=Solverptr->g;
    dt=Solverptr->Tstep;
    dl=Parptr->dx ;  	
    dtol=Solverptr->DepthThresh; // same as LF-FP 

    // DEMs related
	z0=Arrptr->DEM[p0];
    z1=Arrptr->DEM[p1];

	// Slope (why not declared at the beginning ?)
    // Roe solver State-variables	
	//Index equivalence
	//ir=il+1;
	//hl=h[il][jl];
	//hr=h[ir][jl];

	h0=hl=Arrptr->H[p0];
  	h1=hr=Arrptr->H[p1];

	hul=Arrptr->HU[p0] ;
	hur=Arrptr->HU[p1] ;

	hvl=Arrptr->HV[p0] ;
	hvr=Arrptr->HV[p1] ;
	
    /* Index of variables to be defined in both cells and in the intercell 
   	// Intercell -flows
	Arrptr->FH[pq0]= ;
 	Arrptr->FHU[pq0]= ;
 	Arrptr->FHV[pq0]= ;
	// Momentum sources 	
	Arrptr->RSHU[p0]= ;	
	Arrptr->LSHU[p1]= ;
	Arrptr->RSHV[p0]= ;	
	Arrptr->LSHV[p1]= ;
    */

    // Type of flow cases:
    //----------------------------------------------------------------------------------------	

	if (hl >= dtol && hr >= dtol){	/* both normal */

		       	shr= pow(hr,0.5);
				shl= pow(hl,0.5);
				aux1= shr + shl ;

				ul= hul/hl;
				ur= hur/hr;
				vl= hvl/hl;
				vr= hvr/hr;

				ubarra= (shr*ur+shl*ul)/aux1;
				vbarra= (shr*vr+shl*vl)/aux1;
				cbarra= pow(0.5*g*(hl+hr),0.5);

				a1=ubarra+cbarra ;
				a2=ubarra ;
				a3=ubarra-cbarra ;

				a1m=fabs(a1) ;
				a2m=fabs(a2) ;
				a3m=fabs(a3) ;
						
				a1l=ul+pow(g*hl,0.5) ;
				a1r=ur+pow(g*hr,0.5) ;

				epsilon=maximum(0.,(a1-a1l),(a1r-a1)) ;
				a1m=(a1m >= epsilon) ? a1m : epsilon  ;
						
				a3l=ul-pow(g*hl,0.5) ; 
				a3r=ur-pow(g*hr,0.5) ;

				epsilon=maximum(0.,(a3-a3l),(a3r-a3)) ;
				a3m=(a3m >= epsilon) ? a3m : epsilon  ;
						
				dh=hr-hl ;
				dhu=hur-hul ;
				dhv=hvr-hvl ;

				alfa1=0.5*dh+0.5*(dhu-ubarra*dh)/cbarra ;
				alfa2=(dhv-vbarra*dh)/cbarra ;
				alfa3=0.5*dh-0.5*(dhu-ubarra*dh)/cbarra ;

				e11=1. ;
				e12=ubarra+cbarra;
				e13=vbarra ;

				e21=0. ;
				e22=0. ;
				e23=cbarra;

				e31=1. ;
				e32=ubarra-cbarra;
				e33=vbarra;

				/* flux variables */ 

				f1lp=hul ;
				f2lp=hul*ul+0.5*g*hl*hl ;
				f3lp=hul*vl ;
						
				f1rp=hur ;
				f2rp=hur*ur+0.5*g*hr*hr ;
				f3rp=hur*vr ;		
							
				Arrptr->FHx[pq0]=0.5*(f1rp+f1lp-a1m*alfa1*e11-a2m*alfa2*e21-a3m*alfa3*e31) ;
				Arrptr->FHUx[pq0]=0.5*(f2rp+f2lp-a1m*alfa1*e12-a2m*alfa2*e22-a3m*alfa3*e32) ;
				Arrptr->FHVx[pq0]=0.5*(f3rp+f3lp-a1m*alfa1*e13-a2m*alfa2*e23-a3m*alfa3*e33) ;
						
				s0=-(Arrptr->DEM[p1]-Arrptr->DEM[p0])/dl;

				Arrptr->RSHU[p0]=0.5*g*hl*s0 ;	
				Arrptr->LSHU[p1]=0.5*g*hr*s0 ;
	
				Arrptr->RSHV[p0]=0. ;	
				Arrptr->LSHV[p1]=0. ;								
															
			}
			else if(hl < dtol && hr > dtol){  /* left dry  */	
				
				fctm=1.-pow(1.,0.5)*fabs(hur)/(pow(g,0.5)*pow(hr,1.5)) ;
				
				if (fctm < 0.) fctm=0.;
				
				if (hr*fctm+Arrptr->DEM[p1] > Arrptr->DEM[p0]  && (hur/hr-pow(g*hr,0.5)) <= 0. ){	/* overtopping */
					
					wsurface=(Arrptr->DEM[p0]-hr-Arrptr->DEM[p1])/dl ;

					Arrptr->FHx[pq0]=hur ;
					Arrptr->FHUx[pq0]=hur*hur/hr ;
					Arrptr->FHVx[pq0]=hvr*hur/hr ;

					Arrptr->H[p0]=dtol;
					
					Arrptr->RSHU[p0]=g*hr*wsurface ;	
					Arrptr->LSHU[p1]=g*hr*wsurface ;
					
					Arrptr->RSHV[p0]=0. ;	
					Arrptr->LSHV[p1]=0. ;	
					
				}
				else{					/* wall */
						// To verify that left cell is not updated
						// WDC modifies Arrptr->FHx[pq0]=0 ;
/*
						if (hur <= 0.){
							Arrptr->FHUx[pq0]=hur*hur/hr ;  // remind (-) sign while upgrade
							Arrptr->FHVx[pq0]=0. ;
						}
						else
						{
							Arrptr->FHUx[pq0]=0. ;
							Arrptr->FHVx[pq0]=0. ;
						}

						//Raw formulation from Brufau (Dont add RSH[p0] for avoiding UpdateRoeH)


							ubarra= ur;
							cbarra= pow(0.5*g*hr,0.5);	
								
							if (z0 < z1 + hr*(1.+fabs(ubarra)/cbarra) && fabs(ubarra) < cbarra){
								Arrptr->FHx[pq0]=0.5*(cbarra*(z0-z1)-hr*(fabs(ubarra)+cbarra));
							}
							else{
								Arrptr->FHx[pq0]=0 ;
							}
							// If z1 >> z0  -> great losses
							// If ubarra=0  what happens 
*/

					if (Statesptr->Roe_slow==ON)
					{
							hl = hr;
							hul = -hur;
							hvl = hvr;

							//s0=-(Arrptr->DEM[p1+1]-Arrptr->DEM[p1])/dl;
							//Arrptr->RSHU[p0]=0.0 ;	
							//if (s0 >= 0) Arrptr->LSHU[p1]=0.0;
							//else if (Arrptr->H[p1+1] < dtol) Arrptr->LSHU[p1]=0.0;
							//else Arrptr->LSHU[p1]=0.5*g*hr*s0;
					
							Arrptr->RSHV[p0]=0.0 ;	
							Arrptr->LSHV[p1]=0.0 ;

							//Arrptr->LSHU[p1]=0.0 ;
							//Arrptr->LSHV[p1]=0.0 ;
							Arrptr->RSHU[p0]=0.0 ;//GUSTAVO
							Arrptr->LSHU[p1]=0.0 ;//GUSTAVO

							shr= pow(hr,0.5); shl= pow(hl,0.5); aux1= shr + shl ;
							ul= hul/hl; ur= hur/hr; vl= hvl/hl; vr= hvr/hr;
							ubarra= (shr*ur+shl*ul)/aux1; vbarra= (shr*vr+shl*vl)/aux1; cbarra= pow(0.5*g*(hl+hr),0.5);
							a1=ubarra+cbarra ; a2=ubarra ; a3=ubarra-cbarra ; 
							a1m=fabs(a1) ; a2m=fabs(a2) ; a3m=fabs(a3) ; 
							a1l=ul+pow(g*hl,0.5) ; a1r=ur+pow(g*hr,0.5) ;
							epsilon=maximum(0.,(a1-a1l),(a1r-a1)) ;
							a1m=(a1m >= epsilon) ? a1m : epsilon  ;
							a3l=ul-pow(g*hl,0.5) ; a3r=ur-pow(g*hr,0.5) ;
							epsilon=maximum(0.,(a3-a3l),(a3r-a3)) ; a3m=(a3m >= epsilon) ? a3m : epsilon  ;
							dh=hr-hl ; dhu=hur-hul ; dhv=hvr-hvl ;
							alfa1=0.5*dh+0.5*(dhu-ubarra*dh)/cbarra ;
							alfa2=(dhv-vbarra*dh)/cbarra ;
							alfa3=0.5*dh-0.5*(dhu-ubarra*dh)/cbarra ;
							e11=1. ; e12=ubarra+cbarra; e13=vbarra ;
							e21=0. ; e22=0. ; e23=cbarra;
							e31=1. ; e32=ubarra-cbarra; e33=vbarra;
							/* flux variables */ 
							f1lp=hul ; f2lp=hul*ul+0.5*g*hl*hl ; f3lp=hul*vl ;
							f1rp=hur ; f2rp=hur*ur+0.5*g*hr*hr ; f3rp=hur*vr ;
							Arrptr->FHx[pq0]=0.5*(f1rp+f1lp-a1m*alfa1*e11-a2m*alfa2*e21-a3m*alfa3*e31) ; 
							Arrptr->FHUx[pq0]=0.5*(f2rp+f2lp-a1m*alfa1*e12-a2m*alfa2*e22-a3m*alfa3*e32) ;
							Arrptr->FHVx[pq0]=0.5*(f3rp+f3lp-a1m*alfa1*e13-a2m*alfa2*e23-a3m*alfa3*e33) ;
					}
					else
					{
						                        // FLAT SOLVER

						//if (hur <= 0.){
							Arrptr->FHUx[pq0]=hur*(hur/hr-pow(g*hr,0.5))+0.5*g*hr*hr ;  // remind (-) sign while upgrade
							Arrptr->FHVx[pq0]=0. ;
							Arrptr->FHx[pq0]=0.;//GUSTAVO

							Arrptr->RSHV[p0]=0.0 ;//GUSTAVO
							Arrptr->LSHV[p1]=0.0 ;//GUSTAVO
							Arrptr->RSHU[p0]=0.0 ;//GUSTAVO
							Arrptr->LSHU[p1]=0.0 ;//GUSTAVO
						//}
						//else
						//{
						//	Arrptr->FHUx[pq0]=0. ;
						//	Arrptr->FHVx[pq0]=0. ;
						//}
					}
				}		
			}
			else if(hr < dtol && hl > dtol){	/* right dry */
			
				fctm=1.-pow(1.,0.5)*fabs(hul)/(pow(g,0.5)*pow(hl,1.5)) ;
					
				if (fctm < 0.) fctm=0.;
					
				if (hl*fctm+Arrptr->DEM[p0]  > Arrptr->DEM[p1]  && (hul/hl+pow(g*hl,0.5)) >= 0.){	/* overtopping */
				
					wsurface=(hl+Arrptr->DEM[p0]-Arrptr->DEM[p1])/dl;
					
					Arrptr->FHx[pq0]=hul ;
					Arrptr->FHUx[pq0]=hul*hul/hl ;
					Arrptr->FHVx[pq0]=hvl*hul/hl ;
					
					Arrptr->H[p1]=dtol;
					
					Arrptr->RSHU[p0]=g*hl*wsurface ;	
					Arrptr->LSHU[p1]=g*hl*wsurface ;
					
					Arrptr->RSHV[p0]=0. ;	
					Arrptr->LSHV[p1]=0. ;
													
				}
				else
				{   				/* wall  */  
					
					/* wall at j+1 */
					
						// WDC modifies Arrptr->FHx[pq0]=0 ;		
					if (Statesptr->Roe_slow==ON)
					{
					hr = hl;
					hur = -hul;
					hvr = hvl;

					//Arrptr->RSHU[p0]=0.0 ;
					//if (Arrptr->LSHU[p0] > 0) Arrptr->RSHU[p0]; 
					//else Arrptr->RSHU[p0]= -Arrptr->LSHU[p0] ;
					
					//s0=-(Arrptr->DEM[p1+1]-Arrptr->DEM[p1])/dl;
					//if (Arrptr->LSHU[p0] <= 0) Arrptr->RSHU[p0] = 0.0;
					//else if (Arrptr->H[p0-1] < dtol) Arrptr->RSHU[p0] = 0.0; // if other side is also dry assume zero
					//else Arrptr->RSHU[p0]=Arrptr->LSHU[p0] ;	
					//Arrptr->LSHU[p1]=0.0 ; //0.5*g*hr*s0;
					
					Arrptr->RSHV[p0]=0.0 ;	
					Arrptr->LSHV[p1]=0.0 ;

					//Arrptr->RSHU[p0]=0.0 ;
					//Arrptr->RSHV[p0]=0.0 ;
					Arrptr->RSHU[p0]=0.0 ;//GUSTAVO
					Arrptr->LSHU[p1]=0.0 ;//GUSTAVO

					shr= pow(hr,0.5); shl= pow(hl,0.5); aux1= shr + shl ;
					ul= hul/hl; ur= hur/hr; vl= hvl/hl; vr= hvr/hr;
					ubarra= (shr*ur+shl*ul)/aux1; vbarra= (shr*vr+shl*vl)/aux1; cbarra= pow(0.5*g*(hl+hr),0.5);
					a1=ubarra+cbarra ; a2=ubarra ; a3=ubarra-cbarra ; 
					a1m=fabs(a1) ; a2m=fabs(a2) ; a3m=fabs(a3) ; 
					a1l=ul+pow(g*hl,0.5) ; a1r=ur+pow(g*hr,0.5) ;
					epsilon=maximum(0.,(a1-a1l),(a1r-a1)) ;
					a1m=(a1m >= epsilon) ? a1m : epsilon  ;
					a3l=ul-pow(g*hl,0.5) ; a3r=ur-pow(g*hr,0.5) ;
					epsilon=maximum(0.,(a3-a3l),(a3r-a3)) ; a3m=(a3m >= epsilon) ? a3m : epsilon  ;
					dh=hr-hl ; dhu=hur-hul ; dhv=hvr-hvl ;
					alfa1=0.5*dh+0.5*(dhu-ubarra*dh)/cbarra ;
					alfa2=(dhv-vbarra*dh)/cbarra ;
					alfa3=0.5*dh-0.5*(dhu-ubarra*dh)/cbarra ;
					e11=1. ; e12=ubarra+cbarra; e13=vbarra ;
					e21=0. ; e22=0. ; e23=cbarra;
					e31=1. ; e32=ubarra-cbarra; e33=vbarra;
					/* flux variables */ 
					f1lp=hul ; f2lp=hul*ul+0.5*g*hl*hl ; f3lp=hul*vl ;
					f1rp=hur ; f2rp=hur*ur+0.5*g*hr*hr ; f3rp=hur*vr ;
					Arrptr->FHx[pq0]=0.5*(f1rp+f1lp-a1m*alfa1*e11-a2m*alfa2*e21-a3m*alfa3*e31) ; 
					Arrptr->FHUx[pq0]=0.5*(f2rp+f2lp-a1m*alfa1*e12-a2m*alfa2*e22-a3m*alfa3*e32) ;
					Arrptr->FHVx[pq0]=0.5*(f3rp+f3lp-a1m*alfa1*e13-a2m*alfa2*e23-a3m*alfa3*e33) ;
					}
					else
					{
						/* wall at j+1 */
					
						// WDC modifies Arrptr->FHx[pq0]=0 ;		
						
						//if (hul >= 0.){
						
							Arrptr->FHUx[pq0]=hul*(hul/hl+pow(g*hl,0.5))+0.5*g*hl*hl;
							Arrptr->FHVx[pq0]=0.;
							
							Arrptr->FHx[pq0]=0.;//GUSTAVO

							Arrptr->RSHV[p0]=0.0 ;//GUSTAVO
							Arrptr->LSHV[p1]=0.0 ;//GUSTAVO
							Arrptr->RSHU[p0]=0.0 ;//GUSTAVO
							Arrptr->LSHU[p1]=0.0 ;//GUSTAVO
						//}
						//else{
						
						//	Arrptr->FHUx[pq0]=0.;
						//	Arrptr->FHVx[pq0]=0.;
						
						//}
					}
					/*
						if (hul >= 0.){
						
							Arrptr->FHUx[pq0]=hul*hul/hl;
							Arrptr->FHVx[pq0]=0.;
						}
						else{
						
							Arrptr->FHUx[pq0]=0.;
							Arrptr->FHVx[pq0]=0.;
						
						}
							
						// To verify that right cell is not updated
						
						Arrptr->RSHU[p0]=0 ;
						Arrptr->RSHV[p0]=0 ; 


						// WDC from Brufau et Al    

							
							//Raw formulation from Brufau (Dont add RSH[p0] for avoiding UpdateRoeH)


							ubarra= ul;
							cbarra= pow(0.5*g*hl,0.5);	
								
							if (z1 < z0 + hl*(1.+fabs(ubarra)/cbarra) && fabs(ubarra) < cbarra){

								Arrptr->FHx[pq0]=-0.5*(cbarra*(z1-z0)-hl*(fabs(ubarra)+cbarra));
							}
							else{
								Arrptr->FHx[pq0]=0.;
							}

							// If z1 >> z0  -> great losses
							// If ubarra=0  what happens 


						// 2nd Add all betas to Arrptr->RSHU[p0], Arrptr->RSHV[p0];


							beta1=0.;
							beta2=0.;
							beta3=0.;	*/					
				}	
			}
    //---------------------------------------------------------------------------------------- 	
  Q=Arrptr->FHx[pq0]*dl ; // origional version included addition *dl but removed for compatability with UpdateH JCN
  // option to save V's
  if (Statesptr->voutput==ON)
  {
	  if (Q!=0.0)
	  {
	    hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
        hflow=getmax(hflow,0);
        hflow=getmin(hflow,Solverptr->MaxHflow);
	    // calc Vx
	    Arrptr->Vx[pq0]=Arrptr->FHx[pq0]/hflow;
	    // get max Vx
	    Arrptr->maxVx[pq0]=getmax(Arrptr->maxVx[pq0],fabs(Arrptr->Vx[pq0]));
	  }
	  else Arrptr->Vx[pq0]=0.0;
  }
  return(Q);      	
}

//-----------------------------------------------------------------------------------
// CALCULATE FLOODPLAIN (ie NON WIER) FLOW BETWEEN A POINT AND ITS S NEIGHBOUR USING ACCELERATION FORMULATION
double CalcFPQyRoe(int i,int j,States *Statesptr,Pars *Parptr,Solver *Solverptr,Arrays *Arrptr)
{
  double z0,z1,h0,h1,hflow,dh=0.0,Q=0.0,g,dt;
  int p0,p1,pq0;

  double dl,hl,hr,s0,hul,hur,hvl,hvr;
  double dtol;
  double shl,shr,aux1,ul,ur,vl,vr,ubarra,vbarra,cbarra,a1,a2,a3,a1m,a2m,a3m,a1l,a1r,a3l,a3r,epsilon,dhu,dhv; 	
  double alfa1,alfa2,alfa3,e11,e12,e13,e21,e22,e23,e31,e32,e33,g1lp,g2lp,g3lp,g1rp,g2rp,g3rp,fctm,wsurface;
  double maximum(double a,double b,double c);

  // Indexing	
	
  p0=i+j*Parptr->xsz;
  p1=i+(j+1)*Parptr->xsz;

  pq0=i+(j+1)*(Parptr->xsz+1);

  // Predefined-variables	

  g=Solverptr->g;
  dt=Solverptr->Tstep;
  dl=Parptr->dx ;  	

  dtol=Solverptr->DepthThresh ; // same as LF-FP ?	

    // DEMs related
  z0=Arrptr->DEM[p0];
  z1=Arrptr->DEM[p1];

	// Slope (why not declared at the beginning ?)

    // Roe solver State-variables	

	//Index equivalence
	//jr=jl+1;
	//hl=h[il][jl];
	//hr=h[il][jr];

	h0=hl=Arrptr->H[p0];
  	h1=hr=Arrptr->H[p1];

	hul=Arrptr->HU[p0] ;
	hur=Arrptr->HU[p1] ;

	hvl=Arrptr->HV[p0] ;
	hvr=Arrptr->HV[p1] ;

    /* Index of variables to be defined in both cells and in the intercell 
	
   	// Intercell -flows

	Arrptr->FH[pq0]= ;
 	Arrptr->FHU[pq0]= ;
 	Arrptr->FHV[pq0]= ;

	// Momentum sources 	
 
	Arrptr->BSHU[p0]= ;	
	Arrptr->TSHU[p1]= ;

	Arrptr->BSHV[p0]= ;	
	Arrptr->TSHV[p1]= ;
	
    */

    // Type of flow cases:

    //----------------------------------------------------------------------------------------	

    		if (hl >= dtol && hr >= dtol){	/* both normal */
			       				
				shr= pow(hr,0.5);
				shl= pow(hl,0.5);
				aux1= shr + shl ;

				ul= hul/hl;
				ur= hur/hr;
				vl= hvl/hl;
				vr= hvr/hr;

				ubarra= (shr*ur+shl*ul)/aux1;
				vbarra= (shr*vr+shl*vl)/aux1;
				cbarra= pow(0.5*g*(hl+hr),0.5);

				a1=vbarra+cbarra ;
				a2=vbarra ;
				a3=vbarra-cbarra ;

				a1m=fabs(a1) ;
				a2m=fabs(a2) ;
				a3m=fabs(a3) ;
						
				a1l=vl+pow(g*hl,0.5) ;
				a1r=vr+pow(g*hr,0.5) ;

				epsilon=maximum(0.,(a1-a1l),(a1r-a1)) ;
				a1m=(a1m >= epsilon) ? a1m : epsilon  ;

				a3l=vl-pow(g*hl,0.5) ;
				a3r=vr-pow(g*hr,0.5) ;

				epsilon=maximum(0.,(a3-a3l),(a3r-a3)) ;
				a3m=(a3m >= epsilon) ? a3m : epsilon  ;
						
				dh=hr-hl ;
				dhu=hur-hul ;
				dhv=hvr-hvl ;

				alfa1=0.5*dh+0.5*(dhv-vbarra*dh)/cbarra ;
				alfa2=-(dhu-ubarra*dh)/cbarra ;
				alfa3=0.5*dh-0.5*(dhv-vbarra*dh)/cbarra ;

				e11=1. ;
				e12=ubarra;
				e13=vbarra+cbarra ;

				e21=0. ;
				e22=-cbarra ;
				e23=0.;

				e31=1. ;
				e32=ubarra ;
				e33=vbarra-cbarra ;

				/* flux variables */ 
				g1lp=hvl ;
				g2lp=hul*vl ;
				g3lp=hvl*vl+0.5*g*hl*hl ;

				g1rp=hvr ;
				g2rp=hur*vr ;
				g3rp=hvr*vr+0.5*g*hr*hr ;
	
				Arrptr->FHy[pq0]=0.5*(g1rp+g1lp-a1m*alfa1*e11-a2m*alfa2*e21-a3m*alfa3*e31) ;
				Arrptr->FHUy[pq0]=0.5*(g2rp+g2lp-a1m*alfa1*e12-a2m*alfa2*e22-a3m*alfa3*e32) ;
				Arrptr->FHVy[pq0]=0.5*(g3rp+g3lp-a1m*alfa1*e13-a2m*alfa2*e23-a3m*alfa3*e33) ;

				s0=-(Arrptr->DEM[p1]-Arrptr->DEM[p0])/dl;

				Arrptr->BSHU[p0]=0. ;	
				Arrptr->TSHU[p1]=0. ;
	
				Arrptr->BSHV[p0]=0.5*g*hl*s0 ;	
				Arrptr->TSHV[p1]=0.5*g*hr*s0 ;																				
			}
			else if(hl < dtol && hr > dtol){		        /* left dry  */		

				fctm=1.-pow(1.,0.5)*fabs(hvr)/(pow(g,0.5)*pow(hr,1.5)) ;

				if (fctm < 0.) fctm=0.;

				if (hr*fctm+Arrptr->DEM[p1] > Arrptr->DEM[p0] && (hvr/hr-pow(g*hr,0.5))<= 0.){	/* overtopping */
				
					wsurface=(hl+Arrptr->DEM[p0]-hr-Arrptr->DEM[p1])/dl;	
								
					Arrptr->FHy[pq0]=hvr ;
					Arrptr->FHUy[pq0]=hvr*hur/hr ;
					Arrptr->FHVy[pq0]=hvr*hvr/hr ;
													
					Arrptr->H[p0]=dtol;

					Arrptr->BSHU[p0]=0. ;	
					Arrptr->TSHU[p1]=0. ;

					Arrptr->BSHV[p0]=g*hr*wsurface ;	
					Arrptr->TSHV[p1]=g*hr*wsurface ;

				}
				else
				{					/* wall */

						// To verify that left cell is not updated

						//Arrptr->FHy[pq0]=0 ;
					if (Statesptr->Roe_slow==ON)
					{
					hl = hr;
					hul = hur;
					hvl = -hvr;

					//s0=-(Arrptr->DEM[p1]-Arrptr->DEM[p1+Parptr->xsz])/dl; //work out slope on other side
					//if (s0 < 0) Arrptr->TSHV[p1]=0.0; // if away from boundary use zero
					//else Arrptr->TSHV[p1]=-0.5*g*hr*s0 ; // if towrads boundary balance woth ghost cell
					//s0=-(Arrptr->DEM[p1+Parptr->xsz]-Arrptr->DEM[p1])/dl;

					Arrptr->BSHU[p0]=0.0 ;	
					//Arrptr->TSHU[p1]=0.0 ;
	
					Arrptr->BSHV[p0]=0.0 ;	
					//if (s0 >= 0) Arrptr->TSHV[p1]=0.0; 
					//else if (Arrptr->H[p1+Parptr->xsz] < dtol) Arrptr->TSHV[p1]=0.0;
					//else Arrptr->TSHV[p1]=0.5*g*hr*s0 ;	
					
					Arrptr->TSHU[p1]=0.0 ;
					Arrptr->TSHV[p1]=0.0 ;

					// implement solver for dry edge
					shr= pow(hr,0.5); shl= pow(hl,0.5); aux1= shr + shl ;
					ul= hul/hl; ur= hur/hr; vl= hvl/hl; vr= hvr/hr;
					ubarra= (shr*ur+shl*ul)/aux1; vbarra= (shr*vr+shl*vl)/aux1;	cbarra= pow(0.5*g*(hl+hr),0.5); 
					a1=vbarra+cbarra ; a2=vbarra ; a3=vbarra-cbarra ; 
					a1m=fabs(a1) ; a2m=fabs(a2) ; a3m=fabs(a3) ;
					a1l=vl+pow(g*hl,0.5) ; a1r=vr+pow(g*hr,0.5) ; 
					epsilon=maximum(0.,(a1-a1l),(a1r-a1)) ; a1m=(a1m >= epsilon) ? a1m : epsilon  ; 
					a3l=vl-pow(g*hl,0.5) ; a3r=vr-pow(g*hr,0.5) ;
					epsilon=maximum(0.,(a3-a3l),(a3r-a3)) ;	a3m=(a3m >= epsilon) ? a3m : epsilon  ;
					dh=hr-hl ; dhu=hur-hul ; dhv=hvr-hvl ;

					alfa1=0.5*dh+0.5*(dhv-vbarra*dh)/cbarra ; alfa2=-(dhu-ubarra*dh)/cbarra ; alfa3=0.5*dh-0.5*(dhv-vbarra*dh)/cbarra ; 		
					e11=1. ; e12=ubarra; e13=vbarra+cbarra ; 
					e21=0. ; e22=-cbarra ; e23=0.; 		
					e31=1. ; e32=ubarra ; e33=vbarra-cbarra ; 
					/* flux variables */ 
					g1lp=hvl ; g2lp=hul*vl ; g3lp=hvl*vl+0.5*g*hl*hl ; 
					g1rp=hvr ; g2rp=hur*vr ; g3rp=hvr*vr+0.5*g*hr*hr ;
					Arrptr->FHy[pq0]=0.5*(g1rp+g1lp-a1m*alfa1*e11-a2m*alfa2*e21-a3m*alfa3*e31) ;
					Arrptr->FHUy[pq0]=0.5*(g2rp+g2lp-a1m*alfa1*e12-a2m*alfa2*e22-a3m*alfa3*e32) ;
					Arrptr->FHVy[pq0]=0.5*(g3rp+g3lp-a1m*alfa1*e13-a2m*alfa2*e23-a3m*alfa3*e33) ;
					}
					else
					{
						// To verify that left cell is not updated

						//Arrptr->FHy[pq0]=0 ;

						//if (hvr <= 0){
						
							Arrptr->FHUy[pq0]=0.;
							Arrptr->FHVy[pq0]=hvr*(hvr/hr-pow(g*hr,0.5))+0.5*g*hr*hr; // remind (-) sign while upgrade
							Arrptr->FHy[pq0]=0.;//GUSTAVO

							Arrptr->BSHU[p0]=0.0 ;//GUSTAVO
							Arrptr->BSHV[p0]=0.0 ;//GUSTAVO
							Arrptr->TSHU[p1]=0.0 ;//GUSTAVO
							Arrptr->TSHV[p1]=0.0 ;//GUSTAVO

						//}
						//else{
						
						//	Arrptr->FHUy[pq0]=0.;
						//	Arrptr->FHVy[pq0]=0.;
						
						//}
					}

					/*
						if (hvr <= 0){
						
							Arrptr->FHUy[pq0]=0.;
							Arrptr->FHVy[pq0]=hvr*hvr/hr; // remind (-) sign while upgrade
						
						}
						else{
						
							Arrptr->FHUy[pq0]=0.;
							Arrptr->FHVy[pq0]=0.;	
						
						}
	
						Arrptr->TSHU[p1]=0. ;
						Arrptr->TSHV[p1]=0. ;

						//Raw formulation from Brufau (Dont add RSH[p0] for avoiding UpdateRoeH)


							ubarra= ur;
							cbarra= pow(0.5*g*hr,0.5);	
								
							if (z0 < z1 + hr*(1.+fabs(ubarra)/cbarra) && fabs(ubarra) < cbarra){
								Arrptr->FHy[pq0]=0.5*(cbarra*(z0-z1)-hr*(fabs(ubarra)+cbarra));
							}
							else{
								Arrptr->FHy[pq0]=0 ;
							}

							// If z1 >> z0  -> great losses
							// If ubarra=0  what happens 
						
*/
						
				}
			}
			else if(hr < dtol && hl > dtol){			/* right dry */
			
				
				fctm=1.-pow(1.,0.5)*fabs(hvl)/(pow(g,0.5)*pow(hl,1.5)) ;
				
				if (fctm < 0.) fctm=0.;
					
				if (hl*fctm+Arrptr->DEM[p0] > Arrptr->DEM[p1] && (hvl/hl +pow(g*hl,0.5) ) >= 0.){	/* overtopping */
				
					
					wsurface=(hl+Arrptr->DEM[p0]-hr-Arrptr->DEM[p1])/dl;	
					
					Arrptr->FHy[pq0]=hvl ;
 					Arrptr->FHUy[pq0]=hvl*hul/hl ;
 					Arrptr->FHVy[pq0]=hvl*hvl/hl ;
				
					Arrptr->H[p1]=dtol;

					Arrptr->BSHU[p0]=0. ;	
					Arrptr->TSHU[p1]=0. ;
	

					Arrptr->BSHV[p0]=g*hl*wsurface ;	
					Arrptr->TSHV[p1]=g*hl*wsurface ;
													
				}
				else
				{   				/* wall  */  
			
					/* wall at j+1 */
					
						// To verify that right cell is not updated
					if (Statesptr->Roe_slow==ON)
					{
					hr = hl;
					hur = hul;
					hvr = -hvl;

					//Arrptr->BSHU[p0]=0.0;
					//if (Arrptr->TSHV[p0] > 0) Arrptr->BSHV[p0] = 0.0; 	
					//else Arrptr->BSHV[p0]= -Arrptr->TSHV[p0];

					//s0=-(Arrptr->DEM[p1+Arrptr->xsz]-Arrptr->DEM[p1])/dl;

					//Arrptr->BSHU[p0]=0.0 ;	
					Arrptr->TSHU[p1]=0.0 ;
	
					//if (Arrptr->TSHV[p0] <= 0) Arrptr->BSHV[p0]=0.0;
					//else Arrptr->BSHV[p0]=Arrptr->TSHV[p0] ;	
					Arrptr->TSHV[p1]=0.0;	

					Arrptr->BSHU[p0]=0.0;
					Arrptr->BSHV[p0]=0.0;

					// implement solver for dry edge
					shr= pow(hr,0.5); shl= pow(hl,0.5); aux1= shr + shl ;
					ul= hul/hl; ur= hur/hr; vl= hvl/hl; vr= hvr/hr;
					ubarra= (shr*ur+shl*ul)/aux1; vbarra= (shr*vr+shl*vl)/aux1;	cbarra= pow(0.5*g*(hl+hr),0.5); 
					a1=vbarra+cbarra ; a2=vbarra ; a3=vbarra-cbarra ; 
					a1m=fabs(a1) ; a2m=fabs(a2) ; a3m=fabs(a3) ;
					a1l=vl+pow(g*hl,0.5) ; a1r=vr+pow(g*hr,0.5) ; 
					epsilon=maximum(0.,(a1-a1l),(a1r-a1)) ; a1m=(a1m >= epsilon) ? a1m : epsilon  ; 
					a3l=vl-pow(g*hl,0.5) ; a3r=vr-pow(g*hr,0.5) ;
					epsilon=maximum(0.,(a3-a3l),(a3r-a3)) ;	a3m=(a3m >= epsilon) ? a3m : epsilon  ;
					dh=hr-hl ; dhu=hur-hul ; dhv=hvr-hvl ; 

					alfa1=0.5*dh+0.5*(dhv-vbarra*dh)/cbarra ; alfa2=-(dhu-ubarra*dh)/cbarra ; alfa3=0.5*dh-0.5*(dhv-vbarra*dh)/cbarra ; 		
					e11=1. ; e12=ubarra; e13=vbarra+cbarra ; 
					e21=0. ; e22=-cbarra ; e23=0.; 		
					e31=1. ; e32=ubarra ; e33=vbarra-cbarra ; 
					/* flux variables */ 
					g1lp=hvl ; g2lp=hul*vl ; g3lp=hvl*vl+0.5*g*hl*hl ; 
					g1rp=hvr ; g2rp=hur*vr ; g3rp=hvr*vr+0.5*g*hr*hr ;
					Arrptr->FHy[pq0]=0.5*(g1rp+g1lp-a1m*alfa1*e11-a2m*alfa2*e21-a3m*alfa3*e31) ;
					Arrptr->FHUy[pq0]=0.5*(g2rp+g2lp-a1m*alfa1*e12-a2m*alfa2*e22-a3m*alfa3*e32) ;
					Arrptr->FHVy[pq0]=0.5*(g3rp+g3lp-a1m*alfa1*e13-a2m*alfa2*e23-a3m*alfa3*e33) ;
					}
					else
					{
						/* wall at j+1 */
					
						// To verify that right cell is not updated
						
						//Arrptr->FHy[pq0]=0 ;

						//if (hvl >= 0){
						
							Arrptr->FHUy[pq0]=0.;
							Arrptr->FHVy[pq0]=hvl*(hvl/hl+pow(g*hl,0.5))+0.5*g*hl*hl;
							Arrptr->FHy[pq0]=0.;//GUSTAVO

							Arrptr->TSHU[p1]=0.0 ;//GUSTAVO
							Arrptr->TSHV[p1]=0.0;//GUSTAVO
							Arrptr->BSHU[p0]=0.0;//GUSTAVO
							Arrptr->BSHV[p0]=0.0;//GUSTAVO
                                                       
						//}
						//else{
						
						//	Arrptr->FHUy[pq0]=0.;
						//	Arrptr->FHVy[pq0]=0.;
						
						//}
					}
		/*				
						//Arrptr->FHy[pq0]=0 ;

						if (hvl >= 0){
						
							Arrptr->FHUy[pq0]=0.;
							Arrptr->FHVy[pq0]=hvl*hvl/hl;
						}
						else{
						
							Arrptr->FHUy[pq0]=0.;
							Arrptr->FHVy[pq0]=0.;	
						
						}
						
						Arrptr->BSHU[p0]=0.;
						Arrptr->BSHV[p0]=0.;


						// WDC from Brufau et Al    

							
							//Raw formulation from Brufau (Dont add RSH[p0] for avoiding UpdateRoeH)


							ubarra= ul;
							cbarra= pow(0.5*g*hl,0.5);	
								
							if (z1 < z0 + hl*(1.+fabs(ubarra)/cbarra) && fabs(ubarra) < cbarra){

								Arrptr->FHy[pq0]=-0.5*(cbarra*(z1-z0)-hl*(fabs(ubarra)+cbarra));
							}
							else{
								Arrptr->FHy[pq0]=0.;
							}
		*/
				}				
		
			}

    //----------------------------------------------------------------------------------------

  Q=Arrptr->FHy[pq0]*dl  ;	// origional version included addition *dl but removed for compatability with UpdateH JCN
  // option to save V's
  if (Statesptr->voutput==ON)
  {
	  if (Q!=0)
	  {
	    hflow=getmax(z0+h0,z1+h1)-getmax(z0,z1);
        hflow=getmax(hflow,0);
        hflow=getmin(hflow,Solverptr->MaxHflow);
	    // calc Vy
	    Arrptr->Vy[pq0]=Arrptr->FHy[pq0]/hflow;
	    // get max Vy
	    Arrptr->maxVy[pq0]=getmax(Arrptr->maxVy[pq0],fabs(Arrptr->Vy[pq0]));
	  }
	  else Arrptr->Vy[pq0]=0.0;
  }
  return(Q);
}
////-----------------------------------------------------------------------------------
//// Calculate max(H) for each timestep
/*double CalcMaxHRoe(Pars *Parptr, Arrays *Arrptr) ... duplicate of CalcMaxH removed JCN
{
	int i,j,p0;
	double h0, Hmax=0.0;

// Loop through to calculate maximum water depth
//#pragma omp parallel for private(i,p0,h0) reduction(MAX:Hmax) // JN: omp MAX only supported in FORTRAN
	for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
	{
		p0=i+j*Parptr->xsz;
		h0=Arrptr->H[p0];
		if(MaskTestRoe(Arrptr->ChanMask[p0])) 
		{
			Hmax=getmax(h0,Hmax);
		}
	}

return(Hmax);
}*/
//-----------------------------------------------------------------------------------
// Calculate timestep for acceleration version based on 2D shallow water CFL condition
void CalcTRoe(Pars *Parptr, Solver *Solverptr, Arrays *Arrptr)
{
	int i,j,p0,pqx,pqy;
	double cfl=0.7, g;
	//double MH=0.0;
	double locH=0.0,locV=0.0,locU=0.0;
	double locT,locTx,locTy;

	g=Solverptr->g;
	cfl=Solverptr->cfl;
	
	// Calculate maximum water depth every timestep
	//MH=CalcMaxHRoe(Parptr, Arrptr);
	// Apply local 2D CFL condition and use global minimum timestep to ensure global stability
	//if(MH >	Solverptr->DepthThresh) //Don't apply where h=0 as h appears in equation denominator
	//{
	// set time step to initial time step
	Solverptr->Tstep=Solverptr->InitTstep;

	for(j=0;j<Parptr->ysz;j++) for(i=0;i<Parptr->xsz;i++)
	{
		p0=i+j*Parptr->xsz;	
		locH=Arrptr->H[p0];
			
		if(locH > Solverptr->DepthThresh)
		{
			pqx=i+j*(Parptr->xsz+1)+1;
			pqy=i+(j+1)*(Parptr->xsz+1);
			
			locV=Arrptr->HV[p0]/locH;
			locU=Arrptr->HU[p0]/locH;
		
			locTx=cfl*(Parptr->dx/(fabs(locU)+(sqrt(g*locH))));
			locTy=cfl*(Parptr->dx/(fabs(locV)+(sqrt(g*locH))));
			locT=(locTx < locTy) ? locTx : locTy ;
			
			Solverptr->Tstep=getmin(Solverptr->Tstep,locT);
		}	
	}
	//}
	//else // Set to initial timestep if h=0
	//{
	//	Solverptr->Tstep=Solverptr->InitTstep;
	//}
	
return;
}
//----------------------------------------------------------------------------
// MaskTest for channel for MaxH calculation ... duplicate of MaskTest removed JCN
/*int MaskTestRoe(int m1)
{
  if(m1==-1) return(1);
  
  return(0);
}
*/
//----------------------------------------------------------------------------
// Update all the old qs with new qs for acceleration version
void UpdateQsRoe(Pars *Parptr, Solver *Solverptr, Arrays *Arrptr)
{
	int i,j;

	int pc,pxl,pxr,pyl,pyr;
	double fn,sf,norm,manning;

// Loop over all  cells to update HU and HV
   #pragma omp parallel for private(i,pc,pxl,pxr,pyl,pyr,fn,sf,norm,manning) // parellisation added JCN
	for(j=0;j<Parptr->ysz;j++)
    {
      for(i=0;i<Parptr->xsz;i++)
      {			
          pc=i+j*Parptr->xsz;
       
		  if(Arrptr->H[pc] > Solverptr->MomentumThresh)
		  {
			pxl=i+j*(Parptr->xsz+1);	
			pxr=i+1+j*(Parptr->xsz+1);
			pyl=i+j*(Parptr->xsz+1);	
			pyr=i+(j+1)*(Parptr->xsz+1);		

			Arrptr->HU[pc]+=-(Solverptr->Tstep/Parptr->dx)*(Arrptr->FHUx[pxr]-Arrptr->FHUx[pxl]+Arrptr->FHUy[pyr]-Arrptr->FHUy[pyl])+Solverptr->Tstep*(Arrptr->BSHU[pc]+Arrptr->TSHU[pc]+Arrptr->LSHU[pc]+Arrptr->RSHU[pc]) ;
			Arrptr->HV[pc]+=-(Solverptr->Tstep/Parptr->dx)*(Arrptr->FHVx[pxr]-Arrptr->FHVx[pxl]+Arrptr->FHVy[pyr]-Arrptr->FHVy[pyl])+Solverptr->Tstep*(Arrptr->BSHV[pc]+Arrptr->TSHV[pc]+Arrptr->LSHV[pc]+Arrptr->RSHV[pc]) ;
	        /* Implicit friction */	
			norm=pow((Arrptr->HU[pc]/Arrptr->H[pc])*(Arrptr->HU[pc]/Arrptr->H[pc])+(Arrptr->HV[pc]/Arrptr->H[pc])*(Arrptr->HV[pc]/Arrptr->H[pc]),0.5) ;		
			if(Arrptr->Manningsn!=NULL)
			{
			  fn=Arrptr->Manningsn[pc];   // centered
			}		
		    else
		    {
			  fn=Parptr->FPn;
		    }	
		    manning=pow((pow(fn,1.5)/Arrptr->H[pc]+0.*pow(fn,1.5)/Parptr->dx),(4./3.)) ;
		    sf=manning*norm  ;
		    Arrptr->HU[pc]*=(1./(1.+Solverptr->g*Solverptr->Tstep*sf)) ;
		    Arrptr->HV[pc]*=(1./(1.+Solverptr->g*Solverptr->Tstep*sf)) ;
		  }
		  else{
			Arrptr->HU[pc]=0.;
			Arrptr->HV[pc]=0.;
		  }
	  }
	}	
	return;
}
//----------------------------------------------------------------------------


double maximum(double a,double b,double c) 
{
double max ;

max=(a>b) ? a : b ;
max=(max>c) ? max : c ;
return(max) ;
}
