/*
#####################################################################################
LISFLOOD-FP flood inundation model
#####################################################################################

@ copyright Bristol University Hydrology Research Group 2012

webpage -	http://www.ggy.bris.ac.uk/research/hydrology/models/lisflood
contact -	Professor Paul Bates, email: paul.bates@Bristol.ac.uk,
Tel: +44-117-928-9108, Fax: +44-117-928-7878

*/

// version setting *** make sure you update before submitting to subversion ***
// don't forget to update checkpointing version if it is affected
#define VersionMajor 5
#define VersionMinor 99
#define VersionInc   0
#define CheckVersion 0

/*
########################### Version History #########################################
Reserved for major changes, use list below this one and subversion for full details

Ver  Date         Details
---  ----         -------------------------------------------------------------------------
5.99  Oct 07 2016  Re-structured code to be executable via BMI-adapter (Jannis Hoch)
5.9  Aug 05 2013  Support for LatLong added (CCS)
5.8  Feb 15 2013  Regression model wetted perimiter added, bug fix for SGC wider than cell (JCN)
5.7  Sep 26 2012  Fully tested and bug fixed bridge implementation in subgrid (MT)
5.6  Aug 06 2012  Alternative sub-grid channel geometries added (JCN)
5.5	 Jun 27 2012  2D version of friction (x-y coupled) implemented. Old (1D) version still available using the "1Dfriction" keyword in the .par file.(Gustavo A.M. de Almeida) 
5.4	 Mar 21 2012  Rainfall routing added to allow rainfall to be simulated over complex terrain (Chris Sampson)
5.3  Sep 26 2011  q-centred numerical scheme implemented and tested for the solution of the simplified St. Venant eq. (Gustavo A.M. de Almeida) 
5.2  Jul 22 2011  Weir flows implemented in sub-grid channel and floodplain models (Jeff Neal)
5.1  May 31 2011  Initial sub-grid channel implementation
5.0  Jan 21 2011  Fully functional Roe solver (Jeff Neal) and multiple river capability (Chris Sampson)
4.4  Feb 01 2010  Tested version of Roe solver. 2D only, point source closed boundary only added by Jeff Neal
4.3  Sep 04 2009  Full dynamic and diffusive initial steady state channel solution added and tested by Tim Fewtrell
4.1  Nov 10 2008  TRENT formulation added. Integrated version tested by Jeff Neal. TRENT solver not tested
4.0  Aug 18 2008  Acceleration formulation fully implemented and tested by Paul Bates and Tim Fewtrell
3.6  Jul 31 2008  Decouple river channel timestep from floodplain timestep by Mark Trigg
3.5  Jun 13 2008  OpenMP version implemented and tested on Bustcot by Jeff Neal
3.4  Apr 21 2008  Double precision version by Mark Trigg
3.3  Jan 11 2008  Diffusive channel solver & Bug fixed branching channels by Mark Trigg
3.1  Oct 08 2007  Fully tested and bug fixed modular code by Mark Trigg
3.0  May 25 2006  Modularised the code and added porosity scaling algorithm by Tim Fewtrell
2.7  Feb 25 2005  Evaporation and Infiltration added by Matt Wilson
2.6  Dec 20 2004  Add more output file and command line options by Matt Wilson
2.5  Nov 25 2004  Checkpointing functionality added by Matt Wilson
2.0  Jun 08 2004  Adaptive timestep implemented by Neil Hunter
1.0  2003         First public release version by Matt Horritt
0.9  2003         Increase output file and command line options by Matt Wilson
0.8  2001         Prototype C++ version created by Matt Horritt
0.5  2001         Original version created by Paul Bates and Ad De Roo




Date         Version  Details
----         -------  -----------------------------------------------------------------------
May 02 2014  5.9.11   Corrected point HVAR and HFIX boundary (JCN)
Feb 04 2014  5.9.9    Correction to rain-(evap+infil) calculation (JCN)
Jan 28 2014  5.9.8    Added point free boundary to subgrid model (JCN)
Jan 24 2014  5.9.7    Added SGC velocity output and new user manual (JCN)
Oct 23 2013  5.9.6    Added SGC bank full depth mode to SGC model and distributed velocity to routing model, correction to flood area OMP calc (JCN)
Oct 06 2013  5.9.5    Correction to FloodedArea calc for latlong mode (JCN)
Sep 30 2013  5.9.4    Correction to SGC max elevation calculation (.mxe file) (JCN)
Aug 20 2013  5.9.3    Edit to .pram file read to support \n, \r and EOF line endings (JCN) (Removed!)
Aug 13 2013  5.9.2    Change to SGCgroup file read to make it more robust to a mixture of float and integer data (JCN)
Aug 06 2013  5.9.1    Change to channel length calculation to support different dx and dy (JCN)
Aug 05 2013  5.9.0    Support for LatLong added (CCS), merged with code and revisions 5.8.7-12 (JCN)
Aug 02 2013  5.8.12   Bug fux in CalcSGCz introduced by update 5.8.11 (JCN)
Jul 16 2013  5.8.11   Channel accumulation area added to SGC channel (JCN)
Jul 09 2013  5.8.10   Meander coefficient added to sub-grid model (JCN)
Jul 02 2013  5.8.9    Increase max array size from 10,000 to 20,000 to allow more point sources for user (MT).
Jun 24 2013  5.8.8    Corrected weir free/drowned flow condition (JCN)
Jun 14 2013  5.8.7    Added hazard calculation to subgrid model and changed hazard equation to DEFRA 2006 (ALD)
May 13 2013	 5.8.6	  Routing scheme enabled for subgrid floodplain.  Routing can now be used in place of inertial eqn. for
					  steep water surface gradients using an assumed fixed velocity.  Routing can be used without rainfall (CCS)
Apr 16 2013  5.8.5    Added flow width and channel Q output to sub-grid model, added variable gravity for mars experiment (JCN) 
Mar 25 2013  5.8.4    Voutput added to sub-grid model (JCN)
Mar 22 2013  5.8.3    Fix to ascii_write SGC wdfp output file (JCN)
Mar 06 2013  5.8.2    New gamma parametrs for SGC power shaped channel, bug fix for SGC channel wider than cell width (JCN)
Mar 04 2013  5.8.1    Bug fix to weir file to prevent write to QxSGold and QySGold during non SGC runs (JCN)
Feb 15 2013  5.8.0    Implementation of regression model for SGC exponential channel wetted perimiter estimate (JCN)
Feb 13 2013  5.7.8	  Bug fix for HVAR boundary condition when used within acceleration mode (CAB)
Jan 31 2013  5.7.7    Minor change to reading of banks, width and changroup to make them robust to no_data_values (JCN)
Jan 21 2013  5.7.6    Minor correction to .discharge file output to round up rather than down to nearest cell width (JCN)
Jan 21 2013  5.7.5    Minor correction to SGC channel Mannings to initalise array as NULL (JCN)
Jan 14 2013  5.7.4    SGC fix for channels with widths greater than a cell, distributed SGC channel Mannings added. (JCN)
Dec 17 2012  5.7.3    New hydraulic radius calculation added for SGC, independent SGC parameters and regions added (JCN)
Oct 23 2012  5.7.2    Significant change to SGC updateH to a calculation based on volume. Function CalcSGC_UpV added. (JCN)
Oct 09 2012  5.7.1    Upgrade to sub-grid model to account for tributary volume in cell. Output to command line added for SGC channel parameters.
                      Fix to CalcFPQxSGC and CalcFPQySGC bug introduced in version 5.6.10 (JCN)
Sep 25 2012  5.7.0    Fully tested and bug fixed bridge implimention in subgrid with new test cases added to test directory (MT)
Sep 25 2012  5.6.10   New sub-grid channel type implemented. Edit to pars.cpp to force acceleration=ON with SGC model (JCN)
Sep 18 2012  5.6.9    Another SGC FREE boundary bug fix edit to por_flow to prevent compiler double to int warnings (JCN)
Sep 17 2012  5.6.8    Fix to SGC FREE boundary bug (JCN)
Sep 14 2012  5.6.7    UpdateH() & SGC_UpdateH() Moved point source BC dH calc to before dV check,  
					  in case flow out of cell results in negative H reset and the point inflow would prevent this (MT)
Sep 12 2012  5.6.6    Included subgrid arrays in debug output. (MT)
Sep 11 2012  5.6.5    Included energy gradient elevation in bridge orifice flow calculation - important for high approach velocities. (MT)
Sep 03 2012  5.6.4    Added function for bridge to SGC model. Only operates when bridge soffit is exceeded. (JCN)
Aug 08 2012  5.6.3    Added functions for spatially distributed SGC model parameters (JCN)
<<<<<<< .mine
Aug 07 2012	 5.6.2	Minor change to allow output directory inclusion on MAC OSX systems (GAMA).
=======
Aug 07 2012	 5.6.1    Minor correction to SGC parabolic channel update h (JCN).
>>>>>>> .r229
Aug 06 2012	 5.6.0    Alternative sub-grid channel geometries added - 1.Rectangular, 2.Triangular (right angle), 
					  3. Triangular (equalateral), 4.Parabolic, 5.Trapazoidal, 6.No banks (JCN)
Jul 22 2012  5.4.5    Roe solver HVAR and QVAR subcritical boundary conditions implemented using ghost cells (GAMA)
Jul 17 2012  5.4.4    HVAR boundary condition implemented for Roe solver (GAMA)
Jul 13 2012  5.4.3    QVAR boundary condition implemented for Roe solver (GAMA)
Jun 22 2012  5.4.2    Minor correction to SGC evap routine mass balance calculation to account for bank transitions (JCN) 
Mar 27 2012	 5.4.1	  Added flow limiter to rainfall routing routine based on surface water heights (CCS)
Mar 21 2012  5.4.0    Rainfall routing added to allow rainfall to be simulated over complex terrain (CCS)
Feb 21 2012  5.3.8    Corrections to the Roe solver (wet/dry fronts) which significantly reduce mass errors (GAMA). 
Feb 21 2012  5.3.7    Change in the Roe_slow solver. The variables RSHU and LSHU should be set to zero at internal walls (GAMA).
Feb 06 2012  5.3.6    SGC_Evaporation parallelised (JCN)
Feb 03 2012  5.3.5    Improvements to efficiency of SGC update H and change to R calculations in rectangular SG channel (JCN)
					  Correction to load evap and load rain to prevent overwrite of data if an empty line is present at the end of the file see note in code.
Jan 16 2012  5.3.4    Minor edit to -T.op output to correct filename (JCN).
Jan 13 2012  5.3.3    SGC rectangular channel changed to use R and A. New function to export floodplain water levels from SGC model. 
					  inclomplete SGC trapazoidal channel momentum code added. Minor edit to output.cpp save more than 9999 files.
					  Edit to SGC model to allow channel widths greater than cell width. (JCN)
Nov 15 2011  5.3.2    Minor edit to output.cpp to allow users to save more than 9999 output files (GAMA)
Nov 08 2011  5.3.1    Minor edit to fp_acc.cpp because the copysign function isn't available in visual studio. (JCN)
Sep 26 2011  5.3.0    q-centred numerical scheme implemented and tested for the solution of the simplified St. Venant eq. (Gustavo A.M. de Almeida)
Sep 07 2011  5.2.7    Load weir amended to check weir crest is above DEM elevation warning message added. In the case of the sub-grid model only 
					  a weir below DEM or SGCz will have its crest height increased to SGCz or DEM to prevent model instability (JCN) 
Aug 17 2011  5.2.6    Correction to SGC hotstart (JCN)
Aug 09 2011  5.2.5    Minor change to stage time series output when t != 0 (JCN) 
Aug 04 2011  5.2.4    startq implemented for SGC model to initialise q using Manning's equation. tested on SGC straight channel. (JCN)
Aug 03 2011  5.2.3    Functionality added to read water surface eleveation files for initial condition. Correction of binary water surface
					  elevation output. Update to overpass profiles output to only output when profiles is ON. Update to allow water surface 
					  elevation output for non-SGC models. Functionality added to export binary max water surface elevation. (e.g. if you 
					  ask for binary_out all raster are now binary files rather than a mixture). tstart command added to .par file for starting
					  at any time in bdyfile. (JCN)
Jul 28 2011  5.2.2    Functionality added to read binary start files (JCN)
Jul 26 2011  5.2.1    Virtual discharge sections added (JCN)
Jul 22 2011  5.2.0    Weir flows implemented in sub-grid channel and floodplain models (JCN)
Jun 29 2011  5.1.9    Extentsion of HVAR correction in 5.1.7 to HFIX. Correction of openMP private variable list is SGC updateH (JCN)
Jun 29 2011  5.1.8    Minor change to inflow.cpp boundary conditions to prevent overwrite of point sources by boundary conditions (JCN)
Jun 29 2011  5.1.7    Minor changes in HVAR BC: i) change water slope sign ii) set Q to zero when both depths are bellow threshold (GAMA)
Jun 24 2011  5.1.6    Minor change in boundary condition to allow HVAR BC calculation when depth in the first cell is below threshold (GAMA)
Jun 23 2011  5.1.5    Upload of windows executable and associated .dll files. Minor fix to acceleration boundaries from 5.1.3 (JCN)
Jun 22 2011  5.1.4    Minor bug fix to SGC boundary conditions (JCN).
Jun 20 2011  5.1.3    Implementation of SGC HVAR and HFIX boundaries + minor correction to acceleration HVAR and QVAR boundaries. (JCN)
Jun 10 2011  5.1.2    keyword 'binary_out' added to export rasters as binary files. Further updates to SGC HVAR  (JCN).
Jun 07 2011  5.1.1    Update to SGC point QVAR to prevent mass errors during bank overtopping at inflow (JCN).
May 31 2011  5.1.0    Initial implementation of sub-grid channels sgc.cpp added (JCN)
Feb 03 2011  5.0.2    max velocity calculations changed to calculate maximum absolute velocity, unnecessary recording of hflow removed. (JCN)
Jan 26 2011  5.0.1	  ChannelQ and ChannelQ_Diff now parallelised in OpenMP for multiple rivers. Fixed profiling bug when using multiple 
					  rivers. (CCS)
Jan 21 2011	 5.0.0	  Restructuring of channel model to allow multiple unconnected river networks using a .rivers file (CCS)
Nov 30 2010  4.4.13   Addition of Roe_slow keyword to chose between Ghost cell and non-ghost cell wet/dry method in Roe solver(JN/IV)
Nov 25 2010  4.4.12   Updated Roe solver internal wet/dry boundaries to use ghost cells, bug fix on close boundaries. Tested on 
                      Glasgow2m. (JCN) 
Nov 23 2010  4.4.11   Roe closed boundary update (JCN)
Sep 07 2010  4.4.10   Roe dry edge updates (JCN for IV)
May 28 2010  4.4.9    Option to calculate hazard added using parameter file keyword hazard. hfowx and hflowy variables 
                      added to store Hflow. (JCN)
May 27 2010  4.4.8    Implementation of a simple steady-state simulation control, such that a simulation will
                      automatically end when Qout matches Qin (within a certain tolerance). Enabled using the keyword -steady
					  on the command line. Tolerance defaults to 0.0005 (difference between Qout and Qin, in m3s-1), with the
					  keyword -steadytol (followed by value) allowing user control of this (MDW).
May 25 2010  4.4.7    Update to Roe model. MomentumThresh and FREE boundary added. DepthThresh and MomentumThresh can now
                      be changes from the parameter file using depththresh and momentumthresh. (JCN for IV).
Mar 26 2010	 4.4.6	  Bug fix for voutput and stage output files so that stage output can occur independently of the location
					  specific voutput. (TJF)
Mar 09 2010	 4.4.5    OpenMP critical section bug fixes for time step update to acceleration, adaptive and Roe. A bunch
					  of bug fixes to rainfall routine and startfile/DEM loading for comparing no_data values.
					  Windows and Linux executables updated. (TJF)
Feb 13 2010  4.4.4    Updateds to Roe solver. keywork voutput will now export point velocities if stage file present (JCN).
Feb 04 2010  4.4.3    Functions to output velocity (keyword "voutput") added to diffusive, intertial and Roe solvers (JCN).
Feb 02 2010  4.4.2    Roe solver update. 2D only FREE or closed boundary. Tested on Carlisle. (JCN)
Feb 01 2010  4.4.1    Tested version of Roe solver. 2D only, point source closed boundary only. OMP implemented with Roe (JCN).
Jan 08 2010  4.3.9    Time varying, spatially uniform rainfall module incoporated. Rainfall option activated using
					  "rainfall" followed by filename in the parameter file. Format of rainfall file is the same as the
					  evaporation file format. Mass balance upkeeping updated to include new input. Units in mm t-1 (TJF)
Dec 17 2009  4.3.8    Bug fix for startq fabs() instead of abs() in exit criteria (JCN).
Dec 08 2009  4.3.7    FREE boundary update for adaptive/qlim version to allow user specified slope. General tidy
                      up of BCs to prevent acceleratin verson calulating both adaptive and acceleration boundaries (JCN).
Nov 25 2009  4.3.6    FREE boundary update for acceleration verion to alloww flow in all directions. (JCN)
Nov 19 2009  4.3.5    Mass error calculation fix. Verror now calculated. Inf and Evap combined. (JCN)
Sep 29 2009	 4.3.4    FREE boundary update for acceleration version implemented and tested by Jeff Neal.
					  To use a user specified slope at the boundary add the slope to the .bci file after FREE keyword.
					  Option added to parameter file to turn off drycheck function. Use "drycheckoff" keyword and only
					  use with adaptive or acceleration versions (JCN).
Sep 29 2009	 4.3.3	  Option to overwrite dhlin added for command line and parameter file. Use "-dhlin" and value on
					  command line and "dhlin" keyword and value in parameter file. Value to read in is the dh threshold
					  for linearisation. Default value is now set as a function of dx and a gradient of 0.0002
					  (from Cunge et al., 1980) to improve stability at high resolutions (TJF).
Sep 15 2009  4.3.2    FREE boundary implementation for diffusive channel solver implemented and tested by Tim Fewtrell
					  Last line of function vector and Jacobian matrices edited to substitute Manning's into diffusive
					  form of momentum equations. New exit criteria created for diffusive solver using the norm or max
					  of x or f(x) where x is the solution to the matrix. Uses norm(f(x)) as default and hardcoded. See
					  manual for detailed explanation of implementation. Diffusive solver commented in more detail (TJF).
Sep 10 2009	 4.3.1    Bug fix for uninitialised maximumH for acceleration time step calculation (should lead to
					  longer initial time steps) (TJF)
Sep 04 2009  4.3.0    Full dynamic and diffusive initial steady state channel solution added and tested by Tim Fewtrell
					  Turned on using "startq" in parameter file. Diffusive steady state used by default but can be changed
					  to full dynamic using "ch_dynamic" in parameter file or "-dynsw" on the command line. Kinematic still
					  uses kinematic initial solution. (TJF)
Jun 25 2009  4.2.3    Force flush of the stage file to ensure it stays in sync when using checkpointing in particular and
					  infomration added to stage file for a checkpoint restart. (TJF)
Jun 24 2009  4.2.2    CFL value for acceleration version is now adjustable in the parameter file and on the command line.
					  -cfl flag followed by the value for the command line and cfl option in parameter file followed by the
					  value. Guidance suggests values between 0.2 and 0.7. Default value is 0.7. (TJF)
Jun 18 2009	 4.2.1    Secondary bug fix for FREE boundary calculation - unlikely to have a major effect but is now coded
					  correctly. (TJF)
Jun 16 2009  4.2.0	  Bug fix for boundary flow calculation. Previously, HFIX or HVAR boundaries on the floodplain were
					  forced based on dhlin resulting in small fluxes and water depths. Update corrects this and thus allows
					  downstream boundary conditions to be implemented correctly. (TJF)
Jun 12 2009	 4.1.7	  Bug fix in FloodArea calculation when using OpenMP which should speed things up. Also implemented
					  more elegant coding for -kill option to reduce duplicate code - now uses a break command and finished
					  final output. (TJF)
Jun 05 2009  4.1.6    Added a -kill option to the command line: this is followed by the desired maximum computation time
                      in hours. With this option, the model will exit after this length of runtime. This is useful where
					  a cluster has a maximum runtime in a queue and model output isn't returned to killed jobs. This allows
					  the model to exit gracefully before it is kicked off the queue. (MDW)
May 22 2009  4.1.5	  Tstep calculation updated in acceleration version. OMP critical section added to MaxH calculation for
					  acceleration and OMP functionality added to UpdateQs. (TJF)
Mar 27 2009  4.1.4	  Adjusted Tstep update to include omp critical section around Tstep if its being changed. This prevents
					  a race occurring when more than one tread attempts to update Tstep... will significantly increase speedup
					  for some models. (JN)
Mar 12 2009  4.1.3    Fixed bug in mass calc for infiltration and ET. Added new ET and Infiltration
                      test to test case directory T012_ETInfTest. (MT)
Mar 11 2009  4.1.2    Added new parameter "qlimfact" to allow user to specify a relaxation factor for Qlimit.
                      Qlimit applies in non adaptive mode to ensure solution remains stable, unfortunately
                      it also artificially modifies the progress of the floodwave. Under some hydraulic conditions,
                      the qlimit threshold at which flow from one cell to the next is limited can be relaxed
                      without destabilising the solution. A very high factor will effectively switch off Qlimit.
                      Be careful using this option and make sure your solution is stable and believable.
                      None adaptive mode and qlimfact should only be used where run times are excessive and you
                      understand the consequences of using this mode. (MT)
Feb 18 2009	 4.1.1	  Checkpointing fixed for all x86 (i.e. little endian) architectures. Check for channel,
                      checkpoint version and LISFLOOD-FP version added. Checkpoint version moved to VersionHistory.h and
                      should ONLY be changed if checkpointing is changed. If either versions are different,
                      checkpoint file is ignored and starts from scratch. (TJF)
Nov 10 2008  4.1.0    Initial version of Ignacio's TRENT formulation in LISFLOOD. New file fp_trent.cpp
                      added. Specify "Roe" in parameter file to use new version (JN/IV)
Sep 2  2008  4.0.1    iterateq.cpp updated to allow calculation of maxH, maxHtm, totalHtm and initHtm at
                      the mass interval if mint_hk is included in the .par file. Corrected #pragma
                      statement in fp_acc.cpp (JCN). Testcases run against 4.0.0 with ipcp.
Aug 18 2008  4.0.0	  Acceleration version fully implemented and tested. Use keyword "acceleration"
                      in the parameter file or command line and ensure "adaptoff" not specified. Choice
                      of initial timestep now important as this is the upper limit of the timestep. (PDB/TJF)
Aug 13 2008  3.6.3    Fixed minor bug in new channel chainage code bug which produced #IO for
                      the vol in mass file. Related to river files that ovelap dem significantly.
Aug 11 2008  3.6.2    Added some code to make river channel chainage independent of cell size.
                      This uses the straightline chainage between the entered cross-sections,
                      rather than that derived from the cell dx dimension. It can be switched
                      off to allow backwards compatibility by entering the keyword chainageoff
                      in the parameter file. The points defined in the river file will now
                      control the chainage, so these should define the river with enough detail
                      that you get the 'true' centreline chainage. (MT)
Jul 31 2008  3.6.1    Fixed minor bug in BCs  relating to recording timestep. (MT)
Jul 31 2008  3.6.0    Decouple river channel timestep from floodplain timestep. Lisflood will
                      now default to a x1 multiple, or you can specify the multiple to use in
                      the parameter file, eg "ts_multiple 10". Note that the multiple has to be an
                      integer. It is defined as a multiple rather than as a channel timestep
                      so it can follow the fp adaptive timestep. Tests show up to x10 gives
                      almost identical results to x1. You can go higher than x10, but the error
                      increases and the speed gain is marginal. If you do use multiples of
                      x100 or more, check the sensitivity of your results to this.(MT)
July 17 2008 3.5.3    Increase precision of stage file output from 2 to 4 decimal places. (MT)
July 14 2008 3.5.2    Fixed bug in stage file read - displaced location South by one cell. (MT)
July 7 2008  3.5.1    Fixed bug in infevap (incorrect array indexing). (JCN, TJF)
June 13 2008 3.5.0	  OpenMP version implemented and tested on Bustcot. Accessed using the makefile
                      by changing omp variable to "1". Compiles using ICPC as a static exectuable. (JCN)
June 12 2008 3.4.2	  Kinematic channel code bug fix. Remove double counting of bank flows for
                      junction node between channels. (MT)
May 20 2008	 3.4.1	  Time reporting with model/comp ratio make as optional with "comp_out" keyword (TJF)
Apr 21 2008  3.4.0    All code converted to Double Precision from Single Precison Floats. Reduces
                      mass error on smaller timesteps and smaller mesh sizes. No speed penalty as
                      compilers and cpu native precision is now double anyway. Fixed a minor bug
                      related to transfer of the d/s main channel water level to trib BC.
                      Add new "-log" option to command line - redirects screen output to a log
                      file. River channel profile output now very detailed including flows and bank levels.
                      htol is now an optional parameter in the par file, but defaults to 1m. The
                      checkpoint function no longer overwrites previous results files when restarting and
                      new mass file output lines are appended with a clear break point. (MT)
Jan 14 2008  3.3.2    Add new function "SetChannelStartH()" to set starting water depth in channel.
                      User can leave as default 2m water depth for all channels or set the level
                      with "ch_start_h" and a number in parameter file (as described in 3.3.0).
                      Or user can use "startq" in parameter file which means lisflood will use
                      the starting inflow to each channel and assume this is the starting flow
                      down the whole channel and calculate a starting h from this. Also reapply
                      bug fix from 3.2.2 which got left out. (MT)
Jan 14 2008  3.3.1    Add time reporting with model/comp ratio and estimated time to complete run. (MT)
Jan 11 2008  3.3.0    Fully tested and working diffusive channel solver now implemented - derived
                      from Matt Horritt's original single channel diffusive work. New downstream boundary
                      conditions for channel HFIX, HVAR and FREE (bed slope or forced slope) implemented for
                      diffusive (defaults to FREE). Debug parameter and command-line now outputs 3 files;
                      the final dem after burning in the channel and bank mods (*.dem), the channel mask (*.chmask) and
                      the channel segment mask (*.segmask). Numerous bug fixes in branched channel code.
                      Key ones: Hds in mass file is now from main channel, not last one calculated.
                      Now adds only main channel output to boundary flux instead of all channel outputs.
                      Trib connections were out by one cell resulted in incorrect slope. Where
                      more than one trib on one channel, outflows were assigned to wrong locations.
                      Fixed strange slope mod for kinematic solver which set -ve slopes to very small slope resulting
                      in silly water levels - replaced with simple positive slope of same magnitude giving
                      reasonable results (also -ve slope warning to user). Add new parameter ch_start_h.
                      This allows the user to specify the starting water depth in the channel.
                      The default is 2m (as in previous versions) if not explicitly specified. (MT)
Dec 13 2007	 3.2.3	  Fixed minor bug in FloodedArea calculation introduced by porosity scaling. (TJF)
Dec 12 2007	 3.2.2    Fixed porosity technique implementation to scale both volume and flux to maintain
                      mass conservation. Test against fully blocked cell provides identical solutions.
                      Fixed output of Qx and Qy to output at cell boundaries and offset asciiheaders to
                      enable easy integration with GIS. Bug fixed in the storage and output of Qlim
                      and Trec arrays. Note: Qlim output may not work in corner cells. (TJF)
Nov 28 2007  3.2.1    Fix boundary condition bug which stops water exiting north and west under FREE
                      conditions and bug which displaces south and west FREE boundaries by one cell. (MT)
Nov 20 2007  3.2.0    Diffusive channel solver added. Set using parameter "diffusive" in parameter file.
                      Not fully working yet - d/s boundary is currently hardcoded and multichannel
                      not full functional.(MT)
Nov 20 2007  3.1.4    Fixed bug in boundary flux calc for north boundary (added wrong way) (MT)
Nov 19 2007  3.1.3    Fixed bug in counter increment for profiles output.
                      Counter was incrementing before profiles output completed. (TJF)
Oct 31 2007  3.1.2    Bug fix for domains without a channel. Now correctly assigns outflows when ChannelPresent=OFF. (TJF)
Oct 23 2007	 3.1.1    Bug fix in channel code that resulted in the last crossection cell being excluded
                      from the channel mask/calculation and resulted in river results being one cell out.
                      Also modified slope calculation to use the slope at the crossection not the slope
                      of the next segment. Checked results against MacDonald 1997 analytical test 1.
                      Results show a 3 to 12% error on water level for the channel calc.
                      Revert to previous version numbering to avoid random subversion number updates.
                      Avoids meaningless subversion number updates as well as needing to know the number
                      that would be issued by subversion.
Oct 08 2007	 3.1.Sv56	Full option testing carried out on the modularised code. Buscot test case
                      gives numerically identical answers in both adaptive and non-adaptive time step modes.
                      A number of bug fixes added - Smooth banks and checkpointing and parameter file closing. (MT)
May 25 2006	 3.0.Sv55	Modularised the code into separate files. Defined structures to variables in
                      lisflood.h. Full documentation of structures and pointers in manual. Added to Subversion.
                      (Feb 06 2006)	 Added in simple porosity scaling algorithm. Porosity maps derived
                      from porosity.cpp. Keyword in parameter file is "porfile" followed by the filename. (TJF)
Aug 30 2007  2.7.Sv55	Apportioned blame for previous code versions in header file (PB)
Aug 30 2007  2.7.Sv54	Tidied up the history list. (PB)
Aug 30 2007  2.7.Sv53	Tidied up revision numbering (again) (PB)
Aug 30 2007  2.7.Sv52	Tidied up revision numbering (PB)
Aug 30 2007  2.7.Sv51	Edited header information; updated user manual; changed output of maximum elevation file
                      so null data cells are defined as NULLVAL rather than as -99.
Jan 31 2007  2.7.Sv49	Standardised file naming so all files have sensible extensions.
                      All ascii output now through one function with flag to select options - reduces lots of duplication.
                      Unnested output options that were incorrectly nested within an if(save_depth==ON).
                      Fixed minor bug that stopped checkpointing working when using an multiple overpass file.  (MT)
Jan 25 2007  2.7.Sv48	Add optional raster output for Timesteps and Qlimits. New flags in parameter file,
                      "toutput" & "qloutput" will cause the program to output rasters with file extensions
                      .Tx & .Ty, .QLx, .QLy at the normal file output interval.
                      Internal program flags are save_Ts and save_QLs and by default are OFF.
                      Data is store in TRecX,TRecy, LimQz, LimQy. 2 new functions, write_Ts() and write_QLs() output the files.
                      Note write_Qs() has also been updated.
                      Also the ascii null value of -9999.0 has been replaced throughout with the constant NULLVAL.  (MT)
Jan 25 2007  2.7.Sv47	Tidy up file header to make it clearer what the file is.
                      Update versioning so that we have a major.minor number +(subversion number) eg for this version 2.7 (47).
                      Added date/time stamp at start and finish of run - enables tracking of runs better on condor
                      - especially when they are evicted.  (MT)
Jan 25 2007  2.7.Sv46	Clean up compiler warnings. Added two #pragma statements to instruct VC++ and Intel compilers to
                      ignore depreciation warnings. Note - other compilers will ignore these statements.
                      These warnings relate to the fact that we do not use (and probably don't need to) the newer more
                      secure replacements for the ISO standard c++ libraries. This should mean we do not need to suppress
                      warnings (bad practice) when compiling and therefore spot real errors more easily.
                      Also changed two occurrences (lines 670, 1722) of the equivalence operator == where it should have
                      been an assignment operator = Example: changed reset_timeinit==OFF; to reset_timeinit=OFF;  (MT)
Jan 24 2007  2.7.Sv45 Tidy up use of pow() function where the signature of the function call is not clear to the compiler.
                      eg changed pow(hflow,5./3.) to pow(hflow,(double)(5./3.)). Tested - no change to Buscot results.  (MT)
Jan 19 2007  2.7.Sv44	Fixed bug in final call to WriteCheckpoint  (PB)
Jan 09 2007  2.7.Sv43	Corrected some typos (PB)
Dec 18 2006  2.7.Sv42	Changed version number to subversion revision number (MT)
Dec 15 2006  2.7.Sv41	Fixed bug in generation of elev files so these are only created when H<=0.01 (not H<=0).
                      Without this you get significant lateral curvature caused by thin films of water ponding
                      in drying elements. (PB)
############ pre subversion version management ##########################
Feb 21 2006  Version 2.7.5 imported into the Subversion revision control server
Mar 11 2005  Version 2.7.5: Minor fix to force gzip to overwrite files that already
             exist rather than prompting the user. MDW
Mar 11 2005  Version 2.7.4: Profiles of the channel are output at the same time as
             the regular depths - this uses the keyword "profiles" in the parameter file.
             Channel segment numbers have also been added to the profile file. MDW
Mar 11 2005  Version 2.7.3: The time of initial inundation counter can be reset at
             a specified time by using "resettimeinit" in the parameter file, followed
             by the time in seconds. This is useful for calculating time of inundation
             within a flood cycle - e.g. the Amazon has an annual flood cycle and it
             takes a couple of years for the model to stabalise. MDW
Mar 11 2005  Version 2.7.2: Outputs .inittm, .totaltm, .maxtm, .max are now written
             at time of regular output - this is to allow this information to be
             extracted from the middle of a long simulation. MDW
Feb 28 2005  Version 2.7.1: InfilTotalLoss and EvapTotalLoss added to checkpoint
             file - checkpoint version is now 3. MDW
Feb 25 2005  Version 2.7.0: Time variable evaporation added - keyword in paramater
             file is "evaporation" followed by the filename of the evaporation data.
             The file format is similar to the bdy file: line 1 - comment, line 2 - number
             of data then units (seconds, hours or days), line 3+ evaporation (mm/day) and
             time. Note - no account is made of vegetation (i.e. evapotranspiration) yet. MDW
Feb 23 2005  Version 2.6.4: Units for inflow data (.bdy file) can now also be given
             in number of days. MDW
Feb 22 2005  Version 2.6.3: Bug fix in total time of inundation. MDW
Feb 16 2005  Version 2.6.2: Bug fix in checkpointing. MSH
Dec 20 2004  Version 2.6.1: Version number of the executable can be checked without
             running the model. Use tag -version on the command line. MDW
Dec 17 2004  Version 2.6.0: Total time of inundation per cell added - will output a .totaltm file. This
             will be important for assessing things like the Amazon wetlands... MDW
Nov 25 2004  Version 2.5.1: Sim_Time can be specified on the command line - could be useful to
             extend simulations after checkpointing. Use -simtime followed by the
             duration in seconds. MDW
Nov 25 2004  Version 2.5.0: Checkpointing functionality added. See chkpnt.cpp for details. MDW
Sep 10 2004  Version 2.4.0: An alternative ASCII header can be specified in the parameter file
             with the keyword "ascheader" followed by the filename. This allows
             output with geographic coordinates (e.g. lat/lon) rather than model
             coordintates (in metres). This is a rough and ready fix - input
             still needs to be referenced in metres. Ideally I'd like to write
             something that converts geographic input, but... MDW
Aug 25 2004  Version 2.3.0: Output can now be compressed on the fly by using -gzip on the command line.
             Note: this causes a system command to be issued to run gzip for each overpass.
             Obviously gzip must be installed and on the system path. MDW
Jul 29 2004  Version 2.2.0: Static floodplain infiltration added - keyword "infiltration" in parameter file followed by value
             in ms-1. Command line over-ride is -inf. MDW
Jul 01 2004  Version 2.1.0: "qoutput" tag in param file will output 2 Q arrays in addition to depths. Can be used
             in conjuntion with the quiver function in Matlab. MDW
Jun 22 2004: Version 2.0.2: "Tribsoff" tag in param file disabled. Instead, the keyword "Tribs" is added before
             the number of segments in the river file. If this keyword is not present, the number
             of segments is set at 1. MDW
Jun 22 2004: Version 2.0.1: Adaptive time-stepping can be turned off using keyword "adaptoff" in param file.
             In this case, the old QLim version is used. MDW
Jun 08 2004: Version 2.0.0: Converted to adaptive timestepping version with generic linear BC interpolation. NMH, MSH & MDW.
Apr 08 2004: Weir file can now be given on command line using -weir tag.
Apr 05 2004: Water stage over time can now be output for any number of given points.
             Use the keyword "stagefile" in the parameter file, followed by the filename.
             This file should start with the number of stages required, followed by
             the x/y coordinates. Output frequency is the same as the mass balance file. MDW.
Mar 18 2004: Command line can now override parameter file for three variables:
             static floodplain friction, -nfp, static channel friction, -nch, and output
             directory, -dir. MDW.
Feb 13 2004: Fixed direction weirs (culverts) added - use tags NF, EF, SF, WF in weir file. MDW.
Jan 12 2004, Feb 3 2004 Weir equation implementation altered by Rich Dawson
Aug 28 2003: Multiple overpasses now possible - scanned in from file. Keyword "overpassfile"
             followed by the filename in parameter file. File includes the number of overpasses
             on the first line, then a list of the iterations. MDW
Aug 28 2003: Total iteration time added : MDW
Aug 28 2003: Time of initial flood inundation calculated and file written (.inittm) : MDW
Aug 28 2003: Time of maximum depth calculated and file written (.maxtm) : MDW
Apr 01 2003: Point sources in .bci file ("P") with HFIX,HVAR,QFIX or
       QVAR options. Flows given per unit width (m2s-1).
Mar 31 2003: Depth/Elevation output switched off using depthoff/elevoff
                    in .par file.
Mar 31 2003: Maximum water elevation file  (.mxe) written
Mar 31 2003: Fluxes in .bci file for domain boundary (NOT channel)
       are now	given in m2s-1 (flow per unit width) rather than
       m3s-2 (flow per pixel).
Mar 31 2003: Weir locations now entered in Eastings/Northings
Revised 11:21 Feb 6 2003 ----> mdsflood.c.310303_1516
*/


/*
Commands to create the user manual and documentation directly from source code
files. Only update if new files (.html) are added to the user manual. For clarity
the .html files are section2.html where the new page to create is called section.
*/

/*!
\page intro Introduction
\htmlinclude intro2.html

\page install Installation guide
\htmlinclude install2.html

\page input Data requirements, input files and file formats
\htmlinclude input2.html

\page simsetup Setting up a simulation
\htmlinclude simsetup2.html

\page simrun Running a simulation
\htmlinclude simrun2.html

\page ref References and bibliography
\htmlinclude ref2.html

\page examples Example Applications
\htmlinclude examples2.html

\page technote Technical Note
\htmlinclude technote2.html


Important documentation for options added to comments at location in files. Clear documentation implemented
through dOxygen - stored in manual/html/index.html.

*/
