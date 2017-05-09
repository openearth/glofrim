#ifndef GLOBAL_H
#define GLOBAL_H

#include "lisflood.h"

// Variables that are shared between functions and possible exposed using the library interface
extern char t1[80];
extern char tmp_sys_com[255]; // temporary string to hold system command

extern SGCprams *SGCptr;
extern Fnames ParFp;
extern Fnames *Fnameptr;
extern Files Fps;
extern States *Statesptr;
extern Pars *Parptr;
extern Solver ParSolver;
extern States SimStates;
extern Solver *Solverptr;
extern BoundCs *BCptr;
extern Stage *Stageptr;
extern ChannelSegmentType *CSTypePtr;
extern Arrays *Arrptr;
extern vector<int> *RiversIndexVecPtr;
extern int *RiversIndexPtr;
extern vector<ChannelSegmentType> *ChannelSegmentsVecPtr;
extern int *verbose;

extern ChannelSegmentType *ChannelSegments;
extern double Previous_t;      // previous time channel was calculated

extern int tstep_counter;
extern double tstep_channel; // channel timestep
extern double FloodArea;
extern int steadyCount;
/* double loss; //temp variable to keep track of losses since last mass interval */
/* double Comp_time, Model_Comp_Ratio, Model_time_left, Est_Time_Tot, Est_Time_Fin; */

extern Files *Fptr;
extern Stage *Locptr;

#endif /* GLOBAL_H */
