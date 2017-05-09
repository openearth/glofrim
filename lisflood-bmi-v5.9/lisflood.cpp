/*
#####################################################################################
LISFLOOD-FP flood inundation model
#####################################################################################

© copyright Bristol University Hydrology Research Group 2008

webpage -	http://www.ggy.bris.ac.uk/research/hydrology/models/lisflood
contact -	Professor Paul Bates, email: paul.bates@Bristol.ac.uk,
Tel: +44-117-928-9108, Fax: +44-117-928-7878

*/


#include "lisflood.h"
#include "VersionHistory.h"
#include "initialize.h"
#include "finalize.h"
#include "global.h"
//---------------------------------------------------------------------------

extern "C" int init(int, char* []);
extern "C" int init_iterateq();


int main(int argc, char *argv[])
{
  int result = 0;
  // initialize all global variables (defined in global.h)
  result = init(argc, argv);
  if (result) {
    // if result is not 0, exit (don't exit in library code)
    exit(result);
  }
  init_iterateq();
  IterateQ(); // Fnameptr,&Fps,Statesptr,Parptr,Solverptr,BCptr,Stageptr,CSTypePtr,Arrptr,RiversIndexVecPtr, RiversIndexPtr, ChannelSegmentsVecPtr, verbose);
  final_iterateq();
  final();
  return result;
}
//---------------------------------------------------------------------------
