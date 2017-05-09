// SampleCppLibrary.cpp : Defines the exported functions for the DLL application.
//

#include <cstdio>
#include <string>
#include <sstream>
#include "bmi.h"
#include "lisflood.h"
#include "global.h"
#include "initialize.h"
#include "update.h"
#include "finalize.h"

// define some arrays for exchange
// list of variables to be potentially exposed via BMI can be found in lisflood.h
enum VARIABLE_LABEL {
  H,
  DEM,
  Qx,
  Qy,
  rain,
  dx,
  dx_sqrt,
  dy,
  dA,
  SGCwidth,
  SGCQin,
  blx,
  bly,
  BC_numPS,
  BC_xpi,
  BC_ypi,
  // if variable is not found
  NOVAR
};

VARIABLE_LABEL variable2label (std::string const& variable) {
  if (variable == "H") return H;
  if (variable == "DEM") return DEM;
  if (variable == "Qx") return Qx;
  if (variable == "Qy") return Qy;
  if (variable == "rain") return rain;
  if (variable == "dx") return dx;
  if (variable == "dx_sqrt") return dx_sqrt;
  if (variable == "dy") return dy;
  if (variable == "dA") return dA;
  if (variable == "SGCwidth") return SGCwidth;
  if (variable == "SGCQin") return SGCQin;
  if (variable == "blx") return blx;
  if (variable == "bly") return bly;
  if (variable == "BC.numPS") return BC_numPS;
  if (variable == "BC.xpi") return BC_xpi;
  if (variable == "BC.ypi") return BC_ypi;
  //if (variable == "FArea") return FArea;
  else return NOVAR;
};

/* Store callback */
Logger logger = NULL;

/* Logger function */
void _log(Level level, std::string msg);

extern "C" int init(int, const char*[]);
extern "C" int init_iterateq();

extern "C" {
  BMI_API int initialize(const char *config_file)
  {
    const char * argv[] = {
      "bmi",
      "-v",
      config_file
    };
    int result = init(3, argv);
    init_iterateq();
    return result;
  }

  BMI_API int update(double dt)
  {
    // set timestep and do:
    double oldStep = Solverptr->Tstep;
    if (dt != -1) {
      Solverptr->Tstep = dt;
    }
    iterateq_step();
    // restore dt to default
    Solverptr->Tstep = oldStep;
    return 0;
  }

  BMI_API int finalize()
  {
    final_iterateq();
    final();
    return 0;
  }

  BMI_API void get_start_time(double *t)
  {
    *t = 0;
  }

  BMI_API void get_end_time(double *t)
  {
    *t = Solverptr->Sim_Time;
  }

  BMI_API void get_current_time(double *t)
  {
    *t = Solverptr->t;
  }

  BMI_API void get_time_step(double *dt)
  {
    *dt = Solverptr->Tstep;
  }

  BMI_API void get_var(const char *name, void **ptr)
  {
    VARIABLE_LABEL label = variable2label(name);
    switch (label) {
    case H: {
      *ptr = (void*)(Arrptr->H);
      break;
    }
    case DEM: {
      *ptr = (void*)(Arrptr->DEM);
      break;
    }
    case Qx: {
      *ptr = (void*)(Arrptr->Qx);
      break;
    }
    case Qy: {
      *ptr = (void*)(Arrptr->Qy);
      break;
    }
    case rain: {
      *ptr = (void*)(Arrptr->rain);
      break;
    }
    case dx: {
      *ptr = &(Parptr->dx);
      break;
    }
    case dx_sqrt: {
      *ptr = &(Parptr->dx_sqrt);
      break;
    }
    case dy: {
      *ptr = &(Parptr->dy);
      break;
    }
    case dA: {
      *ptr = (void*)(Arrptr->dA);
      break;
    }
    case SGCwidth: {
      *ptr = (void*)(Arrptr->SGCwidth);
      break;
    }
    case SGCQin: {
      *ptr = (void*)(Arrptr->SGCQin);
      break;
    }
    case blx: {
      *ptr = &(Parptr->blx);
      break;
    }
    case bly: {
      *ptr = &(Parptr->bly);
      break;
    }
    case BC_numPS: {
      *ptr = &(BCptr->numPS);
      break;
    }
    case BC_xpi: {
      *ptr = BCptr->xpi;
      break;
    }
    case BC_ypi: {
      *ptr = BCptr->ypi;
      break;
    }
    case NOVAR: {
      break;
    }
    default:
      break;
    }
  }
  BMI_API void get_var_rank(const char *name, int *rank)
  {
    VARIABLE_LABEL label = variable2label(name);
    switch (label) {
    case H:
    case DEM:
    case Qx:
    case Qy:
    case dA:
    case SGCwidth:
    case SGCQin:
    case rain:
    {
      *rank = 2;
      break;
    }
    case BC_xpi:
    case BC_ypi:
    {
      *rank = 1;
      break;
    }
    case dx:
    case dx_sqrt:
    case dy:
    case blx:
    case bly:
    case BC_numPS:
    case NOVAR: {
      *rank = 0;
      break;
    }
    default:
      *rank = 0;
      break;
    }
  }

  BMI_API void get_var_shape(const char *name, int *shape)
  {
    VARIABLE_LABEL label = variable2label(name);
    switch (label) {
    case H:
    case DEM:
    case dA:
    case SGCwidth:
    case SGCQin:
    case Qx:
    case Qy:
    case rain:
    //case FArea:
    {
      // TODO: check/fortran/c memory
      shape[1] = Parptr->xsz;
      shape[0] = Parptr->ysz;
      break;
    }
    case BC_xpi:
    case BC_ypi:
    {
      shape[0] = BCptr->numPS;
      break;
    }
    case dx:
    case dx_sqrt:
    case dy:
    case blx:
    case bly:
    case BC_numPS:
    case NOVAR: {
      break;
    }
    default:
      break;
    }
  }
  BMI_API void get_var_type(const char *name, char *type)
  {
    VARIABLE_LABEL label = variable2label(name);
    switch (label) {
    case H:
    case DEM:
    case Qx:
    case Qy:
    case rain:
    case dx:
    case dx_sqrt:
    case dy:
    case dA:
    case SGCwidth:
    case SGCQin:
    case blx:
    case bly:
    {
      strcpy(type, "double");
      break;
    }
    case BC_xpi:
    case BC_ypi:
    case BC_numPS:
    {
      strcpy(type, "int");
      break;
    }
    case NOVAR: {
      break;
    }
    default:
      break;
    }

  }

  BMI_API void set_logger(Logger callback)
  {
    Level level = INFO;
    std::string msg = "Logging attached to cxx model";
    logger = callback;
    logger(level, msg.c_str());
  }
}

void _log(Level level, std::string msg) {
  if (logger != NULL) {
    logger(level, msg.c_str());
  }
}

// placeholder function, all dll's need a main.. in windows only
#if defined _WIN32
void main()
{
}
#endif
