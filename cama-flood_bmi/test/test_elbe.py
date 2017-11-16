import logging
import numpy as np
import bmi.wrapper
import ctypes
import sys

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG) # showing messages

if __name__ == "__main__":
    model = bmi.wrapper.BMIWrapper('../src/libcama.so')
    model.set_logger(logger)
    model.initialize("../../test_Elbe/CFM_Elbe/")

    #TODO: UPDATE RUNOFF HERE, later obtained from PCR
    sample_runoff = np.fromfile('../../sample_runoff/ELSE_GPCC/Roff/Roff____19900101.one', 'f').reshape(7,10)
#    print(sample_runoff)
    model.set_var('runoff', sample_runoff)

# Check updated runoff
#    runoff = model.get_var('runoff').reshape(7,10)
#    print(runoff)
#    print(model.get_var_shape('runoff'))
#    diff = sample_runoff - runoff
#    print(diff.max(), diff.min())
#    sys.exit()

    start_time = model.get_start_time()
    current_time = model.get_current_time()
    end_time = model.get_end_time()
    logger.info(
        "start_time: %s, current_time %s, end_time: %s",
        start_time,
        current_time,
        end_time
    )
    logger.info("timestep 1")

    model.update(10)
    current_time = model.get_current_time()
    time_step = model.get_time_step()
    logger.info(
        "start_time: %s, current_time %s, timestep %s",
        start_time,
        current_time,
        time_step
    )

# For the model runs, input runoff should be updated from PCR by set_var function
#    while model.get_current_time() < model.get_end_time():
#        model.update(model.get_time_step())
    model.finalize()
