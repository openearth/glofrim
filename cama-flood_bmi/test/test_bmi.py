import logging
import numpy as np
import bmi.wrapper
import ctypes
import sys

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG) # showing messages

def test_library():
    model = bmi.wrapper.BMIWrapper('src/libcama.so')
    model.set_logger(logger)
    #TODO: UPDATE RUNOFF HERE
    sample_runoff = np.fromfile('../sample_runoff/ELSE_GPCC/Roff/Roff____19900101.one', 'f').reshape(7,10)
#    model.set_var('runoff', sample_runoff)
    model.initialize("../test_Elbe/CFM_Elbe/")
    runoff = model.get_var('runoff')
    print runoff
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

    logger.info("type: %s", model.get_var_type('sfcelv'))
    logger.info("rank: %s", model.get_var_rank('sfcelv'))
    logger.info("shape: %s", model.get_var_shape('sfcelv'))
    logger.info("sfcelv: %s", model.get_var('sfcelv'))

    names =[
        "rivinf", "rivdph", "rivvel", "fldinf",
        "flddph", "fldfrc", "fldare", "pthout",
        "pthinf", "sfcelv", "outflw", "storge",
        "rivout_avg", "outflw_avg", "fldout_avg",
        "rivvel_avg", "pthout_avg", "pthflw_avg"
    ]
    for name in names[:2]:
        arr = model.get_var(name)
        logger.info("variable %s:\n%s", name, arr)

    model.update(10)
    before = model.get_var('sfcelv')[0, 0]
    arr = model.get_var('sfcelv').copy()
    arr[0] += 1
    model.set_var('sfcelv', arr)
    after = model.get_var('sfcelv')[0, 0]
    logger.info("before %s, after %s, ok: %s", before, after, (after - before) == 1)
    # logger.info("sfcelv: %s", model.get_var('sfcelv'))

#    while model.get_current_time() < model.get_end_time():
#        model.update(model.get_time_step())
    model.finalize()
