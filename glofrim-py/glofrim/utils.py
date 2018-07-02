import sys
import logging
import logging.handlers

def setlogger(logfilename=None, loggername='glofrim', thelevel=logging.INFO):
    """
    Set-up the logging system and return a logger object. Exit if this fails
    """

    try:    
        #create logger
        logger = logging.getLogger(loggername)
        if not isinstance(thelevel, int):
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(thelevel)
        #create formatter
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(module)s - %(levelname)s - %(message)s")
        # console logger
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        console.setFormatter(formatter)
        logger.addHandler(console)
        #add filehadnler to logger
        if logfilename is not None:
            ch = logging.FileHandler(logfilename,mode='w')
            ch.setFormatter(formatter)
            ch.setLevel(logging.DEBUG)
            logger.addHandler(ch)
            logger.debug("File logging to {}".format(logfilename))
        return logger
    except IOError:
        print("ERROR: Failed to initialize logger with logfile: " + logfilename)
        sys.exit(2)

def closeLogger(logger, ch):
    logger.removeHandler(ch)
    ch.flush()
    ch.close()
    return logger, ch

def close_with_error(logger, ch, msg):
    logger.error(msg)
    logger, ch = closeLogger(logger, ch)
    del logger, ch
    sys.exit(1)