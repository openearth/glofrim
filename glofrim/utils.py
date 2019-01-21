import sys
import logging
import logging.handlers

_formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

def setlogger(logfilename=None, loggername='glofrim', thelevel=logging.INFO, show_in_console=True, formatter=_formatter):
    """
    Set-up the logging system and return a logger object. Exit if this fails
    """

    try:    
        #create logger
        logger = logging.getLogger(loggername)
        # remove any existing handlers
        while logger.handlers:
            i = logger.handlers[0]
            logger.removeHandler(i)
            i.flush()
            i.close()
        logger.setLevel(thelevel)
        # console logger
        if show_in_console:
            console = logging.StreamHandler()
            console.setLevel(logger.level)
            console.setFormatter(formatter)
            logger.addHandler(console)
        #add filehandler to logger
        if logfilename is not None:
            add_file_handler(logger, logfilename)
        return logger
    except IOError:
        print("ERROR: Failed to initialize logger with logfile: " + logfilename)
        sys.exit(2)

def add_file_handler(logger, logfilename, formatter=_formatter):
    ch = logging.FileHandler(logfilename, mode='w')
    ch.setFormatter(formatter)
    ch.setLevel(logger.level)
    logger.addHandler(ch)
    logger.debug("File logging to {}".format(logfilename))

def closelogger(logger):
    while logger.handlers:
        i = logger.handlers[0]
        logger.removeHandler(i)
        i.flush()
        i.close()
    return logger

def close_with_error(logger, msg):
    logger.error(msg)
    closelogger(logger)
    sys.exit(1)