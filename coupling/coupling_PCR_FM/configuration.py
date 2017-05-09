# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 13:32:04 2016

This set of functions executes the actual model functions to manipulate model
arrays and to compute the required volumes and depths, as well as updating
water heights accordingly. 

Please note: currently only to be used for 1way-coupling!!!

@author: Jannis M.Hoch (j.m.hoch@uu.nl), Department of Physical Geography, Utrecht University
"""

import ConfigParser
import datetime
import os
import optparse

class Configuration(object):
    
    def __init__(self):
        
        object.__init__(self)
        
    def parse_configuration_file(self, modelFileName):
        
        config = ConfigParser.ConfigParser()
        config.optionxform = str
        config.read(modelFileName)

        # all sections provided in the configuration/ini file
        self.allSections  = config.sections()

        # read all sections 
        for sec in self.allSections:
            vars(self)[sec] = {}                               # example: to instantiate self.globalOptions 
            options = config.options(sec)                      # example: logFileDir
            for opt in options:
                val = config.get(sec, opt)                     # value defined in every option 
                self.__getattribute__(sec)[opt] = val          # example: self.globalOptions['logFileDir'] = val
