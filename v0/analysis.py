#!/usr/bin/env python

import logging
import os
import random
import sys

import marcoporoversion

_processname = 'analyse'
_jobprefix = 'mpn'

def Prerequisites(args, P, mylogger, myhandler, processname):
    'Exit program if some prerequisites are not met.'
    return 0

def Analyse(args, P, mylogger, myhandler, processname, exptidL, E):
    return 0

def Process(args, P, mylogger, myhandler, processname):
    'Run aggregateone for all experiments in experiments.txt file.'
    exptidL, E = P.expt_read(args.experiments)
    Analyse(args, P, mylogger, myhandler, processname, exptidL, E)
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    'Execute this sub-tool.'
    mylogger.info('Started')
    args.execjobs = P.str_2bool(args.execjobs)
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname)
    mylogger.info('Finished')
    return 0
