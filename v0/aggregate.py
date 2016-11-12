#!/usr/bin/env python

import logging
import os
import random
import sys

import marcoporoversion

_processname = 'aggregate'
_jobprefix = 'mpa'

def Prerequisites(args, P, mylogger, myhandler, processname):
    'Exit program if some prerequisites are not met.'
    return 0

def AggregateOne(args, P, mylogger, exptid, jobprefix):
    'Aggregate ont stats table metrics into timebucket windows for one experiment.'
    mylogger.info('Processing experiment {0}'.format(exptid))
    cmdpath = os.path.join(args.outdir, '{exptid}_aggregate.sh'.format(exptid=exptid))
    logpath = os.path.join(args.outdir, '{exptid}_aggregate.sh.log'.format(exptid=exptid))
    jobidpath = os.path.join(args.outdir, '{exptid}_aggregate.sh.jobid'.format(exptid=exptid))
    runwithqsub = P.option['aggregate']['resourcestype'] == 'sge'
    if runwithqsub:
        jobname = '{jobprefix}{digits:04d}'.format(jobprefix=jobprefix, digits=random.randint(1, 9999))
        qsubparamL = P.cmdfile_qsubparams('aggregate', jobname, logpath)
    else:
        qsubparamL = []
    cmdL = [
        '{marcoporo} aggregateone'.format(marcoporo=P.option['program']['marcoporo']),
        ' -bin {dirpath}'.format(dirpath=args.bin),
        ' -profile {filepath}'.format(filepath=args.profile),
        ' -config {filepath}'.format(filepath=args.config),
        ' -exptid {exptid}'.format(exptid=exptid),
        ' -indir {dirpath}'.format(dirpath=args.indir),
        ' -timebucket {hours}'.format(hours=args.timebucket),
        ' -maxrunlen {hours}'.format(hours=args.maxrunlen),
        ' -outdir {dirpath}'.format(dirpath=args.outdir)]
    cmd = ' \\\n'.join(cmdL)
    rv  = P.cmdfile_create(qsubparamL, cmd, cmdpath)
    if rv != 0:
        mylogger.warning('Failed to create command file ({0})'.format(cmdpath))
        return 1
    if args.execjobs:
        rv, errmsg = P.cmdfile_run(cmdpath, jobidpath, runwithqsub)
        if rv != 0:
            mylogger.warning('Job failed ({0}), errmsg={1}'.format(cmdpath, rrmsg))
            return 1
    else:
        mylogger.info('Not executing job {0}'.format(cmdpath))
    return 0

def Aggregate(args, P, mylogger, myhandler, processname, exptidL, E):
    '''
    Create one shell script per exptid containing qsub options and
    aggregateone command-line, schedule by executing with qsub and
    save the jobid to a file.
    '''
    for exptid in exptidL:
        AggregateOne(args, P, mylogger, exptid, _jobprefix)
    return 0

def Process(args, P, mylogger, myhandler, processname):
    'Run aggregateone for all experiments in experiments.txt file.'
    exptidL, E = P.expt_read(args.experiments)
    Aggregate(args, P, mylogger, myhandler, processname, exptidL, E)
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    'Execute this sub-tool.'
    mylogger.info('Started')
    args.execjobs = P.str_2bool(args.execjobs)
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname)
    mylogger.info('Finished')
    return 0
