#!/usr/bin/env python

import logging
import numpy as np
import os
import sys
import time

import marcoporoversion

_processname = 'marginalign'
MAXALIGNMENTLENGTHTOSAMPLE=10000000

def Prerequisites(args, P, mylogger, myhandler, processname):
    'Exit program if some prerequisites are not met.'
    return 0

def MarginAlign_AlignEMNo(args, P, mylogger, myhandler, processname, exptid, E, fastqpath, exptjobtreedir):

    sampath = os.path.join(args.outdir, os.path.basename(fastqpath).replace('.fastq', '')+'_BWANoRealign.sam')
    if not args.overwrite and os.path.exists(sampath) and os.path.getsize(sampath) > 0:
        mylogger.info('Skipping marginalign BWANoRealign {fastqpath}'.format(fastqpath=fastqpath))
        return sampath
    cmdjobtreedir = os.path.join(exptjobtreedir, os.path.basename(fastqpath).replace('.fastq', '')+'_BWANoRealign')
    if os.path.exists(cmdjobtreedir):
        # Could have done a recursive removal of this directory, but was worried it could be error-prone and delete data
        mylogger.error('Skipping marginalign BWANoRealign - jobtree dir already exists and must be manually removed before re-running {jobtreedir}'.format(jobtreedir=cmdjobtreedir))
        return None
    cmd = '{marginalign} {infastq} {inreffasta} {outsam} --jobTree {jobtreedir} --bwa --noRealign'.format(
        marginalign=P.option['program']['marginalign'],
        infastq=fastqpath,
        inreffasta=P.option['program']['refpath'],
        outsam=sampath,
        jobtreedir=cmdjobtreedir
    )
    mylogger.info('Running marginalign BWANoRealign - {cmd}'.format(cmd=cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed marginalign BWANoRealign: returnval={0}, stdout={1}, stderr={2}'.format(rv, ro, re))
        return None
    return sampath

def MarginAlign_AlignEMYes(args, P, mylogger, myhandler, processname, exptid, E, fastqpath, exptjobtreedir):
    sampath = os.path.join(args.outdir, os.path.basename(fastqpath).replace('.fastq', '')+'_BWAMEM10M.sam')
    if not args.overwrite and os.path.exists(sampath) and os.path.getsize(sampath) > 0:
        mylogger.info('Skipping marginalign BWAMEM10M {fastqpath}'.format(fastqpath=fastqpath))
        return sampath
    cmdjobtreedir = os.path.join(exptjobtreedir, os.path.basename(fastqpath).replace('.fastq', '')+'_BWAMEM10M')
    if os.path.exists(cmdjobtreedir):
        # Could have done a recursive removal of this directory, but was worried it could be error-prone and delete data
        mylogger.error('Skipping marginalign BWAMEM10M - jobtree dir already exists and must be manually removed before re-running {jobtreedir}'.format(jobtreedir=cmdjobtreedir))
        return None
    cmd = '{marginalign} {infastq} {inreffasta} {outsam} --jobTree {jobtreedir} --bwa --em --maxAlignmentLengthToSample={maxalignlen}'.format(
        marginalign=P.option['program']['marginalign'],
        infastq=fastqpath,
        inreffasta=P.option['program']['refpath'],
        outsam=sampath,
        jobtreedir=cmdjobtreedir,
        maxalignlen=MAXALIGNMENTLENGTHTOSAMPLE
    )
    mylogger.info('Running marginalign BWAMEM10M - {cmd}'.format(cmd=cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed marginalign BWAMEM10M: returnval={0}, stdout={1}, stderr={2}'.format(rv, ro, re))
        return None
    return sampath

def MarginAlign_Stats(args, P, mylogger, myhandler, processname, exptid, E, fastqpath, sampath):

    if not os.path.exists(sampath):
        mylogger.error('Skipping marginStats - no input file {sampath}'.format(sampath=sampath))
        return False
    statspath = os.path.join(args.outdir, os.path.basename(sampath).replace('.sam', '.stats'))
    if not args.overwrite and os.path.exists(statspath) and os.path.getsize(statspath) > 0:
        mylogger.info('Skipping marginStats {sampath}'.format(sampath=sampath))
        return True

    cmd = '{marginstats} {sampath} {fastqpath} {inreffasta} --readIdentity --readCoverage --mismatchesPerAlignedBase --deletionsPerReadBase --insertionsPerReadBase --localAlignment --readLength --printValuePerReadAlignment'.format(
        marginstats=P.option['program']['marginstats'],
        sampath=sampath,
        fastqpath=fastqpath,
        inreffasta=P.option['program']['refpath']
    )
    mylogger.info('Running marginalign BWAMEM10M - {cmd}'.format(cmd=cmd))
    print('Debug: {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed marginalign BWAMEM10M: returnval={0}, stdout={1}, stderr={2}'.format(rv, ro, re))
        return False
    with open(statspath, 'w') as out_fp:
        out_fp.write('{0}\n'.format(ro))
    return True

def Process_Experiment(args, P, mylogger, myhandler, processname, exptid, E, exptjobtreedir):
    'Run marginalign align+stats withEM+withoutEM for each FASTQ file for the experiment.'
    for fastqfile in [f for f in os.listdir(args.extractdir) if f.startswith(exptid) and f.endswith('fastq')]:
        fastqpath = os.path.join(args.extractdir, fastqfile)
        if not os.path.exists(exptjobtreedir):
            os.makedirs(exptjobtreedir)
        sam_norealign_path = MarginAlign_AlignEMNo(args, P, mylogger, myhandler, processname, exptid, E, fastqpath, exptjobtreedir)
        sam_bwaem10m_path = MarginAlign_AlignEMYes(args, P, mylogger, myhandler, processname, exptid, E, fastqpath, exptjobtreedir)
        if sam_norealign_path is not None:
            MarginAlign_Stats(args, P, mylogger, myhandler, processname, exptid, E, fastqpath, sam_norealign_path)
        if sam_bwaem10m_path is not None:
            MarginAlign_Stats(args, P, mylogger, myhandler, processname, exptid, E, fastqpath, sam_bwaem10m_path)
    return 0

def Process(args, P, mylogger, myhandler, processname, jobtreedir):
    'Run nanook extract, align and analyse for each experiment in the experiment file.'
    exptidL, E = P.expt_read(args.experiments)
    for exptid in exptidL:
        exptjobtreedir = os.path.join(jobtreedir, exptid)
        Process_Experiment(args, P, mylogger, myhandler, processname, exptid, E, exptjobtreedir)

def run(parser, args, P, mylogger, myhandler, argv):
    'Execute this sub-tool.'
    mylogger.info('Started')
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    jobtreedir = os.path.join(args.outdir, 'jobtree')
    if not os.path.exists(jobtreedir):
        os.makedirs(jobtreedir)
    args.overwrite = P.str_2bool(args.overwrite)
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname, jobtreedir)
    mylogger.info('Finished')
    return 0
