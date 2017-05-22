#!/usr/bin/env python

import logging
import numpy as np
import os
import subprocess as sp
import sys
import time

import marcoporoversion

_processname = 'mapwithbwa'

def Prerequisites(args, P, mylogger, myhandler, processname):
    'Exit program if some prerequisites are not met.'
    return 0

def MapWithBwaMem(args, P, mylogger, myhandler, processname, exptid, E, fastqpath):
    'Run bwa and samtools to get coord-sorted BAM file, return bampath (or None on error).'
    # Variables
    bampath = os.path.join(args.bwamemdir, os.path.basename(fastqpath).replace('.fastq', '.bam'))
    baipath = bampath + '.bai'
    outbampathstem = bampath.replace('.bam', '')
    outbamerrpath = bampath.replace('.bam', '.err')
    # Return if nothing to be done
    if not args.overwrite \
    and os.path.exists(bampath) and os.path.getsize(bampath) > 0 \
    and os.path.exists(baipath) and os.path.getsize(baipath) > 0:
        #mylogger.info('Skipping mapwithbwa - BAM exists ({0})'.format(bampath))
        return bampath, False
    # Generated sorted BAM and the BAI file
    #mylogger.info('Running mapwithbwa - BAM exists ({0})'.format(bampath))
    cmd1 = '{bwa_prog} mem -x ont2d -M -t {threads} {reffasta} {fastqpath}'.format(
        bwa_prog=P.option['program']['bwa'],
        threads=P.option['program']['server_threads'],
        reffasta=P.option['program']['refpath'],
        fastqpath=fastqpath)
    cmd2 = '{samtools_prog} view -b -S -'.format(
        samtools_prog=P.option['program']['samtools'])
    cmd3 = '{samtools_prog} sort - {outbampathstem}'.format(
        samtools_prog=P.option['program']['samtools'],
        outbampathstem=outbampathstem)
    p1 = sp.Popen(cmd1.split(), stdout=sp.PIPE, stderr=sp.PIPE)
    p2 = sp.Popen(cmd2.split(), stdin=p1.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    p3 = sp.Popen(cmd3.split(), stdin=p2.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
    p1.stdout.close() # mimic pipe
    p2.stdout.close() # mimic pipe
    ro, re = p3.communicate()	# Does not return returncode, only stdout and stderr text as a tuple
    if not os.path.exists(bampath) and os.path.getsize(bampath) > 0:
        mylogger.error('Failed to map reads with bwa - file missing or empty ({0})'.format(bampath))
        return None, False
    with open(outbamerrpath, 'w') as out_fp:
        out_fp.write('{0}\n'.format(re))
    # Always re-generate the BAI file if you generated the BAM file.
    cmd = '{samtools_prog} index {bampath}'.format(
        samtools_prog=P.option['program']['samtools'],
        bampath=bampath)
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('Failed to create BAM index file ({0})'.format(baipath))
        return None, False
    return bampath, True

_readtypeL = ['1T', '1C', '2D']
_readclassL = ['pass', 'fail']

def RunPoremapstats(args, P, mylogger, myhandler, processname, exptid, E, bampath):
    'Run poremapstats on the BAM file.'
    for readtype in _readtypeL:
        for readclass in _readclassL:
            #outprefix = exptid + '_' + readtype + '_' + readclass
            #set bampath
            outprefix = os.path.basename(bampath).replace('.bam', '')
            if not os.path.exists(bampath):
                continue
            initstatspath = os.path.join(args.bwamemdir, outprefix+'_initstats.txt')
            readstatspath = os.path.join(args.bwamemdir, outprefix+'_readstats.txt')
            runstatspath = os.path.join(args.bwamemdir, outprefix+'_runstats.txt')
            if not args.overwrite \
            and os.path.exists(initstatspath) and os.path.getsize(initstatspath) > 0 \
            and os.path.exists(readstatspath) and os.path.getsize(readstatspath) > 0 \
            and os.path.exists(runstatspath) and os.path.getsize(runstatspath) > 0:
                mylogger.info('Not running poremapstats - files already exist ({0}, {1}, {2})'.format(initstatspath, readstatspath, runstatspath))
                return
            logpath = os.path.join(args.bwamemdir, outprefix + '_poremapstats.log') 
            if readtype == '1T':
                readtypelongfmt = 'temp'
            elif readtype == '1C':
                readtypelongfmt = 'comp'
            elif readtype == '2D':
                readtypelongfmt = '2d'
            else:
                readtypelongfmt = 'unknown'
            cmd = '{poremapstats_prog}'.format(poremapstats_prog=P.option['program']['poremapstats'])
            cmd += ' --bindir {bindir}'.format(bindir=args.bin)
            cmd += ' --profilepath None'
            cmd += ' --runid {runid}'.format(runid=exptid)
            cmd += ' --readtype {readtype}'.format(readtype=readtypelongfmt)
            cmd += ' --readclass {readclass}'.format(readclass=readclass)
            cmd += ' --datatype minion'
            cmd += ' --mapprog bwa'
            cmd += ' --mapparams "-x ont2d -M"'
            cmd += ' --alignclasspath None'
            cmd += ' --readsbam {bampath}'.format(bampath=bampath)
            cmd += ' --targetrefpath {targetfasta}'.format(targetfasta=P.option['program']['targetpath'])
            cmd += ' --controlrefpath {controlfasta}'.format(controlfasta=P.option['program']['controlpath'])
            cmd += ' --outdir {outdir}'.format(outdir=args.bwamemdir)
            cmd += ' --outprefix {outprefix}'.format(outprefix=outprefix)
            cmd += ' --savealignments False'
            cmd += ' --fastalinewidth {fastalinewidth}'.format(fastalinewidth=P.option['program']['fastalinewidth'])
            cmd += ' --overwrite {overwrite}'.format(overwrite=args.overwrite)
            mylogger.info('Running {0}'.format(cmd))
            rv, ro, re = P.sys_exec(cmd)
            if rv != 0:
                mylogger.error('Failed to run poremapstats: returnval={0}, stdout={1}, stderr={2}'.format(rv, ro, re))
                return False
            with open(logpath, 'w') as out_fp:
                if len(ro):
                    out_fp.write(ro)
                if len(re):
                    out_fp.write(re) 
    return True

def MapReads(args, P, mylogger, myhandler, processname, exptid, E):
    'Run bwa mem to get sam file, then use samtools to sort by coordinate and get bam index file.'
    for fastqfile in [f for f in os.listdir(args.extractdir) if f.startswith(exptid) and f.endswith('fastq')]:
        fastqpath = os.path.join(args.extractdir, fastqfile)
        bampath, needtorunstats = MapWithBwaMem(args, P, mylogger, myhandler, processname, exptid, E, fastqpath)
        if bampath is not None and needtorunstats:
            #mylogger.info('Running poremapstats for exptid {0}'.format(exptid))
            RunPoremapstats(args, P, mylogger, myhandler, processname, exptid, E, bampath)
        else:
            pass
            #mylogger.info('Skipping poremapstats for exptid {0}'.format(exptid))
    return 0

def Process(args, P, mylogger, myhandler, processname):
    'Run nanook extract, align and analyse for each experiment in the experiment file.'
    exptidL, E = P.expt_read(args.experiments)
    for exptid in exptidL:
        mylogger.info('Processing experiment {0}'.format(exptid))
        MapReads(args, P, mylogger, myhandler, processname, exptid, E)
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    'Execute this sub-tool.'
    mylogger.info('Started')
    if not os.path.exists(args.bwamemdir):
        os.makedirs(args.bwamemdir)
    args.overwrite = P.str_2bool(args.overwrite)
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname)
    mylogger.info('Finished')
    return 0
