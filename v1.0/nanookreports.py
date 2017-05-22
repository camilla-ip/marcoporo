#!/usr/bin/env python

import logging
import numpy as np
import os
import sys
import time

import marcoporoversion

_processname = 'nanookreports'
_NANOOKEXTRACTMINPCT = 99.0

def Prerequisites(args, P, mylogger, myhandler, processname):
    'Exit program if some prerequisites are not met.'
    return 0

def Nanook_SetupOutdirs(args, P, mylogger, myhandler, processname, exptid, E):
    'To prevent writing of output to the sampledir, create outdir with symlinks to the input data.'
    # Variables
    exptoutdir = os.path.join(args.outdir, exptid)
    # Create output directory for experiment
    try:
        if not os.path.exists(exptoutdir):
            os.makedirs(exptoutdir)
    except:
        mylogger.error('Failed to create nanook outdir ({0})'.format(exptoutdir))
        return None
    # Create symbolic link to reads for this experiment
    srcdir = os.path.join(E[exptid]['dirpath'], 'reads')
    dstdir = os.path.join(exptoutdir, 'reads')
    try:
        if not os.path.exists(dstdir):
            os.symlink(srcdir, dstdir)
    except:
        mylogger.error('Failed to create symbolic link to source data ({0} -> {1})'.format(dstdir, dstdir))
        return None
    return exptoutdir

def Nanook_ExtractFasta(args, P, mylogger, myhandler, processname, exptid, E, exptindir):
    'Run nanook extract to get one FASTA from each FAST5 file in the experiment.'
    # Output already exists if subdirs required for the 1D or 2D library exist,
    # as well as at least 99% of the FAST5 reads have been extracted
    fast5keyD = {}
    indir = os.path.join(E[exptid]['dirpath'], 'reads', 'downloads')
    for readclass in [x for x in os.listdir(indir) if os.path.isdir(os.path.join(indir, x))]:
        for fast5file in [f for f in os.listdir(os.path.join(indir, readclass)) if f.endswith('fast5')]:
            fast5keyD[fast5file] = 1
    fast5cnt = len(fast5keyD.keys())
    fastakeyD = {}
    indir = os.path.join(exptindir, 'fasta')
    if os.path.exists(indir):
        for readclass in [x for x in os.listdir(indir) if os.path.isdir(os.path.join(indir, x))]:
            dir = os.path.join(indir, readclass)
            for readtype in [x for x in os.listdir(dir) if os.path.isdir(os.path.join(dir, x))]:
                fastathingL = [f for f in os.listdir(os.path.join(dir, readtype)) if f.endswith('fasta')] # Fmt: READINFO.fast5_BaseCalled_Template.fasta
                for fastathing in fastathingL:
                    S = fastathing.split('.')
                    fastafile = '.'.join(S[:-2]) + '.' + S[-1].split('_')[0]
                    fastakeyD[fastafile] = 1
    fastacnt = len(fastakeyD.keys())
    pctextracted = round(fastacnt / float(fast5cnt) * 100.0 if fastacnt else 0, 2)
    if not args.overwrite:
        if pctextracted >= _NANOOKEXTRACTMINPCT:
            mylogger.info('Skipping nanook extract FASTA for exptid {0} - {1}% of files already processed'.format(exptid, pctextracted))
            return
    else:
        if fastacnt > 0:
            mylogger.info('Overwriting nanook extract FASTA output files for exptid {0} - {1} fasta output files'.format(exptid, fastacnt))
    # Run the extraction command
    cmd = '{nanook} extract -s {sampledir} -f {readsdir} -a -basecallindex {instanceN} -printpath -t {threads}'.format(
        nanook=P.option['program']['nanook'],
        sampledir=exptindir,
        readsdir=os.path.join('reads', 'downloads'),
        instanceN=E[exptid]['instanceN'],
        threads=P.option['program']['server_threads']
    )
    mylogger.info('Running {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('nanook extract FASTA failed: stdout={0}, stderr={1}'.format(ro, re))
        return False
    return True

def Nanook_ExtractFastq(args, P, mylogger, myhandler, processname, exptid, E, exptindir):
    'Run nanook extract to get one FASTQ from each FAST5 file in the experiment.'
    # Output already exists if subdirs required for the 1D or 2D library exist,
    # as well as at least 99% of the FAST5 reads have been extracted
    fast5keyD = {}
    indir = os.path.join(E[exptid]['dirpath'], 'reads', 'downloads')
    for readclass in [x for x in os.listdir(indir) if os.path.isdir(os.path.join(indir, x))]:
        for fast5file in [f for f in os.listdir(os.path.join(indir, readclass)) if f.endswith('fast5')]:
            fast5keyD[fast5file] = 1
    fast5cnt = len(fast5keyD.keys())
    fastqkeyD = {}
    indir = os.path.join(exptindir, 'fastq')
    if os.path.exists(indir):
        for readclass in [x for x in os.listdir(indir) if os.path.isdir(os.path.join(indir, x))]:
            dir = os.path.join(indir, readclass)
            for readtype in [x for x in os.listdir(dir) if os.path.isdir(os.path.join(dir, x))]:
                fastqthingL = [f for f in os.listdir(os.path.join(dir, readtype)) if f.endswith('fastq')] # Fmt: READINFO.fast5_BaseCalled_Template.fastq
                for fastqthing in fastqthingL:
                    S = fastqthing.split('.')
                    fastqfile = '.'.join(S[:-2]) + '.' + S[-1].split('_')[0]
                    fastqkeyD[fastqfile] = 1
    fastqcnt = len(fastqkeyD.keys())
    pctextracted = round(fastqcnt / float(fast5cnt) * 100.0 if fastqcnt else 0, 2)
    if not args.overwrite:
        if pctextracted >= _NANOOKEXTRACTMINPCT:
            mylogger.info('Skipping nanook extract FASTQ for exptid {0} - {1}% of files already processed'.format(exptid, pctextracted))
            return
    else:
        if fastqcnt > 0:
            mylogger.info('Overwriting nanook extract FASTQ output files for exptid {0} - {1} fastq output files'.format(exptid, fastqcnt))
    # Run the extraction command
    cmd = '{nanook} extract -s {sampledir} -f {readsdir} -q -basecallindex {instanceN} -printpath -t {threads}'.format(
        nanook=P.option['program']['nanook'],
        sampledir=exptindir,
        readsdir=os.path.join('reads', 'downloads'),
        instanceN=E[exptid]['instanceN'],
        threads=P.option['program']['server_threads']
    )
    mylogger.info('Running {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('nanook extract FASTQ failed: stdout={0}, stderr={1}'.format(ro, re))
        return False
    return True

def Nanook_Align(args, P, mylogger, myhandler, processname, exptid, E, exptindir):
    'Run nanook align to use bwa mem -x ont2d -M to align reads.'
    # Count number of missing FASTA.sam files. If zero are missing and overwrite is False, skip this step.

    cmd = '{nanook} align -s {sampledir} -r {referencefasta} -aligner bwa -alignerparams "-x ont2d -M" -showaligns -t {threads}'.format(
        nanook=P.option['program']['nanook'],
        sampledir=exptindir,
        referencefasta=P.option['program']['refpath'],
        threads=P.option['program']['server_threads']
    )
    mylogger.info('Running {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('nanook align failed: stdout={0}, stderr={1}'.format(ro, re))
        return False
    # Remove any .sam.last files - not sure why they are here at all if requesting alignment with BWA only.
    fastadir = os.path.join(exptindir, 'fasta')
    if os.path.exists(fastadir):
        for readclass in [x for x in os.listdir(fastadir) if os.path.isdir(os.path.join(fastadir, x))]:
            readclassdir = os.path.join(fastadir, readclass)
            for readtype in [x for x in os.listdir(readclassdir) if os.path.isdir(os.path.join(readclassdir, x))]:
                readtypedir = os.path.join(fastadir, readclass, readtype)
                for fastafile in [f for f in os.listdir(readtypedir) if f.endswith('fasta')]:
                    samlastpath = os.path.join(exptindir, 'bwa', readclass, readtype, fastafile) + '.sam.last'
                    if os.path.exists(samlastpath):
                        try:
                            #mylogger.info('Removing {0}'.format(samlastpath))
                            os.remove(samlastpath)
                        except:
                            mylogger.warning('Failed to remove {0}'.format(samlastpath))
    return True

def Nanook_Analyse(args, P, mylogger, myhandler, processname, exptid, E, exptindir):
    'Run nanook analyse to generate statistics on the experiment.'
    # Generate PDF report for pass+fail reads
    cmd = '{nanook} analyse -s {sampledir} -r {referencefasta} -aligner bwa -coveragebin 100 -t {threads}'.format(
        nanook=P.option['program']['nanook'],
        sampledir=exptindir,
        referencefasta=P.option['program']['refpath'],
        threads=P.option['program']['server_threads']
    )
    if E[exptid]['libtype'] == '1D':
        cmd += ' -templateonly'
    mylogger.info('Running {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('nanook analyse passfail failed: stdout={0}, stderr={1}'.format(ro, re))
        return False
    origpassfailpdf = os.path.join(exptindir, 'latex_bwa_passfail', exptid+'.pdf')
    passfailpdf = os.path.join(exptindir, 'latex_bwa_passfail', exptid+'_passfail.pdf')
    try:
        os.rename(origpassfailpdf, passfailpdf)
    except:
        mylogger.warning('Failed to rename nanook pdf {0} -> {1}'.format(origfailonlypdf, failonlypdf))
    # Generate PDF report for pass reads only
    cmd = '{nanook} analyse -s {sampledir} -r {referencefasta} -aligner bwa -coveragebin 100 -t {threads}'.format(
        nanook=P.option['program']['nanook'],
        sampledir=exptindir,
        referencefasta=P.option['program']['refpath'],
        threads=P.option['program']['server_threads']
    )
    if E[exptid]['libtype'] == '1D':
        cmd += ' -templateonly'
    cmd += ' -passonly'
    mylogger.info('Running {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('nanook analyse passonly failed: stdout={0}, stderr={1}'.format(ro, re))
        return False
    origpassonlypdf = os.path.join(exptindir, 'latex_bwa_passonly', exptid+'.pdf')
    passonlypdf = os.path.join(exptindir, 'latex_bwa_passonly', exptid+'_passonly.pdf')
    try:
        os.rename(origpassonlypdf, passonlypdf)
    except:
        mylogger.warning('Failed to rename nanook pdf {0} -> {1}'.format(origfailonlypdf, failonlypdf))
    # Generate PDF report for fail reads only
    cmd = '{nanook} analyse -s {sampledir} -r {referencefasta} -aligner bwa -coveragebin 100 -t {threads}'.format(
        nanook=P.option['program']['nanook'],
        sampledir=exptindir,
        referencefasta=P.option['program']['refpath'],
        threads=P.option['program']['server_threads']
    )
    if E[exptid]['libtype'] == '1D':
        cmd += ' -templateonly'
    cmd += ' -failonly'
    mylogger.info('Running {0}'.format(cmd))
    rv, ro, re = P.sys_exec(cmd)
    if rv != 0:
        mylogger.error('nanook analyse failonly failed: stdout={0}, stderr={1}'.format(ro, re))
        return False
    origfailonlypdf = os.path.join(exptindir, 'latex_bwa_failonly', exptid+'.pdf')
    failonlypdf = os.path.join(exptindir, 'latex_bwa_failonly', exptid+'_failonly.pdf')
    try:
        os.rename(origfailonlypdf, failonlypdf)
    except:
        mylogger.warning('Failed to rename nanook pdf {0} -> {1}'.format(origfailonlypdf, failonlypdf))

    # If any of the output files are missing, print a warning message.
    outpathL = [passfailpdf, passonlypdf, failonlypdf]
    for outpath in outpathL:
        if not os.path.exists(outpath):
            mylogger.warning('nanook analyse failed to output {0}'.format(outpath))
    return True

def Process(args, P, mylogger, myhandler, processname):
    'Run nanook extract, align and analyse for each experiment in the experiment file.'
    exptidL, E = P.expt_read(args.experiments)
    for exptid in exptidL:
        exptindir = Nanook_SetupOutdirs(args, P, mylogger, myhandler, processname, exptid, E)
        if exptindir is not None:
            #Nanook_ExtractFastq(args, P, mylogger, myhandler, processname, exptid, E, exptindir)
            Nanook_ExtractFasta(args, P, mylogger, myhandler, processname, exptid, E, exptindir)
            Nanook_Align(args, P, mylogger, myhandler, processname, exptid, E, exptindir)
            Nanook_Analyse(args, P, mylogger, myhandler, processname, exptid, E, exptindir)
            #sys.exit(99)
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    'Execute this sub-tool.'
    mylogger.info('Started')
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    args.overwrite = P.str_2bool(args.overwrite)
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname)
    mylogger.info('Finished')
    return 0

