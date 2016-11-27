#!/usr/bin/env python

import logging
import numpy as np
import os
import sys
import time

import marcoporoversion

_processname = 'aggregateone'

def Prerequisites(args, P, mylogger, myhandler, processname):
    'Exit program if some prerequisites are not met.'
  # Extractone output stats files must already exist, even if they are empty
    for suffix in P.extractstatfilesuffix:
        inpath = os.path.join(args.indir, args.exptid+suffix)
        if not os.path.exists(inpath):
            mylogger.error('Input file does not exist {0}'.format(inpath))
            sys.exit(P.err_code('ErrorFileMissing'))
  # If all empty, then probably something wrong
    allhavedata = True
    numfiles = len(P.extractstatfilesuffix)
    numempty = 0
    for suffix in P.extractstatfilesuffix:
        inpath = os.path.join(args.indir, args.exptid+suffix)
        if not os.path.getsize(inpath):
            numempty += 1
    if numempty == numfiles:
        mylogger.error('All input files in indir are empty, nothing to do ({0})'.format(args.indir))
        sys.exit(P.err_code('ErrorEmptyFile'))
    return 0

def binmean(valA, timeA, binA):
    indices = np.digitize(timeA, binA)
    binmeans = [np.nanmean([x for x in valA[indices == i] if x != -1]) for i in range(0, len(binA))]
    return binmeans

def Aggregate_readevent(args, P, mylogger, myhandler, processname, exptid):
  # Output file
    outpath = os.path.join(args.outdir, args.exptid+'_aggregate_readevent.txt')
  # Time bins, every timebucket hours between 0 and maxrunlen hours
    binA = np.arange(0, args.maxrunlen+args.timebucket, args.timebucket)
    readeventpath = os.path.join(args.indir, exptid+'_readeventstats.txt')
    readevent = np.genfromtxt(readeventpath, skiprows=1, delimiter='\t', dtype=P.ontreadeventstatsH, missing_values="NA")
    endtimehrs = readevent[:]['eventendtimesec']/60.0/60.0
  # readevent metrics aggregated by time bins
    readevent_meandurationsec = binmean(readevent[:]['eventdurationsec'], endtimehrs, binA)
    readevent_meaneventcount = binmean(readevent[:]['eventcount'], endtimehrs, binA)
    readevent_meaneventspersec = binmean(readevent[:]['eventspersec'], endtimehrs, binA)
  # Create final 2D matrix (rows=timebuckets, columns=variables) and save to file
    H = ['exptid', 'timehr',
         'durationsec', 'eventcount', 'eventspersec']
    exptidA = np.array([exptid]*len(binA))
    A = np.column_stack((
        exptidA,
        binA,
        readevent_meandurationsec, readevent_meaneventcount, readevent_meaneventspersec))
    np.savetxt(outpath, A, fmt='%s', delimiter='\t', newline='\n', comments='', header='\t'.join(H))
    return 0

def Aggregate_read1d(args, P, mylogger, myhandler, processname, exptid):
  # Output file
    outpath = os.path.join(args.outdir, args.exptid+'_aggregate_read1d.txt')
  # Time bins, every timebucket hours between 0 and maxrunlen hours
    binA = np.arange(0, args.maxrunlen+args.timebucket, args.timebucket)
    read1dpath = os.path.join(args.indir, exptid+'_read1dstats.txt')
    read1d = np.genfromtxt(read1dpath, skiprows=1, delimiter='\t', dtype=P.ontread1dstatsH, missing_values="NA")
    endtimehrs1T = read1d[read1d[:]['readtype'] == '1T']['strandendtimesec']/60.0/60.0
    endtimehrs1C = read1d[read1d[:]['readtype'] == '1C']['strandendtimesec']/60.0/60.0
  # read1d_1T metrics aggregated by time bins
    read1d_1T_stranddurationsec = binmean(read1d[read1d[:]['readtype'] == '1T']['stranddurationsec'], endtimehrs1T, binA)
    read1d_1T_meanqscore = binmean(read1d[read1d[:]['readtype'] == '1T']['meanqscore'], endtimehrs1T, binA)
    read1d_1T_meanseqlen = binmean(read1d[read1d[:]['readtype'] == '1T']['seqlen'], endtimehrs1T, binA)
    read1d_1T_meanbq = binmean(read1d[read1d[:]['readtype'] == '1T']['bqmean'], endtimehrs1T, binA)
    read1d_1T_meangcpct = binmean(read1d[read1d[:]['readtype'] == '1T']['gcpct'], endtimehrs1T, binA)
    read1d_1T_meanbps = binmean(read1d[read1d[:]['readtype'] == '1T']['basespersecond'], endtimehrs1T, binA)
  # read1d_1C metrics aggregated by time bins
    read1d_1C_stranddurationsec = binmean(read1d[read1d[:]['readtype'] == '1C']['stranddurationsec'], endtimehrs1C, binA)
    read1d_1C_meanqscore = binmean(read1d[read1d[:]['readtype'] == '1C']['meanqscore'], endtimehrs1C, binA)
    read1d_1C_meanseqlen = binmean(read1d[read1d[:]['readtype'] == '1C']['seqlen'], endtimehrs1C, binA)
    read1d_1C_meanbq = binmean(read1d[read1d[:]['readtype'] == '1C']['bqmean'], endtimehrs1C, binA)
    read1d_1C_meangcpct = binmean(read1d[read1d[:]['readtype'] == '1C']['gcpct'], endtimehrs1C, binA)
    read1d_1C_meanbps = binmean(read1d[read1d[:]['readtype'] == '1C']['basespersecond'], endtimehrs1C, binA)
  # Create final 2D matrix (rows=timebuckets, columns=variables) and save to file
    H = ['exptid', 'timehr',
         'durationsec1T', 'qscore1T', 'seqlen1T', 'bq1T', 'gcpct1T', 'basesps1T',
         'durationsec1C', 'qscore1C', 'seqlen1C', 'bq1C', 'gcpct1C', 'basesps1C']
    exptidA = np.array([exptid]*len(binA))
    A = np.column_stack((
        exptidA,
        binA,
        read1d_1T_stranddurationsec, read1d_1T_meanqscore, read1d_1T_meanseqlen, read1d_1T_meanbq, read1d_1T_meangcpct, read1d_1T_meanbps,
        read1d_1C_stranddurationsec, read1d_1C_meanqscore, read1d_1C_meanseqlen, read1d_1C_meanbq, read1d_1C_meangcpct, read1d_1C_meanbps))
    np.savetxt(outpath, A, fmt='%s', delimiter='\t', newline='\n', comments='', header='\t'.join(H))
    return 0

def Aggregate_read2d(args, P, mylogger, myhandler, processname, exptid):
  # Output file
    outpath = os.path.join(args.outdir, args.exptid+'_aggregate_read2d.txt')
  # Time bins, every timebucket hours between 0 and maxrunlen hours
    binA = np.arange(0, args.maxrunlen+args.timebucket, args.timebucket)
    read2dpath = os.path.join(args.indir, exptid+'_read2dstats.txt')
    read2d = np.genfromtxt(read2dpath, skiprows=1, delimiter='\t', dtype=P.ontread2dstatsH, missing_values="NA")
    endtimehrs = read2d[:]['endtimesec']/60.0/60.0
  # read2d metrics aggregated by time bins
    read2d_meandurationsec = binmean(read2d[:]['durationsec'], endtimehrs, binA)
    read2d_meanqscore = binmean(read2d[:]['meanqscore'], endtimehrs, binA)
    read2d_meanseqlen = binmean(read2d[:]['seqlen'], endtimehrs, binA)
    read2d_meanbq = binmean(read2d[:]['bqmean'], endtimehrs, binA)
    read2d_meangcpct = binmean(read2d[:]['gcpct'], endtimehrs, binA)
    read2d_meanbps = binmean(read2d[:]['basespersecond'], endtimehrs, binA)
  # Create final 2D matrix (rows=timebuckets, columns=variables) and save to file
    H = ['exptid', 'timehr',
         'durationsec2D', 'qscore2D', 'seqlen2D', 'bq2D', 'gcpct2D', 'basesps2D']
    exptidA = np.array([exptid]*len(binA))
    A = np.column_stack((
        exptidA,
        binA,
        read2d_meandurationsec, read2d_meanqscore, read2d_meanseqlen, read2d_meanbq, read2d_meangcpct,
        read2d_meanbps))
    np.savetxt(outpath, A, fmt='%s', delimiter='\t', newline='\n', comments='', header='\t'.join(H))
    return 0

def Process(args, P, mylogger, myhandler, processname, exptid):
    'Aggregate per-read ont*stats.txt metrics into windows of X hours.'
    Aggregate_readevent(args, P, mylogger, myhandler, processname, exptid)
    Aggregate_read1d(args, P, mylogger, myhandler, processname, exptid)
    Aggregate_read2d(args, P, mylogger, myhandler, processname, exptid)
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    'Execute this sub-tool.'
    mylogger.info('Started')
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname, args.exptid)
    mylogger.info('Finished')
    return 0
