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
    binmeans = [valA[indices == i].mean() for i in range(0, len(binA))]
    return binmeans

def Process(args, P, mylogger, myhandler, processname, exptid):
    'Read file of constant field names produced by marcoporo exptconstants.'
    # Output file
    outpath = os.path.join(args.outdir, args.exptid+'_aggregate.txt')
    # Time bins, every 0.25 hours between 0 and 48 hours
    binA = np.arange(0,2*24+0.25,0.25)
    read1dpath = os.path.join(args.indir, exptid+'_read1dstats.txt')
    read1d = np.genfromtxt(read1dpath, skiprows=1, delimiter='\t', dtype=P.ontread1dstatsH, missing_values="NA")
    endtimehrs = read1d[read1d[:]['readtype'] == '1T']['strandendtimesec']/60.0/60.0
    # read1d_1T metrics aggregated by time bins
    read1d_1T_stranddurationsec = binmean(read1d[read1d[:]['readtype'] == '1T']['stranddurationsec'], endtimehrs, binA)
    read1d_1T_meanqscore = binmean(read1d[read1d[:]['readtype'] == '1T']['meanqscore'], endtimehrs, binA)
    read1d_1T_meanseqlen = binmean(read1d[read1d[:]['readtype'] == '1T']['seqlen'], endtimehrs, binA)
    read1d_1T_meanbq = binmean(read1d[read1d[:]['readtype'] == '1T']['bqmean'], endtimehrs, binA)
    read1d_1T_meangcpct = binmean(read1d[read1d[:]['readtype'] == '1T']['gcpct'], endtimehrs, binA)
    read1d_1T_meanbps = binmean(read1d[read1d[:]['readtype'] == '1T']['basespersecond'], endtimehrs, binA)
    # read1d_1C metrics aggregated by time bins
    read1d_1C_stranddurationsec = binmean(read1d[read1d[:]['readtype'] == '1C']['stranddurationsec'], endtimehrs, binA)
    read1d_1C_meanqscore = binmean(read1d[read1d[:]['readtype'] == '1C']['meanqscore'], endtimehrs, binA)
    read1d_1C_meanseqlen = binmean(read1d[read1d[:]['readtype'] == '1C']['seqlen'], endtimehrs, binA)
    read1d_1C_meanbq = binmean(read1d[read1d[:]['readtype'] == '1C']['bqmean'], endtimehrs, binA)
    read1d_1C_meangcpct = binmean(read1d[read1d[:]['readtype'] == '1C']['gcpct'], endtimehrs, binA)
    read1d_1C_meanbps = binmean(read1d[read1d[:]['readtype'] == '1C']['basespersecond'], endtimehrs, binA)
    # Create final 2D matrix (X=variables, Y=timebuckets)
    H = ['exptid', 'timehr',
         '1Tdurationsec1T', '1Tqscore', '1Tseqlen', '1Tbq', '1Tgcpct', '1Tbasesps',
         '1Cdurationsec1C', '1Cqscore', '1Cseqlen', '1Cbq', '1Cgcpct', '1Cbasesps']
    #A = np.column_stack((
    #    read1d_exptid, binA,
    #    read1d_1T_stranddurationsec, read1d_1T_meanqscore, read1d_1T_meanseqlen, read1d_1T_meanbq, read1d_1T_meangcpct, read1d_1T_meanbps,
    #    read1d_1C_stranddurationsec, read1d_1C_meanqscore, read1d_1C_meanseqlen, read1d_1C_meanbq, read1d_1C_meangcpct, read1d_1C_meanbps))
    exptidA = np.array([exptid]*len(binA))
    A = np.column_stack((
        exptidA,
        binA,
        read1d_1T_stranddurationsec, read1d_1T_meanqscore, read1d_1T_meanseqlen, read1d_1T_meanbq, read1d_1T_meangcpct, read1d_1T_meanbps,
        read1d_1C_stranddurationsec, read1d_1C_meanqscore, read1d_1C_meanseqlen, read1d_1C_meanbq, read1d_1C_meangcpct, read1d_1C_meanbps))
    np.savetxt(outpath, A, fmt='%s', delimiter='\t', newline='\n', comments='', header='\t'.join(H))
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    'Execute this sub-tool.'
    mylogger.info('Started')
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname, args.exptid)
    mylogger.info('Finished')
    return 0
