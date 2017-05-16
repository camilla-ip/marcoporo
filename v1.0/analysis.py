#!/usr/bin/env python

import logging
import numpy as np
import os
import random
import sys

import marcoporoversion

_processname = 'analyse'
_jobprefix = 'mpn'

def Prerequisites(args, P, mylogger, myhandler, processname):
    'Exit program if some prerequisites are not met.'
    if not os.path.exists(args.outdir):
        try:
            os.makedirs(args.outdir)
        except:
            mylogger.error('Failed to create outdir ({0})'.format(args.outdir))
    return 0

def ReadClass(topdir, exptidL):
    '''
    Find the readclass of all readids from the readclass field of the
    _readstats.txt files, store in pass[(exptid, readid)] = 1 and
    fail[(exptid, readid)] = 1.
    '''
    passD = {}
    failD = {}
    rcnkD = {}
    for exptid in exptidL:
        inpath = os.path.join(topdir, '03-extract', exptid+'_readstats.txt')
        with open(inpath, 'r') as in_fp:
            linecnt = 0
            for line in in_fp:
                linecnt += 1
                if linecnt == 1:
                    continue
                L = line.rstrip('\n').split('\t')
                readid = L[3]
                readclass = L[8]
                if readclass == 'pass':
                    passD[(exptid, readid)] = True
                elif readclass == 'fail':
                    failD[(exptid, readid)] = True
                else:
                    rcnkD[(exptid, readid)] = True
    return passD, failD, rcnkD

def BaseYield_GB(P, topdir, exptid, libtype, readclass):
    'Compute the base yield from the seqlen field of the read1dstats file (rows with readtype 1T or 1C) or read2dstats seqlen field.'
  # Store readclass of each read in a dictionary
    #passD, failD, rcnkD = ReadClass(topdir, exptidL)
  # Find the total
    if libtype == '1D':
        inpath = os.path.join(topdir, '08-analysis', exptid+'_read1dstats.txt')
        read1d = np.genfromtxt(inpath, skiprows=1, delimiter='\t', dtype=P.ontread1dstatsbH[:-2], missing_values="NA")
        n = sum(read1d[np.array(read1d[:]['readtype'] == '1T') & np.array(read1d[:]['readclass'] == readclass)]['seqlen']/1000000000.0)
    elif libtype == '2D':
        inpath = os.path.join(topdir, '08-analysis', exptid+'_read2dstats.txt')
        read2d = np.genfromtxt(inpath, skiprows=1, delimiter='\t', dtype=P.ontread2dstatsbH[:-2], missing_values="NA")        
        n = sum(read2d[read2d[:]['readclass'] == readclass]['seqlen']/1000000000.0)
    else:
        n = None
    return n

def Table_YieldAndQuality(args, P, mylogger, myhandler, processname, exptidL, E):
    'Produce X-Y matrix of experiments-statistics on experiment yield and quality.'
    readclassL = ['pass', 'fail']  
  # Collate matrix of data
    T = []
  # L1 : experiment - phase-lab-replicate
    T.append(['Statistic / Experiment'] + exptidL)
  # L2 : libtype - 1D, 2D
    T.append(['library & basecall type'] + [E[e]['libtype'] for e in exptidL])
  # L3 : readclass - pass, fail
    T.append(['read class'] + readclassL*len(exptidL))
  # L4: Read count (K) - split by readclass
    row = ['Read count (K)']
    for exptid in exptidL:
        for readclass in readclassL:
            dir = os.path.join(args.topdir, '01-fast5', exptid, 'reads', 'downloads', readclass)
            n = str(round(len([x for x in os.listdir(dir) if x.endswith('.fast5')])/1000.0, 1))
            row.append(n)
    T.append(row)
  # L5 : Read count (K) - total of readclasses
    row = ['total']
    for exptid in exptidL:
        n = 0
        for readclass in readclassL:
            dir = os.path.join(args.topdir, '01-fast5', exptid, 'reads', 'downloads', readclass)
            n += len([x for x in os.listdir(dir) if x.endswith('.fast5')])/1000.0
        row.append(str(round(n, 1)))
        row.append('')
    T.append(row)
  # L6 : Base yield (G) - split by readclass
    row = ['Base yield (GB)']
    for exptid in exptidL:
        for readclass in readclassL:
            n = str(round(BaseYield_GB(P, args.topdir, exptid, E[exptid]['libtype'], readclass), 3))
            row.append(n)
    T.append(row)
  # L7 : Base yield (G) - total of readclasses

  # Save matrix to file
    outpath = os.path.join(args.outdir, 'table_yieldandquality.txt')
    with open(outpath, 'w') as out_fp:
        for row in T:
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in row])))

def AddReadClassColumnToStatsTables(args, P, mylogger, myhandler, processname, exptidL, E):
    'Add readclass column from _readstats to _readevents, _read1d and _read2d tables in the 08-analysis dir.'

    passD, failD, rcnkD = ReadClass(args.topdir, exptidL)
    for exptid in exptidL:
        mylogger.info('Adding readclass column for {0}'.format(exptid))
      # _readeventstats.txt
        inpath = os.path.join(args.topdir, '03-extract', exptid+'_readeventstats.txt')
        outpath = os.path.join(args.outdir, exptid+'_readeventstats.txt')
        data = np.genfromtxt(inpath, skiprows=1, delimiter='\t', dtype=P.ontreadeventstatsH, missing_values="NA")
        n = 6
        tablename = 'ontreadeventstats'
        with open(outpath, 'w') as out_fp:
            H = P.fast5_headernames(tablename)
            L = H[:n] + ['readclass'] + H[n:]
            out_fp.write('{0}\n'.format('\t'.join(L)))
            for row in data:
                readid = row['readid']
                readclass = 'pass' if passD.has_key((exptid, readid)) else ('fail' if failD.has_key((exptid, readid)) else 'NK')
                L = [row[x] for x in H[:n]] + [readclass] + [row[x] for x in H[n:]]
                out_fp.write('{0}\n'.format('\t'.join([str(x) for x in L])))
      # _read1dstats.txt
        inpath = os.path.join(args.topdir, '03-extract', exptid+'_read1dstats.txt')
        outpath = os.path.join(args.outdir, exptid+'_read1dstats.txt')
        data = np.genfromtxt(inpath, skiprows=1, delimiter='\t', dtype=P.ontread1dstatsH, missing_values="NA")
        n = 7
        tablename = 'ontread1dstats'
        H = P.fast5_headernames(tablename)
        with open(outpath, 'w') as out_fp:
            L = H[:n] + ['readclass'] + H[n:]
            out_fp.write('{0}\n'.format('\t'.join(L)))
            for row in data:
                readid = row['readid']
                readclass = 'pass' if passD.has_key((exptid, readid)) else ('fail' if failD.has_key((exptid, readid)) else 'NK')
                L = [row[x] for x in H[:n]] + [readclass] + [row[x] for x in H[n:]]
                out_fp.write('{0}\n'.format('\t'.join([str(x) for x in L])))
      # _read2dstats.txt
        inpath = os.path.join(args.topdir, '03-extract', exptid+'_read2dstats.txt')
        outpath = os.path.join(args.outdir, exptid+'_read2dstats.txt')
        data = np.genfromtxt(inpath, skiprows=1, delimiter='\t', dtype=P.ontread2dstatsH, missing_values="NA")
        n = 6
        tablename = 'ontread2dstats'
        H = P.fast5_headernames(tablename)
        with open(outpath, 'w') as out_fp:
            H = P.fast5_headernames(tablename)
            L = H[:n] + ['readclass'] + H[n:]
            out_fp.write('{0}\n'.format('\t'.join(L)))
            for row in data:
                readid = row['readid']
                readclass = 'pass' if passD.has_key((exptid, readid)) else ('fail' if failD.has_key((exptid, readid)) else 'NK')
                L = [row[x] for x in H[:n]] + [readclass] + [row[x] for x in H[n:]]
                out_fp.write('{0}\n'.format('\t'.join([str(x) for x in L])))
    return 0

def Analyse(args, P, mylogger, myhandler, processname, exptidL, E):
    Table_YieldAndQuality(args, P, mylogger, myhandler, processname, exptidL, E)
    return 0

def Process(args, P, mylogger, myhandler, processname):
    'Run aggregateone for all experiments in experiments.txt file.'
    exptidL, E = P.expt_read(args.experiments)
    #AddReadClassColumnToStatsTables(args, P, mylogger, myhandler, processname, exptidL, E)
    Analyse(args, P, mylogger, myhandler, processname, exptidL, E)
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    'Execute this sub-tool.'
    mylogger.info('Started')
    args.overwrite = P.str_2bool(args.overwrite)
    args.execjobs = P.str_2bool(args.execjobs)
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname)
    mylogger.info('Finished')
    return 0
