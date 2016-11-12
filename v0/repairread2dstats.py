#!/usr/bin/env python

import os
import sys

if len(sys.argv) != 4:
    print('Usage: repairread2dstats.py read1dfile read2dfile newread2dfile')
    print('       Use the information in the read1d file to add 5 extra time columns to the read2dfile.')
    sys.exit(1)
in1dpath, in2dpath, out2dpath = sys.argv[1:]

def Process(in1dpath, in2dpath, out2dpath):

  # Read in the 1dpath data
  # Input has columns: exptid batchid runid readid bc1dinstanceN readtype returnstatus numevents numskip numstays numcalled strandstarttimesec strandendtimesec stranddurationsec strandstarttimeiso strandendtimeiso meanqscore strandscore seqlen bqlen bqmean bqmedian gcpct basespersecond
    data1d = {}
    with open(in1dpath, 'r') as in_fp:
        linecnt = 0
        for line in in_fp:
            L = line.rstrip('\n').split('\t')
            readid = L[3]
            readtype = L[5]
            strandstarttimesec = L[11]
            strandendtimesec = L[12]
            stranddurationsec = L[13]
            strandstarttimeiso = L[14]
            strandendtimeiso = L[15]
            key = (readid, readtype)
            data1d[key] = {
                'strandstarttimesec': strandstarttimesec,
                'strandendtimesec': strandendtimesec,
                'stranddurationsec': stranddurationsec,
                'strandstarttimeiso': strandstarttimeiso,
                'strandendtimeiso': strandendtimeiso
            }

  # Read in the 2dpath data
  # Input has columns: exptid batchid runid readid bc2instanceN returnstatus meanqscore seqlen bqlen bqmean bqmedian gcpct basespersecond
  # Output has columns: exptid batchid runid readid bc2instanceN returnstatus STARTTIMESEC ENDTIMESEC DURATIONSEC STARTTIMEISO ENDTIMEISO meanqscore seqlen bqlen bqmean bqmedian gcpct basespersecond
    data2d = {}
    data2dreadidL = []
    with open(in2dpath, 'r') as in_fp:
        linecnt = 0
        for line in in_fp:
            linecnt += 1
            if linecnt == 1:
                continue
            L = line.rstrip('\n').split('\t')
            readid = L[3]
            data2dreadidL.append(readid)
            data2d[readid] = L

  # Output a read2d file with 2 extra columns (endtimehrs and durationhrs)
    H = ['exptid', 'batchid', 'runid', 'readid', 'bc2instanceN',
        'returnstatus', 'starttimesec', 'endtimesec', 'durationsec', 'starttimeiso',
        'endtimeiso', 'meanqscore', 'seqlen', 'bqlen', 'bqmean',
        'bqmedian', 'gcpct', 'basespersecond']
    with open(out2dpath, 'w') as out_fp:
        out_fp.write('{0}\n'.format('\t'.join(H)))
        for readid in data2dreadidL:
            has1T = data1d.has_key((readid, '1T')) and data1d[(readid, '1T')]['strandstarttimesec'] != 'NA'
            has1C = data1d.has_key((readid, '1C')) and data1d[(readid, '1C')]['strandstarttimesec'] != 'NA'
            if has1T and has1C:
                starttimesec = data1d[(readid, '1T')]['strandstarttimesec']
                endtimesec = data1d[(readid, '1C')]['strandendtimesec']
                durationsec = str(round(float(endtimesec) - float(starttimesec), 3))
                starttimeiso = data1d[(readid, '1T')]['strandstarttimeiso']
                endtimeiso = data1d[(readid, '1C')]['strandendtimeiso']
            elif has1T and not has1C:
                starttimesec = data1d[(readid, '1T')]['strandstarttimesec']
                endtimesec = data1d[(readid, '1T')]['strandendtimesec']
                durationsec = data1d[(readid, '1T')]['stranddurationsec']
                starttimeiso = data1d[(readid, '1T')]['strandstarttimeiso']
                endtimeiso = data1d[(readid, '1T')]['strandendtimeiso']
            elif not has1T and has1C:
                starttimesec = data1d[(readid, '1C')]['strandstarttimesec']
                endtimesec = data1d[(readid, '1C')]['strandendtimesec']
                durationsec = data1d[(readid, '1C')]['stranddurationsec']
                starttimeiso = data1d[(readid, '1C')]['strandstarttimeiso']
                endtimeiso = data1d[(readid, '1C')]['strandendtimeiso']
            else:
                starttimesec = 'NA'
                endtimesec = 'NA'
                durationsec = 'NA'
                starttimeiso = 'NA'
                endtimeiso = 'NA'
            newrow = L[:6] + [starttimesec, endtimesec, durationsec, starttimeiso, endtimeiso] + L[6:]
            out_fp.write('{0}\n'.format('\t'.join(newrow)))

Process(in1dpath, in2dpath, out2dpath)

