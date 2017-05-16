#!/usr/bin/env python

import logging
import numpy as np
import os
import sys
import time

import marcoporoversion

_processname = 'aggregate'

def Prerequisites(args, P, mylogger, myhandler, processname, exptidL, E):
    'Exit program if some prerequisites are not met.'
  # marcoporo extract output stats files must already exist, even if they are empty
    for suffix in P.extractstatfilesuffix:
        for exptid in exptidL:
            inpath = os.path.join(args.extractdir, exptid+suffix)
            if not os.path.exists(inpath):
                mylogger.error('Input file does not exist {0}'.format(inpath))
                sys.exit(P.err_code('ErrorFileMissing'))
  # If all empty, then probably something wrong
    allhavedata = True
    numfiles = len(P.extractstatfilesuffix)
    numempty = 0
    for suffix in P.extractstatfilesuffix:
        for exptid in exptidL:
            inpath = os.path.join(args.extractdir, exptid+suffix)
            if not os.path.getsize(inpath):
                numempty += 1
    if numempty == numfiles:
        mylogger.error('All input files in extractdir are empty, nothing to do ({0})'.format(args.extractdir))
        sys.exit(P.err_code('ErrorEmptyFile'))
    return 0

def binmean(valA, timeA, binA):
    if not len(valA) and not len(timeA):
        return [0.0] * len(binA), [0.0] * len(binA)
    indexL = np.digitize(timeA, binA)
    binmeans = [np.nanmean([x for x in valA[indexL == i] if x != -1]) for i in range(0, len(binA))]
    bincount = [len([x for x in valA[indexL == i] if x != -1]) for i in range(0, len(binA))]
    return binmeans, bincount

def Aggregate_readevent(args, P, mylogger, myhandler, processname, exptid):
  # Output file
    outpath = os.path.join(args.outdir, exptid+'_aggregate_readevent.txt')
  # Time bins, every timebucket hours between 0 and maxrunlen hours
    timeh = np.arange(0, args.maxrunlen+args.timebucket, args.timebucket)
    readeventpath = os.path.join(args.extractdir, exptid+'_readeventstats.txt')
    readevent = np.genfromtxt(readeventpath, skiprows=1, delimiter='\t', dtype=P.ontreadeventstatsH, missing_values="NA")
    endtimehrs = readevent[:]['eventendtimesec']/60.0/60.0
  # readevent metrics aggregated by time bins
    readevent_meandurationsec, readevent_meandurationsec_bincount = binmean(readevent[:]['eventdurationsec'], endtimehrs, timeh)
    readevent_meaneventcount, readevent_meaneventcount_bincount = binmean(readevent[:]['eventcount'], endtimehrs, timeh)
    readevent_meaneventspersec, readevent_meaneventspersec_bincount = binmean(readevent[:]['eventspersec'], endtimehrs, timeh)
  # Create final 2D matrix (rows=timebuckets, columns=variables) and save to file
    H = ['exptid', 'timehr',
         'durationsec', 'eventcount', 'eventspersec']
    exptidA = np.array([exptid]*len(timeh))
    A = np.column_stack((
        exptidA,
        timeh,
        readevent_meandurationsec, readevent_meaneventcount, readevent_meaneventspersec))
    np.savetxt(outpath, A, fmt='%s', delimiter='\t', newline='\n', comments='', header='\t'.join(H))
    return 0

def Parse_poremapstats_read1dstats(passpath, failpath):
    'Return var with mapstats[\'header\'] = L and mapstats[\'data\'][(exptid, readid, readtype)] = L.'
    stats = {}
    stats['header'] = []
    stats['data'] = {}
    if os.path.exists(passpath):
        with open(passpath, 'r') as in_fp:
            linecnt = 0
            for line in in_fp:
                linecnt += 1
                L = line.rstrip('\n').split('\t')
                if linecnt == 1:
                    if not len(stats['header']):
                        stats['header'] = L
                else:
                    exptid = L[0].split('_')[0]
                    readid = L[1].split('_')[0]
                    readtype = '1T' if L[2] == 'temp' else ('1C' if L[2] == 'comp' else '2D')
                    readclass = L[3]
                    key = (exptid, readid, readtype)
                    stats['data'][key] = L
    if os.path.exists(failpath):
        with open(failpath, 'r') as in_fp:
            linecnt = 0
            for line in in_fp:
                linecnt += 1
                L = line.rstrip('\n').split('\t')
                if linecnt == 1:
                    if not len(stats['header']):
                        stats['header'] = L
                else:
                    exptid = L[0].split('_')[0]
                    readid = L[1].split('_')[0]
                    readtype = '1T' if L[2] == 'temp' else ('1C' if L[2] == 'comp' else '2D')
                    readclass = L[3]
                    key = (exptid, readid, readtype)
                    stats['data'][key] = L
    return stats

def Parse_poremapstats_read2dstats(passpath, failpath):
    'Return var with mapstats[\'header\'] = L and mapstats[\'data\'][(exptid, readid, readtype)] = L.'
    stats = {}
    stats['header'] = []
    stats['data'] = {}
    if os.path.exists(passpath):
        with open(passpath, 'r') as in_fp:
            linecnt = 0
            for line in in_fp:
                linecnt += 1
                L = line.rstrip('\n').split('\t')
                if linecnt == 1:
                    if not len(stats['header']):
                        stats['header'] = L
                else:
                    exptid = L[0].split('_')[0]
                    readid = L[1].split('_')[0]
                    readtype = '2D'
                    readclass = L[3]
                    key = (exptid, readid, readtype)
                    stats['data'][key] = L
    if os.path.exists(failpath):
        with open(failpath, 'r') as in_fp:
            linecnt = 0
            for line in in_fp:
                linecnt += 1
                L = line.rstrip('\n').split('\t')
                if linecnt == 1:
                    if not len(stats['header']):
                        stats['header'] = L
                else:
                    exptid = L[0].split('_')[0]
                    readid = L[1].split('_')[0]
                    readtype = '2D'
                    readclass = L[3]
                    key = (exptid, readid, readtype)
                    stats['data'][key] = L
    return stats


def Aggregate_merge1d_headerL(statname):
    H = ['exptid', 'timeh']
    for readtype in ['1T', '1C']:
        H.append('durationsec_'+readtype)
        for readclass in ['passfail', 'passonly', 'failonly']:
            for maptype in ['mapa', 'mapy', 'mapn']:
                s = statname+'_'+readtype+'_'+readclass+'_'+maptype+'_mean'
                H.append(s)
            for maptype in ['mapa', 'mapy', 'mapn']:
                s = statname+'_'+readtype+'_'+readclass+'_'+maptype+'_count'
                H.append(s)
    return H

def Aggregate_merge2d_headerL(statname):
    H = ['exptid', 'timeh']
    for readtype in ['2D']:
        H.append('durationsec_'+readtype)
        for readclass in ['passfail', 'passonly', 'failonly']:
            for maptype in ['mapa', 'mapy', 'mapn']:
                s = statname+'_'+readtype+'_'+readclass+'_'+maptype+'_mean'
                H.append(s)
            for maptype in ['mapa', 'mapy', 'mapn']:
                s = statname+'_'+readtype+'_'+readclass+'_'+maptype+'_count'
                H.append(s)
    return H

def Create_merg1d(args, P, exptid, merg1dpath, read1dpath, read1d, endtimehrs1T, endtimehrs1C):
  # Read in the poremapstats 'readstats.txt' output files for 1D reads
    stats1tpasspath = os.path.join(args.bwamemdir,
        #'{exptid}_{readtype}_{readclass}'.format(exptid=exptid, readtype='1T', readclass='pass'),
        '{exptid}_{readtype}_{readclass}_readstats.txt'.format(exptid=exptid, readtype='1T', readclass='pass'))
    stats1tfailpath = os.path.join(args.bwamemdir,
        #'{exptid}_{readtype}_{readclass}'.format(exptid=exptid, readtype='1T', readclass='fail'),
        '{exptid}_{readtype}_{readclass}_readstats.txt'.format(exptid=exptid, readtype='1T', readclass='fail'))
    stats1cpasspath = os.path.join(args.bwamemdir,
        #'{exptid}_{readtype}_{readclass}'.format(exptid=exptid, readtype='1C', readclass='pass'),
        '{exptid}_{readtype}_{readclass}_readstats.txt'.format(exptid=exptid, readtype='1C', readclass='pass'))
    stats1cfailpath = os.path.join(args.bwamemdir,
        #'{exptid}_{readtype}_{readclass}'.format(exptid=exptid, readtype='1C', readclass='fail'),
        '{exptid}_{readtype}_{readclass}_readstats.txt'.format(exptid=exptid, readtype='1C', readclass='fail'))
    stats1t = Parse_poremapstats_read1dstats(stats1tpasspath, stats1tfailpath)
    stats1c = Parse_poremapstats_read1dstats(stats1cpasspath, stats1cfailpath)
  # Save merged tables to a file called exptid_read1dmerged.txt with one line for 1T and one line for 1C stats
    merg1dpath = os.path.join(args.outdir, exptid+'_merged1dstats.txt')
    with open(merg1dpath, 'w') as out_fp:
      # Print header for merged table
        row = P.fast5_headernames('ontread1dstats') + stats1t['header'][5:]
        out_fp.write('{0}\n'.format('\t'.join([str(x) for x in row])))
        out_fp.flush()
      # Print extended rows
        for rowstt in read1d:
            exptid = rowstt[0]
            readid = rowstt[3]
            readtype = rowstt[5]
            key = (exptid, readid, readtype)
            if stats1t['data'].has_key(key):
                rowend = stats1t['data'][key][5:]
            elif stats1c['data'].has_key(key):
                rowend = stats1c['data'][key][5:]
            else:
                rowend = [''] * len(stats1t['header'][5:])
            row = list(rowstt)[:-2] + rowend
            row = ['NA' if (str(x)=='nan' or str(x) == '-1' or not len(str(x))) else str(x) for x in row]
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in row])))
            out_fp.flush()
    return 0

def Create_merg2d(args, P, exptid, merg2dpath, read2dpath, read2d, endtimehrs2D):
  # Read in the poremapstats 'readstats.txt' output files for 2D reads
    stats2dpasspath = os.path.join(args.bwamemdir,
        '{exptid}_{readtype}_{readclass}'.format(exptid=exptid, readtype='2D', readclass='pass'),
        '{exptid}_{readtype}_{readclass}_readstats.txt'.format(exptid=exptid, readtype='2D', readclass='pass'))
    stats2dfailpath = os.path.join(args.bwamemdir,
        '{exptid}_{readtype}_{readclass}'.format(exptid=exptid, readtype='2D', readclass='fail'),
        '{exptid}_{readtype}_{readclass}_readstats.txt'.format(exptid=exptid, readtype='2D', readclass='fail'))
    if not os.path.exists(stats2dpasspath) or not os.path.exists(stats2dfailpath):
        return 1
    stats2d = Parse_poremapstats_read2dstats(stats2dpasspath, stats2dfailpath)
  # Save merged tables to a file called exptid_read2dmerged.txt with one line for 2D stats
    merg2dpath = os.path.join(args.outdir, exptid+'_merged2dstats.txt')
    with open(merg2dpath, 'w') as out_fp:
      # Print header for merged table
        row = P.fast5_headernames('ontread2dstats') + stats2d['header'][5:]
        out_fp.write('{0}\n'.format('\t'.join([str(x) for x in row])))
        out_fp.flush()
      # Print extended rows
        for rowstt in read2d:
            exptid = rowstt[0]
            readid = rowstt[3]
            readtype = '2D'
            key = (exptid, readid, readtype)
            if stats2d['data'].has_key(key):
                rowend = stats2d['data'][key][5:]
            else:
                rowend = [''] * len(stats2d['header'][5:])
            row = list(rowstt)[:-2] + rowend
            row = ['NA' if (str(x)=='nan' or str(x) == '-1' or not len(str(x))) else str(x) for x in row]
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in row])))
            out_fp.flush()
    return 0

def Print_Aggregate_Read1D_Statistics_File(merg1d, outdir, exptid, timeh, stranddurationsec1T, stranddurationsec1C, mask, statname):

    colmeans = { '1T':{}, '1C':{} }
    colcount = { '1T':{}, '1C':{} }
    for key in mask['1T'].keys():
        colmeans['1T'][key], colcount['1T'][key] = binmean(merg1d[mask['1T'][key]][statname], merg1d[mask['1T'][key]]['strandendtimesec']/60.0/60.0, timeh)
    for key in mask['1C'].keys():
        colmeans['1C'][key], colcount['1C'][key] = binmean(merg1d[mask['1C'][key]][statname], merg1d[mask['1C'][key]]['strandendtimesec']/60.0/60.0, timeh)
    H = Aggregate_merge1d_headerL(statname)
    exptidA = np.array([exptid]*len(timeh))
    A = np.column_stack((
        exptidA,		# col 1-2
        timeh,
        stranddurationsec1T,	# col 3-22
        colmeans['1T']['passfail_mapa'], colmeans['1T']['passfail_mapy'], colmeans['1T']['passfail_mapn'],
        colcount['1T']['passfail_mapa'], colcount['1T']['passfail_mapy'], colcount['1T']['passfail_mapn'],
        colmeans['1T']['passonly_mapa'], colmeans['1T']['passonly_mapy'], colmeans['1T']['passonly_mapn'],
        colcount['1T']['passonly_mapa'], colcount['1T']['passonly_mapy'], colcount['1T']['passonly_mapn'],
        colmeans['1T']['failonly_mapa'], colmeans['1T']['failonly_mapy'], colmeans['1T']['failonly_mapn'],
        colcount['1T']['failonly_mapa'], colcount['1T']['failonly_mapy'], colcount['1T']['failonly_mapn'],
        stranddurationsec1C,	# col 23-41
        colmeans['1C']['passfail_mapa'], colmeans['1C']['passfail_mapy'], colmeans['1C']['passfail_mapn'],
        colcount['1C']['passfail_mapa'], colcount['1C']['passfail_mapy'], colcount['1C']['passfail_mapn'],
        colmeans['1C']['passonly_mapa'], colmeans['1C']['passonly_mapy'], colmeans['1C']['passonly_mapn'],
        colcount['1C']['passonly_mapa'], colcount['1C']['passonly_mapy'], colcount['1C']['passonly_mapn'],
        colmeans['1C']['failonly_mapa'], colmeans['1C']['failonly_mapy'], colmeans['1C']['failonly_mapn'],
        colcount['1C']['failonly_mapa'], colcount['1C']['failonly_mapy'], colcount['1C']['failonly_mapn']
        ))
    AL = [[x if str(x)!='nan' else 'NA' for x in row] for row in A.tolist()]
    outpath = os.path.join(outdir, exptid+'_aggregate_read1d_'+statname+'.txt')
    np.savetxt(outpath, AL, fmt='%s', delimiter='\t', newline='\n', comments='', header='\t'.join(H))
    return 0

def Print_Aggregate_Read2D_Statistics_File(merg2d, outdir, exptid, timeh, stranddurationsec2D, mask, statname):
    colmeans = { '2D':{} }
    colcount = { '2D':{} }
    for key in mask['2D'].keys():
        colmeans['2D'][key], colcount['2D'][key] = binmean(merg2d[mask['2D'][key]][statname], merg2d[mask['2D'][key]]['endtimesec']/60.0/60.0, timeh)
    H = Aggregate_merge2d_headerL(statname)
    exptidA = np.array([exptid]*len(timeh))
    A = np.column_stack((
        exptidA,                # col 1-2
        timeh,
        stranddurationsec2D,    # col 3-22
        colmeans['2D']['passfail_mapa'], colmeans['2D']['passfail_mapy'], colmeans['2D']['passfail_mapn'],
        colcount['2D']['passfail_mapa'], colcount['2D']['passfail_mapy'], colcount['2D']['passfail_mapn'],
        colmeans['2D']['passonly_mapa'], colmeans['2D']['passonly_mapy'], colmeans['2D']['passonly_mapn'],
        colcount['2D']['passonly_mapa'], colcount['2D']['passonly_mapy'], colcount['2D']['passonly_mapn'],
        colmeans['2D']['failonly_mapa'], colmeans['2D']['failonly_mapy'], colmeans['2D']['failonly_mapn'],
        colcount['2D']['failonly_mapa'], colcount['2D']['failonly_mapy'], colcount['2D']['failonly_mapn']
        ))
    AL = [[x if str(x)!='nan' else 'NA' for x in row] for row in A.tolist()]
    outpath = os.path.join(outdir, exptid+'_aggregate_read2d_'+statname+'.txt')
    np.savetxt(outpath, AL, fmt='%s', delimiter='\t', newline='\n', comments='', header='\t'.join(H))
    return 0

def Aggregate_read1d(args, P, mylogger, myhandler, processname, exptid):
  # Create the merg1d file containing the extract and poremapstats data, if it doesn't already exist
    merg1dpath = os.path.join(args.outdir, exptid+'_merged1dstats.txt')
    read1dpath = os.path.join(args.extractdir, exptid+'_read1dstats.txt')
    read1d = np.genfromtxt(read1dpath, skiprows=1, delimiter='\t', dtype=P.ontread1dstatsH, missing_values="NA")
    endtimehrs1T = read1d[read1d[:]['readtype'] == '1T']['strandendtimesec']/60.0/60.0
    endtimehrs1C = read1d[read1d[:]['readtype'] == '1C']['strandendtimesec']/60.0/60.0
    if not os.path.exists(merg1dpath) or os.path.getsize(merg1dpath) == 0:
        retval = Create_merg1d(args, P, exptid, merg1dpath, read1dpath, read1d, endtimehrs1T, endtimehrs1C)
  # Read in the merged table file
    merg1d = np.genfromtxt(merg1dpath, skiprows=1, delimiter='\t', dtype=P.ontmerg1dstatsH, missing_values='NA')
  # Set up time bins, every timebucket hours between 0 and maxrunlen hours
    timeh = np.arange(0, args.maxrunlen+args.timebucket, args.timebucket)
  # Set up various masks for the 1T and 1C read sets
  # 0-based-index: 5=readtype, 28=ismapped, 6=returnstatus
    mask = { '1T':{}, '1C':{} }

    mask['1T']['passfail_mapa'] = merg1d[:]['readtype'] == '1T'
    mask['1T']['passfail_mapy'] = np.logical_and(merg1d[:]['readtype'] == '1T', merg1d[:]['ismapped'] == 1)
    mask['1T']['passfail_mapt'] = np.logical_and(
        np.logical_and(merg1d[:]['readtype'] == '1T', merg1d[:]['ismapped'] == 1),
        ((merg1d[:]['alltargetalignbp']/merg1d[:]['seqlen'].astype(float)>=0.75)))
    mask['1T']['passfail_mapn'] = np.logical_and(merg1d[:]['readtype'] == '1T', merg1d[:]['ismapped'] != 1)

    mask['1T']['passonly_mapa'] = np.logical_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'pass')
    mask['1T']['passonly_mapy'] = np.logical_and(
        np.logical_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'pass'), merg1d[:]['ismapped'] == 1)
    mask['1T']['passonly_mapt'] = np.logical_and(
        np.logical_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'pass'),
        np.logical_and(merg1d[:]['ismapped'] == 1, (merg1d[:]['alltargetalignbp']/merg1d[:]['seqlen'].astype(float)) >= 0.75))
    mask['1T']['passonly_mapn'] = np.logical_and(
        np.logical_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'pass'), merg1d[:]['ismapped'] != 1)

    mask['1T']['failonly_mapa'] = np.logical_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'fail')
    mask['1T']['failonly_mapy'] = np.logical_and(
        np.logical_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'fail'), merg1d[:]['ismapped'] == 1)
    mask['1T']['failonly_mapt'] = np.logical_and(
        np.logical_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'fail'),
        np.logical_and(merg1d[:]['ismapped'] == 1, ((merg1d[:]['alltargetalignbp']/merg1d[:]['seqlen'].astype(float))>=0.75)))
    mask['1T']['failonly_mapn'] = np.logical_and(
        np.logical_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'fail'), merg1d[:]['ismapped'] != 1)

    mask['1C']['passfail_mapa'] = merg1d[:]['readtype'] == '1C'
    mask['1C']['passfail_mapy'] = np.logical_and(merg1d[:]['readtype'] == '1C', merg1d[:]['ismapped'] == 1)
    mask['1C']['passfail_mapt'] = np.logical_and(
        np.logical_and(merg1d[:]['readtype'] == '1C', merg1d[:]['ismapped'] == 1),
        ((merg1d[:]['alltargetalignbp']/merg1d[:]['seqlen'].astype(float))>=0.75))
    mask['1C']['passfail_mapn'] = np.logical_and(merg1d[:]['readtype'] == '1C', merg1d[:]['ismapped'] != 1)

    mask['1C']['passonly_mapa'] = np.logical_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'pass')
    mask['1C']['passonly_mapy'] = np.logical_and(
        np.logical_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'pass'), merg1d[:]['ismapped'] == 1)
    mask['1C']['passonly_mapt'] = np.logical_and(
        np.logical_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'pass'),
        np.logical_and(merg1d[:]['ismapped'] == 1, ((merg1d[:]['alltargetalignbp']/merg1d[:]['seqlen'].astype(float))>=0.75)))
    mask['1C']['passonly_mapn'] = np.logical_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'pass', merg1d[:]['ismapped'] != 1)

    mask['1C']['failonly_mapa'] = np.logical_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'fail')
    mask['1C']['failonly_mapy'] = np.logical_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'fail', merg1d[:]['ismapped'] == 1)
    mask['1C']['failonly_mapt'] = np.logical_and(
        np.logical_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'fail'),
        np.logical_and(merg1d[:]['ismapped'] == 1, ((merg1d[:]['alltargetalignbp']/merg1d[:]['seqlen'].astype(float))>=0.75)))
    mask['1C']['failonly_mapn'] = np.logical_and(
        np.logical_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'fail'),
        merg1d[:]['ismapped'] != 1)
  # Compute the read durations (in seconds) for 1T and 1C components
    stranddurationsec1T, stranddurationsec1T_bincount = binmean(merg1d[mask['1T']['passfail_mapa']]['stranddurationsec'], merg1d[mask['1T']['passfail_mapa']]['strandendtimesec']/60.0/60.0, timeh)
    stranddurationsec1C, stranddurationsec1C_bincount = binmean(merg1d[mask['1C']['passfail_mapa']]['stranddurationsec'], merg1d[mask['1C']['passfail_mapa']]['strandendtimesec']/60.0/60.0, timeh)
  # Collate the aggregate statistics files for each statistic
    Print_Aggregate_Read1D_Statistics_File(merg1d, args.outdir, exptid, timeh, stranddurationsec1T, stranddurationsec1C, mask, 'seqlen')
    Print_Aggregate_Read1D_Statistics_File(merg1d, args.outdir, exptid, timeh, stranddurationsec1T, stranddurationsec1C, mask, 'basespersecond')
    Print_Aggregate_Read1D_Statistics_File(merg1d, args.outdir, exptid, timeh, stranddurationsec1T, stranddurationsec1C, mask, 'meanqscore')
    Print_Aggregate_Read1D_Statistics_File(merg1d, args.outdir, exptid, timeh, stranddurationsec1T, stranddurationsec1C, mask, 'bqmean')
    Print_Aggregate_Read1D_Statistics_File(merg1d, args.outdir, exptid, timeh, stranddurationsec1T, stranddurationsec1C, mask, 'gcpct')
    return 0

def Aggregate_read2d(args, P, mylogger, myhandler, processname, exptid):
  # Create the merg2d file containing the extract and poremapstats data, if it doesn't already exist
    merg2dpath = os.path.join(args.outdir, exptid+'_merged2dstats.txt')
    read2dpath = os.path.join(args.extractdir, exptid+'_read2dstats.txt')
    read2d = np.genfromtxt(read2dpath, skiprows=1, delimiter='\t', dtype=P.ontread2dstatsH, missing_values="NA")
    endtimehrs2D = read2d[:]['endtimesec']/60.0/60.0
    if not os.path.exists(merg2dpath) or os.path.getsize(merg2dpath) == 0:
        retval = Create_merg2d(args, P, exptid, merg2dpath, read2dpath, read2d, endtimehrs2D)
        if retval == 1:
            return 0
  # Read in the merged table file
    merg2d = np.genfromtxt(merg2dpath, skiprows=1, delimiter='\t', dtype=P.ontmerg2dstatsH, missing_values='NA')
  # Set up time bins, every timebucket hours between 0 and maxrunlen hours
    timeh = np.arange(0, args.maxrunlen+args.timebucket, args.timebucket)
  # Set up various masks for the 1T and 1C read sets
  # 0-based-index: 5=readtype, 28=ismapped, 6=returnstatus
    mask = { '2D':{} }

    mask['2D']['passfail_mapa'] = merg2d[:]['exptid'] == exptid
    mask['2D']['passfail_mapy'] = merg2d[:]['ismapped'] == 1
    mask['2D']['passfail_mapt'] = np.logical_and(merg2d[:]['ismapped'] == 1, ((merg1d[:]['alltargetalignbp']/merg1d[:]['seqlen'].astype(float))>=0.75))
    mask['2D']['passfail_mapn'] = merg2d[:]['ismapped'] != 1

    mask['2D']['passonly_mapa'] = merg2d[:]['readclass'] == 'pass'
    mask['2D']['passonly_mapy'] = np.logical_and(merg2d[:]['readclass'] == 'pass', merg2d[:]['ismapped'] == 1)
    mask['2D']['passonly_mapt'] = np.logical_and(
        np.logical_and(merg2d[:]['readclass'] == 'pass', merg2d[:]['ismapped'] == 1),
        ((merg1d[:]['alltargetalignbp']/merg1d[:]['seqlen'].astype(float))>=0.75))
    mask['2D']['passonly_mapn'] = np.logical_and(merg2d[:]['readclass'] == 'pass', merg2d[:]['ismapped'] != 1)

    mask['2D']['failonly_mapa'] = merg2d[:]['readclass'] == 'fail'
    mask['2D']['failonly_mapy'] = np.logical_and(merg2d[:]['readclass'] == 'fail', merg2d[:]['ismapped'] == 1)
    mask['2D']['failonly_mapt'] = np.logical_and(
        np.logical_and(merg2d[:]['readclass'] == 'fail', merg2d[:]['ismapped'] == 1),
        ((merg1d[:]['alltargetalignbp']/merg1d[:]['seqlen'].astype(float))>=0.75))
    mask['2D']['failonly_mapn'] = np.logical_and(merg2d[:]['readclass'] == 'fail', merg2d[:]['ismapped'] != 1)

  # Compute the read durations (in seconds) for 2D components
    durationsec2D, durationsec2D_bincount = binmean(merg2d[mask['2D']['passfail_mapa']]['durationsec'], merg2d[mask['2D']['passfail_mapa']]['endtimesec']/60.0/60.0, timeh)
  # Collate the aggregate statistics files for each statistic
    Print_Aggregate_Read2D_Statistics_File(merg2d, args.outdir, exptid, timeh, durationsec2D, mask, 'seqlen')
    Print_Aggregate_Read2D_Statistics_File(merg2d, args.outdir, exptid, timeh, durationsec2D, mask, 'basespersecond')
    Print_Aggregate_Read2D_Statistics_File(merg2d, args.outdir, exptid, timeh, durationsec2D, mask, 'meanqscore')
    Print_Aggregate_Read2D_Statistics_File(merg2d, args.outdir, exptid, timeh, durationsec2D, mask, 'bqmean')
    Print_Aggregate_Read2D_Statistics_File(merg2d, args.outdir, exptid, timeh, durationsec2D, mask, 'gcpct')
    return 0

def Process(args, P, mylogger, myhandler, processname, exptidL, E):
    'Aggregate per-read ont*stats.txt metrics into windows of X hours.'
    for exptid in exptidL:
        Aggregate_readevent(args, P, mylogger, myhandler, processname, exptid)
        Aggregate_read1d(args, P, mylogger, myhandler, processname, exptid)
        Aggregate_read2d(args, P, mylogger, myhandler, processname, exptid)
    done
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    'Execute this sub-tool.'
    mylogger.info('Started')
    args.overwrite = P.str_2bool(args.overwrite)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    exptidL, E = P.expt_read(args.experiments)
    Prerequisites(args, P, mylogger, myhandler, _processname, exptidL, E)
    Process(args, P, mylogger, myhandler, _processname, exptidL, E)
    mylogger.info('Finished')
    return 0
