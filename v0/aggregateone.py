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
        inpath = os.path.join(args.extractdir, args.exptid+suffix)
        if not os.path.exists(inpath):
            mylogger.error('Input file does not exist {0}'.format(inpath))
            sys.exit(P.err_code('ErrorFileMissing'))
  # If all empty, then probably something wrong
    allhavedata = True
    numfiles = len(P.extractstatfilesuffix)
    numempty = 0
    for suffix in P.extractstatfilesuffix:
        inpath = os.path.join(args.extractdir, args.exptid+suffix)
        if not os.path.getsize(inpath):
            numempty += 1
    if numempty == numfiles:
        mylogger.error('All input files in extractdir are empty, nothing to do ({0})'.format(args.extractdir))
        sys.exit(P.err_code('ErrorEmptyFile'))
    return 0

def binmean(valA, timeA, binA):
    if not len(valA) and not len(timeA):
        return [0.0] * len(binA)
    indexL = np.digitize(timeA, binA)
    binmeans = [np.nanmean([x for x in valA[indexL == i] if x != -1]) for i in range(0, len(binA))]
    return binmeans

def Aggregate_readevent(args, P, mylogger, myhandler, processname, exptid):
  # Output file
    outpath = os.path.join(args.outdir, args.exptid+'_aggregate_readevent.txt')
  # Time bins, every timebucket hours between 0 and maxrunlen hours
    timeh = np.arange(0, args.maxrunlen+args.timebucket, args.timebucket)
    readeventpath = os.path.join(args.extractdir, exptid+'_readeventstats.txt')
    readevent = np.genfromtxt(readeventpath, skiprows=1, delimiter='\t', dtype=P.ontreadeventstatsH, missing_values="NA")
    endtimehrs = readevent[:]['eventendtimesec']/60.0/60.0
  # readevent metrics aggregated by time bins
    readevent_meandurationsec = binmean(readevent[:]['eventdurationsec'], endtimehrs, timeh)
    readevent_meaneventcount = binmean(readevent[:]['eventcount'], endtimehrs, timeh)
    readevent_meaneventspersec = binmean(readevent[:]['eventspersec'], endtimehrs, timeh)
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

def Parse_poremapstats_readstats(passpath, failpath):
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
    # Previously, tried the following.
    #stats1tpass = np.genfromtxt(stats1tpasspath, skiprows=0, delimiter='\t', dtype='str', missing_values='NA')
    #stats1tfail = np.genfromtxt(stats1tfailpath, skiprows=0, delimiter='\t', dtype='str', missing_values='NA')
    #stats1cpass = np.genfromtxt(stats1cpasspath, skiprows=0, delimiter='\t', dtype='str', missing_values='NA')
    #stats1cfail = np.genfromtxt(stats1cfailpath, skiprows=0, delimiter='\t', dtype='str', missing_values='NA')

def Aggregate_merge1d_headerL(statname):
    H = ['exptid', 'timeh']
    for readtype in ['1T', '1C']:
        H.append('durationsec_'+readtype)
        for readclass in ['passfail', 'passonly', 'failonly']:
            for maptype in ['mapa', 'mapy', 'mapn']:
                s = statname+'_'+readtype+'_'+readclass+'_'+maptype
                H.append(s)
    return H

def Create_merg1d(args, P, exptid, merg1dpath, read1dpath, read1d, endtimehrs1T, endtimehrs1C):
  # Read in the 'marcoporo extract' output file
    #read1dpath = os.path.join(args.extractdir, exptid+'_read1dstats.txt')
    #read1d = np.genfromtxt(read1dpath, skiprows=1, delimiter='\t', dtype=P.ontread1dstatsH, missing_values="NA")
    #endtimehrs1T = read1d[read1d[:]['readtype'] == '1T']['strandendtimesec']/60.0/60.0
    #endtimehrs1C = read1d[read1d[:]['readtype'] == '1C']['strandendtimesec']/60.0/60.0
  # Read in the poremapstats 'readstats.txt' output files for 1D reads
    stats1tpasspath = os.path.join(args.bwamemdir,
        '{exptid}_{readtype}_{readclass}'.format(exptid=exptid, readtype='1T', readclass='pass'),
        '{exptid}_{readtype}_{readclass}_readstats.txt'.format(exptid=exptid, readtype='1T', readclass='pass'))
    stats1tfailpath = os.path.join(args.bwamemdir,
        '{exptid}_{readtype}_{readclass}'.format(exptid=exptid, readtype='1T', readclass='fail'),
        '{exptid}_{readtype}_{readclass}_readstats.txt'.format(exptid=exptid, readtype='1T', readclass='fail'))
    stats1cpasspath = os.path.join(args.bwamemdir,
        '{exptid}_{readtype}_{readclass}'.format(exptid=exptid, readtype='1C', readclass='pass'),
        '{exptid}_{readtype}_{readclass}_readstats.txt'.format(exptid=exptid, readtype='1C', readclass='pass'))
    stats1cfailpath = os.path.join(args.bwamemdir,
        '{exptid}_{readtype}_{readclass}'.format(exptid=exptid, readtype='1C', readclass='fail'),
        '{exptid}_{readtype}_{readclass}_readstats.txt'.format(exptid=exptid, readtype='1C', readclass='fail'))
    stats1t = Parse_poremapstats_readstats(stats1tpasspath, stats1tfailpath)
    stats1c = Parse_poremapstats_readstats(stats1cpasspath, stats1cfailpath)
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
            #rowend = stats1t['data'][key][5:] if readtype == '1T' else stats1c['data'][key][5:]
            row = list(rowstt)[:-2] + rowend
            row = ['NA' if (str(x)=='nan' or str(x) == '-1' or not len(str(x))) else str(x) for x in row]
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in row])))
            out_fp.flush()
    return 0

def Print_Aggregate_Statistics_File(merg1d, outdir, exptid, timeh, stranddurationsec1T, stranddurationsec1C, mask, statname):

    coldata = { '1T':{}, '1C':{} }
    for key in mask['1T'].keys():
        coldata['1T'][key] = binmean(merg1d[mask['1T'][key]][statname], merg1d[mask['1T'][key]]['strandendtimesec']/60.0/60.0, timeh)
    for key in mask['1C'].keys():
        coldata['1C'][key] = binmean(merg1d[mask['1C'][key]][statname], merg1d[mask['1C'][key]]['strandendtimesec']/60.0/60.0, timeh)
    H = Aggregate_merge1d_headerL(statname)
    exptidA = np.array([exptid]*len(timeh))
    A = np.column_stack((
        exptidA,		# col 1-2
        timeh,
        stranddurationsec1T,	# col 3-12
        coldata['1T']['passfail_mapa'], coldata['1T']['passfail_mapy'], coldata['1T']['passfail_mapn'],
        coldata['1T']['passonly_mapa'], coldata['1T']['passonly_mapy'], coldata['1T']['passonly_mapn'],
        coldata['1T']['failonly_mapa'], coldata['1T']['failonly_mapy'], coldata['1T']['failonly_mapn'],
        stranddurationsec1C,	# col 13-22
        coldata['1C']['passfail_mapa'], coldata['1C']['passfail_mapy'], coldata['1C']['passfail_mapn'],
        coldata['1C']['passonly_mapa'], coldata['1C']['passonly_mapy'], coldata['1C']['passonly_mapn'],
        coldata['1C']['failonly_mapa'], coldata['1C']['failonly_mapy'], coldata['1C']['failonly_mapn']
        ))
    AL = [[x if str(x)!='nan' else 'NA' for x in row] for row in A.tolist()]
    outpath = os.path.join(outdir, exptid+'_aggregate_read1d_'+statname+'.txt')
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
    #merg1d = np.genfromtxt(merg1dpath, skiprows=0, delimiter='\t', dtype='str', missing_values='NA')
    #merg1dpath = os.path.join(args.outdir, exptid+'_merged1dstats_100rows.txt')
    merg1d = np.genfromtxt(merg1dpath, skiprows=1, delimiter='\t', dtype=P.ontmerg1dstatsH, missing_values='NA')

  # Set up time bins, every timebucket hours between 0 and maxrunlen hours
    timeh = np.arange(0, args.maxrunlen+args.timebucket, args.timebucket)
  # Set up various masks for the 1T and 1C read sets
  # 0-based-index: 5=readtype, 28=ismapped, 6=returnstatus
    mask = { '1T':{}, '1C':{} }
    mask['1T']['passfail_mapa'] = merg1d[:]['readtype'] == '1T'
    mask['1T']['passfail_mapy'] = np.bitwise_and(merg1d[:]['readtype'] == '1T', merg1d[:]['ismapped'] == 1)
    mask['1T']['passfail_mapn'] = np.bitwise_and(merg1d[:]['readtype'] == '1T', merg1d[:]['ismapped'] != 1)
    mask['1T']['passonly_mapa'] = np.bitwise_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'pass')
    mask['1T']['passonly_mapy'] = np.bitwise_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'pass', merg1d[:]['ismapped'] == 1)
    mask['1T']['passonly_mapn'] = np.bitwise_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'pass', merg1d[:]['ismapped'] != 1)
    mask['1T']['failonly_mapa'] = np.bitwise_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'fail')
    mask['1T']['failonly_mapy'] = np.bitwise_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'fail', merg1d[:]['ismapped'] == 1)
    mask['1T']['failonly_mapn'] = np.bitwise_and(merg1d[:]['readtype'] == '1T', merg1d[:]['returnstatus'] == 'fail', merg1d[:]['ismapped'] != 1)
    mask['1C']['passfail_mapa'] = merg1d[:]['readtype'] == '1C'
    mask['1C']['passfail_mapy'] = np.bitwise_and(merg1d[:]['readtype'] == '1C', merg1d[:]['ismapped'] == 1)
    mask['1C']['passfail_mapn'] = np.bitwise_and(merg1d[:]['readtype'] == '1C', merg1d[:]['ismapped'] != 1)
    mask['1C']['passonly_mapa'] = np.bitwise_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'pass')
    mask['1C']['passonly_mapy'] = np.bitwise_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'pass', merg1d[:]['ismapped'] == 1)
    mask['1C']['passonly_mapn'] = np.bitwise_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'pass', merg1d[:]['ismapped'] != 1)
    mask['1C']['failonly_mapa'] = np.bitwise_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'fail')
    mask['1C']['failonly_mapy'] = np.bitwise_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'fail', merg1d[:]['ismapped'] == 1)
    mask['1C']['failonly_mapn'] = np.bitwise_and(merg1d[:]['readtype'] == '1C', merg1d[:]['returnstatus'] == 'fail', merg1d[:]['ismapped'] != 1)
  # Compute the read durations (in seconds) for 1T and 1C components
    #stranddurationsec1T = binmean(merg1d[mask['1T']['passfail_mapa']]['stranddurationsec'], endtimehrs1T, timeh)
    #stranddurationsec1C = binmean(merg1d[mask['1C']['passfail_mapa']]['stranddurationsec'], endtimehrs1C, timeh)
    stranddurationsec1T = binmean(merg1d[mask['1T']['passfail_mapa']]['stranddurationsec'], merg1d[mask['1T']['passfail_mapa']]['strandendtimesec']/60.0/60.0, timeh)
    stranddurationsec1C = binmean(merg1d[mask['1C']['passfail_mapa']]['stranddurationsec'], merg1d[mask['1C']['passfail_mapa']]['strandendtimesec']/60.0/60.0, timeh)
  # Collate the aggregate statistics files for each statistic
    Print_Aggregate_Statistics_File(merg1d, args.outdir, exptid, timeh, stranddurationsec1T, stranddurationsec1C, mask, 'seqlen')
    Print_Aggregate_Statistics_File(merg1d, args.outdir, exptid, timeh, stranddurationsec1T, stranddurationsec1C, mask, 'basespersecond')
    Print_Aggregate_Statistics_File(merg1d, args.outdir, exptid, timeh, stranddurationsec1T, stranddurationsec1C, mask, 'meanqscore')
    Print_Aggregate_Statistics_File(merg1d, args.outdir, exptid, timeh, stranddurationsec1T, stranddurationsec1C, mask, 'bqmean')
    Print_Aggregate_Statistics_File(merg1d, args.outdir, exptid, timeh, stranddurationsec1T, stranddurationsec1C, mask, 'gcpct')

    return 0

  # Aggregated gcpct
    merg1d_1T_gcpct_passfail_mapa = binmean(merg1d[merg1t_passfail_mapa_mask]['gcpct'], endtimehrs1T, timeh)
    merg1d_1T_gcpct_passfail_mapy = binmean(merg1d[merg1t_passfail_mapy_mask]['gcpct'], endtimehrs1T, timeh)
    merg1d_1T_gcpct_passfail_mapn = binmean(merg1d[merg1t_passfail_mapn_mask]['gcpct'], endtimehrs1T, timeh)
    merg1d_1T_gcpct_passonly_mapa = binmean(merg1d[merg1t_passonly_mapa_mask]['gcpct'], endtimehrs1T, timeh)
    merg1d_1T_gcpct_passonly_mapy = binmean(merg1d[merg1t_passonly_mapy_mask]['gcpct'], endtimehrs1T, timeh)
    merg1d_1T_gcpct_passonly_mapn = binmean(merg1d[merg1t_passonly_mapn_mask]['gcpct'], endtimehrs1T, timeh)
    merg1d_1T_gcpct_failonly_mapa = binmean(merg1d[merg1t_failonly_mapa_mask]['gcpct'], endtimehrs1T, timeh)
    merg1d_1T_gcpct_failonly_mapy = binmean(merg1d[merg1t_failonly_mapy_mask]['gcpct'], endtimehrs1T, timeh)
    merg1d_1T_gcpct_failonly_mapn = binmean(merg1d[merg1t_failonly_mapn_mask]['gcpct'], endtimehrs1T, timeh)
    merg1d_1C_gcpct_passfail_mapa = binmean(merg1d[merg1c_passfail_mapa_mask]['gcpct'], endtimehrs1C, timeh)
    merg1d_1C_gcpct_passfail_mapy = binmean(merg1d[merg1c_passfail_mapy_mask]['gcpct'], endtimehrs1C, timeh)
    merg1d_1C_gcpct_passfail_mapn = binmean(merg1d[merg1c_passfail_mapn_mask]['gcpct'], endtimehrs1C, timeh)
    merg1d_1C_gcpct_passonly_mapa = binmean(merg1d[merg1c_passonly_mapa_mask]['gcpct'], endtimehrs1C, timeh)
    merg1d_1C_gcpct_passonly_mapy = binmean(merg1d[merg1c_passonly_mapy_mask]['gcpct'], endtimehrs1C, timeh)
    merg1d_1C_gcpct_passonly_mapn = binmean(merg1d[merg1c_passonly_mapn_mask]['gcpct'], endtimehrs1C, timeh)
    merg1d_1C_gcpct_failonly_mapa = binmean(merg1d[merg1c_failonly_mapa_mask]['gcpct'], endtimehrs1C, timeh)
    merg1d_1C_gcpct_failonly_mapy = binmean(merg1d[merg1c_failonly_mapy_mask]['gcpct'], endtimehrs1C, timeh)
    merg1d_1C_gcpct_failonly_mapn = binmean(merg1d[merg1c_failonly_mapn_mask]['gcpct'], endtimehrs1C, timeh)
    H = Aggregate_merge1d_headerL(gcpct)
    exptidA = np.array([exptid]*len(timeh))
    A = np.column_stack((
        exptidA,
        timeh,
        merg1d_1T_stranddurationsec,
        merg1d_1T_gcpct_passfail_mapa, merg1d_1T_gcpct_passfail_mapy, merg1d_1T_gcpct_passfail_mapn,
        merg1d_1T_gcpct_passonly_mapa, merg1d_1T_gcpct_passonly_mapy, merg1d_1T_gcpct_passonly_mapn,
        merg1d_1T_gcpct_failonly_mapa, merg1d_1T_gcpct_failonly_mapy, merg1d_1T_gcpct_failonly_mapn,
        merg1d_1C_stranddurationsec,
        merg1d_1C_gcpct_passfail_mapa, merg1d_1C_gcpct_passfail_mapy, merg1d_1C_gcpct_passfail_mapn,
        merg1d_1C_gcpct_passonly_mapa, merg1d_1C_gcpct_passonly_mapy, merg1d_1C_gcpct_passonly_mapn,
        merg1d_1C_gcpct_failonly_mapa, merg1d_1C_gcpct_failonly_mapy, merg1d_1C_gcpct_failonly_mapn
        ))
    outpath = os.path.join(args.outdir, args.exptid+'_aggregate_read1d_gcpct.txt')
    np.savetxt(outpath, A, fmt='%s', delimiter='\t', newline='\n', comments='', header='\t'.join(H))
    return 0

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
    #H = ['exptid', 'timehr',
    #     'durationsec1T', 'qscore1T', 'seqlen1T', 'bq1T', 'gcpct1T', 'basesps1T',
    #     'durationsec1C', 'qscore1C', 'seqlen1C', 'bq1C', 'gcpct1C', 'basesps1C']

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
    read2dpath = os.path.join(args.extractdir, exptid+'_read2dstats.txt')
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
    # XXXX - REMOVED TO DEBUG
    #Aggregate_readevent(args, P, mylogger, myhandler, processname, exptid)
    Aggregate_read1d(args, P, mylogger, myhandler, processname, exptid)
    #Aggregate_read2d(args, P, mylogger, myhandler, processname, exptid)
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    'Execute this sub-tool.'
    mylogger.info('Started')
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname, args.exptid)
    mylogger.info('Finished')
    return 0
