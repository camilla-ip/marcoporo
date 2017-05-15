#!/usr/bin/env python

#
# Known issues:
# 1. The FAST5 from ONT contains a FASTQ record like:
#    @fd0f474b-5264-4a67-9513-c55b62c6f29b_Basecall_2D_complemena
#    To be unique, it should have the instanceN as part of the name, for example:
#    @fd0f474b-5264-4a67-9513-c55b62c6f29b_Basecall_2D_000_complement
#

import h5py
import itertools
import logging
import numpy as np
import os
import random
import re
import sys
import time
from Bio.SeqUtils import GC

import marcoporoversion

_processname = 'extract'
#_messageinterval = 1000

def Attr(attrD, defaultvalue, pathL):
    'Return the first retrievable value in pathL.'
    result = defaultvalue
    for path in pathL:
        if attrD.has_key(path):
            result = attrD[path][1]
            break
    return result

def Get_Batchid(fast5path):
    '''
    For inferring batchids, can cope with fast5 files like:
    R9: MinION2_20160802_FNFAD22824_MN16454_sequencing_run_Chip93_MARC_R9_1D_UBC_77825_ch100_read249_strand.fast5 -> 77825
    R7: makeson_PC_MA_286_R7.3_MARC_K12_Ib_04_16_15_5048_1_ch299_file184_strand.fast5 -> 5048_1
    '''
    batchid = None
    file = os.path.basename(fast5path)
    filestem = '.'.join(file.split('.')[0:-1])
    L = filestem.split('_')
    if L[-1].startswith('strand') and L[-2].startswith('read') and L[-3].startswith('ch'):
        batchid = L[-4]
    elif L[-1].startswith('strand') and L[-2].startswith('file') and L[-3].startswith('ch'):
        batchid = '_'.join(L[-5:-3])
    return batchid

def Extract_Pairs(P, constD, exptid, batchid, instanceN, attrD, fpD):
    attrL = attrD.keys()
    attrL.sort()
    for attr in attrL:
        var = attr
        val = attrD[var][1]
        newvar = P.fast5_attribute_to_instanceN(attr)
        if constD.has_key(newvar):
            row = [exptid, batchid, instanceN, newvar, val]
            fpD['exptpairs'].write('{0}\n'.format('\t'.join(row)))
        else:
            try:
                readid = attrD['Analyses/EventDetection_{0}/Configuration/general/uuid'.format(instanceN)][1]
            except:
                readid = 'NK'       # Need to change this to something based on the filename (chN and readN)
            row = [exptid, batchid, readid, instanceN, newvar, val]
            fpD['readpairs'].write('{0}\n'.format('\t'.join(row)))
    fpD['exptpairs'].flush()
    fpD['readpairs'].flush()
    return 0

def Extract_Fastq(mylogger, readclass, fastqD, fpD):
    if readclass == 'pass':
        if fastqD['1T'] is not None:
            fpD['fq1Tpass'].write('{0}\n'.format('\n'.join(fastqD['1T'])))
            fpD['fq1Tpass'].flush()
        if fastqD['1C'] is not None:
            fpD['fq1Cpass'].write('{0}\n'.format('\n'.join(fastqD['1C'])))
            fpD['fq1Cpass'].flush()
        if fastqD['2D'] is not None:
            fpD['fq2Dpass'].write('{0}\n'.format('\n'.join(fastqD['2D'])))
            fpD['fq2Dpass'].flush()
    elif readclass == 'fail':
        if fastqD['1T'] is not None:
            fpD['fq1Tfail'].write('{0}\n'.format('\n'.join(fastqD['1T'])))
            fpD['fq1Tfail'].flush()
        if fastqD['1C'] is not None:
            fpD['fq1Cfail'].write('{0}\n'.format('\n'.join(fastqD['1C'])))
            fpD['fq1Cfail'].flush()
        if fastqD['2D'] is not None:
            fpD['fq2Dfail'].write('{0}\n'.format('\n'.join(fastqD['2D'])))
            fpD['fq2Dfail'].flush()
    else:
        mylogger.error('Unrecognised readclass {0} ({1})'.format(readclass, fast5path))
    return 0

def Print_ontexptstats(P, exptid, batchid, instanceN, attrD, fpD):
    'Only one row per experiment required.'
    runid = attrD['UniqueGlobalKey/tracking_id/run_id'][1]
    samplingrate = float(attrD['UniqueGlobalKey/channel_id/sampling_rate'][1])
    asicid = attrD['UniqueGlobalKey/tracking_id/asic_id'][1]
    deviceid = attrD['UniqueGlobalKey/tracking_id/device_id'][1]
    flowcellid = attrD['UniqueGlobalKey/tracking_id/flow_cell_id'][1]
    scriptname = attrD['UniqueGlobalKey/tracking_id/exp_script_name'][1]
    scriptpurpose = attrD['UniqueGlobalKey/tracking_id/exp_script_purpose'][1]
    expstarttime = int(attrD['UniqueGlobalKey/tracking_id/exp_start_time'][1])
    expstarttimeiso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(expstarttime))
    try:
        versionname = attrD['UniqueGlobalKey/tracking_id/version_name'][1]
    except:
        versionname = ''
    version = attrD['UniqueGlobalKey/tracking_id/version'][1]
    workflowfullname = Attr(attrD, 'NK',
        ['Analyses/EventDetection_{0}/Configuration/general/workflow_name'.format(instanceN),
         'Analyses/Basecall_2D_{0}/Configuration/recipes/basecall'.format(instanceN)])
    comment = ''
    dbuserid = ''
    rowNP = np.array(
        (exptid, batchid, runid, samplingrate, asicid,
        deviceid, flowcellid, scriptname, scriptpurpose, expstarttime,
        expstarttimeiso, versionname, version, workflowfullname,
        comment, dbuserid),
        dtype=P.ontexptstatsH)
    rowL = list(np.atleast_1d(rowNP).tolist()[0])
    rowL = [x if x != -1 and x != '-1' and x != '-1.0' else 'NA' for x in rowL]
    fpD['exptstats'].write('{0}\n'.format('\t'.join([str(x) for x in rowL])))
    return 0

def Print_ontreadstats(P, exptid, batchid, readclass, instanceN, attrD, fpD, fast5file):
  # Intermediate
    #readnumberS = attrD['Analyses/Hairpin_Split_{0}/Configuration/general/read_id'.format(instanceN)][1]
    readnumberS = Attr(attrD, '-1',
        ['Analyses/Hairpin_Split_{0}/Configuration/general/read_id'.format(instanceN),
         'Analyses/Basecall_2D_{0}/Configuration/general/read_id'.format(instanceN),
         'Analyses/Segment_Linear_{0}/Configuration/general/read_id'.format(instanceN)])
    if readnumberS == '-1':
        return 1
    try:
        exp_start_time = float(attrD['UniqueGlobalKey/tracking_id/exp_start_time'][1])
    except:
        return 1
    try:
        samplingrate = float(attrD['UniqueGlobalKey/channel_id/sampling_rate'][1])
    except:
        return 1
  # Returned
    try:
        runid = attrD['UniqueGlobalKey/tracking_id/run_id'][1]
    except:
        return 1
    #readid = attrD['Raw/Reads/Read_{0}/read_id'.format(readnumberS)][1]
    readid = Attr(attrD, 'NK',
        ['Raw/Reads/Read_{0}/read_id'.format(readnumberS),
         'Analyses/EventDetection_{0}/Reads/Read_{1}/read_id'.format(instanceN, readnumberS)])
    filename = re.sub('.fast5$', '', os.path.basename(fast5file))
    channelnumber = int(attrD['UniqueGlobalKey/channel_id/channel_number'][1])
    readnumber = int(readnumberS)
    #filenumber = int(attrD['Analyses/EventDetection_{0}/Configuration/general/file_number'.format(instanceN)][1])
    filenumber = Attr(attrD, '-1',
        ['Analyses/EventDetection_{0}/Configuration/general/file_number'.format(instanceN),
         'Analyses/Basecall_2D_{0}/Configuration/general/file_number'.format(instanceN)])
    #readclass ='X'
    asictemp = float(attrD['UniqueGlobalKey/tracking_id/asic_temp'][1])
    heatsinktemp = float(attrD['UniqueGlobalKey/tracking_id/heatsink_temp'][1])
    #readstarttime = float(attrD['Raw/Reads/Read_{0}/start_time'.format(readnumberS)][1])
    readstarttime = float(Attr(attrD, '-1',
        ['Raw/Reads/Read_{0}/start_time'.format(readnumberS),
         'Analyses/EventDetection_{0}/Reads/Read_{1}/start_time'.format(instanceN, readnumberS)]))
    #readduration = float(attrD['Raw/Reads/Read_{0}/duration'.format(readnumberS)][1])
    readduration = float(Attr(attrD, '0',
        ['Raw/Reads/Read_{0}/duration'.format(readnumberS),
         'Analyses/EventDetection_{0}/Reads/Read_{1}/duration'.format(instanceN, readnumberS)]))
    #readstarttimesec = exp_start_time + readstarttime / samplingrate
    #readendtimesec = exp_start_time + readstarttime / samplingrate + readduration / samplingrate
    readstarttimesec = (readstarttime / samplingrate) if samplingrate else 0.0
    readendtimesec = (readstarttimesec + readduration / samplingrate) if samplingrate else 0.0
    readdurationsec = (readduration / samplingrate) if samplingrate else 0.0
    readstarttimeiso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(exp_start_time + readstarttimesec))
    readendtimeiso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(exp_start_time + readendtimesec))
    comment = ''
    dbuserid = ''
    rowNP = np.array(
        (exptid, batchid, runid, readid, filename,
         channelnumber, readnumber, filenumber, readclass, asictemp,
         heatsinktemp, readstarttime, readduration, readstarttimesec, readendtimesec,
         readdurationsec, readstarttimeiso, readendtimeiso, comment, dbuserid),
        dtype=P.ontreadstatsH)
    rowL = list(np.atleast_1d(rowNP).tolist()[0])
    rowL = [x if x != -1 and x != '-1' and x != '-1.0' else 'NA' for x in rowL]
    fpD['readstats'].write('{0}\n'.format('\t'.join([str(x) for x in rowL])))
    fpD['readstats'].flush()
    return 0

def Print_ontreadeventstats(P, exptid, batchid, readclass, instanceN, attrD, fpD):
  # Intermediate
    #readnumberS = attrD['Analyses/Hairpin_Split_{0}/Configuration/general/read_id'.format(instanceN)][1]
    readnumberS = Attr(attrD, '-1',
        ['Analyses/Hairpin_Split_{0}/Configuration/general/read_id'.format(instanceN),
         'Analyses/Basecall_2D_{0}/Configuration/general/read_id'.format(instanceN),
         'Analyses/Segment_Linear_{0}/Configuration/general/read_id'.format(instanceN)])
    if readnumberS == '-1':
        return 1
    exp_start_time = float(attrD['UniqueGlobalKey/tracking_id/exp_start_time'][1])
    samplingrate = float(attrD['UniqueGlobalKey/channel_id/sampling_rate'][1])
  # Returned
    runid = attrD['UniqueGlobalKey/tracking_id/run_id'][1]
    #readid = attrD['Raw/Reads/Read_{0}/read_id'.format(readnumberS)][1]
    readid = Attr(attrD, '0',
        ['Raw/Reads/Read_{0}/read_id'.format(readnumberS),
         'Analyses/EventDetection_{0}/Reads/Read_{1}/read_id'.format(instanceN, readnumberS)])
    eventinstanceN = instanceN
    returnstatus = Attr(attrD, 'NR',
        ['Analyses/EventDetection_{0}/Summary/return_status'.format(eventinstanceN)])
    eventstarttime = float(attrD['Analyses/EventDetection_{0}/Reads/Read_{1}/start_time'.format(eventinstanceN, readnumberS)][1])
    eventduration = float(attrD['Analyses/EventDetection_{0}/Reads/Read_{1}/duration'.format(eventinstanceN, readnumberS)][1])
    #eventstarttimesec = exp_start_time + eventstarttime / samplingrate
    #eventendtimesec = exp_start_time + eventstarttime / samplingrate + eventduration / samplingrate
    eventstarttimesec = eventstarttime / samplingrate if samplingrate else -1
    eventendtimesec = (eventstarttimesec + eventduration / samplingrate) if samplingrate else -1
    eventdurationsec = (eventduration / samplingrate) if samplingrate else -1
    eventstarttimeiso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(exp_start_time + eventstarttimesec))
    eventendtimeiso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(exp_start_time + eventendtimesec))
    eventcount = int(Attr(attrD, '-1',
        ['Analyses/EventDetection_{0}/Summary/event_detection/num_events'.format(eventinstanceN),
         'event_count']))
    #if attrD.has_key('Analyses/EventDetection_{0}/Summary/event_detection/num_events'.format(eventinstanceN)):
    #    eventcount = int(attrD['Analyses/EventDetection_{0}/Summary/event_detection/num_events'.format(eventinstanceN)][1])
    #elif attrD.has_key('Analyses/EventDetection_{0}/Reads/Read_{0}/Events'.format(eventinstanceN, readnumberS)):
    #    key = 'Analyses/EventDetection_{0}/Reads/Read_{0}/Events'.format(eventinstanceN, readnumberS)
    #    eventcount = len(hdf[key][()])
    eventspersec = (eventcount / eventdurationsec) if (eventdurationsec and eventcount > 100) else -1
    comment = ''
    dbuserid = ''
    rowNP = np.array(
        (exptid, batchid, runid, readid, eventinstanceN,
        readclass, returnstatus, eventstarttime, eventduration, eventstarttimesec,
        eventendtimesec, eventdurationsec, eventstarttimeiso, eventendtimeiso, eventcount,
        eventspersec, comment, dbuserid),
        dtype=P.ontreadeventstatsH)
    rowL = list(np.atleast_1d(rowNP).tolist()[0])
    rowL = [x if x != -1 and x != '-1' and x != '-1.0' else 'NA' for x in rowL]
    fpD['readeventstats'].write('{0}\n'.format('\t'.join([str(x) for x in rowL])))
    fpD['readeventstats'].flush()
    return 0

def Print_read1dblank(P, exptid, batchid, fpD):
    rowNP = np.array(
        (exptid, batchid, 'NA', 'NA', 'NA',
        'NA', 'NA', -1, -1, -1,
        -1, -1.0, -1.0, -1.0, 'NA',
        'NA', -1.0, -1.0, -1, -1,
        -1.0, -1.0, -1.0, -1.0, '',
        ''),
        dtype=P.ontread1dstatsH)
    rowL = list(np.atleast_1d(rowNP).tolist()[0])
    rowL = [x if x != -1 and x != '-1' and x != '-1.0' else 'NA' for x in rowL]
    fpD['read1dstats'].write('{0}\n'.format('\t'.join([str(x) for x in rowL])))
    fpD['read1dstats'].flush()
    return 0

def Print_ontread1tstats(P, exptid, batchid, readclass, instanceN, readtype, attrD, fastqD, returnstatus, fpD):
  # Intermediate
    #readnumberS = attrD['Analyses/Hairpin_Split_{0}/Configuration/general/read_id'.format(instanceN)][1]
    readnumberS = Attr(attrD, '-1',
        ['Analyses/EventDetection_{0}/Configuration/general/read_id'.format(instanceN),
         'Analyses/Hairpin_Split_{0}/Configuration/general/read_id'.format(instanceN),
         'Analyses/Basecall_2D_{0}/Configuration/general/read_id'.format(instanceN),
         'Analyses/Segment_Linear_{0}/Configuration/general/read_id'.format(instanceN)])
    if readnumberS == '-1':
        Print_read1dblank(P, exptid, batchid, fpD)
        return None
    try:
        exp_start_time = float(attrD['UniqueGlobalKey/tracking_id/exp_start_time'][1])
        samplingrate = float(attrD['UniqueGlobalKey/channel_id/sampling_rate'][1])
        hpsinstanceN = instanceN
        bqnumA = np.array([ord(x)-33 for x in fastqD[readtype][3]]) if fastqD[readtype] is not None else None
        bc2dinstanceN = instanceN
      # Returned
        runid = attrD['UniqueGlobalKey/tracking_id/run_id'][1]
        #readid = attrD['Raw/Reads/Read_{0}/read_id'.format(readnumberS)][1]
        readid = Attr(attrD, '-1',
            ['Raw/Reads/Read_{0}/read_id'.format(readnumberS),
             'Analyses/EventDetection_{0}/Reads/Read_{1}/read_id'.format(instanceN, readnumberS)])
        bc1dinstanceN = instanceN
        numevents = int(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/Summary/basecall_1d_template/num_events'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/Summary/basecall_1d_template/num_events'.format(bc2dinstanceN)]))
        numskips = int(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/Summary/basecall_1d_template/num_skips'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/Summary/basecall_1d_template/num_skips'.format(bc2dinstanceN)]))
        numstays = int(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/Summary/basecall_1d_template/num_stays'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/Summary/basecall_1d_template/num_stays'.format(bc2dinstanceN)]))
        numcalled = int(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/Summary/basecall_1d_template/called_events'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/Summary/basecall_1d_template/called_events'.format(bc2dinstanceN)]))
        strandstarttimesec = float(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/BaseCalled_template/Events/start_time'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/BaseCalled_template/Events/start_time'.format(bc2dinstanceN)]))
        stranddurationsec = float(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/BaseCalled_template/Events/duration'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/BaseCalled_template/Events/duration'.format(bc2dinstanceN)]))
        strandendtimesec = (strandstarttimesec + stranddurationsec) if strandstarttimesec != -1 and stranddurationsec != -1 else -1 
        #strandstarttimesec = exp_start_time + strandstarttime / samplingrate
        #strandendtimesec = exp_start_time + strandstarttime / samplingrate + strandduration / samplingrate
        #strandstarttimesec = strandstarttime / samplingrate
        #strandendtimesec = strandstarttimesec + strandduration / samplingrate
        #stranddurationsec = strandduration / samplingrate
        if strandstarttimesec != '-1':
            strandstarttimeiso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(exp_start_time + strandstarttimesec))
            strandendtimeiso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(exp_start_time + strandendtimesec))
        else:
            strandstarttimeiso = 'NA'
            strandendtimeiso = 'NA'
        meanqscore = float(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/Summary/basecall_1d_template/mean_qscore'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/Summary/basecall_1d_template/mean_qscore'.format(bc1dinstanceN)]))
        strandscore = float(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/Summary/basecall_1d_template/strand_score'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/Summary/basecall_1d_template/strand_score'.format(bc1dinstanceN)]))
        if attrD.has_key('Analyses/Basecall_1D_{0}/Summary/basecall_1d_template/sequence_length'.format(bc1dinstanceN)):
            seqlen = int(attrD['Analyses/Basecall_1D_{0}/Summary/basecall_1d_template/sequence_length'.format(bc1dinstanceN)][1])
        else:
            seqlen = len(fastqD[readtype][1]) if fastqD[readtype] is not None else -1
        bqlen = len(fastqD[readtype][3]) if fastqD[readtype] is not None else 0
        bqmean = np.mean(bqnumA) if bqnumA is not None else -1
        bqmedian = np.median(bqnumA) if bqnumA is not None else -1
        gcpct = round(GC(fastqD[readtype][1]), 1) if fastqD[readtype] is not None else -1
        basespersecond = seqlen / stranddurationsec if stranddurationsec != -1 else -1
        comment = ''
        dbuserid = ''
    except:
        Print_read1dblank(P, exptid, batchid, fpD)
        return None
    rowNP = np.array(
        (exptid, batchid, runid, readid, bc1dinstanceN,
        readtype, readclass, returnstatus, numevents, numskips,
        numstays, numcalled, strandstarttimesec, strandendtimesec, stranddurationsec,
        strandstarttimeiso, strandendtimeiso, meanqscore, strandscore, seqlen,
        bqlen, bqmean, bqmedian, gcpct, basespersecond,
        comment, dbuserid),
        dtype=P.ontread1dstatsH)
    rowL = list(np.atleast_1d(rowNP).tolist()[0])
    rowL = [x if x != -1 and x != '-1' and x != '-1.0' else 'NA' for x in rowL]
    fpD['read1dstats'].write('{0}\n'.format('\t'.join([str(x) for x in rowL])))
    fpD['read1dstats'].flush()
    return rowL

def Print_ontread1cstats(P, exptid, batchid, readclass, instanceN, readtype, attrD, fastqD, returnstatus, fpD):
  # Intermediate
    #readnumberS = attrD['Analyses/Hairpin_Split_{0}/Configuration/general/read_id'.format(instanceN)][1]
    readnumberS = Attr(attrD,  '-1',
        ['Analyses/EventDetection_{0}/Configuration/general/read_id'.format(instanceN),
         'Analyses/Hairpin_Split_{0}/Configuration/general/read_id'.format(instanceN),
         'Analyses/Basecall_2D_{0}/Configuration/general/read_id'.format(instanceN),
         'Analyses/Segment_Linear_{0}/Configuration/general/read_id'.format(instanceN)])
    if readnumberS == '-1':
        Print_read1dblank(P, exptid, batchid, fpD)
        return None
    try:
        exp_start_time = float(attrD['UniqueGlobalKey/tracking_id/exp_start_time'][1])
        samplingrate = float(attrD['UniqueGlobalKey/channel_id/sampling_rate'][1])
        hpsinstanceN = instanceN
        bqnumA = np.array([ord(x)-33 for x in fastqD[readtype][3]]) if fastqD[readtype] is not None else None
        bc2dinstanceN = instanceN
      # Returned
        runid = attrD['UniqueGlobalKey/tracking_id/run_id'][1]
        #readid = attrD['Raw/Reads/Read_{0}/read_id'.format(readnumberS)][1]
        readid = Attr(attrD, '-1',
            ['Raw/Reads/Read_{0}/read_id'.format(readnumberS),
             'Analyses/EventDetection_{0}/Reads/Read_{1}/read_id'.format(instanceN, readnumberS)])
        bc1dinstanceN = instanceN
        #numevents = int(attrD['Analyses/Basecall_1D_{0}/Summary/basecall_1d_complement/num_events'.format(bc1dinstanceN)][1])
        numevents = int(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/Summary/basecall_1d_complement/num_events'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/Summary/basecall_1d_complement/num_events'.format(bc2dinstanceN)]))
        numskips = int(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/Summary/basecall_1d_complement/num_skips'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/Summary/basecall_1d_complement/num_skips'.format(bc2dinstanceN)]))
        numstays = int(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/Summary/basecall_1d_complement/num_stays'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/Summary/basecall_1d_complement/num_stays'.format(bc2dinstanceN)]))
        numcalled = int(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/Summary/basecall_1d_complement/called_events'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/Summary/basecall_1d_complement/called_events'.format(bc2dinstanceN)]))
        strandstarttimesec = float(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/BaseCalled_complement/Events/start_time'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/BaseCalled_complement/Events/start_time'.format(bc2dinstanceN)]))
        stranddurationsec = float(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/BaseCalled_complement/Events/duration'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/BaseCalled_complement/Events/duration'.format(bc2dinstanceN)]))
        strandendtimesec = (strandstarttimesec + stranddurationsec) if strandstarttimesec != -1 and stranddurationsec != -1 else -1 
        #strandstarttimesec = exp_start_time + strandstarttime / samplingrate
        #strandendtimesec = exp_start_time + strandstarttime / samplingrate + strandduration / samplingrate
        #strandstarttimesec = strandstarttime / samplingrate
        #strandendtimesec = strandstarttimesec + strandduration / samplingrate
        #stranddurationsec = strandduration / samplingrate
        if strandstarttimesec != -1:
            strandstarttimeiso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(exp_start_time + strandstarttimesec))
            strandendtimeiso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(exp_start_time + strandendtimesec))
        else:
            strandstarttimeiso = 'NA'
            strandendtimeiso = 'NA'
        meanqscore = float(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/Summary/basecall_1d_complement/mean_qscore'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/Summary/basecall_1d_complement/mean_qscore'.format(bc1dinstanceN)]))
        strandscore = float(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/Summary/basecall_1d_complement/strand_score'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/Summary/basecall_1d_complement/strand_score'.format(bc1dinstanceN)]))
        if attrD.has_key('Analyses/Basecall_1D_{0}/Summary/basecall_1d_complement/sequence_length'.format(bc1dinstanceN)):
            seqlen = int(attrD['Analyses/Basecall_1D_{0}/Summary/basecall_1d_complement/sequence_length'.format(bc1dinstanceN)][1])
        else:
            seqlen = len(fastqD[readtype][1]) if fastqD[readtype] is not None else -1
        bqlen = len(fastqD[readtype][3]) if fastqD[readtype] is not None else -1
        bqmean = np.mean(bqnumA) if bqnumA is not None else -1
        bqmedian = np.median(bqnumA) if bqnumA is not None else -1
        gcpct = round(GC(fastqD[readtype][1]), 1) if fastqD[readtype] is not None else -1
        basespersecond = (seqlen / stranddurationsec) if stranddurationsec != -1 else -1
        comment = ''
        dbuserid = ''
    except:
        Print_read1dblank(P, exptid, batchid, fpD)
        return None
#    rowNP = np.array(
#        (exptid, batchid, runid, readid, bc1dinstanceN,
#        readtype, returnstatus, numevents, numskips, numstays,
#        numcalled, strandstarttimesec, strandendtimesec, stranddurationsec, strandstarttimeiso,
#        strandendtimeiso, meanqscore, strandscore, seqlen, bqlen,
#        bqmean, bqmedian, gcpct, basespersecond, comment,
#        dbuserid),
#        dtype=P.ontread1dstatsH)
    rowNP = np.array(
        (exptid, batchid, runid, readid, bc1dinstanceN,
        readtype, readclass, returnstatus, numevents, numskips,
        numstays, numcalled, strandstarttimesec, strandendtimesec, stranddurationsec,
        strandstarttimeiso, strandendtimeiso, meanqscore, strandscore, seqlen,
        bqlen, bqmean, bqmedian, gcpct, basespersecond,
        comment, dbuserid),
        dtype=P.ontread1dstatsH)
    rowL = list(np.atleast_1d(rowNP).tolist()[0])
    rowL = [x if x != -1 and x != '-1' and x != '-1.0' else 'NA' for x in rowL]
    fpD['read1dstats'].write('{0}\n'.format('\t'.join([str(x) for x in rowL])))
    fpD['read1dstats'].flush()
    return rowL

def Print_ontread1dstats(P, exptid, batchid, readclass, instanceN, attrD, fastqD, fpD):
    returnstatus = Attr(attrD, 'NR',
        ['Analyses/Basecall_1D_{0}/Summary/return_status'.format(instanceN)])
    read1tstats = Print_ontread1tstats(P, exptid, batchid, readclass, instanceN, '1T', attrD, fastqD, returnstatus, fpD)
    read1cstats = Print_ontread1cstats(P, exptid, batchid, readclass, instanceN, '1C', attrD, fastqD, returnstatus, fpD)
    return read1tstats, read1cstats

def Print_read2dblank(P, exptid, batchid, fpD):
    rowNP = np.array(
        (exptid, batchid, 'NA', 'NA', 'NA',
        'NA', -1.0, -1.0, -1.0, '-1',
        '-1', -1.0, -1, -1, -1.0,
        -1.0, -1.0, -1.0, '', ''),
        dtype=P.ontread2dstatsH)
    rowL = list(np.atleast_1d(rowNP).tolist()[0])
    rowL = [x if x != -1 and x != '-1' and x != '-1.0' else 'NA' for x in rowL]
    fpD['read2dstats'].write('{0}\n'.format('\t'.join([str(x) for x in rowL])))
    fpD['read2dstats'].flush()
    return 0

def Print_ontread2dstats(P, exptid, batchid, readclass, instanceN, attrD, fastqD, fpD, read1tstats, read1cstats):
  # Intermediate
    readtype = '2D'
    #readnumberS = attrD['Analyses/Hairpin_Split_{0}/Configuration/general/read_id'.format(instanceN)][1]
    readnumberS = Attr(attrD, '-1',
        ['Analyses/Hairpin_Split_{0}/Configuration/general/read_id'.format(instanceN),
         'Analyses/Basecall_2D_{0}/Configuration/general/read_id'.format(instanceN),
         'Analyses/Segment_Linear_{0}/Configuration/general/read_id'.format(instanceN)])
    if readnumberS == '-1':
        Print_read2dblank(P, exptid, batchid, fpD)
        return None
    try:
        bc1dinstanceN = instanceN
        samplingrate = float(attrD['UniqueGlobalKey/channel_id/sampling_rate'][1])
        exp_start_time = float(attrD['UniqueGlobalKey/tracking_id/exp_start_time'][1])
        tempduration = float(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/BaseCalled_template/Events/duration'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/BaseCalled_template/Events/duration'.format(bc1dinstanceN)]))
        compduration = float(Attr(attrD, '-1',
            ['Analyses/Basecall_1D_{0}/BaseCalled_complement/Events/duration'.format(bc1dinstanceN),
             'Analyses/Basecall_2D_{0}/BaseCalled_complement/Events/duration'.format(bc1dinstanceN)]))
        #tempdurationsec = tempduration
        #compdurationsec = compduration
        #meandurationsec = ((tempdurationsec + compdurationsec) / 2.0) if tempdurationsec != -1 and compdurationsec != -1 else -1
        if read1tstats is not None and read1cstats is not None:
            tempdurationsec = read1tstats[13]
            compdurationsec = read1cstats[13]
            meandurationsec = ((tempdurationsec + compdurationsec) / 2.0) \
                if (tempdurationsec != -1 and compdurationsec != -1 and tempdurationsec != 'NA' and compdurationsec != 'NA') else -1
        else:
            meandurationsec = -1
        bqnumA = np.array([ord(x)-33 for x in fastqD[readtype][3]]) if fastqD[readtype] is not None else None
      # Returned
        runid = attrD['UniqueGlobalKey/tracking_id/run_id'][1]
        #readid = attrD['Raw/Reads/Read_{0}/read_id'.format(readnumberS)][1]
        readid = Attr(attrD, '-1',
            ['Raw/Reads/Read_{0}/read_id'.format(readnumberS),
             'Analyses/EventDetection_{0}/Reads/Read_{1}/read_id'.format(instanceN, readnumberS)])
        bc2dinstanceN = instanceN
        returnstatus = Attr(attrD, 'NR',
            ['Analyses/Basecall_2D_{0}/Summary/return_status'.format(bc2dinstanceN)])
        #starttimesec = -1
        #endtimesec = -1
        #durationsec = -1
        #starttimeiso = 'NA'
        #endtimeiso = 'NA'
        if read1tstats is not None and read1tstats[12] != -1 and read1tstats[12] != 'NA':
            starttimesec = read1tstats[12]
        else:
            starttimesec = -1
        if read1cstats is not None and read1cstats[13] != -1 and read1cstats[13] != 'NA':
            endtimesec = read1cstats[13]
        elif read1tstats is not None and read1tstats[13] != -1 and read1tstats[13] != 'NA':
            endtimesec = read1tstats[13]
        else:
            endtimesec = -1
        durationsec = (endtimesec - starttimesec) \
            if (starttimesec != -1 and endtimesec != -1 and starttimesec != 'NA' and endtimesec != 'NA') else -1
        if durationsec < 0 and durationsec != -1.0:
            pass
        if starttimesec != -1 and starttimesec != '-1':
            starttimeiso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(exp_start_time + starttimesec))
            endtimeiso = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(exp_start_time + endtimesec))
        else:
            starttimeiso = 'NA'
            endtimeiso = 'NA'
        try:
            meanqscore = float(attrD['Analyses/Basecall_2D_{0}/Summary/basecall_2d/mean_qscore'.format(bc2dinstanceN)][1])
        except:
            meanqscore = -1
        try:
            seqlen = int(attrD['Analyses/Basecall_2D_{0}/Summary/basecall_2d/sequence_length'.format(bc2dinstanceN)][1])
        except:
            seqlen = len(fastqD[readtype][1]) if fastqD[readtype] is not None else -1
        bqlen = len(fastqD[readtype][3]) if fastqD[readtype] is not None else -1
        bqmean = np.mean(bqnumA) if bqnumA is not None else -1
        bqmedian = np.median(bqnumA) if bqnumA is not None else -1
        gcpct = round(GC(fastqD[readtype][1]), 1) if fastqD[readtype] is not None else -1
        basespersecond = (seqlen / meandurationsec) if meandurationsec != -1 else -1
        comment = ''
        dbuserid = ''
    except:
        Print_read2dblank(P, exptid, batchid, fpD)
        return 1
    rowNP = np.array(
        (exptid, batchid, runid, readid, bc2dinstanceN,
        returnstatus, readclass, starttimesec, endtimesec, durationsec,
        starttimeiso, endtimeiso, meanqscore, seqlen, bqlen,
        bqmean, bqmedian, gcpct, basespersecond, comment,
        dbuserid),
        dtype=P.ontread2dstatsH)
    rowL = list(np.atleast_1d(rowNP).tolist()[0])
    rowL = [x if x != -1 and x != '-1' and x != '-1.0' else 'NA' for x in rowL]
    fpD['read2dstats'].write('{0}\n'.format('\t'.join([str(x) for x in rowL])))
    fpD['read2dstats'].flush()
    return rowL

def Extract_Stats(filecnt, P, exptid, batchid, readclass, instanceN, attrD, fastqD, fpD, fast5file):
    'Print records from current fast5 file to ont[expt|read]stats, ontread[event|1d|2d]stats files.'
    if not len(attrD):
        return 0
    if filecnt == 1:
        Print_ontexptstats(P, exptid, batchid, instanceN, attrD, fpD)
    Print_ontreadstats(P, exptid, batchid, readclass, instanceN, attrD, fpD, fast5file)
    Print_ontreadeventstats(P, exptid, batchid, readclass, instanceN, attrD, fpD)
    read1tstats, read1cstats = Print_ontread1dstats(P, exptid, batchid, readclass, instanceN, attrD, fastqD, fpD)
    read2dstats = Print_ontread2dstats(P, exptid, batchid, readclass, instanceN, attrD, fastqD, fpD, read1tstats, read1cstats)
    return 0

def Extract_Fast5_Data(filecnt, args, P, mylogger, exptid, fast5path, readclass, fpD, constD, instanceN):
    'Open the FAST5 file, extract all requested information, write it to the file pointer.'
    #mylogger.debug('Extract_Fast5_Data : Processing fast5 from experiment {0} readclass {1} {2}\n'.format(exptid, readclass, fast5path))
    attrD, runnumberD, readnumberD, fastqD = P.fast5_extract(fast5path, instanceN, True, True, args.fastq, True, args.fastqheaderformat)
    if not len(attrD.keys()):
        return None
    filteredattrD, filterok = P.fast5_attributes_filter(attrD, instanceN)
    batchid = Get_Batchid(fast5path)
    if args.pairs:
        Extract_Pairs(P, constD, exptid, batchid, instanceN, attrD, fpD)
    if args.fastq:
        Extract_Fastq(mylogger, readclass, fastqD, fpD)
    if args.stats:
        Extract_Stats(filecnt, P, exptid, batchid, readclass, instanceN, attrD, fastqD, fpD, fast5path)
    return batchid

def Extract_Expt_Data(args, P, mylogger, myhandler, processname, exptid, exptdir, exptinstanceN, constD, fp):
    'Iterate through each FAST5 file for this experiment, save metadata to files.'
    mylogger.info('Processing experiment {0}'.format(exptid))
  # Print headers, if required
    fp['batch'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontbatch'))))
    if args.pairs:
        fp['exptpairs'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontexptpairs'))))
        fp['readpairs'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontreadpairs'))))
    if args.stats:
        fp['exptstats'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontexptstats'))))
        fp['readstats'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontreadstats'))))
        fp['readeventstats'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontreadeventstats'))))
        fp['read1dstats'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontread1dstats'))))
        fp['read2dstats'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontread2dstats'))))
  # Start processing
    #mylogger.debug('Extract_Expt_Data : Processing fast5 from experiment {0}\n'.format(exptid))
    passdir = os.path.join(exptdir, 'reads', 'downloads', 'pass')
    faildir = os.path.join(exptdir, 'reads', 'downloads', 'fail')
    fast5L = []
    fast5L += [(passdir, x, 'pass') for x in os.listdir(passdir) if x.endswith('.fast5')]
    fast5L += [(faildir, x, 'fail') for x in os.listdir(faildir) if x.endswith('.fast5')]
    maxfiles = min(args.samplesize, len(fast5L))
    fcnt = 0
    batchD = {}
    for fast5dir, fast5, readclass in fast5L:
        fcnt += 1
        if ((fcnt % P.processmessageinterval) == 0):
            mylogger.info('Processing {0}-th file of experiment {1}'.format(fcnt, exptid))
        batchid = Extract_Fast5_Data(fcnt, args, P, mylogger, exptid, os.path.join(fast5dir, fast5), readclass, fp, constD, exptinstanceN)
        if batchid is not None and not batchD.has_key(batchid):
            batchD[batchid] = [exptid, batchid, '', exptinstanceN]
            fp['batch'].write('{0}\n'.format('\t'.join(batchD[batchid])))
            fp['batch'].flush()
        if fcnt == maxfiles:
            break
    return 0

def Files_List(outdir, exptid):
    'Return list of output files expected for this program.'
    outpathD = {
      # -fastq
        'fq1Tpass' : os.path.join(outdir, '{exptid}_1T_pass.fastq'.format(exptid=exptid)),
        'fq1Tfail' : os.path.join(outdir, '{exptid}_1T_fail.fastq'.format(exptid=exptid)),
        'fq1Cpass' : os.path.join(outdir, '{exptid}_1C_pass.fastq'.format(exptid=exptid)),
        'fq1Cfail' : os.path.join(outdir, '{exptid}_1C_fail.fastq'.format(exptid=exptid)),
        'fq2Dpass' : os.path.join(outdir, '{exptid}_2D_pass.fastq'.format(exptid=exptid)),
        'fq2Dfail' : os.path.join(outdir, '{exptid}_2D_fail.fastq'.format(exptid=exptid)),
      # -pairs
        'batch': os.path.join(outdir, '{exptid}_batch.txt'.format(exptid=exptid)),
        'exptpairs' : os.path.join(outdir, '{exptid}_exptpairs.txt'.format(exptid=exptid)),
        'readpairs' : os.path.join(outdir, '{exptid}_readpairs.txt'.format(exptid=exptid)),
      # -stats
        'exptstats' : os.path.join(outdir, '{exptid}_exptstats.txt'.format(exptid=exptid)),
        'readstats' : os.path.join(outdir, '{exptid}_readstats.txt'.format(exptid=exptid)),
        'readeventstats' : os.path.join(outdir, '{exptid}_readeventstats.txt'.format(exptid=exptid)),
        'read1dstats' : os.path.join(outdir, '{exptid}_read1dstats.txt'.format(exptid=exptid)),
        'read2dstats' : os.path.join(outdir, '{exptid}_read2dstats.txt'.format(exptid=exptid))
    }
    return outpathD

def Files_Open(outpathD, P, outdir, dofastq, dopairs, dostats, exptid, mylogger):
    'Open all the output files requested at the command-line.'
    keyL = []
    keyL += ['batch']
    if dofastq:
        keyL += ['fq1Tpass', 'fq1Tfail', 'fq1Cpass', 'fq1Cfail', 'fq2Dpass', 'fq2Dfail']
    if dopairs:
        keyL += ['exptpairs', 'readpairs']
    if dostats:
        keyL += ['exptstats', 'readstats', 'readeventstats', 'read1dstats', 'read2dstats']
    failed = False
    fp = {}
    for key in keyL:
        try:
            fp[key] = open(outpathD[key], 'w')
        except:
            mylogger.error('Failed to open file for writing ({0})'.format(outpathD[key]))
            failed = True
    if failed:
        for key in fp.keys():
            try:
                fp[key].close()
            except:
                pass
        sys.exit(P.err_code('ErrorFileOpen'))
    return fp

def Files_Close(fp):
    for key in fp.keys():
        try:
            fp[key].close()
        except:
            pass
    return 0

def Prerequisites(args, P, mylogger, myhandler, processname):
    'Exit program if some prerequisites are not met.'
    exptconstants_path = args.exptconstants
    if not os.path.exists(exptconstants_path):
        mylogger.error('Input file missing ({0}). Please run \'marcoporo.py exptconstants\' now'.format(exptconstants_path))
        sys.exit(P.err_code('ErrorFileMissing'))
    if not os.path.exists(args.extractdir):
        os.makedirs(args.extractdir)
    return 0

def Files_NeedGenerating(fileL, overwrite):
    'Return True if .'
    if overwrite:
        return True
    allok = True
    for outpath in fileL:
        if os.path.exists(outpath) and os.path.getsize(outpath) > 0:
            mylogger.error('Output file already exists and overwrite is false ({0})'.format(outpath))
            allok = False
    if not allok:
        return False
    return True

def Process(args, P, mylogger, myhandler, processname):
    'Read file of constant field names produced by marcoporo exptconstants.'
    exptconstants_path = args.exptconstants
    constL = open(args.exptconstantfields, 'r').read().strip().split('\n')
    constD = dict(itertools.izip(constL, len(constL)*[None]))
    exptidL, E = P.expt_read(args.experiments)
    if args.fastq or args.pairs or args.stats:

        for exptid in exptidL:

            outpathD = Files_List(args.extractdir, exptid)
            proceed = Files_NeedGenerating(outpathD, args.overwrite)
            if proceed:
                fp = Files_Open(outpathD, P, args.extractdir, args.fastq, args.pairs, args.stats, exptid, mylogger)
                Extract_Expt_Data(args, P, mylogger, myhandler, processname, exptid, E[exptid]['dirpath'], E[exptid]['instanceN'], constD, fp)
                Files_Close(fp)
        else:
            mylogger.info('Not processing experiment {0} - no data requested'.format(exptid))
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    'Execute this sub-tool.'
    mylogger.info('Started')
    args.fastq = P.str_2bool(args.fastq)
    args.pairs = P.str_2bool(args.pairs)
    args.stats = P.str_2bool(args.stats)
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname)
    mylogger.info('Finished')
    return 0
