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
import sys

import marcoporoversion

_processname = 'extractone'
_fast5samplesize = 3

def Get_Batchid(fast5path):
    '''
    For inferring batchids, can cope with fast5 files like:
    MinION2_20160802_FNFAD22824_MN16454_sequencing_run_Chip93_MARC_R9_1D_UBC_77825_ch100_read249_strand.fast5 -> 77825
    '''
    batchid = None
    file = os.path.basename(fast5path)
    filestem = '.'.join(file.split('.')[0:-1])
    L = filestem.split('_')
    if L[-1].startswith('strand') and L[-2].startswith('read') and L[-3].startswith('ch'):
        batchid = L[-4]
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

def Extract_Stats():
    pass

def Extract_Fast5_Data(args, P, mylogger, exptid, fast5path, readclass, fpD, constD, instanceN):
    'Open the FAST5 file, extract all requested information, write it to the file pointer.'
    mylogger.debug('Extract_Fast5_Data : Processing fast5 from experiment {0} readclass {1}\n'.format(exptid, readclass))
    attrD, runnumberD, readnumberD, fastqD = P.fast5_extract(fast5path, instanceN, args.pairs, True, args.fastq, True, args.fastqheaderformat)
    filteredattrD, filterok = P.fast5_attributes_filter(attrD, instanceN)
    batchid = Get_Batchid(fast5path)
    if args.pairs:
        Extract_Pairs(P, constD, exptid, batchid, instanceN, attrD, fpD)
    if args.fastq:
        Extract_Fastq(mylogger, readclass, fastqD, fpD)
    if args.stats:
        Extract_Stats()
    return batchid

def Extract_Expt_Data(args, P, mylogger, myhandler, processname, exptid, exptdir, exptinstanceN, constD, fp):
    'Iterate through each FAST5 file for this experiment, save metadata to files.'
    mylogger.info('Processing experiment {0}'.format(exptid))
  # Print headers, if required
    if args.pairs:
        fp['batch'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontbatch'))))
        fp['exptpairs'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontexptpairs'))))
        fp['readpairs'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontreadpairs'))))
    if args.stats:
        fp['exptstats'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontexptstats'))))
        fp['readstats'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontreadstats'))))
        fp['readeventstats'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontreadeventstats'))))
        fp['read1dstats'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontread1dstats'))))
        fp['read2dstats'].write('{0}\n'.format('\t'.join(P.fast5_headernames('ontread2dstats'))))
  # Start processing
    mylogger.debug('Extract_Expt_Data : Processing fast5 from experiment {0}\n'.format(exptid))
    passdir = os.path.join(exptdir, 'reads', 'downloads', 'pass')
    faildir = os.path.join(exptdir, 'reads', 'downloads', 'fail')
    fast5L = [(passdir, x, 'pass') for x in os.listdir(passdir) if x.endswith('.fast5')]
    fast5L += [(faildir, x, 'fail') for x in os.listdir(faildir) if x.endswith('.fast5')]
    maxfiles = min(args.samplesize, len(fast5L))
    fcnt = 0
    batchD = {}
    for fast5dir, fast5, readclass in fast5L:
        fcnt += 1
        if fcnt == maxfiles:
            break
        batchid = Extract_Fast5_Data(args, P, mylogger, exptid, os.path.join(fast5dir, fast5), readclass, fp, constD, exptinstanceN)
        if batchid is not None and not batchD.has_key(batchid):
            batchD[batchid] = [exptid, batchid, '', args.instanceN]
            fp['batch'].write('{0}\n'.format('\t'.join(batchD[batchid])))
            fp['batch'].flush()
    return 0

def Files_Open(outdir, dofastq, dopairs, dostats, exptid, mylogger):
    'Open all the output files requested at the command-line.'
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
    keyL = []
    if dofastq:
        keyL += ['fq1Tpass', 'fq1Tfail', 'fq1Cpass', 'fq1Cfail', 'fq2Dpass', 'fq2Dfail']
    if dopairs:
        keyL += ['batch', 'exptpairs', 'readpairs']
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
    exptconstants_path = os.path.join(args.outdir, P.file_exptconstants())
    if not os.path.exists(exptconstants_path):
        mylogger.error('Input file missing ({0}). Please run \'marcoporo.py exptconstants\' now'.format(exptconstants_path))
        sys.exit(P.err_code('ErrorFileMissing'))

def Process(args, P, mylogger, myhandler, processname, exptid):
    'Read file of constant field names produced by marcoporo exptconstants.'
    exptconstants_path = os.path.join(args.outdir, P.file_exptconstants())
    constL = open(os.path.join(args.outdir, P.file_exptconstantfields()), 'r').read().strip().split('\n')
    constD = dict(itertools.izip(constL, len(constL)*[None]))
    if args.fastq or args.model or args.pairs or args.stats:
        fp = Files_Open(args.outdir, args.fastq, args.pairs, args.stats, args.exptid, mylogger)
        Extract_Expt_Data(args, P, mylogger, myhandler, processname, args.exptid, args.indir, args.instanceN, constD, fp)
        Files_Close(fp)
    else:
        mylogger.info('Not processing experiment {0} - no data requested'.format(exptid))
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    'Execute this sub-tool.'
    mylogger.info('Started')
    args.fastq = P.str_2bool(args.fastq)
    args.model = P.str_2bool(args.model)
    args.pairs = P.str_2bool(args.pairs)
    args.stats = P.str_2bool(args.stats)
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname, args.exptid)
    mylogger.info('Finished')
