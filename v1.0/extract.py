#!/usr/bin/env python

import h5py
import itertools
import logging
import numpy as np
import os
import random
import sys

import marcoporoversion

_mylogger = None
_processname = 'extract'
_fast5samplesize = 3
_ontbatchH = ['exptid', 'batchid', 'batchds', 'bestnnn']
_ontexptpairH = ['exptid', 'batchid', 'basecallN', 'var', 'val']
_ontreadpairH = ['exptid', 'batchid', 'basecallN', 'var', 'val']

def Get_BasecallN(attr):
    '''
    Assumes attributes like Analyses/Basecall_1D_000/Configuration/general/sampling_rate
    where basecallN is always the last token in the second part of the attribute path.
    '''
    basecallN = 'NA'
    keypartL = attr.split('/')
    if len(keypartL) >= 2: 
        lastpart = keypartL[1].split('_')[-1]
        if lastpart.isdigit():
            basecallN = lastpart
    return basecallN

def Get_Batchid(fast5path):
  # For inferring batchids, can cope with fast5 files like:
  # MinION2_20160802_FNFAD22824_MN16454_sequencing_run_Chip93_MARC_R9_1D_UBC_77825_ch100_read249_strand.fast5
    batchid = ''
    file = os.path.basename(fast5path)
    filestem = '.'.join(file.split('.')[0:-1])
    L = filestem.split('_')
    if L[-1] == 'strand' and L[-2].startswith('read') and L[-3].startswith('ch'):
        batchid = L[-4]
    return batchid

def Extract_Fast5_Data(args, P, mylogger, exptid, fast5path, readclass, fpD, constD, basecallN):
    'Open the FAST5 file, extract all requested information, write it to the file pointer.'
    mylogger.debug('Extract_Fast5_Data : Processing fast5 from experiment {0} readclass {1}\n'.format(exptid, readclass))
    attrD,runnumberD, readnumberD  = P.fast5_attributes(fast5path, 'all')
    filteredattrD, filterok = P.fast5_attributes_filter(attrD, basecallN)
    batchid = Get_Batchid(fast5path)
    if args.pairs:
        attrL = attrD.keys()
        attrL.sort()
        for attr in attrL:
            #basecallN = Get_BasecallN(attr)
            var = attr
            val = attrD[var][1]
            newvar = P.fast5_attribute_to_NNN(attr)
            if constD.has_key(newvar):
                row = [exptid, batchid, basecallN, newvar, val]
                fpD['ontexptpair'].write('{0}\n'.format('\t'.join(row)))
            else:
                try:
                    readid = attrD['Analyses/EventDetection_{0}/Configuration/general/uuid'.format(basecallN)][1]
                except:
                    readid = 'NK'	# Need to change this to something based on the filename (chN and readN)
                row = [exptid, batchid, readid, basecallN, newvar, val]
                fpD['ontreadpair'].write('{0}\n'.format('\t'.join(row)))
    fpD['ontexptpair'].flush()
    fpD['ontreadpair'].flush()
    return 0

def Extract_Expt_Data(args, P, mylogger, myhandler, processname, exptid, exptdir, exptbasecallN, constD, fp):
    'Iterate through each FAST5 file for this experiment, save metadata to files.'

    mylogger.info('Processing experiment {0}'.format(exptid))

  # Write the header rows
    if args.pairs:
        fp['ontexptpair'].write('{0}\n'.format('\t'.join(_ontexptpairH)))
        fp['ontreadpair'].write('{0}\n'.format('\t'.join(_ontreadpairH)))

  # Iterate across all FAST5 files
    mylogger.debug('Extract_Expt_Data : Processing fast5 from experiment {0}\n'.format(exptid))
    passdir = os.path.join(exptdir, 'reads', 'downloads', 'pass')
    faildir = os.path.join(exptdir, 'reads', 'downloads', 'fail')
    passL = [x for x in os.listdir(passdir) if x.endswith('.fast5')]
    failL = [x for x in os.listdir(faildir) if x.endswith('.fast5')]
    maxfiles = min(args.samplesize, len(passL))
    fcnt = 0
    for fast5 in passL:
        fcnt += 1
        if fcnt == maxfiles:
            break
        Extract_Fast5_Data(args, P, mylogger, exptid, os.path.join(passdir, fast5), 'pass', fp, constD, exptbasecallN)
    maxfiles = min(args.samplesize, len(failL))
    fcnt = 0
    for fast5 in failL:
        fcnt += 1
        if fcnt == maxfiles:
            break
        Extract_Fast5_Data(args, P, mylogger, exptid, os.path.join(faildir, fast5), 'fail', fp, constD, exptbasecallN)
        
    #mylogger.info('Experiment {0} : Processing finished'.format(exptid))

def Files_Open(outpathD):
    fp = {}
    keyL = outpathD.keys()
    keyL.sort()
    failed = False
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

def PairFiles_Open(outdir):
    outpathD = {
        'ontbatch': os.path.join(outdir, 'ontbatch.txt'),
        'ontexptpair' : os.path.join(outdir, 'ontexptpair.txt'),
        'ontreadpair' : os.path.join(outdir, 'ontreadpair.txt')
    }
    fp = Files_Open(outpathD)
    return fp

def FastqFiles_Open(outdir, exptid, passdir, faildir):
    #passdir = os.path.join(exptdir, 'reads', 'downloads', 'pass')
    #faildir = os.path.join(exptdir, 'reads', 'downloads', 'fail')
    outpathD = {
      # -fastq
        'fq1Tpass' : os.path.join(outdir, '{exptid}_1T_pass.fastq'.format(exptid=exptid)),
        'fq1Tfail' : os.path.join(outdir, '{exptid}_1T_fail.fastq'.format(exptid=exptid)),
        'fq1Cpass' : os.path.join(outdir, '{exptid}_1C_pass.fastq'.format(exptid=exptid)),
        'fq1Cfail' : os.path.join(outdir, '{exptid}_1C_fail.fastq'.format(exptid=exptid)),
        'fq2Dpass' : os.path.join(outdir, '{exptid}_2D_pass.fastq'.format(exptid=exptid)),
        'fq2Dfail' : os.path.join(outdir, '{exptid}_2D_fail.fastq'.format(exptid=exptid))
    }
    fp = Files_Open(outpathD)
    return fp

def ModelFiles_Open(outdir):
    outpathD = {
        'model' : os.path.join(args.outdir, '{exptid}_model.txt'.format(exptid=exptid))
    }
    fp = Files_Open(outpathD)
    return fp

def StatsFiles_Open(outdir):
    outpathD = {
        'ontexptstats' : os.path.join(args.outdir, 'ontexptstat.txt'),
        'ontreadstats' : os.path.join(args.outdir, 'ontreadstat.txt')
    }
    fp = Files_Open(outpathD)
    return fp

def Iterate_Fast5(args, P, mylogger, myhandler, processname, exptidL, E, constD, fp):
    'Iterate through each FAST5 file and store the info specified in the command-line flags.'
    mylogger.info('Start iterating across all FAST5 files.')
    for exptid in exptidL:
        Extract_Expt_Data(args, P, mylogger, myhandler, processname, exptid, E[exptid]['dirpath'], E[exptid]['basecallN'], constD, fp)
    mylogger.info('Finished iterating across all FAST5 files.')
    return 0

def Prerequisites(args, P, mylogger, myhandler, processname):
    'Exit program if some prerequisites are not met.'
    exptconstants_path = os.path.join(args.outdir, P.file_exptconstants())
    if not os.path.exists(exptconstants_path):
        mylogger.error('Input file missing ({0}). Please run \'marcoporo.py exptconstants\' now'.format(exptconstants_path))
        sys.exit(P.err_code('ErrorFileMissing'))

def Process(args, P, mylogger, myhandler, processname):
  # Read file of constant field names produced by 'marcoporo exptconstants'
    exptconstants_path = os.path.join(args.outdir, P.file_exptconstants())
    exptidL, E = P.expt_read(args.experiments)
    constL = open(os.path.join(args.outdir, P.file_exptconstantfields()), 'r').read().strip().split('\n')
    constD = dict(itertools.izip(constL, len(constL)*[None]))
    if args.fastq or args.model or args.pairs or args.stats:
        fp = PairFiles_Open(args.outdir)
        Iterate_Fast5(args, P, mylogger, myhandler, processname, exptidL, E, constD, fp)
        Files_Close(fp)
    else:
        mylogger.info('{0}: Not iterating over all FAST5 in all experiments'.format(processname))
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    global _mylogger
    _mylogger = mylogger
    mylogger.info('{0}: Started'.format(_processname))
    args.fastq = P.str_2bool(args.fastq)
    args.model = P.str_2bool(args.model)
    args.pairs = P.str_2bool(args.pairs)
    args.stats = P.str_2bool(args.stats)
    Prerequisites(args, P, mylogger, myhandler, _processname)
    Process(args, P, mylogger, myhandler, _processname)
    mylogger.info('{0}: Finished'.format(_processname))
