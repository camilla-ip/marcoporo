#!/usr/bin/env python

import h5py
import itertools
import logging
import numpy as np
import os
import random

import marcoporoversion

_mylogger = None
_processname = 'extract'
_fast5samplesize = 3
_ontbatchH = ['exptid', 'batchid', 'batchds', 'bestnnn']
_ontexptpairH = ['exptid', 'batchid', 'NNN', 'var', 'val']
_ontreadpairH = ['exptid', 'batchid', 'NNN', 'var', 'val']

def Get_NNN(attr):
    '''
    Assumes attributes like Analyses/Basecall_1D_000/Configuration/general/sampling_rate
    where NNN is always the last token in the second part of the attribute path.
    '''
    NNN = 'NA'
    keypartL = attr.split('/')
    if len(keypartL) >= 2: 
        lastpart = keypartL[1].split('_')[-1]
        if lastpart.isdigit():
            NNN = lastpart
    return NNN

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
    attrD,runnumberD, readnumberD  = P.fast5_attributes(fast5path, 'all')
    filteredattrD, filterok = P.fast5_attributes_filter(attrD, basecallN)
    batchid = Get_Batchid(fast5path)
    if args.pairs:
        attrL = attrD.keys()
        attrL.sort()
        for attr in attrL:
            NNN = Get_NNN(attr)
            var = attr
            val = attrD[var][1]
            row = [exptid, batchid, NNN, var, val]
            if constD.has_key(attr):
                fpD['ontexptpair'].write('{0}\n'.format('\t'.join(row)))
            else:
                fpD['ontreadpair'].write('{0}\n'.format('\t'.join(row)))
        pass
    pass

def Extract_Expt_Data(args, P, mylogger, myhandler, processname, exptid, exptdir, exptbasecallN, constD):
    'Iterate through each FAST5 file for this experiment, save metadata to files.'

    mylogger.info('Started processing experiment {0}'.format(exptid))

  # Set up variables
    passdir = os.path.join(exptdir, 'reads', 'downloads', 'pass')
    faildir = os.path.join(exptdir, 'reads', 'downloads', 'fail')
    passL = [x for x in os.listdir(passdir) if x.endswith('.fast5')]
    failL = [x for x in os.listdir(faildir) if x.endswith('.fast5')]
    outpath = {
      # -fastq
        'fq1Tpass' : os.path.join(args.outdir, '{exptid}_1T_pass.fastq'.format(exptid=exptid)),
        'fq1Tfail' : os.path.join(args.outdir, '{exptid}_1T_fail.fastq'.format(exptid=exptid)),
        'fq1Cpass' : os.path.join(args.outdir, '{exptid}_1C_pass.fastq'.format(exptid=exptid)),
        'fq1Cfail' : os.path.join(args.outdir, '{exptid}_1C_fail.fastq'.format(exptid=exptid)),
        'fq2Dpass' : os.path.join(args.outdir, '{exptid}_2D_pass.fastq'.format(exptid=exptid)),
        'fq2Dfail' : os.path.join(args.outdir, '{exptid}_2D_fail.fastq'.format(exptid=exptid)),
      # _model
        'model' : os.path.join(args.outdir, '{exptid}_model.txt'.format(exptid=exptid)),
      # -pairs
        'ontbatch': os.path.join(args.outdir, 'ontbatch.txt'),
        'ontexptpair' : os.path.join(args.outdir, 'ontexptpair.txt'),
        'ontreadpair' : os.path.join(args.outdir, 'ontreadpair.txt'),
      # -stats
        'ontexptstats' : os.path.join(args.outdir, 'ontexptstat.txt'),
        'ontreadstats' : os.path.join(args.outdir, 'ontreadstat.txt')
    }

  # Open file pointers
    fp = {}
    keyL = outpath.keys()
    keyL.sort()
    failed = False
    for key in keyL:
        try:
            fp[key] = open(outpath[key], 'w')
        except:
            mylogger.error('Failed to open file for writing ({0})'.format(outpath[key]))
            failed = True
    if failed:
        for key in fp.keys():
            try:
                fp[key].close()
            except:
                pass
        sys.exit(P.err_code('ErrorFileOpen'))

  # Write the header rows
    if args.pairs:
        fp['ontexptpair'].write('{0}\n'.format('\t'.join(_ontexptpairH)))
        fp['ontreadpair'].write('{0}\n'.format('\t'.join(_ontreadpairH)))

  # Iterate across all FAST5 files
    maxfiles = min(10, len(passL))
    fcnt = 0
    for fast5 in passL:
        fcnt += 1
        if fcnt == maxfiles:
            break
        Extract_Fast5_Data(args, P, mylogger, exptid, os.path.join(passdir, fast5), 'pass', fp, constD, exptbasecallN)
    maxfiles = min(10, len(failL))
    fcnt = 0
    for fast5 in failL:
        fcnt += 1
        if fcnt == maxfiles:
            break
        Extract_Fast5_Data(args, P, mylogger, exptid, os.path.join(faildir, fast5), 'fail', fp, constD, exptbasecallN)
        
  # Close file pointers
    for key in fp.keys():
        try:
            fp[key].close()
        except:
            pass

    mylogger.info('Finished processing experiment {0}'.format(exptid))

def Iterate_Fast5(args, P, mylogger, myhandler, processname, exptidL, E, constD):
    'Iterate through each FAST5 file and store the info specified in the command-line flags.'
    mylogger.info('Start iterating across all FAST5 files.')
    for exptid in exptidL:
        Extract_Expt_Data(args, P, mylogger, myhandler, processname, exptid, E[exptid]['dirpath'], E[exptid]['basecallN'], constD)
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
    constL = open(os.path.join(args.outdir, P.file_exptconstants()), 'r').read().strip().split('\n')
    constD = dict(itertools.izip(constL, len(constL)*[None]))
    if args.fastq or args.model or args.pairs or args.stats:
        Iterate_Fast5(args, P, mylogger, myhandler, processname, exptidL, E, constD)
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
