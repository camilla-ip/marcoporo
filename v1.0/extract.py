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

def Merge_DictPairs(D):
    'Given D[name]={var:val, ...} return elt[var]=val with same pairs in all name keys.'
    elt = {}
    nameL = D.keys()
    namesz = len(nameL)
    if not len(nameL):
        return elt
    #varsetL = [set(D[name].keys()) for name in D.keys()]
    varsetL = [set(D[name].keys()) for name in D.keys()]
    sharedvarL = set.intersection(*varsetL)
    for var in sharedvarL:
        foundL = [D[name][var][1] for name in nameL if D.has_key(name) and D[name].has_key(var)]
        if len(foundL) == namesz and len(set(foundL)) == 1:
            elt[var] = foundL[0]
    return elt

def Extract_Attributes(fast5_path):
    'Return the attributes in the fast5 file as a dict of {var:value, ...} pairs.'
    attribute = {}
    hdf = h5py.File(fast5_path, 'r')
    list_of_names = []  # Names of all groups and subgroups in the file
    hdf.visit(list_of_names.append)
    for name in list_of_names:
        itemL = hdf[name].attrs.items() # attribute name and value pairs
        for item in itemL:
            attr, val = item
            if type(hdf[name].attrs[attr]) == np.ndarray:
                val = ''.join(hdf[name].attrs[attr])
            val = str(val).replace('\n', '')
            attribute[name+'/'+attr] = val
    hdf.close()
    return attribute

def Compute_Share(args, P, mylogger, myhandler, processname, exptidL, E):
    'Compute table of softed experiment-level attributes for each experiment.'

    mylogger.info('Start computing table of fields shared across each experiment.')

  # For each experiment
    C = {}
    #exptidL = [elt['exptid'] for elt in E]
    for exptid in exptidL:
      # Pick N=10 random reads from the 'pass' folder of each experiment
        indir = os.path.join(E[exptid]['dirpath'], 'reads', 'downloads', 'pass')
        infileL = os.listdir(indir)
        N = min(len(infileL), args.samplesize)
        fast5L = [os.path.join(indir, x) for x in random.sample(infileL, N)]
      # Extract all metadata
        M = {}
        for fast5 in fast5L:
            #M[fast5] = Extract_Attributes(fast5)
             attributeD, runnumberD, readnumberD = P.fast5_attributes(fast5, 'all')
             M[fast5] = attributeD
      # Keep the fields that are constant in this experiment C[exptid] = {var:val, ...}
        C[exptid] = Merge_DictPairs(M)

  # Collate a table of fields vs exptid of all common fields, with
  # fields sorted alphabetically and exptid in the same order as in experiments.txt file
    resultD = {}
    varsetL = [set(C[exptid].keys()) for exptid in C.keys()]
    uniquevarL = list(set.union(*varsetL))
    uniquevarL.sort()
    #resultL.append(['Parameter'] + exptidL)
    for var in uniquevarL:
        #row = [var] + [C[exptid][var] if C[exptid].has_key(var) else '' for exptid in exptidL]
        resultD[var] = [C[exptid][var] if C[exptid].has_key(var) else '' for exptid in exptidL]
  # Create outdir
    if not os.path.exists(args.outdir):
        try:
            os.makedirs(args.outdir)
        except:
            mylogger.error('Failed to create outdir ({0})'.format(args.outdir))
            sys.exit(P.err_code('ErrorDirCreate'))
  # Save table of shared fields and values, and field names one per line
    varL = resultD.keys()
    varL.sort()
    outpath = os.path.join(args.outdir, 'expt_sharedattributes.txt')
    with open(outpath, 'w') as out_fp:
        out_fp.write('{0}\n'.format('\t'.join(['Attribute'] + exptidL)))
        for var in varL:
            out_fp.write('{0}\n'.format('\t'.join([var] + [str(x) for x in resultD[var]])))

  # Save fields that are constant across the run
    outpath = os.path.join(args.outdir, 'expt_sharedattrnames.txt')
    with open(outpath, 'w') as out_fp:
        for var in varL:
            out_fp.write('{0}\n'.format(var))

  # Reformat varL as a dictionary for faster lookup later
    S = dict(itertools.izip(varL, len(varL)*[0]))

    mylogger.info('Finished computing table of fields shared across each experiment (expt_share_[values|fields].txt).')
    return S

def Extract_Fast5_Data(args, P, mylogger, fast5path, readclass, fpD):
    'Open the FAST5 file, extract all requested information, write it to the file pointer.'

    attrD = P.fast5_attributes(fast5path, 'all')
    if args.pairs:
        pass
    pass

def Extract_Expt_Data(args, P, mylogger, myhandler, processname, exptid, exptdir):
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

  # Iterate across all FAST5 files
    for fast5 in passL:
        Extract_Fast5_Data(args, P, mylogger, os.path.join(passdir, fast5), 'pass', fp)
    for fast5 in failL:
        Extract_Fast5_Data(args, P, mylogger, os.path.join(faildir, fast5), 'fail', fp)
        
  # Close file pointers
    for key in fp.keys():
        try:
            fp[key].close()
        except:
            pass

    mylogger.info('Finished processing experiment {0}'.format(exptid))

def Iterate_Fast5(args, P, mylogger, myhandler, processname, exptidL, E):
    'Iterate through each FAST5 file and store the info specified in the command-line flags.'

    mylogger.info('Start iterating across all FAST5 files.')
    #exptidL = [elt['exptid'] for elt in E]
    for exptid in exptidL:
        Extract_Expt_Data(args, P, mylogger, myhandler, processname, exptid, E[exptid]['dirpath'])
    mylogger.info('Finished iterating across all FAST5 files.')
    return 0

def Prerequisites(args, P, mylogger, myhandler, processname):
    

def Process(args, P, mylogger, myhandler, processname):
    exptidL, E = P.expt_read(args.experiments)
    if args.share or args.fastq or args.model or args.pairs or args.stats:
        S = Compute_Share(args, P, mylogger, myhandler, processname, exptidL, E)
    else:
        mylogger.info('{0}: Not computing table of shared attributes and values'.format(processname))
    if args.fastq or args.model or args.pairs or args.stats:
        Iterate_Fast5(args, P, mylogger, myhandler, processname, exptidL, E, S)
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
