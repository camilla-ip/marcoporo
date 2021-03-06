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
_processname = 'exptconstants'
_fast5samplesize = 3
_ontbatchH = ['exptid', 'batchid', 'batchds', 'bestnnn']
_ontexptpairH = ['exptid', 'batchid', 'instanceN', 'var', 'val']
_ontreadpairH = ['exptid', 'batchid', 'instanceN', 'var', 'val']

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

#def Filter_Attributes(attrD, bestN):
#    'Return copy of attrD minus elements of non-bestN attributes.'
#    attrDbestN = {}
#    for attr in attrD.keys():
#        K = attr.split('/')
#        NNN = K[1].split('_')[-1] if len(K) >= 2 else ''
#        matches = (len(NNN)!=3) or (NNN.isdigit() and NNN == bestN)
#        if matches:
#            L = K[1].split('_')
#            L[-1] = 'NNN'
#            K[1] = '_'.join(L)
#            attrcleaned = '/'.join(K)
#            attrDbestN[attrcleaned] = attrD[attr]
#        else:
#            attrDbestN[attr] = attrD[attr]
#    return attrDbestN

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
             attributeD, runnumberD, readnumberD, fastqD = P.fast5_extract(fast5, E[exptid]['instanceN'], True, True, False, True, 'concise')
             #attrDbestN = Filter_Attributes(attributeD)
             attrDbestN, filterok = P.fast5_attributes_filter(attributeD, E[exptid]['instanceN'])
             if not filterok:
                 mylogger.error('Failed to filter instanceN in attributes, skipping file ({0})'.format(fast5))
                 continue
             M[fast5] = attrDbestN
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
    outpath = os.path.join(args.outdir, P.file_exptconstants())
    with open(outpath, 'w') as out_fp:
        out_fp.write('{0}\n'.format('\t'.join(['Attribute'] + exptidL)))
        for var in varL:
            out_fp.write('{0}\n'.format('\t'.join([var] + [str(x) for x in resultD[var]])))

  # Save fields, with each _??? replaced by _NNN, that are constant across the run
    outpath = os.path.join(args.outdir, P.file_exptconstantfields())
    with open(outpath, 'w') as out_fp:
        for var in varL:
            out_fp.write('{0}\n'.format(var))

  # Reformat varL as a dictionary for faster lookup later
    S = dict(itertools.izip(varL, len(varL)*[0]))

    mylogger.info('Finished computing table of fields shared across each experiment (expt_share_[values|fields].txt).')
    return S

def ContinueProcessing(args, P, mylogger):
    'Return without doing anything if args.overwrite is True and output files already exist and do not have zero length.'

    L = [P.file_exptconstants(), P.file_exptconstantfields()]

    if args.overwrite:
        for outfile in L:
            outpath = os.path.join(args.outdir, outfile)
            if os.path.exists(outpath) and os.path.getsize(outpath) > 0:
                mylogger.warning('Overwriting existing non-empty file ({0})'.format(outpath))
        return True
    else:
        allok = True
        for outfile in L:
            outpath = os.path.join(args.outdir, outfile)
            if os.path.exists(outpath) and os.path.getsize(outpath) > 0:
                mylogger.info('Skipping as output file already exists ({0})'.format(outpath))
                allok = False
        if not allok:
            return False
            #sys.exit(P.err_code('ErrorOutfileExists'))
    return True

def Process(args, P, mylogger, myhandler, processname):
    exptidL, E = P.expt_read(args.experiments)
    S = Compute_Share(args, P, mylogger, myhandler, processname, exptidL, E)
    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    global _mylogger
    _mylogger = mylogger
    args.overwrite = P.str_2bool(args.overwrite)
    mylogger.info('{0}: Started'.format(_processname))
    if ContinueProcessing(args, P, mylogger):
        Process(args, P, mylogger, myhandler, _processname)
    mylogger.info('{0}: Finished'.format(_processname))
    return 0
