#!/usr/bin/env python

import h5py
import logging
import numpy as np
import os
import random

import marcoporoversion

_mylogger = None
_processname = 'seqparams'
_fast5samplesize = 3

def Merge_DictPairs(D):
    'Given D[name]={var:val, ...} return elt[var]=val with same pairs in all name keys.'
    elt = {}
    nameL = D.keys()
    namesz = len(nameL)
    if not len(nameL):
        return elt
    varsetL = [set(D[name].keys()) for name in D.keys()]
    sharedvarL = set.intersection(*varsetL)
    for var in sharedvarL:
        foundL = [D[name][var] for name in nameL if D.has_key(name) and D[name].has_key(var)]
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

def Process(args, P, mylogger, myhandler, processname):

  # Read the experiments file
    E = P.expt_read(args.experiments)

  # For each experiment
    C = {}
    exptidL = [elt['exptid'] for elt in E]
    for exptid in exptidL:
      # Pick N=10 random reads from the 'pass' folder of each experiment
        indir = os.path.join([e for e in E if e['exptid'] == exptid][0]['dirpath'],
            'reads', 'downloads', 'pass')
        infileL = os.listdir(indir)
        N = min(len(infileL), args.samplesize)
        fast5L = [os.path.join(indir, x) for x in random.sample(infileL, N)]
      # Extract all metadata
        M = {}
        for fast5 in fast5L:
            M[fast5] = Extract_Attributes(fast5)
      # Keep the fields that are constant in this experiment C[exptid] = {var:val, ...}
        C[exptid] = Merge_DictPairs(M)

  # Collate a table of fields vs exptid of all common fields, with
  # fields sorted alphabetically and exptid in the same order as in experiments.txt file
    resultL = []
    varsetL = [set(C[exptid].keys()) for exptid in C.keys()]
    uniquevarL = list(set.union(*varsetL))
    uniquevarL.sort()
    resultL.append([''] + exptidL)
    for var in uniquevarL:
        row = [var] + [C[exptid][var] if C[exptid].has_key(var) else '' for exptid in exptidL]
        resultL.append(row)
  # Save table as TSV to outdir/seqparams.txt
    if not os.path.exists(args.outdir):
        try:
            os.makedirs(args.outdir)
        except:
            mylogger.error('Failed to create outdir ({0})'.format(args.outdir))
            sys.exit(P.err_code('ErrorDirCreate'))
    outpath = os.path.join(args.outdir, 'seqparams.txt')
    with open(outpath, 'w') as out_fp:
        for row in resultL:
            out_fp.write('{0}\n'.format('\t'.join([str(x) for x in row])))

    return 0

def run(parser, args, P, mylogger, myhandler, argv):
    global _mylogger
    _mylogger = mylogger
    mylogger.info('_processname: Started'.format(_processname))
    Process(args, P, mylogger, myhandler, _processname)
    mylogger.info('_processname: Finished'.format(_processname))
