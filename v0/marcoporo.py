#!/usr/bin/env python

'''
marcoporo nanopore data comparison package
'''

import logging
import os
import sys

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

_errorargsinvalid = 3
if len(sys.argv) > 1:
    _bin = None
    try:
        for i in range(1, len(sys.argv)):
            if sys.argv[i].startswith('-bin'):
                if '=' in sys.argv[i]:
                    _bin = os.path.realpath(os.path.expandvars('='.join(sys.argv[i].split('=')[1:])))
                    break
                else:
                    _bin = os.path.realpath(os.path.expandvars(sys.argv[i+1]))
                    break
    except:
        _bin = os.path.dirname(os.path.realpath(os.path.expandvars(sys.argv[0])))
    if _bin is None:
        _bin = os.path.dirname(os.path.realpath(os.path.expandvars(sys.argv[0])))
    if _bin and os.path.exists(_bin):
        try:
            sys.path.insert(0, _bin)
            import marcoporolib
        except:
            sys.stderr.write('Error: Failed to import marcoporolib module\n')
            sys.exit(_errorargsinvalid)
    else:
        sys.path.insert(0, '/'.join(os.path.dirname(os.path.realpath(sys.argv[0])).split('/')[:-2] + ['bin']))
        sys.path.insert(0, os.path.dirname(os.path.realpath(sys.argv[0])))
        import marcoporolib

    _profile = None
    try:
        for i in range(1, len(sys.argv)):
            if sys.argv[i].startswith('--profile'):
                if '=' in sys.argv[i]:
                    _profile = os.path.realpath(os.path.expandvars('='.join(sys.argv[i].split('=')[1:])))
                else:
                    _profile = os.path.realpath(os.path.expandvars(sys.argv[i+1]))
    except:
        if os.path.exists(_bin):
            _profile = os.path.join(_bin, 'marcoporo.profile')
        else:
            _profile = 'marcoporo.profile'
    if _profile is None:
        _profile = os.path.join(_bin, 'marcoporo.profile')
    if os.path.exists(_profile):
        with open (_profile, 'r') as in_fp:
            for line in in_fp:
                if not line.startswith('export'):
                    continue
                info = line.strip().split(' ')[1].split('=')
                if len(info) == 2:
                    var, val = info
                    sys.stderr.write('Setting {0}={1}\n'.format(var, val))
                    os.environ[var] = val
else:
    sys.path.insert(0, '/'.join(os.path.dirname(os.path.realpath(sys.argv[0])).split('/')[:-2] + ['bin']))
    sys.path.insert(0, os.path.dirname(os.path.realpath(sys.argv[0])))
    import marcoporolib

import marcoporoversion
_P = marcoporolib.marcoporolib()
import argparse
import random
import time

def run_subtool(parser, args, P, mylogger, myhandler):
    if args.command == 'exptconstants':
        import exptconstants as submodule
    elif args.command == 'extractone':
        import extractone as submodule
    elif args.command == 'extract':
        import extract as submodule
    elif args.command == 'aggregateone':
        import aggregateone as submodule
    elif args.command == 'aggregate':
        import aggregate as submodule
    submodule.run(parser, args, P, mylogger, myhandler, sys.argv)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)

def main():

  # Create the top-level parser

    parser = argparse.ArgumentParser(prog='marcoporo',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '-version', help='Print version', action='version',
        version='marcoporo version ' + str(marcoporoversion.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command',
        parser_class=ArgumentParserWithDefaults)

  # Create the individual tool parsers

    p01 = subparsers.add_parser('exptconstants', help='Tabulate constant experimental attributes',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p01.add_argument('-bin', dest='bin', metavar='DIR', required=False, default='./',
        help='marcoporo scripts dir (specify absolute path)')
    p01.add_argument('-profile', dest='profile', metavar='FILE', required=False, default=None,
        help='marcoporo environment statements (specify absolute path)')
    p01.add_argument('-config', dest='config', metavar='FILE', required=True, default='config.txt',
        help='Analysis configuration file')
    p01.add_argument('-experiments', dest='experiments', metavar='FILE', required=True, default=None,
        help='Experiments and analysis parameters')
    p01.add_argument('-samplesize', dest='samplesize', metavar='INT', type=int, required=False, default=250,
        help='Number of FAST5 files to inspect from each expt to infer constant attributes.')
    p01.add_argument('-outdir', dest='outdir', metavar='DIR', required=True, default=None,
        help='Output directory (specify absolute path)')
    p01.set_defaults(func=run_subtool)

    p02 = subparsers.add_parser('extract', help='Extract required data from all experiments',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p02.add_argument('-bin', dest='bin', metavar='DIR', required=False, default='./',
        help='marcoporo scripts dir (specify absolute path)')
    p02.add_argument('-profile', dest='profile', metavar='FILE', required=False, default=None,
        help='marcoporo environment statements (specify absolute path)')
    p02.add_argument('-config', dest='config', metavar='FILE', required=True, default='config.txt',
        help='Analysis configuration file')
    p02.add_argument('-experiments', dest='experiments', metavar='FILE', required=True, default=None,
        help='Experiments and analysis parameters')
    p02.add_argument('-outdir', dest='outdir', metavar='DIR', required=True, default=None,
        help='Output directory (specify absolute path)')
    p02.add_argument('-fastq', dest='fastq', metavar='BOOL', required=False, default='True',
        help='1D and 2D basecalls')
    #p02.add_argument('-model', dest='model', metavar='BOOL', required=False, default='True',
    #    help='Model parameters used in basecalling')
    p02.add_argument('-pairs', dest='pairs', metavar='BOOL', required=False, default='False',
        help='Name-value pairs for each experiment and read attribute')
    p02.add_argument('-stats', dest='stats', metavar='BOOL', required=False, default='True',
        help='Single-row summary stats for each experiment and read')
    p02.add_argument('-samplesize', dest='samplesize', metavar='INT', type=int, required=False, default=10000000,
        help='Number of FAST5 files to inspect from each expt, useful for testing.')
    p02.set_defaults(func=run_subtool)

    p03 = subparsers.add_parser('extractone', help='Extract data from one experiment',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p03.add_argument('-bin', dest='bin', metavar='DIR', required=False, default='./',
        help='marcoporo scripts dir (specify absolute path)')
    p03.add_argument('-profile', dest='profile', metavar='FILE', required=False, default=None,
        help='marcoporo environment statements (specify absolute path)')
    p03.add_argument('-config', dest='config', metavar='FILE', required=True, default='config.txt',
        help='Analysis configuration file')
    p03.add_argument('-exptid', dest='exptid', metavar='STR', required=True, default=None,
        help='Experiment identifier')
    p03.add_argument('-indir', dest='indir', metavar='DIR', required=True, default=None,
        help='Experiment runfolder (specify absolute path)')
    p03.add_argument('-instanceN', dest='instanceN', metavar='STR', required=False, default='000',
        help='The instanceN basecalling instance to extract data from')
    p03.add_argument('-outdir', dest='outdir', metavar='DIR', required=True, default=None,
        help='Output directory (specify absolute path)')
    p03.add_argument('-fastq', dest='fastq', metavar='BOOL', required=False, default='True',
        help='1D and 2D basecalls')
    #p03.add_argument('-model', dest='model', metavar='BOOL', required=False, default='True',
    #    help='Model parameters used in basecalling')
    p03.add_argument('-pairs', dest='pairs', metavar='BOOL', required=False, default='False',
        help='Name-value pairs for each experiment and read attribute')
    p03.add_argument('-stats', dest='stats', metavar='BOOL', required=False, default='True',
        help='Single-row summary stats for each experiment and read')
    p03.add_argument('-samplesize', dest='samplesize', metavar='INT', type=int, required=False, default=10000000,
        help='Number of FAST5 files to inspect from each expt, useful for testing.')
    p03.add_argument('-fastqheaderformat', dest='fastqheaderformat', metavar='FILE', required=False, default='concise',
        help='Output FASTQ header format [fast5|concise|poretools] ')
    p03.set_defaults(func=run_subtool)

    p04 = subparsers.add_parser('aggregateone', help='Aggregate value from one experiment into time buckets',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p04.add_argument('-bin', dest='bin', metavar='DIR', required=False, default='./',
        help='marcoporo scripts dir (specify absolute path)')
    p04.add_argument('-profile', dest='profile', metavar='FILE', required=False, default=None,
        help='marcoporo environment statements (specify absolute path)')
    p04.add_argument('-config', dest='config', metavar='FILE', required=True, default='config.txt',
        help='Analysis configuration file')
    p04.add_argument('-exptid', dest='exptid', metavar='FILE', required=True, default=None,
        help='Experiment identifier')
    p04.add_argument('-indir', dest='indir', metavar='DIR', required=True, default=None,
        help='Experiment runfolder (specify absolute path)')
    p04.add_argument('-maxrunlen', dest='maxrunlen', metavar='FLOAT', type=float, required=False, default=48,
        help='Aggregate metrics for maxrunlen (in hours).')
    p04.add_argument('-timebucket', dest='timebucket', metavar='FLOAT', type=float, required=False, default=0.25,
        help='Time window over which to aggregate metrics (in hours).')
    p04.add_argument('-outdir', dest='outdir', metavar='DIR', required=True, default=None,
        help='Output directory (specify absolute path)')
    p04.add_argument('-execjobs', dest='execjobs', metavar='BOOL', required=False, default='True',
        help='If False, create the job scripts but do not execute.')
    p04.set_defaults(func=run_subtool)

    p05 = subparsers.add_parser('aggregate', help='Aggregate value from one experiment into time buckets',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p05.add_argument('-bin', dest='bin', metavar='DIR', required=False, default='./',
        help='marcoporo scripts dir (specify absolute path)')
    p05.add_argument('-profile', dest='profile', metavar='FILE', required=False, default=None,
        help='marcoporo environment statements (specify absolute path)')
    p05.add_argument('-config', dest='config', metavar='FILE', required=True, default='config.txt',
        help='Analysis configuration file')
    p05.add_argument('-experiments', dest='experiments', metavar='FILE', required=True, default=None,
        help='Experiments and analysis parameters')
    p05.add_argument('-indir', dest='indir', metavar='DIR', required=True, default=None,
        help='Experiment runfolder (specify absolute path)')
    p05.add_argument('-maxrunlen', dest='maxrunlen', metavar='FLOAT', type=float, required=False, default=48,
        help='Aggregate metrics for maxrunlen (in hours).')
    p05.add_argument('-timebucket', dest='timebucket', metavar='FLOAT', type=float, required=False, default=0.25,
        help='Time window over which to aggregate metrics (in hours).')
    p05.add_argument('-outdir', dest='outdir', metavar='DIR', required=True, default=None,
        help='Output directory (specify absolute path)')
    p05.add_argument('-execjobs', dest='execjobs', metavar='BOOL', required=False, default='True',
        help='If False, create the job scripts but do not execute.')
    p05.set_defaults(func=run_subtool)

  # Parse the arguments

    args = parser.parse_args()
    if not args.bin or args.bin is None:
        args.bin = _bin
    if not args.profile or args.profile is None:
        args.profile = _profile

  # Read the -config file

    if not _P.config_read(args.config):
        sys.stderr.write('Error: Failed to read conf file ({0})\n'.format(args.config))
        sys.exit(_P.err_code('ErrorReadingData'))

  # Set logging verbosity and output file

    now = time.localtime()
    logdir = None
    try:
        logdir = _P.option['program']['logdir']
    except:
        sys.stderr.write('Error: Failed to set log dir\n')
        sys.exit(_P.err_code('ErrorDirCreate'))
    #if logdir is None:
    #    try:
    #        logdir = _P.option['program']['logdir']
    #    except:
    #        sys.stderr.write('Error: Failed to set log dir\n')
    #        sys.exit(_P.err_code('ErrorDirCreate'))

    if not os.path.exists(logdir):
        os.makedirs(logdir)
    logfile = '{0}-marcoporo-{1:03d}.log'.format(time.strftime('%Y%m%d-%H%M%S', now), random.randint(1, 100))
    logpath = os.path.join(logdir, logfile)

    handler = logging.FileHandler(logpath)
    handler.setLevel(logging.INFO)
    loggingverbosity = _P.option['program']['loggingverbosity'].upper()
    if loggingverbosity == 'ERROR':
        logger.setLevel(logging.ERROR)
        handler.setLevel(logging.ERROR)
    elif loggingverbosity == 'DEBUG':
        logger.setLevel(logging.DEBUG)
        handler.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s %(levelname)s : %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    logger.info('Program: {0}'.format(marcoporoversion.__program__))
    logger.info('Version: {0}'.format(marcoporoversion.__version__))
    logger.info('Descrip: {0}'.format(marcoporoversion.__description__))
    logger.info('Command: {0}'.format(' '.join(sys.argv)))

  # Call the selected sub-command, ignore SIGPIPEs (32)
    try:
        args.func(parser, args, _P, logger, handler)
    except IOError, e:
        if e.errno != 32:
            logger.info('Received kill signal')
            raise

if __name__ == '__main__':
    main()
