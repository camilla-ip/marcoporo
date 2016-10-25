#!/usr/bin/env python

import h5py
import numpy as np
import os
import psutil
import re
import signal
import sys

class marcoporolib(object):

    def __init__(self):

      # constants

        self.args_calltypeL = ['2D', '1T', '1C', 'mixed', 'unknown']

      # err

        self.err = {
        'SuccessReturn'         : [  0, 'Info', 'Successfully completed' ],
        'SuccessHelp'           : [  1, 'Info', 'Successfully printed help message' ],
        'ErrorArgsMissing'      : [  2, 'Erro', 'Failed to retrieve expected command-line argument(s)' ],
        'ErrorArgsInvalid'      : [  3, 'Erro', 'Failed due to invalid command-line argument(s)' ],
        'ErrorArgsParseFailure' : [  4, 'Erro', 'Failed to parse command-line argument(s)' ],
        'ErrorSagaConnect'      : [  5, 'Erro', 'Failed to connect to SCAMPI db' ],
        'ErrorSagaDisconnect'   : [  6, 'Erro', 'Failed to disconnect from SCAMPI db' ],
        'ErrorFileOpen'         : [  7, 'Erro', 'Failed to open file' ],
        'ErrorSysCall'          : [  8, 'Erro', 'Failed system call' ],
        'ErrorFork1Failure'     : [  9, 'Erro', 'Failed to do first process fork' ],
        'ErrorFork2Failure'     : [ 10, 'Erro', 'Failed to do second process fork' ],
        'SuccessUsage'          : [ 11, 'Info', 'Successfully printed usage message' ],
        'SuccessVersion'        : [ 12, 'Info', 'Successfully printed version message' ],
        'ErrorDirCreate'        : [ 13, 'Erro', 'Failed to create directory' ],
        'ErrorExistingLogFile'  : [ 14, 'Erro', 'Failed because log file already exists' ],
        'ErrorDirDoesNotExist'  : [ 15, 'Erro', 'Failed because output directory does not exist' ],
        'SuccessNothingToDo'    : [ 16, 'Info', 'Successfully terminated because nothing to do' ],
        'ErrorIniFileMissing'   : [ 17, 'Erro', 'Failed to find config.ini file' ],
        'ErrorInvalidData'      : [ 18, 'Erro', 'Failed due to invalid input data' ],
        'ErrorExternalSoftware' : [ 19, 'Erro', 'Failed in call to external software' ],
        'SuccessDryRun'         : [ 20, 'Info', 'Successfully terminated in dryrun mode' ],
        'SuccessKilled'         : [ 21, 'Info', 'Successfully terminated by KILL signal' ],
        'ErrorKilled'           : [ 22, 'Erro', 'Failed because terminated by KILL signal' ],
        'ErrorFileMissing'      : [ 23, 'Erro', 'Failed because file is missing' ],
        'ErrorChmod'            : [ 24, 'Erro', 'Failed to change mode of file or directory' ],
        'ErrorScriptCompute'    : [ 25, 'Erro', 'Failed in spawned script' ],
        'ErrorDeletingFile'     : [ 26, 'Erro', 'Failed to delete file' ],
        'ErrorEmptyFile'        : [ 27, 'Erro', 'Failed due to file empty' ],
        'ErrorOutputIncomplete' : [ 28, 'Erro', 'Failed due to incomplete output data' ],
        'ErrorInvalidRetKey'    : [ 29, 'Erro', 'Failed due to invalid return code key in program' ],
        'ErrorSagaDataGet'      : [ 30, 'Erro', 'Failed to retrieve data from SCAMPI' ],
        'ErrorOutfileExists'    : [ 31, 'Erro', 'Failed because output file already exists' ],
        'ErrorDirRename'        : [ 32, 'Erro', 'Failed to rename dir' ],
        'ErrorFileRename'       : [ 33, 'Erro', 'Failed to rename file' ],
        'ErrorReadingData'      : [ 34, 'Erro', 'Failed to read data' ],
        'ErrorFileCopy'         : [ 35, 'Erro', 'Failed to copy file' ],
        'ErrorFileNotLink'      : [ 36, 'Erro', 'Have a file not a symlink' ],
        'ErrorOutputNotComplete': [ 37, 'Erro', 'Output is incomplete' ]
        }
        self.section = []
        self.option = {}        # self.option[section][var] = val
        self.lasterror = ''

        self.fast5_runnumberpattern = re.compile('.*_(\d*)$')
        self.fast5_readnumberpattern = re.compile('.*[Rr][Ee][Aa][Dd]_(\d*)$')

    # dmn
        self._dmn_sighandlerfn = None

    # err

    def err_dump(self):
        'Return a string containing all the error codes.'
        L = [[self.err[key][1], self.err[key][0], key, self.err[key][2]] for key in self.err.keys()]
        sorted(L, key = lambda x: int(x[0]))
        return '\n'.join(['\t'.join([str(x) for x in elt]) for elt in L])

    def err_code(self, retkey):
        'Return the numerical return code value corresponding to the return code descriptor.'
        return self.err[retkey][0]

    def err_retkey(self, retcode):
        'Return the return code descriptor string corresponding to the numerical return code.'
        return [key for key in self.err.keys() if self.err[key][0] == retcode]

    def err_text(self, retkey):
        'Return a string reporting the return code and the standardised message for that return code.'
        if self.err.has_key(retkey):
            return '{type}: Program finished with marcoporo retcode={code} - {msg}\n'.format(
                type=self.err[retkey][1], code=self.err[retkey][0], msg=self.err[retkey][2])
        return ''

    def err_exit(self, retkey, exitprogram=False):
        'Function to report exit code and exit program.'
        if exitprogram:
            sys.exit(self.err[retkey][0])

    def err_last(self):
        'Return the last error message.'
        return self.lasterror

    # str

    def str_2bool(self, s):
        'Return the Python boolean value True or False depending on what the string is.'
        return (s is not None and s.lower() in ['1', 't', 'true', 'y', 'yes'])

    # fast5

    def fast5_cleantype(self, typestring):
        'Return a type without the extraneous formatting (private function).'
        result = typestring.replace('numpy.', '')
        if result.startswith('S') or result.startswith('|S'):
            result = 'string'
        elif 'string' in result:
            result = 'string'
        else:
            result = result.replace('<type \'', '').replace('\'>', '')
        return result

    def fast5_attributes(self, fast5path, keeptype='all'):
        '''
        Return dictionaries containing the attributes.
        keeptype should be 'all' or 'paramsonly'.
        The code is so ugly because lots of the data access stuff only works in
        particular ONT FAST5 file versions.
        '''
        attributeD = {}
        hdf = h5py.File(fast5path, 'r')
        list_of_names = []
        hdf.visit(list_of_names.append)
      # Create a dictionary that maps all the fields with a run number
      # e.g., runnumberD['Basecall_1D_NNN'] = ['Basecall_1D_000', '000']
        runnumberD = {}
        for name in list_of_names:
            matchL = [x for x in name.split('/') if 'read' not in x.lower() and self.fast5_runnumberpattern.match(x) is not None]
            for match in matchL:
                tokenL = match.split('_')
                if len(tokenL) > 1:
                    n = tokenL[-1]
                    key = '_'.join(tokenL[:-1]) + '_' + 'N'*len(n)
                    if n not in runnumberD.keys():
                        runnumberD[key] = [match, n]
      # Do the same thing to find the read numbers
      # e.g., readnumberD['Read_NN'] = ['Read_24', '24']
        readnumberD = {}
        for name in list_of_names:
            matchL = [x for x in name.split('/') if self.fast5_readnumberpattern.match(x) is not None]
            for match in matchL:
                tokenL = match.split('_')
                if len(tokenL) > 1:
                    n = tokenL[-1]
                    key = '_'.join(tokenL[:-1]) + '_' + 'N'*len(n)
                    if n not in readnumberD.keys():
                        readnumberD[key] = [match, n]
      # Extract all the datasets and attributes
        for name in list_of_names:
            fieldtype = 'unknown'
            if type(hdf[name]) == h5py._hl.dataset.Dataset:
                fieldtype = 'dataset'
                try:
                    if str(hdf[name][()].dtype).startswith('|S'):
                        fieldtype = 'string'
                    else:
                        pass
                except:
                    pass
            else:
                fieldtype = 'unknown'

            if fieldtype == 'dataset':
                typename = 'UNKNOWN'
                try:
                    typename = self.fast5_cleantype(','.join([str(hdf[name][()].dtype[i]) for i in range(0, len(hdf[name][()].dtype.names))]))
                    key = '{0}({1})'.format(name, ','.join([x for x in hdf[name][()].dtype.names]))
                    val = '({0})'.format(','.join([str(x) for x in hdf[name][()][0]]))
                    attributeD[key] = [typename, val]
                except:
                    pass
                if typename == 'UNKNOWN':
                    try:
                        typename = self.fast5_cleantype(type(hdf[name][()]))
                        key = '{0}({1})'.format(name, 'typename')
                        val = '({0})'.format(','.join([str(x) for x in hdf[name][()][0]]))
                        attributeD[key] = [typename, val]
                    except:
                        pass
            elif fieldtype == 'string':
                typename = 'UNKNOWN'
                try:
                    typename = self.fast5_cleantype(str(hdf[name][()].dtype))
                    key = '{0}({1})'.format(name, ','.join([x for x in hdf[name][()].dtype.names]))
                    val = '({0})'.format(','.join([str(x) for x in hdf[name][()][0]]))
                    attributeD[key] = [typename, val]
                except:
                    pass
                if typename == 'UNKNOWN':
                    try:
                        typename = self.fast5_cleantype(type(hdf[name][()]))
                        key = 'name'
                        val = '{0}...{1}'.format(hdf[name][()][0:10], hdf[name][()][-10:]).rstrip()
                        attributeD[key] = [typename, val]
                    except:
                        pass
            try:
                itemL = hdf[name].attrs.items()
                for item in itemL:
                    attr, val = item
                    if keeptype == 'all' or (keeptype == 'paramsonly' and name+'/'+attr in paramL):
                        if type(hdf[name].attrs[attr]) == np.ndarray:
                            val = ''.join(hdf[name].attrs[attr])
                        typename = self.fast5_cleantype(str(type(hdf[name].attrs[attr])))
                        key = name+'/'+attr
                        val = str(val).replace('\n', '')
                        attributeD[key] = [typename, val]
            except:
                pass
        hdf.close()
        return attributeD, runnumberD, readnumberD

    # conf
    def config_read(self, inipath):
        'Read all the data from the inipath into section and option variables; set lasterror if anything goes wrong.'
        line_section = "UNDEFINED"
        optionflat = {}
        with (open(inipath, "r")) as in_fp:
            for line in in_fp:
              # Strip all leading or trailing spaces from the line, replace all tabs by spaces.
              # At this point, there still could be lots of spaces in the 'var' and/or 'val' part.
                line = line.strip().replace("\t", "")
                if (line.startswith(';') or line.startswith('#') or not len(line)):
                    continue
                if (line.startswith('[')):
                    line_section = line.replace("[", "").replace("]", "")
                    self.section.append(line_section)
                    self.option[line_section] = {}
                if ("=" in line):
                    if (line_section == "UNDEFINED"):
                        self.lasterror = "Erro: ini file formatting error: var=val pairs outside of a [section]."
                        in_fp.close()
                        return False
                    L = line.split("=")
                    var = L[0].replace(" ", "")         # Contains no spaces
                    val = "=".join(L[1:]).strip()       # Could contain spaces
                    self.option[line_section][var] = val
                    optionflat[var] = val
            in_fp.close()
      # The section that replaces any instances of {variable} with the value of self.option[section][variable]
      # For every variable, replace {variable} with value in every other value in the dictionary.
        for section, D in self.option.iteritems():
            for var, val in D.iteritems():
                try:
                    localenvvar = re.search("[\w\.]*{(\w*:\w*)}[\w\.]*", val).group(1)
                    localsection, localvar = localenvvar.split(':')
                    oldstring = '{0}{1}{2}'.format('{', localenvvar, '}')
                    newstring = self.option[localsection][localvar]
                    self.option[section][var] = self.option[section][var].replace(oldstring, newstring)
                except:
                    pass
        return True

    # dmn

    def _dmn_set_traps_for_all_signals_except_kill_and_equivalents(self):
        """
        Set trap using the signal handling fn from the calling program for all
        signals except SIGKILL, SIGSTOP and SIG_DFL, which cannot be trapped.
        """
        for i in [x for x in dir(signal) if x.startswith("SIG")]:
            try:
                signum = getattr(signal,i)
                signal.signal(signum, self._dmn_sighandlerfn)
            except:
                pass

    def dmn_daemonize (self, stdin='/dev/null', stdout='/dev/null', stderr='/dev/null'):
        """Do two forks and other operations to change the running process into a daemon."""
        # Based on http://itmanagement.earthweb.com/netsys/article.php/3786866/Creating-a-Daemon-with-Python.htm

        self._dmn_set_traps_for_all_signals_except_kill_and_equivalents()

        # Perform first fork. Exit first parent on failure.
        try:
            pid = os.fork()
            if pid > 0:
                sys.exit(0)
        except OSError, e:
            sys.stderr.write("fork #1 failed: (%d) %s\n" % (e.errno, e.strerror))
            self.err_tidy('ErrorFork1Failure', exitprogram=True)

        # Decouple from parent environment.
        os.chdir("/")
        os.umask(0)
        os.setsid()

        # Perform second fork. Exit second parent on failure.
        try:
            pid = os.fork()
            if pid > 0:
                sys.exit(0)
        except OSError, e:
            sys.stderr.write("fork #2 failed: (%d) %sn" % (e.errno, e.strerror))
            self.err_tidy('ErrorFork2Failure', exitprogram=True)

        # The process is now daemonized and we have an open file descriptor to a log file.

        return

        # Close the 3 standard file descriptors.
        for f in sys.stdout, sys.stderr: f.flush()
        si = file(stdin, 'r')
        so = file(stdout, 'a+')
        se = file(stderr, 'a+', 0)
        os.dup2(si.fileno(), sys.stdin.fileno())
        os.dup2(so.fileno(), sys.stdout.fileno())
        os.dup2(se.fileno(), sys.stderr.fileno())

        # Return to the calling function where you should now start executing
        # the 'payload' of your program.

  # sys

    def sys_pids(self, machineid=None, procname=None, procuser=None, includecpid=False):
        '''
        Return the list of process ids for all processes called procname
        owned by procuser, possibly including the current pid. Does not use
        list comprehensions to find the pid lists because in between getting
        the pid list and testing properties of the pids, the processes may
        terminate and then it is not possible to get the information associated
        with those pids any more.
        ''' 
        pidlist = []
        pidcurr = psutil.get_pid_list()
        for pid in pidcurr:
            try:
                cmdline = ' '.join(psutil.Process(pid).cmdline())
                if procname is None \
                or (len(psutil.Process(pid).cmdline()) >= 2 and psutil.Process(pid).cmdline()[0].lower().startswith('python') and procname in cmdline):
                    if procuser is None or psutil.Process(pid).username() == procuser:
                        if machineid is None or machineid in cmdline:
                            pidlist.append(pid)
            except:
                pass
        cpid = os.getpid()
        if not includecpid and cpid in pidlist:
            pidlist.remove(cpid)
        return pidlist

    def sys_kill(self, pidlist):
        'Try to kill all process ids in pidlist.'
        faillist = []
        for pid in pidlist:
            try:
                psutil.Process(pid).kill()
            except:
                faillist.append(pid)
                return False
        if len(faillist):
            self.lasterr = 'Failed to kill pids: {0}'.format(','.join([str(x) for x in faillist]))
            return False
        return True
