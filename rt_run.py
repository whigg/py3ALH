#!/usr/bin/env python3
#
# Maarten Sneep, 2010
# Purpose: run the radiative transfer module
#

'''
@note: Run DISAMAR using a given configuration file

Created on Feb 25, 2010

@author: Maarten Sneep
@organization: KNMI
'''

__author = "Maarten Sneep"
__date = "2010-02-25"
__version = "0.1a"

import os
import sys
import shutil
import tempfile
import subprocess
import resource

from time import gmtime, time

import rt_cfg
import datastore

class rt_run(object):  #IGNORE:C0103
    '''run the radiative transfer model (and retrieval system) DISAMAR
    '''
    def __init__(self, cfg=None,
                 disamar=None,
                 output=None,
                 spectrum=None,
                 debug=False,
                 quiet=False,
                 tempbase=None):
        '''@type cfg: rt_cfg.RT_configuration object or path to config file.
@type disamar: string (path to DISAMAR executable)
@type output: filename for output.
'''
        if type(cfg) == type(""):
            if os.path.exists(cfg):
                self._cfgfile = os.path.abspath(cfg)
                self.config = rt_cfg.RT_configuration(self._cfgfile)
            else:
                raise ValueError("File {0} not found".format(cfg))
        else:
            self._cfgfile = os.path.abspath(cfg.file())
            self.config = cfg

        self._debug = debug
        self._quiet = quiet

        if not disamar:
            disamar = "/usr/people/sneep/Documents/Disamar/Disamar.exe"
        if not os.path.exists(disamar):
            raise ValueError('Disamar executable not found: {0}'.format(disamar))
        self.disamar = disamar

        if spectrum is not None and os.path.exists(spectrum):
            self.spectrum = os.path.abspath(spectrum)
        else:
            self.spectrum = None

        if output:
            self.output = os.path.abspath(output)
        else:
            f = self._cfgfile
            self.output = os.path.join(os.path.dirname(f), os.path.splitext(os.path.basename(f))[0] + '.hdf5')

        self.attributes = None

        fmt = "{time[0]:04d}{time[1]:02d}{time[2]:02d}_{time[3]:02d}{time[4]:02d}{time[5]:02d}_{pid:05d}"
        date_string = fmt.format(time=gmtime()[0:6], pid=os.getpid())

        if tempbase is not None and os.path.exists(tempbase):
            tempfile.tempdir = os.path.abspath(tempbase)

        self.tmploc = os.path.join(tempfile.gettempdir(), 'DISAMAR.' + date_string)
        os.mkdir(self.tmploc)
        self.make_copy()

    def __del__(self):
        if hasattr(self, 'tmploc'):
            if hasattr(self, '_debug') and not self._debug:
                shutil.rmtree(self.tmploc, True)

    def save_config(self):
        '''Save the configuration object in file 'fname'. '''
        self.copy_spectrum()
        self.config.write(filename=os.path.join(self.tmploc, 'Config.in'))

    def copy_spectrum(self):
        if self.spectrum and os.path.exists(self.spectrum):
            if self._debug:
                print('using spectrum from {0}'.format(self.spectrum))
            filename = os.path.basename(self.spectrum)
            self.config['GENERAL','overall','useReflFromFile'].setvalue(1)
            self.config['GENERAL','fileNames','externalReflectanceFileName'].setvalue(filename)
            os.symlink(self.spectrum, os.path.join(self.tmploc, filename))
        else:
            if self._debug:
                print('Simulate spectrum')
            self.config['GENERAL','overall','useReflFromFile'].setvalue(0)
            self.config['GENERAL','fileNames','externalReflectanceFileName'].setvalue('disamar.sim')
            self.config['ADDITIONAL_OUTPUT','additional','refl_Instr_gridSim'].setvalue(1)

    def __call__(self, *extrafiles, **kwargs):
        self.save_config()
        exename = os.path.basename(self.disamar)
        if os.path.dirname(self.output) != "":
            destwd = os.path.abspath(os.path.dirname(self.output))
        else:
            destwd = os.path.abspath(os.getcwd())
        oldcwd = os.getcwd()
        os.chdir(self.tmploc)
        stdoutf = tempfile.TemporaryFile(mode='w+')
        if self._debug:
            print("TEMPORARY_PATH: {0}".format(self.tmploc))
            for item in os.listdir(self.tmploc):
                print("LISTING --- {0}".format(item))
        try:
            starttime = time()
            failure = subprocess.call(os.path.join(self.tmploc, exename),
                                      shell=True,
                                      stdout=stdoutf,
                                      stderr=subprocess.STDOUT)
            timelapse = time() - starttime

            stdoutf.seek(0)
            for line in stdoutf:
                if 'stopped' in line.lower():
                     raise RuntimeError("Disamar stopped with an error")

            stdoutf.seek(0)
            if failure < 0:
                sys.stderr.write("DISAMAR was terminated by signal {0:d}\n".format(-failure))
                for line in stdoutf:
                    sys.stderr.write(line)
                raise RuntimeError("DISAMAR was terminated by signal {0:d}\n".format(-failure))
            elif failure != 0:
                sys.stderr.write("DISAMAR returned an error ({0:d})\n".format(failure))
                for line in stdoutf:
                    sys.stderr.write(line)
                raise RuntimeError("DISAMAR returned an error ({0:d})\n".format(failure))
        except RuntimeError as e:
            sys.stderr.write("Execution of DISAMAR failed: {0}\nDisamar wrote to the output:\n".format(e))
            stdoutf.seek(0)
            for line in stdoutf:
                sys.stderr.write("\t" + line)
            raise

        os.chdir(destwd)

        attrs = kwargs.copy()
        attrs['uid'] = os.getuid()
        fmt = "{time[0]:04d}{time[1]:02d}{time[2]:02d}T{time[3]:02d}:{time[4]:02d}:{time[5]:02d}"
        attrs['date'] = fmt.format(time=gmtime()[0:6])
        attrs['user'] = os.environ['LOGNAME']
        if 'HOSTNAME' in list(os.environ.keys()):
            attrs['host'] = os.environ['HOSTNAME']
        else:
            attrs['host'] = '<Unknown>'

        attrs['time'] = timelapse

        stdoutf.seek(0)
        saw_error = False
        attrs['Error'] = ''
        attrs['ErrorMsg'] = ''

        for line in stdoutf:
            if not self._quiet:
                print(line[0:-1])
            if line.startswith('time for '):
                label, value = line.split('=')
                try:
                    val = float(value)
                except ValueError:
                    val = float('nan')
                attrs["{0}_{1}".format(label[9:-7].replace(' ', '_'), "time")] = val
            if saw_error:
                saw_error = False
                attrs['ErrorMsg'] = line.strip()
            if 'ERROR' in line:
                saw_error = True
                attrs['Error'] = line.strip()
        stdoutf.close()

        extrafiles = list(extrafiles)
        extrafiles.append('disamar.sim')

        self.attributes = attrs

        try:
            if self._debug:
                print("temploc ({0}) listing:".format(self.tmploc))
                for item in os.listdir(self.tmploc):
                    print("\t{0}".format(item))

            datastore.translate(os.path.join(self.tmploc, 'disamar.asciiHDF'),
                                self.output,
                                attachment=os.path.join(self.tmploc, 'Config.in'),
                                attributes=attrs, debug=self._debug)

            for name in extrafiles:
                if name is not None:
                    if os.path.exists(os.path.join(self.tmploc, os.path.basename(name))):
                        ext = os.path.splitext(name)[1]
                        shutil.copy(os.path.join(self.tmploc, os.path.basename(name)),
                                os.path.join(destwd,
                                             os.path.splitext(os.path.basename(self.output))[0] + ext))
                    else:
                        if os.path.basename(name) != "disamar.sim":
                            print("File {0} not found".format(os.path.basename(name)))
                            print("Temporary location: {0}".format(self.tmploc))
        except (ValueError, IOError) as err:
            print(err)
            for filename in ['disamar.asciiHDF', 'disamar.sim', 'Config.in']:
                ext = os.path.splitext(filename)[1]
                try:
                    shutil.copy(os.path.join(self.tmploc, filename),
                                os.path.join(destwd,
                                  os.path.splitext(os.path.basename(self.output))[0] + ext))
                except:
                    pass
            raise
        os.chdir(oldcwd)

    def make_copy(self):
        '''Make a copy of the DISAMAR executable and support files'''
        disdir = os.path.dirname(os.path.abspath(self.disamar))
        exename = os.path.basename(self.disamar)
        refspecdir = os.path.join(disdir, 'RefSpec')
        scatterdir = os.path.join(disdir, 'expCoefFiles')
        tomsdir = os.path.join(disdir, 'TOMS_v8')

        try:
            os.symlink(refspecdir, os.path.join(self.tmploc, 'RefSpec'))
        except:
            if self._debug:
                sys.write.stderr("Could not create symlink to 'RefSpec'\n")
            raise

        try:
            os.symlink(scatterdir, os.path.join(self.tmploc, 'expCoefFiles'))
        except:
            if self._debug:
                sys.write.stderr("Could not create symlink to 'expCoefFiles'\n"
                                 "Assuming that this version of DISAMAR doesn't need them\n")
        try:
            os.symlink(tomsdir, os.path.join(self.tmploc, 'TOMS_v8'))
        except:
            if self._debug:
                sys.write.stderr("Could not create symlink to ''TOMS_v8''\n"
                                 "Assuming that this version of DISAMAR doesn't need them\n")

        os.symlink(os.path.join(disdir, exename), os.path.join(self.tmploc, exename))



def __main():
    '''Commandline handler'''
    from argparse import ArgumentParser
    usage = """Run Disamar on the given input file."""
    parser = ArgumentParser(description=usage)#, version="%prog {0}".format(__version))

    parser.add_argument("-i", "--input", dest="infile",
                      help="Read Disamar configuration from FILE", metavar='FILE')
    parser.add_argument("-d", "--disamar", dest="disamar",
                      help="Use Disamar located at EXECUTABLE", metavar='EXECUTABLE')
    parser.add_argument("-s", "--spectrum", dest="spectrum_in", metavar='FILE',
                      help="Read input spectrum from FILE", default=None)
    parser.add_argument("-o", "--output", dest="outfile",
                      help="Write output to FILE (hdf-5 format)",
                      metavar="FILE")
    parser.add_argument("-D", "--debug", dest="debug",
                      help="Turn on debugging",
                      action="store_true")
    parser.add_argument("-Q", "--quiet", dest="quiet",
                      help="Be noisy, very noisy",
                      action="store_true")
    parser.add_argument('-t', '--tmp', dest='tempbase', metavar="DIR", default="/tmp",
                        help="temporary directory (base, the run-time directory is created here)")
    parser.add_argument("-I", "--extra", dest="extra_files", nargs="+", metavar="FILE",
                      help="Extra files to copy to the temporary run directory of Disamar")
    #parser.add_argument("input_file", help="Read Disamar configuration from FILE", metavar="FILE")

    args = parser.parse_args()
    if not args.infile:
        args.infile = args.input_file
        
    
    if not os.path.exists(args.outfile):

        runner = rt_run(cfg=args.infile,
                        disamar=args.disamar,
                        output=args.outfile,
                        spectrum=args.spectrum_in,
                        debug=args.debug,
                        quiet=args.quiet,
                        tempbase=args.tempbase)
        runner(args.extra_files)
        
    else:
        
        print('Output file exists')

if __name__ == '__main__':
    __main()
