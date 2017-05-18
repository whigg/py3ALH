'''
Utilities module
================
  
  These are some functions to perform some common tasks are 
  collected for use in the TropomiDataProcessor collection.
  
  @author: Maarten Sneep
  @contact: maarten.sneep@knmi.nl
  @organization: KNMI
  @copyright: Maarten Sneep, KNMI
  @license: Creative Commons Attribution-ShareAlike 3.0 Netherlands License.
'''

from __future__ import division
from __future__ import print_function

import time
import os.path
import subprocess
import sys
import logging
import inspect


def classFromString(className, mod=None):
    """Load a module and get the class from a string.

By default the module name is assumed to be the same as the classname. 
It is possible to set a separate module name through the mod keyword 
parameter.

@param className: The name of the class to load. By default this function 
looks for C{Name.Name} (the class C{Name} in the module C{Name}). 
This can be changed with the C{mod} parameter.
@type className: C{string}
@param mod: The module name, to override the default.
@type mod: C{string}

@return: The class object that is referenced.
@rtype: A class object.

@note: The function returns C{None} if the module can't be loaded. 
ImportError and NameError are handled by the function, other errors are reC{raise}d. 
"""
    if mod is None:
        mod = className
    if className == "NoneType":
        cls = None
    else:
        try:
            __import__(mod, globals(), locals(), [], -1)
            cls = sys.modules[mod].__dict__[className]
        except ImportError:
            try:
                cls = eval("{0}".format(className))
            except NameError:
                print('Class "{0}" from modue "{1}"'
                    ' was not found.'.format(className, mod))
                return
            except:
                print('An unanticipated error occurred '
                      'while trying to find Class "{0}"'
                      ' in module "{1}".'.format(className, mod))
                raise
        except:
            print('Module "{0}" was not found, terminating'.format(mod))
            raise
    return cls

def locateResource(name, loc="tbl", isFile=True, mustExist=True, base=None):
    """Locate a resource within the source distribution.

A resoure can be any file or directory in the TropomiDataProcessor 
distribution. Before looking in the source distribution, we look in the 
current directory. A warning is printed if the file is found there (but 
it will be used).

@param name: The name of the file.
@type name: C{string}
@param loc: The subdirectory in which to look for the file. 
Defaults to C{tbl} for the C{tbl} sibling of the C{src} directory.
@type loc: C{string}
@param isFile: Is the specified resource a file? Defaults to C{True}.
@type isFile: Bool
@param mustExist: Must the specified resource exist? Defaults to C{True}.
@type mustExist: Bool
@param base: the file which is to be used as the starting point for the search.
    The default is to start from C{__file__}.

@return: The fully qualified path to the specified resource.
@rtype: C{string}

@raise ValueError: If the specified resource does not match the 
requirements that are specified with the C{isFile} and C{mustExist} 
parameters.
"""
    if mustExist and isFile and os.path.exists(name):
        if os.path.isabs(name):
            path = name
        else:
            path = os.path.realpath(name)
            sys.stderr.write("""Found file "{0}" in the current directory, not searching in "{1}".\n""".format(name, loc))
    else:
        if base is None:
            base = __file__
        
        path = os.path.join(
                   os.path.dirname(
                       os.path.dirname(
                           os.path.realpath(base))),
                           loc, name)
    
    if mustExist and (not os.path.exists(path)):
        raise ValueError('File "{0}" not found in the distribution'
            ' (in the "{1}" directory).'.format(name, loc))
    if (mustExist) and (isFile) and (not os.path.isfile(path)):
        raise ValueError('Item "{0}" is not a file.'.format(name))
    if (mustExist) and (not isFile) and (not os.path.isdir(path)):
        raise ValueError('Item "{0}" is a file, expected a directory.'.format(name))
    
    return path

_basename = 'knmi'

def logname():
    """
    Raise a fake exception, Find the stack frame of the caller so that we 
    can traverse to the source filename that started the code.  
    From this 'start' code filename the basename is taken without the .py
    
    If the component found is equal to "pdb", it is assumed the process is
    running inside the python debugger and element 5, instead of 0, of the
    stack is examined to find the component.
    
    If the component found is equal to "threading", it is assumed the process is
    running inside a thread object and element 1, instead of 0, of the stack is
    examined to find the component.
    
    Code inspired by (Ian van der Neut).
    """
    global _basename
    
    parent = os.path.splitext(os.path.basename(wheresdaddy()))[0]
    return '.'.join([_basename, os.path.splitext(os.path.basename(sys.argv[0]))[0], parent])

def whoami():
    return inspect.stack()[1][3]
def whosdaddy():
    return inspect.stack()[2][3]
def whereami():
    return inspect.stack()[1][1]
def wheresdaddy():
    return inspect.stack()[2][1]

def indent_string(s, indent="    "):
    return '\n'.join([indent + c for c in s.split('\n')])

class MyFormatter(logging.Formatter):
    """
    Subclass to change the formatting of the message string. 
    
    Indent every line of the message string.
    """
    def format(self, record):
        """
        Format the specified record as text.

        The record's attribute dictionary is used as the operand to a
        string formatting operation which yields the returned string.
        Before formatting the dictionary, a couple of preparatory steps
        are carried out. The message attribute of the record is computed
        using LogRecord.getMessage(). If the formatting string contains
        "%(asctime)", formatTime() is called to format the event time.
        If there is exception information, it is formatted using
        formatException() and appended to the message.
        """
        record.message = indent_string(record.getMessage())
        if "%(asctime)" in self._fmt:
            record.asctime = self.formatTime(record, self.datefmt)
        s = self._fmt % record.__dict__
        if record.exc_info:
            # Cache the traceback text to avoid converting it multiple times
            # (it's constant anyway)
            if not record.exc_text:
                record.exc_text = self.formatException(record.exc_info)
        if record.exc_text:
            if s[-1:] != "\n":
                s = s + "\n"
            s = "{0}    Exception:\n    {1}".format(s, indent_string(record.exc_text))
        return s

def logger():
    name = logname()
    
    # create the logging object
    return logging.getLogger(name)

def setup_logging(name=None, handlers=None, levels=None, filenames=None):
    global _basename
    if name is not None and type(name) == type(""):
        _basename = name
    
    logger = logging.getLogger(_basename)
    
    logger.setLevel(logging.NOTSET)
    
    # create formatter
    formatter = MyFormatter(fmt='%(asctime)s %(name)s %(levelname)s %(module)s %(funcName)s %(lineno)d:\n\t%(message)s',
                            datefmt='%Y-%m-%dT%H:%M:%S')
    formatter.converter = time.gmtime
    
    if handlers is None:
        handlers = [logging.StreamHandler]
    
    if levels is None:
        levels = ["ERROR"]
    
    if filenames is None:
        filenames = [None]
    
    # add all requeste handlers. For now only stream and file handlers are allowed.
    for cls, level, filename in zip(handlers, levels, filenames):
        if cls is logging.StreamHandler:
            h = cls()
        elif cls is logging.FileHandler:
            h = cls(filename, mode='a')
        else:
            continue
        
        numeric_level = getattr(logging, level.upper(), None)
        h.setLevel(numeric_level)
        h.setFormatter(formatter)
        logger.addHandler(h)
    
    # capture warnings.
    logging.captureWarnings(True)

def available_memory(kind=None):
    """Returns the RAM (total, used & free) of a linux system"""
    values = [v.split() for v in subprocess.check_output(["/usr/bin/free", "-m"]).split('\n')[1:] if v.split()]
    
    d = dict([(k[0][0:-1], dict(zip(('total', 'used', 'free'), [int(s) for s in k[1:]]))) 
              for k 
              in [v[0:4] 
                  for v 
                  in values 
                  if v[0] in ('Mem:', 'Swap:')]])
    if kind in ('Mem', 'Swap'):
        return d[kind]
    else:
        return d

