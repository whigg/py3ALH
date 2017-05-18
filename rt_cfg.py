#!/usr/bin/env python
#
# Maarten Sneep, 2009
# Purpose: parser for a configuration file for the radiative transfer module
#

'''
DISAMAR Configuration file parser
=================================

  Read and write a configuration file for Johan's DISIMAR program.

  The format Johan suggested is somewhat convoluted, with sections,
  subsections, and the possibility of repeated keys within a subsection.
  The repeated keys give values for subsequent channels, so order is
  very important. To add to that, a key can contain multiple values
  on a single line. All parts are separated by whitespace.

  Comments are indicated by a # mark at the start of a line, or between ()
  after the last element on a line, and should be preserved as best we can.

  Using '&' at the end of a line, indicates that the next line continues
  the current one.

  In normal operation repeated section, subsections and keys add to,
  rather than replace existing sections. This can be controlled by adding a
  meta-comment. With the following line
  #@ replace_section = True
  existing sections will be replaced from that point in the file onwards.
  The same can be limited to subsections::

    #@ replace_subsection = True

  And finally, keys can be replaced, rather than added to with::

    #@ replace_keys = True

  Note that in the latter mode, the original behavior is I{severely} altered.

  In short: the format is almost too flexible. An example is included at
  the end of this doc string.

  Note that Johan has extensive testing of the contents of the configuration
  file. This class does not attempt to duplicate that behaviour. As long as
  the format of configuration is technically correct, it will be read and
  declared valid.

  Limitations
  -----------

    - Single line # comments are moved to the top of the file on output
    - The order in which key are emitted is random, only the relative order
      of repeated keys is guaranteed.

  @author: Maarten Sneep
  @contact: maarten.sneep@knmi.nl
  @organization: KNMI
  @copyright: Maarten Sneep, KNMI
  @license: Creative Commons Attribution-ShareAlike 3.0 Netherlands License.

  Example file
  ------------

  A short section of a sample configuration file::

    SECTION GENERAL
    subsection fitAltAerCldSurf
    fit_cloud_top         1
    fit_cloud_fraction    0

    subsection filenames
    referenceSpecFileName  refSpecHighLat.dat   (radiance reference spectrum)
    polCorrectionFileName  polCorLat-45.dat     (polarization correction file)

    SECTION INSTRUMENT_SPECIFICATIONS
    subsection FWHM (in nm)
    FWHM_Irradiance_sim    0.45
    FWHM_radiance_sim      0.45d0
    FWHM_Irradiance_retr   0.45
    FWHM_radiance_retr     0.45
    FWHM_Irradiance_sim    0.35
    FWHM_radiance_sim      0.35e0
    FWHM_Irradiance_retr   0.35
    FWHM_radiance_retr     0.35

    subsection slit_index
    # single line comment
    slit_index_Irradiance_sim  1
    slit_index_radiance_sim    1
    slit_index_Irradiance_retr 1
    slit_index_radiance_retr   1

    SECTION SCENE_PARAMETERS
    subsection albedo
    surfacealbedo 0.05 0.04 0.03 0.02 0.02

  Full configuration files are included in the C{InputFiles} dicrecory in the
  DISAMAR distribution.

  Usage
  -----

  A short sample on how to use the code::

    # load the package
    import rt_cfg

    # create configuration object, specifying a file will cause it to be read.
    cfg = rt_cfg.RT_configuration(filename='test.config')

    # select a configuration item
    item = cfg['GENERAL','filenames','referenceSpecFileName']

    # item is now a configuration item object. Its value can be obtained with
    print(item.values()) # which prints 'refSpecHighLat.dat'

    # change the value of an item
    cfg['SCENE_PARAMETERS','albedo','surfacealbedo'].setvalue([0.05, 0.04, 0.04])

    # Add an element. If you want to add a comment, use a tuple.
    cfg['SCENE_PARAMETERS','albedo','cloud'] = (0.8, 'comment')

    # change the file name
    cfg.setfile('new.config')

    # write to the new file
    cfg.write()
'''

# import integer division with sane results


import re, sys, os
import utilities

try:
    from collections import OrderedDict
except ImportError:
    from OrderedDict import OrderedDict

if sys.version_info[0] == 3:
    string_types = (str,)
else:
    string_types = (str,)


def dblformat(item):
    if type(item) == type(1.0):
        s = str(item).lower()
        if 'e' in s:
            rval = s.replace('e', 'D')
        else:
            rval = "{0}D0".format(s.strip())
    else:
        rval = str(item)
    return rval

class RT_dict(OrderedDict):
    '''
An extended dictionary

This extended dictionary includes a comment-field and a comment-block. The
base-class is an ordered dictionary, standard in python 2.7, but defined above
for older python versions. The contents are other L{RT_dict} objects or
L{RT_item} objects.

@ivar _comment: Comment for the section or subsection
@type _comment: C{string}
@ivar _blockcomment: Comment block for the file, section or subsection
@type _blockcomment: C{string}
'''
    def __init__(self, *args, **kwargs):
        '''Initializer

@param args: Positional arguments for the C{dict} constructor
@type args: C{[]}
@param kwargs: Keyword parameters for the C{dict} constructor. A few are
               filtered for use as comments.
@type kwargs: C{{}}

@keyword comment: A comment for the section or subsection, the one between parentheses
@keyword blockcomment: Comments for the block, inserted as #-marked lines.

The specified keyword are removed after they are used. The rest of the arguments
(positional and keywords) is passed on to C{super()}.

@note: This class is now much simpler because the OrderedDict is used from Python 2.7
'''
        if 'comment' in kwargs:
            comment = kwargs['comment']
            del kwargs['comment']
        else:
            comment = ''

        if 'blockcomment' in kwargs:
            blockcomment = kwargs['blockcomment']
            del kwargs['blockcomment']
        else:
            blockcomment = ''

        if len(args) >= 1:
            super(self.__class__, self).__init__(args[0])
        else:
            super(self.__class__, self).__init__(self)

        self._comment = comment
        self._blockcomment = blockcomment

    def setcomment(self, comment):
        '''Change the comment for the object

@param comment: The new comment for the group.
@type comment: C{string}.
'''
        self._comment = comment

    def comment(self):
        '''Obtain the current comment

@return: The current comment.
@rtype: C{string}
'''
        return self._comment

    def setblockcomment(self, comment):
        '''Build a longer comment block

The comments are appended to any existing block comments, separating them with a newline and a hash-mark.

@param comment: The new comment for the group.
@type comment: C{string}.
'''
        if not self._blockcomment:
            self._blockcomment = comment
        else:
            self._blockcomment = self._blockcomment + '\n#' + comment

    def blockcomment(self):
        '''Obtain the current block comment

@return: The current block comment,.
@rtype: C{string} (newline and hash-mark separated).
'''
        return self._blockcomment

    def __str__(self):
        return self.blockcomment() + "\n" + "\n".join([str(v) for v in list(self.values())])

class RT_item(object):
    '''A configuration item.

This object contains the name of the item, the value(s) and the comment string.
This is the finest granularity of a configuration file.

@ivar _value: The container for the values, most likely a list, or even a list of lists.
@ivar _comment: The parenthesized comment string for the item.
@ivar _key: The key to identify the object.
'''

    def __init__(self, key=None, value=None, comment=None):
        '''Initializer

@param value:   can be an existing 'item' object, in which case the remaining
         parameters are ignored. If it is a string, it will be parsed,
         otherwise the value is taken as-is.
@param key:     Required if value is not an 'item' object. String.
@param comment: A comment string for the item.
'''
        self._value = []
        self._comment = []
        if type(value) == type(self):
            self._key = value.key()
            self.setvalue(list(value.values()))
            self._comment = value.comments()
        else:
            if key is None:
                raise ValueError('a key is required')
            self._key = key
            if type(value) == type(''):
                self.appendvalue(value)
            else:
                self.setvalue(value)
            self.appendcomment(comment)

    def convertvalue(self, value): # IGNORE:R0201
        '''Convert a value as read from file to int, float, or string [list].

Always returns a list of item(s), even if there is a single value.
This is convenient when writing the file later. The converter splits the input
string on whitespace. For each item it first tries  to convert the string to an
integer. If that fails, it tries to convert the value to a float value. If that
fails it uses a string value as a fallback option.

@param value: The raw input as read from the file.
@type value: C{string}
'''
        # test if we have a value at all
        if value is None:
            value = ""

        # Split on whitespace
        if isinstance(value, str):
            vallist = value.split()
        elif isinstance(value, list):
            vallist = [str(item) for item in value]
        elif isinstance(value, (float, int)):
            vallist = [str(value)]
        else:
            raise TypeError('Unexpected type for value (got "{0}")'.format(type(value)))

        newvallist = []

        # loop over the elements of the value string
        for elem in vallist:
            # convenience block to try each conversion.
            # Start with the most stringent type (int),
            # and progress towards the ost relaxed (string, no conversion needed)
            while True:
                try:
                    newval = int(elem)
                    break
                except ValueError:
                    pass

                try:
                    # the float interpreter is C-based, and does not like 1d0 as an input
                    newval = float(elem.lower().replace('d', 'e'))
                    break
                except ValueError:
                    pass

                newval = elem
                break
            # Add the converted value to the output list
            newvallist.append(newval)
        return newvallist

    def append(self, key, value, comment):
        '''Append a value and a comment to the item

@param key: The key is included in the call for sanity checking. If the provided
            key is not equal to our stored key, then we fail.
@param value: The new value. This value is added with the L{appendvalue} method.
@param comment: The new or additional comment for the key/value pair. This
            comment is added with the L{appendcomment} method.
'''
        # check that we refer to the same key
        if key != self._key:
            raise ValueError(format('The given key ({0}) does not match the stored key ({1})',
                             key, self.key))
        self.appendvalue(value)
        self.appendcomment(comment)

    def appendvalue(self, value):
        '''Append a value to the item.

The value is first converted with the L{convertvalue} method.

@param value: The raw (string) input.
'''
        newval = self.convertvalue(value)
        if self._value is None:
            self._value = []
        self._value.append(newval)

    def appendcomment(self, comment):
        '''Append a comment to the list of comments.

A value is always appended to the list of comments, even if it is just the empty string.

@param comment: The new comment (can be C{None}).
'''
        if comment is None or comment == 'None':
            comment = ''
        if self._comment is None:
            self._comment = []
        if type(comment) == type(""):
            comment = [comment]

        if type(comment) == type([]):
            self._comment = comment[:]
            new_comment_length = len(self._value)
            old_comment_length = len(self._comment)
            if old_comment_length < new_comment_length:
                for dummy in range(new_comment_length - old_comment_length):
                    self._comment.append('')
            else:
                self._comment = self._comment[0:new_comment_length]

    def comments(self):
        '''Retrieve the comments for this item.

@return: The list of comments (shallow copy of the internal list).
@rtype: C{[]}
'''
        if len(self._comment) > 1:
            return self._comment[:]
        else:
            return self._comment[0]

    def values(self):
        '''Retrieve the values for this item.

The values are always stored as a list of lists. Here we unpack the
items as far as possible.

@return: The value of the item.
@rtype: atom, list, or list of lists.
'''
        if len(self._value) == 1:
            if len(self._value[0]) == 1:
                return self._value[0][0]
            else:
                return self._value[0]
        else:
            return self._value[:]

    value = values

    def rawvalues(self):
        '''get the raw value. Makes it easier to loop over multiple bands.

Returns a shallow copy of the values list.

@return: values.
@rtype: C{[[]]}
'''
        return self._value[:]

    def set_rawvalue(self, newvalue):
        '''Set a raw value.

Ensure that the comments list has the same length as the values list. Otherwise
hardly any checks are performed on the input parameter. Caveat emptor!

@param newvalue: The new (raw) value to use.
@type newvalue: C{[[]]}
'''
        self._value = newvalue
        if len(self._value) != len(self._comment):
            new_comment_length = len(self._value)
            old_comment_length = len(self._comment) if self._comment is not None else 0
            if old_comment_length < new_comment_length:
                for dummy in range(new_comment_length - old_comment_length):
                    self.appendcomment('')
            else:
                self._comment = self._comment[0:new_comment_length]

    def setvalue(self, newvalue):
        '''Change the value for this key.

Ensure that the comments list has the same length as the values list. The values
are checked to ensure that they are packaged eventually as a list of lists of atoms.
A new 'raw' value is produced, and L{set_rawvalue} is used to update the value of C{self}.

@param newvalue: The new (raw) value to use.
@type newvalue: C{[[]]}
'''

        if type(newvalue) != type([]):
            rawvalue = [[newvalue]]
        else:
            if type(newvalue[0]) != type([]):
                rawvalue = [newvalue]
            else:
                rawvalue = newvalue
        self.set_rawvalue(rawvalue)


    def key(self):
        '''Retrieve the key for this item

@return: key
@rtype: C{string}
'''
        return self._key

    def extend(self, valstr, comment):
        '''Continue a configuration line, and extend the last item.

The format of the configuration file allows for continued lines. In that case the
last item of the values list must be extended, rather than a new list appended.

@param valstr: A raw string value.
@type valstr: C{string}
@param comment: Comment for this value
@type comment: C{string}
'''
        values = self.convertvalue(valstr)
        self._value[-1].extend(values)
        if comment:
            if self._comment[-1]:
                self._comment[-1] = self._comment[-1] + '; ' + comment
            else:
                self._comment[-1] = comment

    def __repr__(self):
        '''Create a representation of the values of this item.

The representation is suitable recreating self in Python.
'''
        return '{cls}(key={key!r}, \n\tvalue={val!r}, \n\tcomment={comment!r})'.format(
                cls=self.__class__.__name__,
                key=self.key(),
                val=list(self.values()),
                comment=self.comments())

    def __str__(self):
        '''Produce a string representation in the correct syntax for the configuration file format.

The string output is suitable for generating an item in a configuration file.
'''
        out = []
        for val, comment in zip(self._value, self._comment):
            if comment:
                out.append('{0}  {1}  ({2})'.format(
                           self._key, ' '.join([dblformat(item) for
                           item in val]), comment))
            else:
                out.append('{0}  {1}'.format(self._key, ' '.join([dblformat(item) for
                                             item in val])))
        return '\n'.join(out)


class RT_configuration(object):
    '''Radiative transfer configuration manager and file parser

This is the main entrypoint for this module.'''

    def __init__(self, filename=None, config=None, comment=None, debug=False):
        '''Initializer

Parameters:
        filename:    optional, configuration filename.
        config:  optional, existing RTParser object.
        comment: optional, a comment at the filename-level.'''
        self._file = filename
        if config is not None and type(config) == type(self):
            self.config = config.dictionary()
        else:
            self.config = RT_dict(comment=comment)

        self.logger = utilities.logger()
        if debug:
            self.logger.setLevel(10)
        self.logger.debug("Initializing {0}".format(self.__class__.__name__))

        if not (self._file is None or (config is not None and type(config) == type(self))):
            self.logger.info("Reading config from {0}.".format(self._file))
            self.read()
        self.logger.debug("Created {0}".format(self.__class__.__name__))

    def dictionary(self):
        '''Return a copy of the underlying dictionary'''
        return self.config.copy()

    def file(self):
        '''Retrieve the file that the object is currently linked to.'''
        return self._file

    def setfile(self, filename):
        '''Set the file that will be used for reading or writing.'''
        self._file = filename

    def sections(self):
        '''Obtain a list of sections'''
        sections = list(self.config.keys())
        self.logger.debug("Sections in object: '{0}'.".format("', '".join(sections)))
        return sections

    def subsections(self, section):
        '''Obtain a list of subsections'''
        subsections = list(self[section].keys())
        self.logger.debug("Subsections in section '{1}': '{0}'.".format("', '".join(subsections), section))
        return subsections

    def keys(self, section, subsection):
        '''Obtain a list of keys in a particular subsection'''
        keys = list(self[section, subsection].keys())
        self.logger.debug("Keys in section '{1}', subsection {2}: '{0}'.".format("', '".join(keys), section, subsection))
        return keys

    def items(self, section, subsection):
        '''Obtain a list of key-value pairs for a particular subsection'''
        items = list(self[section, subsection].items())
        self.logger.debug("Items in section '{1}', subsection {2}: {{{0}}}.".format(
                ", ".join(["'{0}': '{1}'".format(k,v) for k,v in items]), section, subsection))
        return items

    def values(self, section, subsection):
        '''Return the values in the subsection'''
        self.logger.debug("Values for SECTION '{0}', subsection '{1}'".format(section, subsection))
        val = []
        for key in list(self[section, subsection].keys()):
            val.append(self[section, subsection, key])
        return val

    def __contains__(self, item):
        rval = 0
        if type(item) == type(''):
            rval = 1 if item in self.config else 0
        elif len(item) == 1:
            rval = 1 if item[0] in self.config else 0
        elif len(item) == 2:
            rval = 1 if item[0] in self.config and item[1] in self.config[item[0]] else 0
        else:
            rval =  1 if ((item[0] in self.config) and
                (item[1] in self.config[item[0]]) and
                (item[2] in self.config[item[0]][item[1]])) else 0
        self.logger.debug('__contains__({0}): {1}'.format(item, 'True' if rval else 'False'))
        return rval

    def __delitem__(self, key):
        '''Remove sections, subsections, and keys'''
        self.logger.debug('__delitem__: {0}'.format(key))
        if type(key) == type(''):
            del self.config[key]
        elif len(key) == 1:
            del self.config[key[0]]
        elif len(key) == 2:
            del self.config[key[0]][key[1]]
        else:
            del self.config[key[0]][key[1]][key[2]]

    def __getitem__(self, key):
        '''Get an item from the configuration dictionary.

Key can be a string or a tuple.
If value is a tuple, then it is assumed that the real value is the first element,
and the second element contains a comment'''
        self.logger.debug('__getitem__: {0}'.format(key))
        if isinstance(key, string_types):
            return self.config[key]
        elif len(key) == 1:
            return self.config[key[0]]
        elif len(key) == 2:
            return self.config[key[0]][key[1]]
        else:
            return self.config[key[0]][key[1]][key[2]]

    def __setitem__(self, key, value):  #IGNORE:R0915 #IGNORE:R0912
        '''Insert an item into the configuration.

Key can be a string (insert a whole section), or a tuple of strings
(insert a subsection, or a key into a subsection). If a section is already
there, then new subsections are merged, same for subsections and keys.

If a key is already there, then the new items extend the existing list,
or the single item is turned into a list, and the new value is appended.'''
        self.logger.debug('__setitem__: [{0}] = {1}'.format(key, value))
        # split the value into a value and a comment
        if type(value) == type(tuple()):
            comment = str(value[1])
            value = value[0]
        else:
            comment = ''

        # handle string or single item tuple values as key.
        # The value must be a dictionary describing the whole section
        if type(key) == type(tuple()) and len(key) == 1:
            key = key[0]

        if type(key) == type(''):
            if type(value) == type({}) or type(value) == type(RT_dict()):
                self.logger.debug("Inserting whole section at '{0}'.".format(key))
                if key in self.config:
                    self.config[key].update(value)
                else:
                    if value:
                        self.config[key] = RT_dict(value, comment=comment)
                    else:
                        self.config[key] = RT_dict(comment=comment)
                return
            else:
                self.logger.error('Value must be a dictionary if used in this way')
                raise ValueError('Value must be a dictionary if used in this way')

        # Key is a tuple with two items. Insert a whole subsection at once.
        if len(key) == 2:
            self.logger.debug("Inserting whole section at '{0[0]}.{0[1]}'.".format(key))
            if type(value) == type({}) or type(value) == type(RT_dict()):
                if key[0] not in self.config:
                    self.config[key[0]] = RT_dict()
                if key[1] in self.config[key[0]]:
                    self.config[key[0]][key[1]].update(value)
                else:
                    if value:
                        self.config[key[0]][key[1]] = RT_dict(value, comment=comment)
                    else:
                        self.config[key[0]][key[1]] = RT_dict(comment=comment)
                return
            else:
                self.logger.error('Value must be a dictionary if used in this way')
                raise ValueError('Value must be a dictionary if used in this way')

        # Key is a tuple with three items. This fully qualifies a single key item
        if len(key) == 3 and (type(value) != type({})
                         and type(value) != type(RT_dict())):
            self.logger.debug("Inserting item at '{0[0]}.{0[1]}.{0[2]}'.".format(key))
            if key[0] not in self.config:
                self.logger.debug("Creating section '{0[0]}'.".format(key))
                self.config[key[0]] = RT_dict()
            if key[1] not in self.config[key[0]]:
                self.logger.debug("Creating subsection '{0[0]}.{0[1]}'.".format(key))
                self.config[key[0]][key[1]] = RT_dict()
            if key[2] in self.config[key[0]][key[1]]:
                self.logger.debug("Appending to value '{0[2]}'.".format(key))
                self.config[key[0]][key[1]][key[2]].append(key[2], value, comment)
            else:
                self.logger.debug("Creating item '{0[0]}.{0[1]}.{0[2]}'.".format(key))
                self.config[key[0]][key[1]][key[2]] = RT_item(key=key[2],
                                                              value=value, comment=comment)
        else:
            self.logger.error('Value can not be a dictionary if used in this way')
            raise ValueError('Value can not be a dictionary if used in this way')

    def __str__(self):
        '''Create a string representation for writing to file.'''
        out = []

        self.logger.debug("Entering {0}.{1}".format(self.__class__.__name__, utilities.whoami()))
        file_comment = self.config.blockcomment()
        if file_comment:
            out.append('# {0}\n'.format(file_comment))

        # obtain the sections
        sections = self.sections()
        self.logger.debug(", ".join(sections))

        # loop over the sections
        for section in sections:
            # write the section header
            sect = self[section]
            comment = sect.comment()
            if comment:
                out.append('\n\nSECTION {0} ({1})\n'.format(section.upper(), comment))
            else:
                out.append('\n\nSECTION {0}\n'.format(section.upper()))

            comment = sect.blockcomment()
            if comment:
                out.append('# {0}\n'.format(comment))

            # obtain the subsections and put them in alphabetical order
            subsections = self.subsections(section)

            # loop over the subsections
            for subsection in subsections:
                # write the subsection header
                subsect = self[section, subsection]
                comment = subsect.comment()
                if comment:
                    out.append('\nsubsection {0} ({1})\n'.format(subsection, comment))
                else:
                    out.append('\nsubsection {0}\n'.format(subsection))
                comment = subsect.blockcomment()
                if comment:
                    out.append('# {0}\n'.format(comment))
                # loop over the keys in this subsection
                for key in list(self[section, subsection].keys()):
                    out.append("{0!s}\n".format(self[section, subsection, key]))
        self.logger.debug("Exit {0}.{1}".format(self.__class__.__name__, utilities.whoami()))
        return ''.join(out)

    def write(self, filename=None):
        '''
        Write self to file.

        If filename is given, then the file attribute is not used.
        If the file attribute of the is not valid and filename is
        not given, write to stdout.
        '''
        if filename is not None:
            fp = open(filename, 'w')
        else:
            if not self._file:
                fp = sys.stdout
            else:
                # Open the file for writing
                fp = open(self._file, 'w')

        fp.write(str(self))
        fp.write('\n')

        if self._file is not None or filename is not None:
            fp.close()

    def read(self, filename=None,
             replace_section=False, #IGNORE:R0915 #IGNORE:R0912
             replace_subsection=False,
             replace_keys=False):
        '''Read configuration from file (either the class attribute or given explicitly) or stdin'''
        if hasattr(self._file, "read"):
            fp = self._file
            lines = fp.readlines()
        elif filename and os.access(filename, os.F_OK) and os.access(filename, os.R_OK):
            # Open the file for reading
            fp = open(filename, 'r')
            lines = fp.readlines()
        elif self._file and os.access(self._file, os.F_OK) and os.access(self._file, os.R_OK):
            fp = open(self._file, 'r')
            lines = fp.readlines()
        elif isinstance(filename, string_types):
            # assume not a file is given but a string
            lines = filename.split('\n')
        else:
            fp = sys.stdin
            lines = fp.readlines()

        # define and compile regular expressions
        re_control = re.compile(
            r'''^#@[ \t]*([a-zA-Z_][a-zA-Z0-9_]+)[ \t]*=[ \t]*(.+?)$''')
        re_section = re.compile(
            r'''^[ \t]*SECTION[ \t]+([a-zA-Z_][a-zA-Z0-9_-]+)[ \t]*\(?(.*?)\)?[ \t]*$''')
        re_subsection = re.compile(
            r'''^[ \t]*subsection[ \t]+([a-zA-Z_][a-zA-Z0-9_-]+)[ \t]*\(?(.*?)\)?[ \t]*$''')
        re_key = re.compile(
            r'''^[ \t]*([a-zA-Z_][a-zA-Z0-9_-]+)[ \t]+([^(&]+)(?:\((.*?)\))?[ \t]*(\&)?[ \t]*$''')  #IGNORE:C0301
        re_continue = re.compile(
            r'''^[ \t]*([^(&]+)(?:\((.*?)\))?[ \t]*(\&)?[ \t]*$''')
        re_comment = re.compile(
            r'''^[ \t]*# ?(.*)$''')
        re_empty = re.compile(
            r'''^[ \t\n\r]*$''')
        re_linebreaks = re.compile(r'''[\r\n]{1,2}''')

        # set some starting values
        section = ''
        subsection = ''
        key = ''
        next_continues = False
        linecount = 0

        # loop over all lines
        for line in lines:
            linecount += 1
            try:
                # strip away linebreaks
                line = re_linebreaks.sub('', line)
                search = re_control.search(line)
                if search:
                    self.logger.debug('Line: {0} ("{1}") matched re_control'.format(linecount, line.replace('\n', '').replace('\r', '')))
                    if search.group(1).lower() == 'replace_section':
                        replace_section = (search.group(2).lower() == 'true')
                    elif search.group(1).lower() == 'replace_subsection':
                        replace_subsection = (search.group(2).lower() == 'true')
                    elif search.group(1).lower() == 'replace_keys':
                        replace_keys = (search.group(2).lower() == 'true')

                    if section:
                        if subsection:
                            self[section, subsection].setblockcomment('@ ' +
                                 search.group(1) + ' = ' + search.group(2))
                        else:
                            self[section].setblockcomment('@ ' +
                                 search.group(1) + ' = ' + search.group(2))
                    else:
                        self.config.setblockcomment('@ ' +
                             search.group(1) + ' = ' + search.group(2))
                    continue
                search = re_empty.search(line)
                if search:
                    self.logger.debug('Line: {0} matched re_empty'.format(linecount))
                    continue
                search = re_section.search(line)
                if search:
                    self.logger.debug('Line: {0} ("{1}") matched re_section'.format(linecount, line))
                    section = search.group(1)
                    subsection = None
                    if replace_section and section in self:
                        del self[section]
                    self[section] = ({}, search.group(2))
                    continue
                search = re_subsection.search(line)
                if search:
                    self.logger.debug('Line: {0} ("{1}") matched re_subsection'.format(linecount, line))
                    subsection = search.group(1)
                    if replace_subsection and (section, subsection) in self:
                        del self[section, subsection]
                    self[section, subsection] = ({}, search.group(2))
                    continue
                search = re_comment.search(line)
                if search:
                    self.logger.debug('Line: {0} ("{1}") matched re_comment'.format(linecount, line))
                    if section:
                        if subsection:
                            self[section, subsection].setblockcomment(' ' + search.group(1))
                        else:
                            self[section].setblockcomment(' ' + search.group(1))
                    else:
                        self.config.setblockcomment(' ' + search.group(1))
                    continue
                search = re_key.search(line)
                if search:
                    self.logger.debug('Line: {0} ("{1}") matched re_key'.format(linecount, line))
                    key = search.group(1)
                    valstr = search.group(2)
                    comment = search.group(3)
                    next_continues = bool(search.group(4))
                    if replace_keys and (section, subsection, key) in self:
                        del self[section, subsection, key]
                    self[section, subsection, key] = (valstr, comment)
                    continue
                search = re_continue.search(line)
                if search and next_continues:
                    self.logger.debug('Line: {0} ("{1}")\nmatched re_continue'.format(linecount, line))
                    valstr = search.group(1)
                    comment = search.group(2)
                    next_continues = bool(search.group(3))
                    self[section, subsection, key].extend(valstr, comment)
                    continue
                raise ValueError('Encountered some unexpected syntax at line # {1}:\n\t{0}'.format(
                                 line, linecount))
            except:
                print("configline ({0:d}):\n{1}".format(linecount, line))
                raise

        if (self._file and os.access(self._file, os.F_OK) and os.access(self._file, os.R_OK) or
            filename and os.access(filename, os.F_OK) and os.access(filename, os.R_OK)):
            fp.close()


def __main():
    try:
        filename=sys.argv[1]
    except:
        sys.stderr.write("Usage: {0} FILE.in\n".format(sys.argv[0]))
        sys.exit(1)

    try:
        cfg = RT_configuration(filename=filename, debug=True)

        print(str(cfg))
    except:
        sys.stderr.write("An exception occurred\n")
        raise

    sys.exit(0)

if __name__ == "__main__":
    __main()
