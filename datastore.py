#!/usr/bin/env python3
'''
A pytables based routine for reading ASCII HDF files,
creating HDF files in the process.

When called as a program it will translate between formats.

@author: Maarten Sneep <maarten.sneep@knmi.nl>
'''

import os
import sys
import json

import tables

import numpy as np

from ahdf import AHDFr

__version = "1.0"
__date = "2017-05-09"
__author = "Maarten Sneep <maarten.sneep@knmi.nl>"

class DatastoreError(Exception):
    '''Error class for all Datastore errors'''
    pass

class DatastoreNode(object):
    '''Base class for the objects in our Datastore hierarchy.

Both groups and leafs inherit from this class, and as a result
all of those can handle attributes.'''
    def __init__(self, name=None, **kwargs):
        '''
        Constructor
        '''
        self._name = name.replace(' ', '_')
        self._parent = None

        if 'attributes' in list(kwargs.keys()):
            self._attributes = kwargs['attributes']
        else:
            self._attributes = {}

        self._debug = 'debug' in list(kwargs.keys()) and kwargs['debug']

    def name(self):
        '''Get the name of the node'''
        return self._name

    def set_name(self, name):
        '''Set the name of the node.'''
        self._name = name

    def set_parent(self, parent):
        '''Set a reference to the parent'''
        self._parent = parent

    def parent(self):
        '''Return a reference to the parent node'''
        return self._parent

    def root(self):
        '''Return a reference to the root object'''
        return self.parent().root()

    def attributes(self):
        '''Return all attributes'''
        return self._attributes

    def set_attributes(self, attrs):
        '''Set attributes to attrs. Replaces all existing attributes'''
        self._attributes = attrs

    def set_attribute(self, key, val):
        '''set a single attribute value'''
        self._attributes[key] = val

    def append_attributes(self, **kwargs):
        '''Add keyword arguments to the attributes'''
        self._attributes.update(kwargs)

    def replace_attributes(self, **kwargs):
        '''Replace existing attributes with keyword values'''
        self.set_attributes(kwargs)

    def attribute(self, key):
        '''Return the value of an attribute'''
        return self._attributes[key]

    def leaf(self):
        '''Determine if the instance is a leaf object'''
        return isinstance(self, DatastoreLeaf)


class DatastoreGroup(DatastoreNode):
    '''Base class for Datastore objects that can contain other groups or leafs'''
    def __init__(self, name=None, **kwargs):
        super(DatastoreGroup, self).__init__(name=name, **kwargs)
        self._children = []
        if 'children' in list(kwargs.keys()):
            for child in kwargs['children']:
                self[child.name()] = child

    def __iter__(self):
        return next(self)

    def __next__(self):
        '''Iterate over all children '''
        for child in self._children:
            yield self.__dict__[child]

    def children(self):
        '''Return a dict with children of this group'''
        children = {}
        for name in self._children:
            children[name] = self.__dict__[name]
        return children

    def add_child(self, child):
        '''Add a child (group or leaf) to the group'''
        self[child.name()] = child

    def items(self):
        '''Return all child items'''
        return self._children[:]

    def keys(self):
        '''return the available child names'''
        return self._children

    def __len__(self):
        return len(self._children)

    def __getitem__(self, key):
        print(key)
        return self.__dict__[key]

    def __setitem__(self, key, val):
        if not issubclass(val.__class__, DatastoreNode):
            raise DatastoreError('The value must be a subclass of DatastoreNode, ' +
                             'instead it is of class {0}'.format(val.__class__.__name__))
        if key != val.name():
            raise KeyError('The key must be equal to the name of the child')
        self.__dict__[key] = val
        self.__dict__[key].set_parent(self)
        self._children.append(key)

    def __contains__(self, item):
        return 1 if item in self._children else 0

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        else:
            raise AttributeError('Attribute {0} not found'.format(name))

    def __setattr__(self, name, val):
        if name in self.__dict__ or name[0] == '_':
            self.__dict__[name] = val
        else:
            self[name] = val

    def __delattr__(self, name):
        if name in self._children:
            del self.__dict__[name]
            self._children.remove(name)
        elif name in self.__dict__:
            del self.__dict__[name]
        else:
            raise AttributeError('Attribute {0} not found'.format(name))


class Datastore(DatastoreGroup):  # IGNORE:R0904
    '''The "root" object. Here the filehandle is maintained'''
    def __init__(self, name=None, **kwargs):
        '''Initialize a datagroup.

The 'name' keyword is used to give a name to the whole dataset. Defaults to 'root'.

If the read keyword is set to a filename, then the data in the group will
be initialized from this file.

If the write keyword is set to a filename, then the data (provided with the 'children'
keyword) are written to that file.

The group is initialized with its children. If the 'children' keyword is set to a
list of DatastoreGroup or DatastoreLeaf objects, then the data object is initialized
with it.

Attributes for the group can be supplied as a dictionary to the 'attributes' keyword.
        '''
        if name is None:
            name = 'root'
        super(Datastore, self).__init__(name=name, **kwargs)
        self._ishdf = None
        self._filename = None
        self._filehandle = None

        keys = list(kwargs.keys())
        if 'read' in keys:
            self.read(kwargs['read'])
        elif 'write' in keys:
            self.write(kwargs['write'])

    def read(self, fname):
        '''Read data from file'''
        self._filename = fname
        if os.path.exists(self._filename):
            self._ishdf = tables.is_hdf5_file(self._filename)
        else:
            raise DatastoreError("File {0} not found".format(self._filename))

        if self._ishdf:
            if self._debug:
                print('Reading from HDF-5', file=sys.stderr)
            self._filehandle = tables.open_file(self._filename, mode='r')
            self.read_hdf()
        else:
            if self._debug:
                print('Reading from ASCII', file=sys.stderr)
            self._filehandle = open(self._filename, 'r')
            self.read_ascii()
        self._filehandle.close()
        self._filehandle = None

    def write(self, fname):
        '''write data to file'''
        self._filename = fname
        if os.path.exists(self._filename):
            self._ishdf = tables.is_hdf5_file(self._filename)
        else:
            self._ishdf = (os.path.splitext(self._filename)[1] in
                           ['.h5', '.hdf5', '.hdf', '.he5'])

        if self._ishdf:
            if self._debug:
                print('Writing to HDF-5', file=sys.stderr)
            filters = tables.Filters(complevel=6, complib='zlib', fletcher32=True)
            self._filehandle = tables.open_file(self._filename, mode='w',
                                               title=self.name(), filters=filters)
            self.write_hdf()
        else:
            if self._debug:
                print('Writing to ASCII', file=sys.stderr)
            self._filehandle = open(self._filename, 'w')
            self.write_ascii()
        self._filehandle.close()
        self._filehandle = None

    def attach_file(self, file_to_attach, attributes=None): #IGNORE:R0912
        '''Attach a file to the hdf-5 file.

The file is attached to the '/attachments' group.'''
        if not self._ishdf:
            raise DatastoreError('Files can only be attached to HDF-5 files')
        if not os.path.exists(file_to_attach):
            raise ValueError('File not found: {0}'.format(file_to_attach))
        if not os.path.exists(self._filename):
            raise DatastoreError('File must exist before ' +
                                 'attaching: {0}'.format(self._filename))

        if self._debug:
            print('Attaching file "{0}" to HDF-5'.format(file_to_attach), file=sys.stderr)

        maxlen = 65535
        filters = tables.Filters(complevel=9, complib='zlib', fletcher32=True)
        atom = tables.StringAtom(itemsize=maxlen, shape=(), dflt='')

        self._filehandle = tables.open_file(self._filename, mode='a')
        try:
            gnode = self._filehandle.get_node('/attachments')
        except  tables.NoSuchNodeError:
            gnode = self._filehandle.create_group('/', 'attachments', title='file attachments')

        with open(file_to_attach, 'r') as fp:
            text = fp.read()

        numattach = 1 + (len(text) // maxlen)
        stop = len(text)

        for ii in range(numattach):
            name = os.path.basename(file_to_attach)
            name = name.replace('.', '_').replace(' ', '_')
            name = name.replace('-', '_').replace('+', '_')
            if numattach > 1:
                name = '{0}_{1:02d}'.format(name, ii + 1)
            carray = self._filehandle.create_carray(gnode, name, atom, (1,), filters=filters)

            # split on whole lines only.
            start = 0 if ii == 0 else stop
            if ii < numattach - 1:
                stop = (ii + 1) * maxlen
                while text[stop - 1] != '\n' and stop > start:
                    stop -= 1
            else:
                stop = len(text)

            carray[:] = bytes(text[start:stop], 'ascii')
            carray.flush()
            if numattach > 1:
                carray._v_attrs['Datastore_part'] = ii + 1

            if attributes:
                for k, v in list(attributes.items()):
                    if type(v) == type(''):
                        v = bytes(v, 'ascii')
                    elif type(v) == type(True):
                        v = 'true' if v else 'false'
                    elif type(v) == type([]) and len(v) > 0:
                        if type(v[0]) == type(''):
                            v = [str(e) for e in v]
                        v = np.asarray(v)
                    carray._v_attrs[k] = v



        self._filehandle.close()
        self._filehandle = None

    def get_attachment(self, attachment_name):
        '''Get an attachment from the hdf-5 file'''
        if not self._ishdf:
            raise DatastoreError('Files can only be attached to HDF files')
        if not os.path.exists(self._filename):
            raise DatastoreError('File must exist before reading ' +
                                 'attachment: {0}'.format(self._filename))

        self._filehandle = tables.open_file(self._filename, mode='r')

        try:
            gnode = self._filehandle.get_node('/attachments')
        except  tables.NoSuchNodeError:
            raise DatastoreError('No attachments were found')

        name = os.path.basename(attachment_name)
        name = name.replace('.', '_').replace(' ', '_')
        name = name.replace('-', '_').replace('+', '_')

        attributes = {}
        text = ''

        for item in gnode:
            if item._v_name.startswith(name):
                if not attributes:
                    for k in item.attrs._f_list():
                        if k != 'Datastore_part':
                            attributes[k] = item.attrs[k]
                text = text + item.read().tostring().strip('\x00')

        self._filehandle.close()
        self._filehandle = None

        return (text, attributes)

    def read_ascii(self):
        '''Read a data-object from ascii file'''
        ascii_read_object = AHDFr(fp=self._filehandle, debug=self._debug)
        dstore = ascii_read_object()
        _translate_old_datastore(dstore, self)

    def write_ascii(self):
        '''Write a data-object to ascii file'''
        _write_ascii(self._filehandle, self)

    def read_hdf(self):
        '''Read a data-object from hdf-5/PyTables file'''
        _translate_pytables(self._filehandle.root, self)

    def write_hdf(self):
        '''Write a data-object to hdf-5/PyTables file'''
        _write_hdf5(self._filehandle, self)

    def __del__(self):
        if self._filehandle is not None:
            self._filehandle.close()

    def root(self):
        '''Return a reference to the root object (endpoint of call chain)'''
        return self


class DatastoreLeaf(DatastoreNode):
    '''Leaf object in our datastore hierarchy'''
    def __init__(self, name=None, **kwargs):
        '''Leaf object.

The name of the object is given in the 'name' keyword.

Data

Attributes for the leaf can be supplied as a dictionary to the 'attributes' keyword.
       '''
        super(DatastoreLeaf, self).__init__(name=name, **kwargs)

        self._data = None
        if 'data' in list(kwargs.keys()):
            self._data = kwargs['data']

    def __getitem__(self, key):
        return self._data.__getitem__(key)

    def __setitem__(self, key, val):
        self._data.__setitem__(key, val)

    def __setattr__(self, name, val):
        if name in self.__dict__ or name[0] == '_':
            self.__dict__[name] = val
        else:
            self._data.__setattr__(name, val)

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        else:
            return self._data.__getattribute__(name)

    def __delattr__(self, name):
        self._data.__delattr__(name)

    def data(self):
        '''Return the data object (numpy array)'''
        return self._data

    def set_data(self, newdata):
        '''Change the data array'''
        self._data = newdata

    def keys(self): #IGNORE:R0201
        '''Instead of letting __getattrs__ fail,
we raise a warning here to indicate that we are a leaf'''
        raise UserWarning("Arrived at leaf of object tree")

def read(filename, debug=False):
    '''Read a structured ASCII or HDF-5 input file'''
    root = Datastore(name='root', read=filename, debug=debug)
    return root

def write(filename, obj):
    '''Write a Datastore object to file'''
    obj.write(filename)

def translate(infile, outfile, attachment=None, attributes=None, debug=False):
    '''Translate a file (asciiHDF or HDF-5) into another format (asciiHDF or HDF-5).

Types are determined from file or filename.'''
    d = read(infile, debug=debug)
    write(outfile, d)
    if attachment is not None and d._ishdf:  # IGNORE:W0212
        if type(attachment) == type(""):
            attachment = [attachment]
        for item in attachment:
            if os.path.exists(item) and os.path.isfile(item):
                d.attach_file(item, attributes=attributes)
            else:
                print("Item '{0}' could not be attached because it doesn't exist or is a directory.".format(item))

def _translate_pytables(node, dest):
    '''Translate a pytables node into a hierarchy of DatastoreNodes'''
    for name in node._v_attrs._f_list(): # IGNORE:W0212
        dest.set_attribute(name, node._v_attrs[name]) # IGNORE:W0212
    for item in node:
        if issubclass(item.__class__, tables.Leaf):
            attrs = {}
            for name in item._v_attrs._f_list(): # IGNORE:W0212
                attrs[name] = item._v_attrs[name] # IGNORE:W0212
            child = DatastoreLeaf(name=item._v_name, # IGNORE:W0212
                                  data=item.read(),
                                  attributes=attrs)
            dest[child.name()] = child
        else:
            child = DatastoreGroup(name=item._v_name) # IGNORE:W0212
            dest[child.name()] = _translate_pytables(item, child)
    return dest

def _translate_old_datastore(dstore, dest):
    '''Translate an old Datastore object into a hierarchy of DatastoreNodes'''
    for key, val in list(dstore.items()):
        dest.set_attribute(key, val)
    for item in dstore.data():
        if item.leaf():
            child = DatastoreLeaf(name=item.name(),
                                  data=item.data(),
                                  attributes=item.copy())
            dest[child.name()] = child
        else:
            child = DatastoreGroup(name=item.name())
            dest[child.name()] = _translate_old_datastore(item, child)
    return dest

class _write_hdf5(object):
    '''Write data to an HDF-5 file'''
    def __init__(self, fp=None, datastore=None):
        self.datastore = datastore
        self.fp = fp
        self.write_group(self.datastore, fp.root)

    def write_group(self, datastore, parent):
        '''Write datastore as a group.'''
        if not datastore.leaf():
            if datastore.name() == 'root':
                group = self.fp.root
            else:
                group = self.fp.create_group(parent, datastore.name(), datastore.name())
            self.write_attributes(datastore, group)
            self.write_data(datastore, group)

    def write_attributes(self, datastore, parent):  #IGNORE:R0201
        '''Write the key-value pairs in datastore.attributes() to an attributes block'''
        attrs = datastore.attributes()
        for k, v in list(attrs.items()):
            # make sure that attribute names do not contain spaces.
            if ' ' in k:
                k = k.replace(' ', '_')

            # do some type cleanup
            if type(v) == type(''):
                v = str(v)
            elif type(v) == type(True):
                v = 'true' if v else 'false'
            elif type(v) == type([]):
                if type(v[0]) == type(''):
                    v = [str(e) for e in v]
                v = np.asarray(v)

            parent._v_attrs[k] = v  #IGNORE:W0212

    def write_data(self, datastore, parent):
        '''Write the data in datastore.

The datastore can be a group (in which case we dispatch to write_group)
or a leaf, and then we write the array.
'''
        if datastore.leaf():
            self.write_array(datastore, parent)
        else:
            for item in datastore:
                if item.leaf():
                    self.write_array(item, parent)
                else:
                    self.write_group(item, parent)

    def write_array(self, datastore, parent):
        '''Write the array to the file'''
        if "U" in str(datastore.data().dtype):
            atom = tables.Atom.from_dtype(np.dtype(str(datastore.data().dtype).replace('U', 'S')))
        else:
            atom = tables.Atom.from_dtype(datastore.data().dtype)
        shape = datastore.data().shape

        if 0 in shape:
            return

        try:
            ca = self.fp.create_carray(parent, datastore.name(), atom, shape)
        except tables.exceptions.NodeError as err:
            print("Warning: {0}".format(err))
            ca = self.fp.get_node(parent, datastore.name())

        ca[...] = datastore.data()

        self.write_attributes(datastore, ca)

class _write_ascii(object):
    '''Write data to an ASCII file'''
    def __init__(self, fp=None, datastore=None):
        self.datastore = datastore
        self.fp = fp
        self.write_group(self.datastore)

    def write_group(self, datastore):
        '''Write datastore as a group.'''
        if not datastore.leaf():
            self.fp.write("BeginGroup({0})\n".format(datastore.name()))
            self.write_attributes(datastore)
            self.write_data(datastore)
            self.fp.write("EndGroup\n")

    def write_attributes(self, datastore):
        '''Write the key-value pairs in datastore to an attributes block'''
        attrs = datastore.attributes()
        if len(attrs) > 0:
            self.fp.write("BeginAttributes\n")
            for key, val in list(attrs.items()):
                self.fp.write("{key} = {value}\n".format(key=key,
                              value=json.dumps(val)))
            self.fp.write("EndAttributes\n")

    def write_data(self, datastore):
        '''Write the data in datastore.

The datastore can be a group (in which case we dispatch to write_group)
or a leaf, and then we write the array.
'''
        if datastore.leaf():
            self.write_array(datastore)
        else:
            for item in datastore:
                if item.leaf():
                    self.write_array(item)
                else:
                    self.write_group(item)

    def write_array(self, data):
        '''Write the array to the file.

The attributes in the datastore are added to the array block. The order (only C),
number of dimensions, and dimensions sizes are determined from the array object.
The numpy tofile method is used to add the data itself.
'''
        if (str(data.data().dtype)[0:2].lower() != '|s'):
            typedescr = str(data.data().dtype)
        else:
            typedescr = 'String'

        self.fp.write("BeginArray({name}, {type})\n".format(name=data.name(),
                                                            type=typedescr))
        self.write_attributes(data)
        self.fp.write("Order = C\n")
        self.fp.write("NumDimensions = {ndims}\n".format(ndims=len(data.data().shape)))
        self.fp.write("Size = {0}\n".format(", ".join([str(d) for d in data.data().shape])))
        if typedescr != 'String':
            data.data().tofile(self.fp, sep=' ')
        else:
            data.data().tofile(self.fp, sep='\n')
        self.fp.write("\nEndArray\n")

def __main():
    '''Commandline handler'''
    from optparse import OptionParser
    usage = """Translate an ascii representation of data into HDF-5 (or vice versa).

Usage: %prog [options] infile outfile"""
    parser = OptionParser(usage=usage, version="%prog {0}".format(__version))

    parser.add_option("-i", "--input", dest="infile",
                      help="Read data from FILE", metavar='FILE')
    parser.add_option("-o", "--output", dest="outfile",
                      help="Write data to FILE", metavar='FILE')
    parser.add_option("-a", "--attach", dest="attach", action="append", default=[],
                      help="Attach FILE to the output (only if output is HDF-5)",
                      metavar="FILE")
    parser.add_option("-A", "--attribute", dest="attrs",
                      help="Add attribute to the attached file",
                      metavar="KEY=VAL", action="append")
    parser.add_option("-D", "--debug", dest="debug",
                      help="Add debug information",
                      action="store_true", default=False)

    (options, args) = parser.parse_args()

    if len(args) == 2:
        infile, outfile = tuple(args)
    else:
        infile, outfile = options.infile, options.outfile

    if not (infile and outfile):
        parser.error("Need to know input and output files.")

    if not os.path.exists(infile):
        parser.error("The input file must exist.")

    attach = None
    attrs = {}
    if options.attach:
        for item in options.attach:
            if not os.path.exists(item):
                parser.error("The file to be attached does not exist.")

        attach = options.attach
        if options.attrs:
            for itm in options.attrs:
                lst = itm.split("=")
                if len(lst) != 2:
                    sys.stderr.write("Ignoring attribute '{0}'\n".format(itm))
                    continue
                key = lst[0].strip()
                val = lst[1].strip()
                try:
                    cval = int(val)
                except ValueError:
                    try:
                        cval = float(val)
                    except ValueError:
                        cval = str(val)
                attrs[key] = cval

    translate(infile, outfile, attachment=attach, attributes=attrs, debug=options.debug)

if __name__ == "__main__":
    __main()
