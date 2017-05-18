#!/usr/bin/env python3 -tt

'''
This module can be used for data storage, and reading and
writing to a simple ascii format.
'''

import re
import sys
import json
import os.path

import numpy as np

__date = "2010-02-26"
__version = "1.0"
__author = "Maarten Sneep"


class _Datastore(dict): # IGNORE:R0904
    '''Data storage object to hold data (groups, arrays) and attributes.

The data is stored as some object:
  - list of _Datastore objects for a group.
  - A numpy array for a leaf object.
the attributes are stored in the inherited dict.
    '''

    def __init__(self, name=None, data=None, attrs=None):
        """Initializer

Copy the data object, and put the attributes in the inherited
dictionary. Since we inherit from dict, all normal methods to handle
key-value pairs are available.
"""
        if not attrs:
            attrs = {}
        dict.__init__(self, attrs)
        self._name = name
        self._data = data
        self._leaf = None
        self.set_leaf()

    def set_leaf(self):
        '''Boolean to indicate whether this is a leaf or a group.

This is automatically determined from the type of data.
- List => group
- Array => leaf
'''
        if self.data is None:
            self._leaf = None
        else:
            if type(self._data) == type([]):
                self._leaf = False
            else:
                self._leaf = True

    def leaf(self):
        '''Is this a leaf or a group?'''
        self.set_leaf()
        return self._leaf

    def attributes(self):
        """Return attributes for the data object (i.e. return self)"""
        return self

    def set_attributes(self, attrs):
        """Replace all attributes"""
        self.clear()
        self.update(attrs)

    def data(self):
        """Return the data object (numpy array or list)"""
        return self._data

    def set_data(self, data):
        """Replace the data object (numpy array or list)"""
        if self._data:
            del self._data
        self._data = data
        self.set_leaf()

    def item(self, aname, avalue):
        """Return the _Datastore objects where attribute 'aname' is set to 'avalue'.

Return value is a tuple of found objects.
"""
        rval = []
        if aname in list(self.keys()) and self[aname] == avalue:
            rval.append(self)

        for obj in self.data():
            if obj.leaf() and aname in list(obj.keys()) and obj[aname] == avalue:
                rval.append(obj)
            else:
                rval.extend(obj.item(aname, avalue))

        return tuple(rval)

    def append(self, data):
        """Add a _Datastore object to the group"""
        if self.leaf():
            raise ValueError("Cannot append to a leaf object")
        if type(data) == type([]):
            self._data.extend(data)
        else:
            self._data.append(data)

    def name(self):
        return self._name

class AHDFr(object):
    """Object to read ASCII HDF representation"""
    def __init__(self, fname=None, fp=None, root=None, debug=False):
        self.re = {}
        self.re['begingroup'] = re.compile(r"^\s*BeginGroup\(([^)]+)\)")
        self.re['endgroup'] = re.compile(r"^\s*EndGroup")
        self.re['beginattrs'] = re.compile(r"^\s*BeginAttributes")
        self.re['attr'] = re.compile(r"^\s*([a-zA-Z][^=]*?)=(.+?)$")
        self.re['endattrs'] = re.compile(r"^\s*EndAttributes")
        self.re['beginarray'] = re.compile(r"^\s*BeginArray\(([^,]+),\s*([^)]+)\)")
        self.re['endarray'] = re.compile(r"^\s*EndArray")
        self.re['comment'] = re.compile(r"^\s*#")
        self.re['order'] = re.compile(r"^\s*Order\s*=\s*(\w+)")
        self.re['ndim'] = re.compile(r"^\s*NumDimensions\s*=\s*(\d+)")
        self.size = lambda n: re.compile(r"^\s*Size\s*=\s*{0}".format(
                                         r"\s*[, ]\s*".join([r"(\d+)"] * n)))
        self.fp = None

        if root is None:
            self.root = _Datastore(name="/", data=[], attrs={})
        else:
            self.root = root

        if fname is not None:
            self.fp = open(fname, "r")
        if fp is not None and self.fp is None:
            self.fp = fp
        self._debug = debug
        self.props = {'continue_scan':True}

    def __del__(self):
        if self.root.name() == "/":
            self.fp.close()

    def __call__(self):
        '''Read a group from the file.

The search routine will reset the continue_scan property when the group is closed.
'''
        self.props['continue_scan'] = True
        while self.props['continue_scan']:
            line = self.fp.readline()
            if not line:
                self.props['continue_scan'] = False
            else:
                for key in ['begingroup', 'endgroup',
                            'beginattrs', 'endattrs',
                            'beginarray', 'endarray',
                            'comment']:
                    search = self.re[key].search(line)
                    if search:
                        if self._debug:
                            print('Found key: {0}'.format(key), file=sys.stderr)
                        f = eval("self." + key)
                        f(search)
                        break

        return self.root

    def begingroup(self, search):
        '''Handle the BeginGroup statement in the data file.

If the name is not "/", then a new AHDF object is created,
and the group is read in this new object.
'''
        self.props['group'] = True

        name = search.group(1)
        if self._debug:
            print('Found group: {0}'.format(name), file=sys.stderr)
        if not (name == "/" or name == ""):
            basename = os.path.basename(name)
            obj = AHDFr(fp=self.fp, root=_Datastore(name=basename, data=list()), debug=self._debug)
            self.root.append(obj())

    def endgroup(self, dummy):
        '''Handle the EndGroup directive

Set continue_scan property to False to end scanning the current group.
'''
        self.props['continue_scan'] = False

    def beginattrs(self, dummy, array=False):
        '''Handle the BeginAttributes directive'''
        attrs = {}
        if self._debug:
            print('Reading {0} attributes'.format('array' if array else 'group'), file=sys.stderr)
        while True:
            before = self.fp.tell()
            line = self.fp.readline()
            if len(line) == 0 or line == '\n' or self.re["endattrs"].search(line):
                self.fp.seek(before)
                break
            else:
                search = self.re["attr"].search(line)
                if search is None and self._debug:
                    print('Attribute search failed: {0}'.format(line[0:-1]))
                try:
                    attrs[search.group(1).strip()] = json.loads(search.group(2).strip())
                except:
                    print('Encountered error while parsing')
                    print(line)
                    raise

        if array:
            self.root.data()[-1].set_attributes(attrs)
        else:
            self.root.update(attrs)

    def attr(self, dummy):# IGNORE:R0201
        '''read attribute. Completely handled in self.beginattrs()'''
        raise ValueError("Unexpected call")

    def endattrs(self, dummy):
        '''End the attributes directive (No-op)'''
        pass

    def beginarray(self, search):
        '''Start reading a new array.

This method handles the array header, and records the name and type of the array.
The data itself is read with the readarraydata() method.
'''
        self.props['name'] = search.group(1).strip()
        self.props['type'] = search.group(2).strip()
        self.props['group'] = False
        if self._debug:
            print('Start reading array "{0}"'.format(search.group(1).strip()), file=sys.stderr)

        self.root.append(_Datastore(name=self.props['name'],
                                    data=np.zeros((0,), dtype='float64')))
        self.readarraydata()

    def endarray(self, dummy):
        '''Remove the properties to prepare the object for the next array'''
        pass

    def readarraydata(self):
        '''Read the data.

Start with the header (Order, NumDimensions, and Size keywords),
and finally read the data itself.
'''

        if self.props['group']:
            if self._debug:
                print('reading array aborted', file=sys.stderr)
            return

        # read array attributes
        line = self.fp.readline()
        search = self.re['beginattrs'].search(line)
        if search:
            self.beginattrs(line, array=True)
            line = self.fp.readline()
            self.endattrs(line)
            line = self.fp.readline()

        # read the order (Fortran or C)
        search = self.re["order"].search(line)
        if not search:
            raise ValueError("Unexpected key:\n\t'{0}'\n\tExpected 'Order'".format(line))
        self.props['reverse'] = (search.group(1).lower() == "fortran")
        if self._debug:
            print('Array order "{0}"'.format(search.group(1)), file=sys.stderr)

        # read the number of dimensions
        line = self.fp.readline()
        search = self.re["ndim"].search(line)
        if not search:
            raise ValueError("Unexpected key:\n\t'{0}'"
                             "\n\tExpected 'NumDimensions'".format(line))
        self.props['ndim'] = int(search.group(1))
        if self._debug:
            print('Number of dimensions "{0}"'.format(int(search.group(1))), file=sys.stderr)

        # read the dimension sizes
        line = self.fp.readline()
        search = self.size(self.props['ndim']).search(line)
        if not search:
            raise ValueError("Unexpected key:\n\t'{0}'\n\tExpected 'Size'".format(line))
        self.props['dimensions'] = np.asarray([int(v) for v in search.groups()])

        dims = self.props['dimensions']
        if self._debug:
            print('Dimensions "{0}"'.format(dims.tolist()), file=sys.stderr)

        # String arrays cannot be read directly by numpy,
        # we have to do this with an intermediate list.
        if self.props['type'].lower() == "string":
            array = []
            for dummy in range(dims.prod()):
                array.append(self.fp.readline().strip())

            if self.props['ndim'] > 1:
                if self.props['reverse']:
                    array = np.asarray(array).reshape(dims[::-1])
                    array.transpose()
                else:
                    array = np.asarray(array).reshape(dims)
            else:
                array = np.asarray(array)
        else:
            array = np.fromfile(file=self.fp, dtype=np.dtype(self.props['type']),
                                    count=dims.prod(),
                                    sep=' ')

            if len(array) != dims.prod():
                if self._debug:
                    print('Array {0} has incorrect size.'.format(self.props['name']), file=sys.stderr)
                    oldsize = len(array)
                    array = np.resize(array, (dims.prod(),))
                    array[oldsize:] = 0.0
                    # finish reading the line
                    line = self.fp.readline()
                    print('Rest of the line:\n"{0}"'.format(line), file=sys.stderr)
                    items = line.split()
                    cnt = -1
                    for item in items[::-1]:
                        try:
                            array[cnt] = float(item)
                        except ValueError:
                            array[cnt] = np.nan
                        cnt -= 1
                else:
                    raise ValueError('Unexpected array size for {0}'.format(self.props['name']))

            if self.props['ndim'] > 1:
                if self.props['reverse']:
                    array = np.reshape(array, dims[::-1])
                    array.transpose()
                else:
                    array = np.reshape(array, dims)

        self.root.data()[-1].set_data(array)

    def comment(self, dummy):
        '''Simply ignore comments'''
        pass

    def order(self, dummy): # IGNORE:R0201
        '''This method should not be called directly'''
        raise ValueError("Unexpected call")

    def ndim(self, dummy): # IGNORE:R0201
        '''This method should not be called directly'''
        raise ValueError("Unexpected call")

class AHDFw(object):
    '''Write data to an ASCII file'''
    def __init__(self, fname=None, datastore=None):
        self.datastore = datastore
        if fname is not None:
            self.fp = open(fname, 'w')
        else:
            self.fp = None

    def __call__(self, datastore=None, fname=None):
        if datastore is not None:
            if self.datastore is not None:
                del self.datastore
            self.datastore = datastore

        if fname is not None:
            if self.fp is not None:
                self.fp.close()
            self.fp = open(fname, 'w')

        if self.datastore.leaf():
            datastore = _Datastore(name="/", data=[self.datastore])
            self.datastore = datastore

        self.write_group(self.datastore)

    def __del__(self):
        self.fp.close()
        self.fp = None

    def write_group(self, datastore):
        '''Write datastore as a group.'''
        if not datastore.leaf():
            self.fp.write("BeginGroup({0})\n".format(datastore.name()))
            self.write_attributes(datastore)
            self.write_data(datastore)
            self.fp.write("EndGroup\n")

    def write_attributes(self, datastore):
        '''Write the key-value pairs in datastore to an attributes block'''
        names = [k for k in list(datastore.keys()) if k[0:2] != '__']
        if len(names) > 0:
            self.fp.write("BeginAttributes\n")
            for k in names:
                self.fp.write("{key} = {value}\n".format(key=k,
                              value=json.dumps(datastore[k])))
            self.fp.write("EndAttributes\n")

    def write_data(self, datastore):
        '''Write the data in datastore.

The datastore can be a group (in which case we dispatch to write_group)
or a leaf, and then we write the array.
'''
        if datastore.leaf():
            self.write_array(datastore)
        else:
            for item in datastore.data():
                if item.leaf():
                    self.write_array(item)
                else:
                    self.write_group(item)

    def write_array(self, data):
        '''Write the array to the file.

The attributes in the datastore are added to the array block. The order (only C),
number of dimensions, and dimensions sizes are determined from the array objet.
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


def read(fname):
    '''Read an ASCII representation, and return the datastore object.'''
    fileobject = AHDFr(fname=fname)
    data_store = fileobject()
    del fileobject
    return data_store

def write(datastore, fname):
    '''Write the datastore object to the file "fname".'''
    fileobject = AHDFw()
    fileobject(datastore, fname)

if __name__ == "__main__":
    root = _Datastore(name="/", data=[])
    group = _Datastore(name='group', data=[], attrs={'group':'subset'})
    labels = _Datastore(name='labels', data=[])

    attrs = {'Synthetic': True,
             'label_array': ['/labels/labels_x', '/labels/labels_y', '/labels/labels_z']}

    data = np.asarray(np.arange(30, 0, -1), dtype='int32').reshape((2, 3, 5))
    root.append(_Datastore(name="arange_data", data=data, attrs=attrs))

    data = 100.0 * np.asarray(np.random.random((2, 3, 5)), dtype='float64')
    root.append(_Datastore(name='random_data', data=data, attrs=attrs))

    data = np.asarray(np.arange(0, 30, 1), dtype='int8').reshape((2, 3, 5))
    group.append(_Datastore(name="arange_data", data=data, attrs=attrs))

    data = 10.0 * np.asarray(np.random.random((2, 3, 5)), dtype='float32')
    group.append(_Datastore(name='random_data', data=data, attrs=attrs))
    root.append(group)

    labels.append(_Datastore(name='labels_x',
                            data=np.asarray(['blub', 'blah']),
                            attrs={'labels_for_axis':1}))
    labels.append(_Datastore(name='labels_y',
                            data=np.asarray(['blub', 'blah', 'foo']),
                            attrs={'labels_for_axis':2}))
    labels.append(_Datastore(name='labels_z',
                            data=np.asarray(['blub', 'blah', 'foo', 'bar', 'baz']),
                            attrs={'labels_for_axis':3}))
    root.append(labels)

    write(root, "datastore.ahdf")
