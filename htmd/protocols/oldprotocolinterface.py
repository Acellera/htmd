# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import abc
from difflib import get_close_matches
import os
import re
from htmd.util import ensurelist

RANGE_ANY = 0
RANGE_POS = 1
RANGE_NEG = 2
RANGE_0POS = 3
RANGE_0NEG = 4

TYPE_INT = 0
TYPE_FLOAT = 1


class ProtocolInterface:
    def __init__(self):
        self._commands = {}
        self._deprecated = {}
        self._numCommands = 0
        self._commandsOrder = {}

    def __setattr__(self, key, value):
        if key[0] != '_':
            self._validate(key, value)
        self.__dict__[key] = value

    def __str__(self):
        sortedkeys = [x[0] for x in sorted(self._commandsOrder.items(), key=lambda x: x[1])]
        s = ''
        for cmd in sortedkeys:
            if isinstance(self.__dict__[cmd], str):
                s += '{} = \'{}\'\n'.format(cmd, self.__dict__[cmd])
            else:
                s += '{} = {}\n'.format(cmd, self.__dict__[cmd])
        return s

    def __misc(self, key, default):
        self._commandsOrder[key] = self._numCommands
        self._numCommands += 1
        self.__dict__[key] = default

    def _cmdFile(self, key, datatype, descr, default, exist=False, writable=False, multiple=False, check=None):
        self._commands[key] = FileValidator(key, datatype, descr, default, exist, writable, multiple, check)
        self.__misc(key, default)

    def _cmdListList(self, key, datatype, descr, length, default=None):
        self._commands[key] = ListListValidator(key, datatype, descr, default, length)
        self.__misc(key, default)

    def _cmdString(self, key, datatype, descr, default, valid_values=None):
        self._commands[key] = StringValidator(key, datatype, descr, default, valid_values)
        self.__misc(key, default)

    def _cmdTimestep(self, key, datatype, descr, default):
        self._commands[key] = TimestepValidator(key, datatype, descr, default)
        self.__misc(key, default)

    def _cmdList(self, key, datatype, descr, default, valid_values=None):
        self._commands[key] = ListValidator(key, datatype, descr, default, valid_values)
        self.__misc(key, default)

    def _cmdDict(self, key, datatype, descr, default, key_type=None, val_type=None):
        self._commands[key] = DictionaryValidator(key, datatype, descr, default, key_type, val_type)
        self.__misc(key, default)

    def _cmdBinary(self, key, datatype, descr, default):
        self._commands[key] = BinaryValidator(key, datatype, descr, default)
        self.__misc(key, default)

    def _cmdValue(self, key, datatype, descr, default, realdatatype, valid_range, list_len=1, units=None):
        self._commands[key] = ValueValidator(key, datatype, descr, default, realdatatype, valid_range, list_len, units)
        self.__misc(key, default)

    def _cmdObject(self, key, datatype, descr, default, classname):
        self._commands[key] = ObjectValidator(key, datatype, descr, default, classname)
        self.__misc(key, default)

    def _cmdClass(self, key, datatype, descr, default, classname):
        self._commands[key] = ClassValidator(key, datatype, descr, default, classname)
        self.__misc(key, default)

    def _cmdBoolean(self, key, datatype, descr, default):
        self._commands[key] = BooleanValidator(key, datatype, descr, default)
        self.__misc(key, default)

    def _cmdFunction(self, key, datatype, descr, default):
        self._commands[key] = FunctionValidator(key, datatype, descr, default)
        self.__misc(key, default)

    def _cmdDeprecated(self, key, newkey=None):
        self._deprecated[key] = newkey

    def _printDocString(self):
        sortedkeys = [x[0] for x in sorted(self._commandsOrder.items(), key=lambda x: x[1])]
        docs = ''
        for k in sortedkeys:
            docs += '{} : {}\n'.format(k, self._commands[k].docString())
        print(docs)

    def _validate(self, key, value, basedir=None):
        if key in self._deprecated:
            newkey = self._deprecated[key]
            if not newkey:
                print(" Attribute '" + key + "'\t is deprecated and is no longer required")
            else:
                print(" Attribute '" + key + "'\t is deprecated and replaced by '" + newkey + "'")
                key = newkey

        if key not in self._commands:
            strerror = "Attribute '" + key + "' not allowed in this class."
            match = get_close_matches(key, self._commands)
            if match:
                strerror += " Try '" + match[0] + "'"
            raise ValueError(strerror)

        cmd = self._commands[key]
        return cmd.validate(value, basedir=basedir)

    def _toArgParse(self, description, parser=None):
        import argparse
        if parser is None:
            parser = argparse.ArgumentParser(description=description)

        sortedkeys = [x[0] for x in sorted(self._commandsOrder.items(), key=lambda x: x[1])]
        for k in sortedkeys:
            parser.add_argument('--{}'.format(k), help=self._commands[k].descr, default=self._commands[k].default)

        return parser

    def _fromArgParse(self, args):
        argsdict = vars(args)
        for k in argsdict:
            self.__dict__[k] = argsdict[k]


class Validator:
    def __init__(self, key, datatype, descr, default):
        self.datatype = datatype
        if not descr:
            desc = os.path.join("help", "en", key)
        if descr and os.path.exists(descr):
            f = open(descr, 'r')
            self.descr = f.read()
            f.close()
        else:
            self.descr = descr
        self.default = default

    def docString(self):
        doc = ''
        if hasattr(self, 'valid_values') and self.valid_values is not None:
            doc += '{}, '.format(self.valid_values)
        doc += '{}'.format(self.datatype)
        # if self.default is not None:
        if isinstance(self.default, str):
            doc += ', default=\'{}\''.format(self.default)
        else:
            doc += ', default={}'.format(self.default)
        if hasattr(self, 'units') and self.units is not None:
            doc += ', units={}'.format(self.units)
        doc += '\n\t{}'.format(self.descr)
        return doc


# TODO: Maybe better to consolidate it into ValueValidator
class BooleanValidator(Validator):
    def __init__(self, key, datatype, descr, default):
        super().__init__(key, datatype, descr, default)

    def args(self):
        return

    def validate(self, object, basedir=None):
        return isinstance(object, bool)


class ObjectValidator(Validator):
    def __init__(self, key, datatype, descr, default, classname):
        super().__init__(key, datatype, descr, default)
        self.classname = classname

    def args(self):
        return

    def validate(self, object, basedir=None):
        classname = self.classname
        object = ensurelist(object)
        classname = ensurelist(classname)

        for obj in object:
            valid = False
            for cl in classname:
                if isinstance(obj, cl):
                    valid = True
                    break
            if not valid:
                raise ValueError('Value must be object of {}'.format(self.classname))

        return object

class ClassValidator(Validator):
    def __init__(self, key, datatype, descr, default, classname):
        super().__init__(key, datatype, descr, default)
        self.classname = classname

    def args(self):
        return

    def validate(self, object, basedir=None):
        classname = self.classname
        object = ensurelist(object)
        classname = ensurelist(classname)

        for obj in object:
            valid = False
            for cl in classname:
                if issubclass(obj, cl):
                    valid = True
                    break
            if not valid:
                raise ValueError('Value must be subclass of {}'.format(self.classname))

        return object


class FunctionValidator(Validator):
    def __init__(self, key, datatype, descr, default):
        super().__init__(key, datatype, descr, default)

    def args(self):
        return

    def validate(self, object, basedir=None):
        if hasattr(object, '__call__'):
            return
        elif isinstance(object, tuple) and hasattr(object[0], '__call__'):
            # Supporting delayed functions which are passed as function/arguments tuples
            return
        else:
            raise ValueError('Value must be a function')


class ListListValidator(Validator):
    def __init__(self, key, datatype, descr, default, length):
        super().__init__(key, datatype, descr, default)
        self.length = length

    def args(self):
        return "A list of lists with length " + self.length

    def validate(self, value_list, basedir=None):
        vv = ""
        for t in value_list:
            vv = vv + " " + t
        vv = vv.strip()
        if vv == "[]":
            self.value = None
            return self.value

        vv = re.sub("\[", "", vv)
        vv = re.sub(",", "", vv)
        vv = re.sub("\]\]$", "", vv)
        vv = re.sub("\],", ",", vv)
        vv = re.sub("'", "", vv)
        vv = vv.split(",")
        value = []
        for t1 in vv:
            phi = []
            for t2 in t1.strip().split():
                phi.append(t2.strip().upper())
            if len(phi) != self.length:
                raise ValueError("list requires " + str(self.length) + " elements")
            value.append(phi)

        if len(value):
            self.value = value

        return self.value


class StringValidator(Validator):
    def __init__(self, key, datatype, descr, default, valid_values=None):
        super().__init__(key, datatype, descr, default)
        self.valid_values = valid_values
        self.units = None

    def args(self):
        if self.valid_values:
            return self.valid_values
        else:
            return "free string"

    def validate(self, value_list, basedir=None):
        import re
        if type(value_list) is not list and type(value_list) is not tuple:
            if value_list is None:
                value_list = 'None'
            value_list = [value_list]
        for v in value_list:
            if self.valid_values and v not in self.valid_values:
                raise ValueError('Value must be one of: {}'.format(self.valid_values))

        value = ' '.join(value_list)
        value = value.strip()
        # Remove starting or ending quotes (single, double)
        value = re.sub('^"', '', value)
        value = re.sub('"$', '', value)
        value = re.sub("^'", '', value)
        value = re.sub("'$", '', value)
        return value


class TimestepValidator(Validator):
    def __init__(self, key, datatype, descr, default):
        super().__init__(key, datatype, descr, default)
        self.units = None

    def args(self):
        return "[ number (ps|ns|us) ]"

    def validate(self, value_list, basedir=None):
        if type(value_list) is not list and type(value_list) is not tuple:
            value_list = [value_list]

        try:
            if len(value_list) == 1:
                value = int(value_list[0])
            elif len(value_list) == 2:
                value = float(value_list[0])
                if value_list[1].lower() == "ps":
                    value = value * 1000 / 4  # TODO FIXME -- note the hard-coded assumption about the timestep
                if value_list[1].lower() == "ns":
                    value = value * 1000000 / 4
                if value_list[1].lower() == "us":
                    value = value * 1000000000 / 4
                value = int(value)
            else:
                raise ValueError()
        except:
            raise ValueError("Value must be >= 0 and may have a unit suffix [ps,ns,us]")
        return value


class ListValidator(Validator):
    def __init__(self, key, datatype, descr, default, valid_values):
        super().__init__(key, datatype, descr, default)
        self.valid_values = valid_values
        self.units = None

    def args(self):
        return str(self.valid_values)

    def validate(self, value_list, basedir=None):
        if type(value_list) is not list and type(value_list) is not tuple:
            value_list = [value_list]
        for value in value_list:
            if self.valid_values is not None and value not in self.valid_values:
                raise ValueError("Value must be one of " + str(self.valid_values))
        return value_list


class DictionaryValidator(Validator):
    def __init__(self, key, datatype, descr, default, key_type=None, val_type=None):
        super().__init__(key, datatype, descr, default)
        self.key_type = key_type
        self.val_type = val_type

    def args(self):
        return

    def validate(self, value, basedir=None):
        if type(value) is not dict:
            raise TypeError('Value is not a dictionary.')
        for k in value.keys():
            if self.key_type is not None and not isinstance(k, self.key_type):
                raise TypeError('Key {} is not of type {}'.format(k, self.key_type))
            if self.val_type is not None and not isinstance(value[k], self.val_type):
                raise TypeError('Value {} of key {} is not of type {}'.format(value[k], k, self.val_type))
        return value


class BinaryValidator(Validator):
    def __init__(self, key, type, descr, default):
        super().__init__(key, type, descr, default)
        self.units = None

    def args(self):
        return "[ on | off ]"

    def validate(self, value_list, basedir=None):
        if type(value_list) is not list and type(value_list) is not tuple:
            value_list = [value_list]
        value = value_list[0]

        try:
            if value:
                value = True
                return value
            if not value or value is None:
                value = False
                return value

            value = value.lower()

            if value == "yes" or value == "on" or value == "1" or value == "all" or value == "true":
                value = True
            elif value == "no" or value == "off" or value == "0" or value == "none" or value == "false":
                value = False
            else:
                raise ValueError("")
        except:
            raise ValueError("Value must be binary ( on|off )")
        return value


class FileValidator(Validator):
    def __init__(self, key, datatype, descr, default, exist=False, writable=False, multiple=False, check=None):
        super().__init__(key, datatype, descr, default)
        self.multiple = multiple
        self.must_exist = exist
        self.writable = writable
        self.check = check

    def args(self):
        dd = "file"
        if self.must_exist:
            dd = "input file"
        if self.writable:
            dd = "output file"
        if self.multiple:
            dd = "list of " + dd + "s"
        return "[ " + dd + " ]"

    def validate(self, value_list, basedir=None):
        if type(value_list) is not list and type(value_list) is not tuple:
            value_list = [value_list]

        l = []
        for value in value_list:
            if self.must_exist:
                found = False
                if basedir and os.path.isfile(os.path.join(basedir, value)):
                    value = os.path.join(basedir, value)
                    found = True
                if not found and os.path.isfile(value):
                    found = True
                if not found:
                    raise NameError("File '" + value + "' does not exist")
                if not os.access(value, os.R_OK):
                    raise NameError("File '" + value + "' cannot be read")
            if self.writable:
                if basedir:
                    value = os.path.join(basedir, value)
                try:
                    f = open(value, "a+")
                    f.close()
                except:
                    raise NameError("File '" + value + "' is not writable")
            l.append(value)

        # Test to see if the file exists
        if len(l) != 1 and not self.multiple:
            raise NameError("Value must be a single filename")
        if not self.multiple:
            value = l[0]
        else:
            value = l
        return value


class ValueValidator(Validator):
    def __init__(self, key, datatype, descr, default, realdatatype, valid_range, list_len=1, units=None):
        super().__init__(key, datatype, descr, default)
        self.valid_range = valid_range
        self.list_len = list_len
        self.units = units
        self.realdatatype = realdatatype

    def validate(self, value_list, basedir=None):
        valid_range = self.valid_range

        if type(value_list) is not list and type(value_list) is not tuple:
            value_list = [value_list]

        l = []
        if self.list_len != len(value_list):
            if self.list_len == 0:
                raise NameError("Value must be a scalar")
            else:
                raise NameError("Value must be a vector with " + str(self.list_len) + " elements")

        for value in value_list:
            if self.realdatatype == TYPE_FLOAT:
                strfudge = "a number"
            elif self.realdatatype == TYPE_INT:
                strfudge = "an integer"

            try:
                try:
                    if self.realdatatype == TYPE_FLOAT:
                        value = float(value)
                    if self.realdatatype == TYPE_INT:
                        value = int(value)
                except:
                    raise ValueError("Value is not " + strfudge)

                if valid_range == RANGE_POS and value <= 0:
                    raise ValueError("Value must be " + strfudge + " > 0")
                if valid_range == RANGE_0POS and value < 0:
                    raise ValueError("Value must be " + strfudge + " >=0")
                if valid_range == RANGE_NEG and value >= 0:
                    raise ValueError("Value must be " + strfudge + " < 0")
                if valid_range == RANGE_0NEG and value > 0:
                    raise ValueError("Value must be " + strfudge + " <=0")
            except NameError as e:
                raise NameError(e)
            if self.list_len == 1:
                l = value
            else:
                l.append(value)
        return l


'''class AttributeParser:
    def __init__(self, parent=None):
        if parent is not None:
            self.attributes = parent.attributes
            self.deprecated = parent.deprecated
            self.numAttributes = parent.numAttributes
            self.attributesOrder = parent.attributesOrder
        else:
            self.attributes = {}
            self.deprecated = {}
            self.numAttributes = 0
            self.attributesOrder = {}

    def addFile(self, key, datatype, descr, default, exist=False, writable=False, multiple=False, check=None):
        self.attributes[key] = FileValidator(datatype, descr, default, exist, writable, multiple, check)
        self.attributesOrder[key] = self.numAttributes
        self.numAttributes += 1

    def addString(self, key, datatype, descr, default, valid_values=None):
        self.attributes[key] = StringValidator(datatype, descr, default, valid_values)
        self.attributesOrder[key] = self.numAttributes
        self.numAttributes += 1

    def addTimestep(self, key, datatype, descr, default):
        self.attributes[key] = TimestepValidator(datatype, descr, default)
        self.attributesOrder[key] = self.numAttributes
        self.numAttributes += 1

    def addList(self, key, datatype, descr, default, valid_values):
        self.attributes[key] = ListValidator(datatype, descr, default, valid_values)
        self.attributesOrder[key] = self.numAttributes
        self.numAttributes += 1

    def addBinary(self, key, datatype, descr, default):
        self.attributes[key] = BinaryValidator(datatype, descr, default)
        self.attributesOrder[key] = self.numAttributes
        self.numAttributes += 1

    def addValue(self, key, datatype, descr, default, realdatatype, valid_range, list_len=1, units=None):
        self.attributes[key] = ValueValidator(datatype, descr, default, realdatatype, valid_range, list_len, units)
        self.attributesOrder[key] = self.numAttributes
        self.numAttributes += 1

    def addObject(self, key, datatype, descr, default, classname):
        self.attributes[key] = ObjectValidator(datatype, descr, default, classname)
        self.attributesOrder[key] = self.numAttributes
        self.numAttributes += 1

    def addBoolean(self, key, datatype, descr, default):
        self.attributes[key] = BooleanValidator(datatype, descr, default)
        self.attributesOrder[key] = self.numAttributes
        self.numAttributes += 1

    def addDeprecated(self, key, newkey=None):
        self.deprecated[key] = newkey

    def printDocString(self):
        docs = ''
        for k in sorted(self.attributesOrder.items(), key=lambda x: x[1]):
            docs += '{} : {}\n'.format(k[0], self.attributes[k[0]].docString())
        print(docs)

    def testDeprecation(self, key):
        if key in self.deprecated:
            newkey = self.deprecated[key]
            if not newkey:
                print(" Attribute '" + key + "'\t is deprecated and is no longer required")
            else:
                print(" Attribute '" + key + "'\t is deprecated and replaced by '" + newkey + "'")
                return newkey
        else:
            return key

    def getDefaultConfiguration(self):
        ret = {}
        for i in self.attributes:
            ret[i] = self.attributes[i].default
        return ret

    def validate(self, key, value, basedir=None):
        key = self.testDeprecation(key)
        if key not in self.attributes:
            strerror = "Attribute '" + key + "' not allowed in this class."
            match = get_close_matches(key, self.attributes)
            if match:
                strerror += " Try '" + match[0] + "'"
            raise ValueError(strerror)

        cmd = self.attributes[key]
        return cmd.validate(value, basedir=basedir)

    def setAttributes(self, obj):
        for a in self.attributes:
            obj.__dict__[a] = self.attributes[a].default'''
