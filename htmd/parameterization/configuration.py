from htmd.parameterization.command import Command
from difflib import get_close_matches
import re
import os
from htmd.protocols.protocolinterface import ProtocolInterface, TYPE_FLOAT, TYPE_INT, RANGE_0POS, RANGE_POS, RANGE_ANY


class ParameterizationConfig(ProtocolInterface):
    def __init__(self, config=None, check=True, **kwargs):
        self._commands = Command(check=check)
        c = self._commands.get_default_configuration(check=check)

        inputfile = None

        self._basedir = os.getcwd()
        if config:
            if os.path.isfile(config):
                self._basedir = os.path.dirname(config)
                inputfile = config
            elif os.path.isdir(config):
                self._basedir = config  # (os.path.abspath(config))
                inputfile = os.path.join(self._basedir, "input")
                #              print(inputfile)
                if not os.path.isfile(inputfile):
                    raise NameError("Input file not found in config directory")
            else:
                print("Configuration is " + config)
                raise NameError("Invalid input. 'config' must be an  input file or a directory containing 'input'")

        for i in c.keys():
            self.__dict__[i] = c[i]

        # TODO: moves to PI
        if inputfile:
            # Read in key-value pairs, excluding any comments beginning #
            f = open(inputfile, "r")
            ll = 1
            for linex in f:
                line = linex
                line = re.sub(r'\t', ' ', line)
                line = re.sub(r'\n', '', line)
                line = re.sub(r'\r', '', line)
                line = re.sub(r'#.*$', '', line, flags=re.M | re.S)
                # Remove any leading or trailing whitespace
                line = line.strip()
                # Remove comment
                if len(line) > 0:
                    key = re.sub('\s.*$', '', line, flags=re.M | re.S)
                    val = re.sub('^[^\s]+\s', '', line, flags=re.M | re.S)
                    key = key.strip()
                    val = val.strip()
                    val = re.split(r'\s+', val)
                    try:
                        self.__setattr__(key, val)
                    except NameError as e:
                        raise NameError(" Line " + str(ll) + "\t: " + linex + "        \t: " + str(e))
                ll += 1
            f.close()

        # Now set any params specified in the constructor
        for key, value in kwargs.items():
            self.__setattr__(key, value)

    # TODO: moves to PI
    def __eq__(self, config):

        for i in self.__dict__:
            if not i.startswith("_"):
                if not (self.__getattr__(i) == config.__getattr__(i)):
                    # print("FAILED 1 " + i + " " + self.__getattr__(i) + " " + config.__getattr(i))
                    return False
        for i in config.__dict__:
            if not i.startswith("_"):
                if not (config.__getattr__(i) == self.__getattr__(i)):
                    # print("FAILED 2 " + i + " " + self.__getattr_(i) + " " + config.__getattr(i))
                    return False

        return True

"""
    def __setattr__(self, key, value):
        if key.startswith("_"):
            self.__dict__[key] = value
            return

            # key = key.replace( '-', '_' )
        #        key = Command.test_deprecation( key, value )
        if key:
            #            print(key)
            value = self._commands._commands[key].validate(value, self._basedir)
            self.__dict__[key] = value
        else:
            key = key.lower()
            #            key = Command.test_deprecation( key, value )
            if key:
                value = self._commands._commands[key].validate(value, self._basedir)
                self.__dict__[key] = value

    def __getattr__(self, key):
        if key.startswith("_"):
            return self.__dict__[key]

        try:
            # key = key.replace( '_', '-' );
            if key in self.__dict__:
                return self.__dict__[key]
            keyx = key.lower()
            return self.__dict__[keyx]
        except:

            errstr = "Command '" + key + "' not found."
            match = get_close_matches(key, self._commands.commands)
            if match:
                errstr = errstr + " Try '" + match[0] + "'"
            raise NameError(errstr)

    def __repr__(self):
        return self.__str__()
"""
    # TODO: moves to PI
    def save(self, filename, fmt="pretty"):
        with open(filename, "w") as fh:
            for cmd in sorted(self.__dict__):
                if self.__dict__[cmd] is not None:
                    if not (cmd.startswith("_")):
                        if fmt == "pretty":
                            print("%25s %s" % (cmd, str(self.__dict__[cmd])), file=fh)
                        elif fmt == "shell":
                            print("export %s=\"%s\"" % (cmd, self.__dict__[cmd]), file=fh)

    def __str__(self):
        s = ""
        for cmd in sorted(self.__dict__):
            if not cmd.startswith("_"):
                s += "%25s %s\n" % (cmd, str(self.__dict__[cmd]))
        return s

    def set(self, key, value):
        self.__setattr__(key, value)

    def help(self, command):
        return self._commands.help(command)

    def stage(self, directory):
        directory = os.path.abspath(directory)
        os.makedirs(directory)
        self.save(os.path.join(directory, "input"))
