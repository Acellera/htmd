# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from difflib import get_close_matches
from htmd.decorators import _Deprecated


def uisetattr(self, key, value, keylist):
    if key.startswith('_'):  # private ok
        self.__dict__[key] = value
        return
    if key in keylist:  # valid key
        if value is not None:  # not None
            self.__dict__[key] = value
        elif key in self.__dict__.keys():
            del(self.__dict__[key])
    else:
        guess = get_close_matches(key, keylist)
        if guess:
            raise NameError('Option "' + key + '" is not recognized. Did you mean "' + str(guess[0]) + '"?')
        else:
            raise NameError('Option "' + key + '" is not recognized.')


@_Deprecated('1.5.15', 'htmd.protocols.protocolinterface.ProtocolInterface')
class UserInterface:
    """
    Creates a class for which only certain public attributes can be created.
    All private attribute can be freely created.
    """

    def __setattr__(self, key, value):
        if '_commands' not in self.__dict__:
            self.__dict__['_commands'] = {}  # create empty dictionary of commands
        uisetattr(self, key, value, self._commands.keys())

    def commands(self):
        cmds = {}
        for f in self.__dict__:
            if not str(f).startswith('_'):  # not private
                cmds[f] = self.__dict__[f]
        return cmds

    def cmdlinelist(self, sep='-'):
        cmds = []
        for f in self.__dict__:
            if not str(f).startswith('_'):  # not private
                cmds.append(sep+f)
                cmds.append(self.__dict__[f])
        return cmds


class Ttt(UserInterface):
    def __init__(self):
        self.__dict__['_commands'] = {'f': None, 'g': None, 'r': None}

#    def __init__(self):
#        self.addCommand('_commands',{'f':None,'g':None,'r':None})


if __name__ == "__main__":
    ui = UserInterface()
    try:
        ui.g = 4
    except:
        print('Erroring as it should. OK')
    test = Ttt()
    test.f = 0
    test.g = 1
    try:
        test.d = 1
    except:
        print('Erroring as it should. OK')
    print(test.commands())
    print(test.cmdlinelist())
