import os
import psutil
from htmd.protocols.protocolinterface import ProtocolInterface, TYPE_FLOAT, TYPE_INT, RANGE_0POS, RANGE_POS, RANGE_ANY


class Command(ProtocolInterface):
    setup = False
    commands = {}

    def __init__(self, check=True):
        super().__init__()

        self._cmdString("Model", "str", None, "Nonpolar", valid_values=["Nonpolar", "Drude", "Match"])
        self._cmdBinary("AsymmetricTorsion", "bool", None, False)
        self._cmdString("ExecutionMode", "str", None, "Inline", valid_values=["Inline", "PBS", "LSF"])
        self._cmdValue('NetCharge', 'int', None, 0, TYPE_INT, RANGE_ANY)
        self._cmdValue("Multiplicity", "int", None, 1, TYPE_INT, RANGE_ANY)
        self._cmdValue("E14Fac", "float", None, 1.0, TYPE_FLOAT, RANGE_0POS)
        self._cmdValue("w_H_Donor_Acceptor", "float", None, 1., TYPE_FLOAT, RANGE_0POS)
        self._cmdValue("w_charge", "float", None, 3., TYPE_FLOAT, RANGE_0POS)
        self._cmdValue("w_water_E_min", "float", None, 0.4, TYPE_FLOAT, RANGE_0POS)
        self._cmdValue("w_water_R_min", "float", None, 8.0, TYPE_FLOAT, RANGE_0POS)
        self._cmdValue("w_alpha", "float", None, 0.4, TYPE_FLOAT, RANGE_0POS)
        self._cmdValue("w_thole", "float", None, 0.2, TYPE_FLOAT, RANGE_0POS)
        self._cmdValue("wd_charge", "float", None, 1., TYPE_FLOAT, RANGE_0POS)
        self._cmdValue("wd_water_E_min", "float", None, 0.4, TYPE_FLOAT, RANGE_0POS)
        self._cmdValue("wd_water_R_min", "float", None, 8.0, TYPE_FLOAT, RANGE_0POS)
        self._cmdValue("wd_alpha", "float", None, 0.4, TYPE_FLOAT, RANGE_0POS)
        self._cmdValue("wd_thole", "float", None, 0.2, TYPE_FLOAT, RANGE_0POS)

        # Mandatory args
        self._cmdFile("FileName", "str", None, None, exist=True, check=check)
        self._cmdFile("Equivalent", "str", None, None, exist=True, check=check)
        self._cmdFile("Neutral", "str", None, None, exist=True, check=check)
        self._cmdFile("FixCharges", "str", None, None, exist=True, check=check)

        self._cmdValue("MaxTorsion", "int", None, 25, TYPE_INT, RANGE_POS)
        # self._cmd[ 'FileName'   ]	= Command.File( None, exist=True, check=check )
        # self._cmd[ 'Torsions'   ]  = Command.Torsionlist( None )
        # self._cmd[ 'JobName'    ]	= Command.Stringx( None )

        self._cmdString("JobName", "str", None, None)
        self._cmdListList("Torsions", "int", None, 4)
        # self._cmd[ 'Equivalent' ]	= Command.File( None, exist=True, check=check )
        # self._cmd[ 'Neutral'    ]	= Command.File( None, exist=True, check=check )
        # self._cmd[ 'FixCharges' ]	= Command.File( None, exist=True, check=check )
        # self._cmd[ 'MaxTorsion' ]    = Command.Value( 25, Command.TYPE_INT, Command.RANGE_POS, 1)

        ncpus = psutil.cpu_count()
        if 'NCPUS' in os.environ:
            try:
                ncpus = int(os.environ['NCPUS'])
            except:
                pass
        mem = psutil.virtual_memory()[0] / 1000000000  # mem in GB
        mem = int(mem)
        mem = int(mem / 2)  # Half the system memory
        if mem < 1:
            mem = 1  # never less than 1GB/core

        if 'MEM' in os.environ:
            try:
                mem = int(os.environ['MEM'])
            except:
                pass

            #        self._cmd[ 'GAUSS_SCRDIR' ]= Command.File( '/tmp', exist=False )
            #        self._cmd[ 'NCORES' ] = Command.Value( ncpus, Command.TYPE_INT, Command.RANGE_POS, 1)
            #        self._cmd[ 'MEMORY' ] = Command.Value( mem, Command.TYPE_INT, Command.RANGE_POS, 1)

            #        self._cmd[ 'Debug'    ]	= Command.Binary( False )

        self._cmdFile("GAUSS_SCRDIR", "str", None, "/tmp", exist=False, check=check)
        self._cmdBinary("Debug", "bool", None, False)
        self._cmdValue("NCORES", "int", None, ncpus, TYPE_INT, RANGE_POS)
        self._cmdValue("MEMORY", "int", None, mem, TYPE_INT, RANGE_POS)

    def get_default_configuration(self):
        ret = dict()
        for i in self._commands.keys():
            ret[i] = self._commands[i].default
        return ret

'''
    @staticmethod
    def get_help_string(cmd):
        libdir = (os.path.dirname(inspect.getfile(Command)))

        libdir = (os.path.join(libdir, "help"))

        libdir = (os.path.join(libdir, "en"))

        helpfile = (os.path.join(libdir, os.path.basename(cmd)))

        helpfile += ".txt"

        if not os.path.isfile(helpfile):
            return "No help found"
        else:
            with open(helpfile) as ff:
                return ff.read().strip()

    @staticmethod
    def pretty_print(cmd):
        print(cmd)
        return

        width = shutil.get_terminal_size(fallback=(80, 25))
        width = width.columns
        if width < 40:
            width = 80
        tok = re.split(r'\s+', cmd)
        l = 4
        print("   ", end="")
        for t in tok:
            if (l + len(t) + 1) <= width:
                print(" " + t, end="")
                l = l + len(t) + 1
            else:
                l = 4
                print("\n    " + t, end="")
                l = l + len(t) + 1

    @staticmethod
    def help(cmd):
        if not cmd:
            print("\n  Valid configuration file commands:\n")
            for a in sorted(self._cmd.keys()):
                print("    " + a)
            print("\n  parameterization --command [command] for detailed help on a specific command\n\n")
            # print sections
            pass
        else:
            match = get_close_matches(cmd, self._cmd)
            if not match or (len(match) < 1):
                print("\nNo matching command for '" + cmd + "' found\n")
                return
            else:
                helpstr = Command.get_help_string(match[0])
                print("\n   " + match[0] + " " + self._cmd[match[0]].args() + "\n")

                Command.pretty_print(helpstr)

                units = self._cmd[match[0]].units
                if not units:
                    units = ""
                print("\n   Default: " + str(self._cmd[match[0]].default) + " " + units + "\n")
'''
# find

# def validate( self,  key, value, basedir=None ):
#
# #       try:
#  #          cmd =self._commands[key]
#   #     except:
#    #        strerror= "Command '" + key + "' not found.";
#            match = get_close_matches( key, self._commands )
#            if match:
#                strerror = strerror + " Try '" + match[0] + "'";
#            raise NameError( strerror )
#
#        return cmd.validate( value, basedir=basedir )
