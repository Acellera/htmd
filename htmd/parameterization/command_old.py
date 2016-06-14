from difflib import get_close_matches
import os
import sys
import re
import shutil
import inspect
import psutil

class Command:
    setup    = False
    commands = {}

    RANGE_ANY  = 0
    RANGE_POS  = 1
    RANGE_NEG  = 2
    RANGE_0POS = 3
    RANGE_0NEG = 4

    TYPE_INT   = 0
    TYPE_FLOAT = 1


    @staticmethod
    def test_deprecation( key, value ):
        if( key in Command.deprecated ):
            newkey = Command.deprecated[key]
            if( not newkey ):
                print( " Command '" + key + "'\t is deprecated and is no longer required" )
            else:
                print( " Command '" + key + "'\t is deprecated and replaced by '" + newkey + "'" )
                return newkey
        else:
            return key

    @staticmethod
    def add_commands( pathcheck=True ):
        if Command.setup:
            return
        Command.deprecated = {
        }

        # Optional Args

#        Command.commands[ 'Ini_Parameters'  ] = Command.List( "Match", [ "Match" ] )
        #Command.commands[ 'TorsionFitting'  ]	= Command.Binary( True )
        Command.commands[ 'AsymmetricTorsion'    ]	= Command.Binary( False )
        Command.commands[ 'Model'    ]	        = Command.List( "Nonpolar", [ "Nonpolar", "Drude", "Match" ] );
        Command.commands[ 'ExecutionMode'    ]	        = Command.List( "Inline", [ "Inline", "PBS" , "LSF" ] );
        Command.commands[ 'NetCharge'       ]	= Command.Value( 0, Command.TYPE_INT, Command.RANGE_ANY, 1  )
        Command.commands[ 'Multiplicity'    ]	= Command.Value( 1, Command.TYPE_INT, Command.RANGE_ANY, 1  )
        Command.commands[ 'E14FAC'          ]	= Command.Value( 1., Command.TYPE_FLOAT, Command.RANGE_0POS, 1  )
        Command.commands[ 'w_H_Donor_Acceptor'          ]	= Command.Value( 1., Command.TYPE_FLOAT, Command.RANGE_0POS, 1  )
        Command.commands[ 'w_charge'          ]	= Command.Value( 3., Command.TYPE_FLOAT, Command.RANGE_0POS, 1  )
        Command.commands[ 'w_water_E_min'     ]	= Command.Value( 0.4, Command.TYPE_FLOAT, Command.RANGE_0POS, 1  )
        Command.commands[ 'w_water_R_min'     ]	= Command.Value( 8.0, Command.TYPE_FLOAT, Command.RANGE_0POS, 1  )
        Command.commands[ 'w_alpha'           ]	= Command.Value( 0.4, Command.TYPE_FLOAT, Command.RANGE_0POS, 1  )
        Command.commands[ 'w_thole'           ]	= Command.Value( 0.2, Command.TYPE_FLOAT, Command.RANGE_0POS, 1  )
        Command.commands[ 'wd_charge'         ]	= Command.Value( 1., Command.TYPE_FLOAT, Command.RANGE_0POS, 1  )
        Command.commands[ 'wd_water_E_min'    ]	= Command.Value( 0.4, Command.TYPE_FLOAT, Command.RANGE_0POS, 1  )
        Command.commands[ 'wd_water_R_min'    ]	= Command.Value( 8.0, Command.TYPE_FLOAT, Command.RANGE_0POS, 1  )
        Command.commands[ 'wd_alpha'          ]	= Command.Value( 0.4, Command.TYPE_FLOAT, Command.RANGE_0POS, 1  )
        Command.commands[ 'wd_thole'          ]	= Command.Value( 0.2, Command.TYPE_FLOAT, Command.RANGE_0POS, 1  )



        # Mandatory args
        Command.commands[ 'FileName'   ]	= Command.File( None, exist=True, pathcheck=pathcheck )
        #Command.commands[ 'FileType'   ]	= Command.List( 'mol2', [ 'mol2', 'pdb' ] )
        Command.commands[ 'Torsions'   ]  = Command.Torsionlist( None )
        Command.commands[ 'JobName'    ]	= Command.Stringx( None )
        Command.commands[ 'Equivalent' ]	= Command.File( None, exist=True, pathcheck=pathcheck )
        Command.commands[ 'Neutral'    ]	= Command.File( None, exist=True, pathcheck=pathcheck )
        Command.commands[ 'FixCharges' ]	= Command.File( None, exist=True, pathcheck=pathcheck )
        Command.commands[ 'MaxTorsion' ]    = Command.Value( 25, Command.TYPE_INT, Command.RANGE_POS, 1)

        ncpus = psutil.cpu_count()
        if 'NCPUS' in os.environ:
           try:
              ncpus=int(os.environ['NCPUS'])
           except:
              pass
        mem   = psutil.virtual_memory()[0] / 1000000000 # mem in GB
        mem   = int(mem)
        mem   = int(mem/2) # Half the system memory
        if mem < 1 : mem = 1 # never less than 1GB/core

        if 'MEM' in os.environ:
           try:
              mem=int(os.environ['MEM'])
           except:
              pass

        Command.commands[ 'GAUSS_SCRDIR' ]= Command.File( '/tmp', exist=False )
        Command.commands[ 'NCORES' ] = Command.Value( ncpus, Command.TYPE_INT, Command.RANGE_POS, 1)
        Command.commands[ 'MEMORY' ] = Command.Value( mem, Command.TYPE_INT, Command.RANGE_POS, 1)

        Command.commands[ 'Debug'    ]	= Command.Binary( False )
 

        Command.setup = True


    @staticmethod
    def get_help_string( cmd ):
        libdir = (os.path.dirname(inspect.getfile(Command)))

        libdir = (os.path.join(libdir, "help"))

        libdir = (os.path.join(libdir, "en"))

        helpfile=  (os.path.join( libdir, os.path.basename( cmd ) ))

        helpfile= helpfile + ".txt"

        if not os.path.isfile( helpfile ):
            return "No help found"
        else:
            with open( helpfile ) as ff:
                return ff.read().strip()

    @staticmethod
    def pretty_print( cmd ):
        print( cmd )
        return

        width = shutil.get_terminal_size( fallback=(80,25) )
        width = width.columns
        if width<40:
            width=80
        tok   = re.split( r'\s+', cmd );
        l=4
        print("   ", end="")
        for t in tok:
            if( l+ len(t)+1<=width):
                print(" "+t, end="");
                l= l + len(t)+1
            else:
                l=4;
                print("\n    "+t, end="");
                l= l + len(t)+1


    @staticmethod
    def help( cmd ):
        Command.add_commands()
        if( not cmd ):
            print("\n  Valid configuration file commands:\n");
            for a in sorted( Command.commands.keys() ):
                print( "    " + a )
            print("\n  parameterization --command [command] for detailed help on a specific command\n\n");
            # print sections
            pass
        else:
            match = get_close_matches( cmd, Command.commands )
            if not match or (len(match)<1):
                print( "\nNo matching command for '" + cmd + "' found\n" )
                return
            else:
                helpstr = Command.get_help_string( match[0] )
                print( "\n   " + match[0] + " "+ Command.commands[ match[0] ].args() + "\n" )

                Command.pretty_print( helpstr );

                units = Command.commands[match[0] ].units
                if not units:
                    units=""
                print( "\n   Default: " + str(Command.commands[ match[0] ].default) +" " + units + "\n" );



            # find
    @staticmethod
    def get_default_configuration( pathcheck=True ):
        Command.add_commands( pathcheck=pathcheck )
        ret = {}
        for i in Command.commands.keys():
            ret[ i ] = Command.commands[i].default
        return ret


    @staticmethod
    def validate(  key, value, basedir=None ):
        try:
            cmd =Command.commands[key]
        except:
            strerror= "Command '" + key + "' not found.";
            match = get_close_matches( key, Command.commands )
            if match:
                strerror = strerror + " Try '" + match[0] + "'";
            raise NameError( strerror )

        return cmd.validate( value, basedir=basedir )


    class Stringx:
        def __init__ (self, default, args=None ):
            self.value = default
            self.default=default
            self.units=None


        def hello(self):
            return "XXX" 

        def args(self):
            return "XXX" 


        def validate(self, value_list, basedir=None ):
            import re
            if type(value_list) is not list and type( value_list) is not tuple:
                value_list = [ value_list ]
            value = ""
            for v in value_list:
                value=value + " " +v
            self.value = value.strip()
            self.value = re.sub( '^"', '', self.value )
            self.value = re.sub( '"$', '', self.value )
            self.value = re.sub( "^'", '', self.value )
            self.value = re.sub( "'$", '', self.value )
            return self.value


    class Timestep:
        def __init__( self,default ):
            self.value   = default
            self.default = default
            self.units   = None

        def args(self):
            return "[ number (ps|ns|us) ]"

        def validate( self, value_list, basedir=None ):
            if type(value_list) is not list and type( value_list) is not tuple:
                value_list = [ value_list ]
            value = value_list[0]
            try:
                if len(value_list) == 1:
                    self.value = int(value_list[0])
                elif len(value_list) == 2:
                    self.value = float(value_list[0])
                    if value_list[1].lower() == "ps":
                        self.value = self.value * 1000/4 # TODO FIXME -- note the hard-coded assumption about the timestep
                    if value_list[1].lower() == "ns":
                        self.value = self.value * 1000000/4
                    if value_list[1].lower() == "us":
                        self.value = self.value * 1000000000/4
                    self.value = int(self.value)
                else:
                    raise ValueError()
            except:
                raise NameError( "Value must be >= 0 and may have a unit suffix [ps,ns,us]")
            return self.value

    class Torsionlist:
        def __init__( self, default ):
                        self.value   = default
                        self.units   = None
                        self.default = []
        def args(self):
            return str(self.sel)

        def validate( self, value_list, basedir=None ):
            vv=""
            for t  in value_list: vv=vv + " " +t 
            vv=vv.strip()
            if( vv == "[]" ): 
                  self.value=None
                  return self.value

            vv=re.sub( "\[", "", vv )
            vv=re.sub( ",", "", vv )
            vv=re.sub( "\]\]$", "", vv )
            vv=re.sub( "\],", ",", vv )
            vv=re.sub( "'", "", vv )
            vv=vv.split(",")
            value=[]
            for t1 in vv:
               phi=[]
               for t2 in t1.strip().split():
                 phi.append(t2.strip().upper())
               if(len(phi) != 4): raise ValueError( "Torsion list requires 4 elements" )
               value.append(phi)

            if(len(value)):
	            self.value = value

            return self.value



    class List:
        def __init__( self, default, sel ):
                        self.value   = default
                        self.default = default
                        self.sel     = sel
                        self.units   = None
        def args(self):
            return str(self.sel)

        def validate( self, value_list, basedir=None ):

            if type(value_list) is not list and type( value_list) is not tuple:
                value_list = [ value_list ]
            value = value_list[0]
            if value in self.sel:
                self.value = value
                return value
            else:
                raise NameError( "Value must be one of " + str( self.sel ) )

    class Binary:
        def __init__( self, default ):
            self.value   = default
            self.default = default
            self.units = None

        def args( self ):
            return "[ on | off ]"

        def validate( self, value_list, basedir=None ):

            if type(value_list) is not list and type( value_list) is not tuple:
                value_list = [ value_list ]
            value = value_list[0]

            try:
                if( value == True ):
                    value = True
                    return value
                if( value == False or value == None ):
                    value = False
                    return value

                value = value.lower()

                if(  value == "yes" or value == "on" or value == "1" or value == "all" or value == "true" ):
                    value = True
                elif(  value == "no" or value == "off" or value == "0" or value == "none" or value == "false" ):
                    value = False
                else:
                    raise ValueError("")
            except:
                raise NameError( "Value must be binary ( on|off )" )
            return value

    class File:

        def __init__( self, default, exist=False, writable=False, multiple=False, check=None, pathcheck=True ):
            self.multiple = multiple
            self.value = default
            self.must_exist = exist
            self.writable   = writable
            self.default = default
            self.check   = check
            self.pathcheck=pathcheck
            self.units = None
           # print("FILE INIT : " + str(pathcheck))
        def args(self):
            dd="file"
            if self.must_exist:
                dd="input file"
            if self.writable:
                dd="output file"
            if self.multiple:
                dd="list of "+dd+"s"
            return "[ " + dd + " ]"

        def validate( self, value_list, basedir=None ):
            if type(value_list) is not list and type( value_list) is not tuple:
                value_list = [ value_list ]

            l = []
            for value in value_list:
                if( self.must_exist  and self.pathcheck ):
                    found = False
                    if basedir and os.path.isfile( os.path.join( basedir, value ) ):
                        value = os.path.join( basedir, value )
                        found = True
                    if not found and  os.path.isfile( value ):
                        found = True
                    if not found:
                        raise NameError( "File '" + value + "' does not exist" )
                    if not os.access( value, os.R_OK ):
                        raise NameError( "File '" + value + "' cannot be read" )
                if( self.writable and self.pathcheck ):
                    if basedir:
                        value = os.path.join( basedir, value )
                    try:
                        f = open( value, "a+" );
                        f.close();
                    except:
                        raise NameError( "File '" + value + "' is not writable" )
                if( self.check ):
                    value = self.check(value, basedir=basedir)

                l.append( value )

            # Test to see if the file exists
            if len(l) != 1 and not self.multiple:
                raise NameError( "Value must be a single filename" )
            if not self.multiple:
                self.value = l[0]
            else:
                self.value = l
            return self.value


    class Value:
        def __init__( self, default, datatype, valid_range, list_len, units=None ):
            self.default     = default
            self.valid_range = valid_range
            self.list_len    = list_len
            self.data_type   = datatype
            self.units       = units

        def validate( self, value_list, basedir=None ):
            valid_range = self.valid_range

            if type(value_list) is not list and type( value_list) is not tuple:
                value_list=[value_list]


            l=[]
            if self.list_len != len( value_list ):
                if( self.list_len == 0 ):
                    raise NameError( "Value must be a scalar" )
                else:
                    raise NameError( "Value must be a vector with " + str( self.list_len ) + " elements" )

            for value in value_list:
                if( self.data_type == Command.TYPE_FLOAT ):
                    strfudge = "a number"
                elif( self.data_type == Command.TYPE_INT ):
                    strfudge = "an integer"

                try:
                    try:
                        if( self.data_type == Command.TYPE_FLOAT ):
                            value = float(value)
                        if( self.data_type == Command.TYPE_INT ):
                            value = int(value)
                    except:
                        raise NameError( "Value is not " + strfudge );

                    if( valid_range == Command.RANGE_POS and value <=0 ):
                        raise NameError( "Value must be " + strfudge + " > 0" )
                    if( valid_range == Command.RANGE_0POS and value <0 ):
                        raise NameError( "Value must be " + strfudge + " >=0" )
                    if( valid_range == Command.RANGE_NEG and value >=0 ):
                        raise NameError( "Value must be " + strfudge + " < 0" )
                    if( valid_range == Command.RANGE_0NEG and value >0 ):
                        raise NameError( "Value must be " + strfudge + " <=0" )

                except NameError as e:
                    raise NameError( e )
                if( self.list_len == 1 ):
                    l = value
                else:
                    l.append( value )

            self.value = l
            return self.value

        def args(self):
                if( self.data_type == Command.TYPE_FLOAT ):
                    strfudge = "number "
                elif( self.data_type == Command.TYPE_INT ):
                    strfudge = "integer "
                if( self.valid_range == Command.RANGE_POS  ):
                    strfudge = strfudge + "> 0"
                if( self.valid_range == Command.RANGE_0POS  ):
                    strfudge = strfudge + ">= 0"
                if( self.valid_range == Command.RANGE_NEG  ):
                    strfudge = strfudge + "< 0"
                if( self.valid_range == Command.RANGE_0NEG  ):
                    strfudge = strfudge + "<= 0"
                if self.list_len > 1:
                    strfudge = str(self.list_len) + " x " + strfudge;

                xx= "[ "+strfudge+" ]"

                return xx
