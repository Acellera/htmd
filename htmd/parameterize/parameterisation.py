import datetime
import time
import os
import logging as log
import tempfile
from htmd.parameterize.step import Step
from htmd.parameterize.nonpolar import NonpolarSet
from htmd.parameterize.drudeset import DrudeSet
from htmd.parameterize.match    import MatchSet
from htmd.parameterize.configuration import Configuration
import shutil
import inspect
import htmd.parameterize
import subprocess
import sys
from htmd.progress.progress import ProgressBar
from htmd import Molecule

log = log.getLogger("parameterize")


class Parameterisation:
    def __init__(self, config=None, name=None, run=True ):
        self.completed=False
        resuming = False



        if config and config.JobName:
            name =os.path.basename( config.JobName ) #sanitise
            wdir = os.path.join(Parameterisation.prefix(), name )
            if os.path.exists(wdir):
                raise ValueError( "Job [" + config.JobName + "] already exists" )

        elif name:
            name =os.path.basename(name) #sanitise
            wdir = os.path.join(Parameterisation.prefix(), name )

            if os.path.exists(wdir):
                config = Configuration( os.path.join( wdir, "configuration"), check=False )
                log.info("Resuming job [" + config.JobName + "]")
                resuming = True
            else:
                raise ValueError("Invalid job name [" + name + "] no such directory [" +  wdir +"]")

        else:
            raise ValueError("Must specify a config (with name) or a name to resume")



        if config.Model == "Drude":
            self.stepset = DrudeSet()
        elif config.Model == "Nonpolar":
            self.stepset = NonpolarSet()
        elif config.Model == "Match":
            self.stepset = MatchSet()
        else:
            raise ValueError("Invalid Model setting")
        try:
            env = self._preflight_test( config )
        except ValueError as e:
            log.info( e )
            raise ValueError( "Prerequisite test failed " + e.__str__() )
        except:
            log.info("Preflight test failed")
            raise ValueError("The prerequisites for running a parameterisation could not be found. Use --verbose to see")


        if not config.JobName:
            # No job name specified -- assume this is going to be a new job
            # and create a new unique ID
            config.JobName = Parameterisation.make_job_name()
            log.info("Created new identifier [" + config.JobName + "]")

        wdir = os.path.join(Parameterisation.prefix(), config.JobName)
        if not os.path.exists(wdir):
            os.mkdir(wdir)

        self.directory=wdir
        self._env=env
        self._config=config
        #print("SET : " + self.directory )
        # even if we are resuming, I might have manually deleted these files while debugging
        #if not os.path.exists( os.path.join( wdir, "routines" ) ):
        self.write_environment_file( os.path.join( wdir, "routines" ), env )
        
        self.copy_qm_file( wdir )





        if( not resuming ):
            self.copy_in_configuration_files( config, wdir )
        if run:
            self.run()

    def run(self):
            self.stepset.run( self.directory, param=self )
            self.completed=True


    def copy_in_configuration_files(self, config, wdir):
        # Copy in the files that form part of the job configuration into the working directory


            # Stage in the various input file
            # 1. The input structure -- mol2 or pdb
            log.info("Staging in input mol2/pdb file")
            indir = os.path.join( wdir, "000-input")
            os.mkdir( indir )
            tmp, suffix = os.path.splitext( config.FileName )
            if not (suffix.lower() in [ ".mol2", ".pdb" ]):
                raise ValueError( "Input file must be in mol2 or PDB format [" + suffix + "]")

            shutil.copyfile(config.FileName, os.path.join( indir, "ligand" + suffix.lower()))
            #2. the atom equivalency file
            if config.Torsions:
                self._map_requested_torsions( config.FileName, config.Torsions, os.path.join( indir, "soft-dih-list.txt" ) )
            else:
                open( os.path.join( indir, "soft-dih-list.txt"), 'a').close()

            if config.Equivalent:
                newpath = os.path.join( indir, "equiv-org.txt" )
                shutil.copyfile(config.Equivalent, newpath )
                config.Equivalent = newpath
            else:
                open( os.path.join( indir, "equiv-org.txt"), 'a').close()

            if config.Neutral:
                newpath = os.path.join( indir, "neutral.txt" )
                shutil.copyfile(config.Neutral, newpath )
                config.Equivalent = newpath
            else:
                open( os.path.join( indir, "neutral.txt"), 'a').close()

            if config.FixCharges:
                newpath = os.path.join( indir, "fixq.txt" )
                shutil.copyfile(config.FixCharges,  newpath )
                config.Equivalent = newpath
            else:
                open( os.path.join( indir, "fixq.txt"), 'a').close()

            # Copy in the input file itself
            config.save( os.path.join( wdir, "configuration" ), fmt="pretty" )
            config.save( os.path.join( wdir, "configuration.sh" ), fmt="shell" )


    def copy_qm_file(self, rootdir):
        import shutil
        import inspect
        d = os.path.dirname(inspect.getfile(Parameterisation))
        d = os.path.join( d, "share" )
        d = os.path.join( d, "scripts" )
        d = os.path.join( d, "QM-para.txt" )
        if not os.path.exists(d):
            raise ValueError( "QM-para.txt not found")
        shutil.copyfile( d, os.path.join( rootdir, "QM-para.txt"))

    @staticmethod

    def prefix():
        import os
        dir = os.path.join(os.path.expanduser("~"), ".htmd")
        try:
            os.mkdir(dir)
        except:
            pass

        dir = os.path.join( dir, "gaamp")
        try:
            os.mkdir(dir)
        except:
            pass


        return dir


    @staticmethod
    def make_job_name():
        d = datetime.date.today()
        d = d.isoformat()
        suffix = 0
        prefix = Parameterisation.prefix()
        name = "%s-%05d" % (d, suffix)
        while not os.path.exists(os.path.join(prefix, name)):
            suffix = suffix + 1
            name = "%s-%05d" % (d, suffix)


    @staticmethod
    def listJobs():
        import os
        for i in os.listdir(Parameterisation.prefix()):
            print("\t" + i)
        pass


    @staticmethod
    def deleteJob(job):
        import shutil
        import os
        jobb = os.path.basename(job)
        dir = os.path.join(Parameterisation.prefix(), jobb)
        if not os.path.exists(dir):
            raise ValueError("[" + job + "] ["+dir+"] not found")
        shutil.rmtree(dir)

        return True


    @staticmethod
    def find_binary( bin, path=None, fatal=True ):
        import shutil

        #print( os.environ['PATH'].split(os.pathsep)  )
        if path:
            for pp in os.environ['PATH'].split(os.pathsep):
                path.append(pp)
        else:
            path = os.environ['PATH']

        pp= ""

        for p in path:
            pp = pp + os.pathsep + p
        #print(bin)
        #print(pp)

        r = shutil.which( bin, mode=os.X_OK, path=pp )
        #print(r)
        if r: return r
        if fatal:
          raise ValueError( "Could not find binary [" + bin + "] in PATH " )
        else:
          return None


    def write_environment_file(self, filename, env ):
        fh = open( filename, "w" )
        for i in env.keys():
            print( "export %s=\"%s\"" % ( i, env[i]), file=fh)

        print( "source ../configuration.sh", file=fh )

        print( """

# Diagnostic calling function. Generated from parameterisation.py

CHECK() {
    if [ ! -x $1 ]; then
        echo "Executable [$1] not found"
        exit 3
    fi
    rm -f error.txt
    b=$(basename "$1")
    STRACE=""
    if [ "$Debug" == "True" ]; then
		STRACE=$(which strace 2>/dev/null)
    fi
    if [ -x "$STRACE" ]; then
       strace -e open -o .$b.strace "$@" > .$b.txt
		else		
		   $@	> .$b.txt
		fi
		ret=$?

    if [ "$ret" != "0" ]; then
        echo "Non-zero exit code from $1: "
        if [ -e .$b.txt ]; then
            cat .$b.txt
        fi
        ret=1
    fi
    if [ -e error.txt ]; then
        echo "Error file created by $1:"
        if [ -e .$b.txt ]; then
            cat error.txt
        fi
        ret=2
    fi

    if [ -e .$b.strace ]; then
        err=$(grep -e 'open("[[:alpha:]]' .$b.strace | grep "RD" | grep ENO )
        if [ "$err" != "" ]; then
            echo "Some input files were not found:"
            echo $err
            ret=3
        fi
    fi

    if [ "$ret" != "0" ]; then
        if [ -e .$b.strace ]; then
             echo "\n-- Files Read --"
             grep -e 'open("[[:alpha:]]' .$b.strace | grep "RD"
             echo "\n-- Files Written --"
             grep -e 'open("[[:alpha:]]' .$b.strace | grep "WR"
        fi
        exit $ret
    fi
}
        """, file=fh )

        fh.close()

    @staticmethod
    def _preflight_test( config ):
        import inspect
        import psutil

        env={}
        gaamp_bin = []
        p = os.path.dirname(inspect.getfile(Parameterisation))
        p = os.path.join( p , "share")
        gaamp_bin.append( os.path.join( p , "bin") )
        p = os.path.join( p , "match")
        gaamp_bin.append( os.path.join( p , "bin") )

        # Locate all of  the third-party binary programs needed

        env['BIN_BABEL']       = Parameterisation.find_binary( "htmd_babel", path=gaamp_bin )
        env['BIN_MATCH']       = Parameterisation.find_binary( "match", path=gaamp_bin )
        #env['BIN_ANTECHAMBER']= Parameterisation.find_binary( "antechamber", path=gaamp_bin )
        has_g09 = False
        has_psi4=False
        try:
           Parameterisation.find_binary( "g09", path=gaamp_bin )
           has_g09=True
        except:
           pass 
        try:
           Parameterisation.find_binary( "psi4", path=gaamp_bin )
           has_psi4=True
        except:
           pass 

        if not has_g09 and not has_psi4:
          raise ValueError( "Neither G09 nor PSI4 found. At least one QM code is required" )          
        env['BIN_G09']         = Parameterisation.find_binary( "g09_wrapper", path=gaamp_bin )
        #env['BIN_GNUPLOT']     = Parameterisation.find_binary( "gnuplot", path=gaamp_bin )

        env['BIN_CGRID']       = Parameterisation.find_binary( "cgrid", path=gaamp_bin )
#        env['BIN_CGRID_DRUDE'] = Parameterisation.find_binary( "drude_cgrid", path=gaamp_bin )
        env['BIN_ADD_TIP3']    = Parameterisation.find_binary( "add_tip3", path=gaamp_bin )
        env['BIN_GEN_XPSF']    = Parameterisation.find_binary( "gen_xpsf", path=gaamp_bin )
#        env['BIN_PDB_TO_CRD']  = Parameterisation.find_binary( "pdb_to_crd", path=gaamp_bin )
        env['BIN_EQUIV_ATOM']  = Parameterisation.find_binary( "equiv_atom", path=gaamp_bin )
        env['BIN_GEN_ESP']     = Parameterisation.find_binary( "gen_esp", path=gaamp_bin )
        env['BIN_CHECK_B0_THETA0'] = Parameterisation.find_binary( "check_b0_theta0", path=gaamp_bin )
        env['BIN_FITCHARGE'] = Parameterisation.find_binary( "fitcharge", path=gaamp_bin )
        env['BIN_UPDATE_XPSF'] = Parameterisation.find_binary( "update_xpsf", path=gaamp_bin )
        env['BIN_UPDATE_TOR_PARA'] = Parameterisation.find_binary( "update_tor_para", path=gaamp_bin )
        env['BIN_ACCEPTOR']    = Parameterisation.find_binary( "acceptor", path=gaamp_bin )
        env['BIN_DONOR']       = Parameterisation.find_binary( "donor"   , path=gaamp_bin )
        env['BIN_DON_ACC_SEPARATE'] = Parameterisation.find_binary( "don_acc_separate"   , path=gaamp_bin )
        env['BIN_EXTRACT_QM_MOL_WATER'] = Parameterisation.find_binary( "extract_qm_mol_water.py", path=gaamp_bin )
        env['BIN_EXTRACT_QM_ROTAMER']= Parameterisation.find_binary( "extract_qm_rotamer.py", path=gaamp_bin )
        env['BIN_GEN_SOFT_LIST']  = Parameterisation.find_binary( "gen_soft_list", path=gaamp_bin )
        env['BIN_FITCHARGE_AGAIN']= Parameterisation.find_binary( "fitcharge_again", path=gaamp_bin )
        env['BIN_MM_PES']         = Parameterisation.find_binary( "mm_pes", path=gaamp_bin )
        env['BIN_MM_PES_LARGE']   = Parameterisation.find_binary( "mm_pes_large", path=gaamp_bin )
        env['BIN_CLUSTERING_PHI'] = Parameterisation.find_binary( "clustering_phi", path=gaamp_bin )
        env['BIN_QM_1D_SCAN_PARA']= Parameterisation.find_binary( "qm_1d_scan_para",path=gaamp_bin )
        env['BIN_QM_1D_SCAN_LARGE_PARA']= Parameterisation.find_binary( "qm_1d_scan_large_para",path=gaamp_bin )
        env['BIN_RUN_QM_JOBS']    = Parameterisation.find_binary( "run_qm_jobs", path=gaamp_bin )
        env['BIN_QM_ROTAMER_SCAN_LARGE']= Parameterisation.find_binary( "qm_rotamer_scan_large", path=gaamp_bin )
        env['BIN_QM_ROTAMER_SCAN']= Parameterisation.find_binary( "qm_rotamer_scan", path=gaamp_bin )
        env['BIN_1D_FITTING']     = Parameterisation.find_binary( "1d_fitting",path=gaamp_bin )
        env['BIN_1D_FITTING_ASYMMETRIC']     = Parameterisation.find_binary( "1d_fitting_asymmetric",path=gaamp_bin )
        env['BIN_1D_ORG']     = Parameterisation.find_binary( "1d_org",path=gaamp_bin )
        env['BIN_1D_ROTAMER_FITTING'] = Parameterisation.find_binary( "1d_rotamer_fitting", path=gaamp_bin )
        env['BIN_PLOT_1D_RESULT_AND_ORG'] = Parameterisation.find_binary( "plot_1d_result_and_org", path=gaamp_bin )
        env['BIN_QSUB'] = Parameterisation.find_binary( "qsub", path=gaamp_bin, fatal=False )
        env['BIN_QSTAT'] = Parameterisation.find_binary( "qstat", path=gaamp_bin, fatal=False )
        env['BIN_BSUB'] = Parameterisation.find_binary( "bsub", path=gaamp_bin, fatal=False )
        env['BIN_BJOBS'] = Parameterisation.find_binary( "bjobs", path=gaamp_bin, fatal=False )

#        print(env)
        if config:
            if ( config.ExecutionMode == "PBS" and ( env['BIN_QSUB'] is None  or  env['BIN_QSTAT'] is None ) ):
                raise ValueError( "Dependencies for PBS execution not found. Ensure QSUB and QSTAT are available" );
            if ( config.ExecutionMode == "LSF" and ( env['BIN_BSUB'] is None  or  env['BIN_BJOBS'] is None ) ):
                raise ValueError( "Dependencies for LSF execution not found. Ensure BSUB and BJOBS are available" );

        return env


    def showMolecule(self):
        m = self.getMolecule()
        return m.view( viewer="ngl", sel="all", style="cpk" )


    @staticmethod
    def listDihedrals( filename ):
           # This routine gets the names of the atoms in the soft dihedrals
           # It's so horrible, since we have to jump through hoops to 
           # set up the input for gen_soft_list
           ll1=[]
           ll2=[]
           env = Parameterisation._preflight_test( None )
           mol = Molecule( filename );
           mol = Parameterisation._rename_mol( mol )
           mol.guessBonds()
        # make a tempdir
           with tempfile.TemporaryDirectory() as td:
             #td=tempfile.mkdtemp()
             #print(td)
             pwd = os.getcwd() 
             os.chdir( td )
             f=open("mol.prm", "w")
             f.close()
             mol.write( "mol.pdb" )
             mol.write( "mol-opt.xyz" )
             mol.write( "mol.xpsf", type="psf" )
             for charge in range(-1,2): 
             	ret =subprocess.check_output( [ env['BIN_MATCH'], "-charge", str(charge), "-forcefield",  "top_all36_cgenff_new",  "mol.pdb" ], stderr = subprocess.STDOUT, shell=False, stdin = None )
             ret =subprocess.check_output( [ env['BIN_GEN_XPSF'], "mol.rtf", "mol.xpsf", "MOL" ], stderr = subprocess.STDOUT, shell=False, stdin = None )
             ret =subprocess.check_output( [ env['BIN_GEN_SOFT_LIST'] ], stderr = subprocess.STDOUT, shell=False, stdin = None )

             f = open( "soft-dih-list.txt", "r" )
             ff = f.readlines()
             for l in ff:
               ss=[]
               tt=[]
               for m in l.split():
                 tt.append( int(m) )
                 ss.append( mol.name[ int(m)-1 ].strip().upper() )
               ll1.append(tt)
               ll2.append(ss)

             f.close()			
             os.chdir( pwd )

           return (ll1, ll2)

    @staticmethod
    def renameStructure( filename ):
        m = Molecule( filename )
        m = Parameterisation._rename_mol( m )
        m.write( filename )

    def _rename_mol(  mol ):
        # This fixes up the atom naming and reside name to be consistent
        # NB this scheme matches what MATCH does. Don't change it
        # Or the naming will be inconsistent with the RTF
        import re

        hh=dict()

        for i in range(len(mol.name)):
            t=re.sub( '[1234567890]*', "", mol.name[i] ).upper()
            idx=0

            if not t in hh:
                hh[t]=idx

            idx = hh[t]+1;
            hh[t]=idx

            t=t + str(idx)
            mol.name[i]=t
            mol.resname[i]="MOL"
        return mol
 
    def getMolecule(self):
        import htmd
        if not self.completed:
            raise ValueError("Parameterisation is not complete")

        # look for the right files  in the output directory
        dir = os.path.join( self.directory, "999-results" )
        dir = os.path.join( dir, "mol.xyz")
        if not os.path.exists(dir):
            raise ValueError("Parameterisation did not produced a minimised structure")
        mol = htmd.Molecule( dir )
        Parameterisation._rename_mol( mol )
        return mol

    def showDihedral(self, index, viewer="ngl" ):
        dih = self.getDihedralIndexes()
        if index<0 or index >= len(dih):
            raise ValueError( "Dihedral index out of range")
        m = self.getMolecule()

        i0=str(dih[index][0]-1)
        i1=str(dih[index][1]-1)
        i2=str(dih[index][2]-1)
        i3=str(dih[index][3]-1)
        m.view( sel="all", style="lines", hold=True )
        m.view( style="cpk", sel="index "+i0, hold=True )
        m.view( style="cpk", sel="index "+i1, hold=True )
        m.view( style="cpk", sel="index "+i2, hold=True )
        return m.view( style="cpk", sel="index "+i3, viewer=viewer )


    def getDihedralIndexes(self, byname=False):
        if not self.completed:
            raise ValueError("Parameterisation is not complete")

        dir = os.path.join( self.directory, "999-results" )

        sdl = os.path.join( dir, "soft-dih-list.txt")

        if not os.path.exists( sdl ):
            raise ValueError( "Parameterisation did not produce any graphs " + sdl )

        if byname:
            m=self.getMolecule()

        f= open( sdl, "r" )
        fl=f.readlines()
        f.close()
        sdl=[]
        for l in fl:
            dl=[]
            l=l.strip().split( )
            for i in range(4):
                if byname:
                    dl.append( m.name[ int(l[i])-1 ])
                else:
                    dl.append( int(l[i]))
            sdl.append(dl)

        return sdl

    def getParameters(self):
        import tempfile
        if not self.completed:
            raise ValueError("Parameterisation is not complete")
        # look for the right files  in the output directory

        dir = os.path.join( self.directory, "999-results" )
        rtf = os.path.join( dir , "mol.rtf" )
        prm = os.path.join( dir , "mol.prm" )
        xyz = os.path.join( dir , "mol.xyz" )
        pdb = os.path.join( dir , "mol.pdb" )

    
        rtf_tmp = tempfile.mkstemp( suffix=".rtf" )
        prm_tmp = tempfile.mkstemp( suffix=".prm" )
        pdb_tmp = tempfile.mkstemp( suffix=".pdb" )
        xyz_tmp = tempfile.mkstemp( suffix=".xyz" )
        os.close(rtf_tmp[0])
        os.close(prm_tmp[0])
        os.close(pdb_tmp[0])
        os.close(xyz_tmp[0])

        shutil.copyfile( rtf, rtf_tmp[1])
        shutil.copyfile( prm, prm_tmp[1])
        shutil.copyfile( xyz, xyz_tmp[1])

        # The Output minimised structure is in XYZ format
        # Need to turn it into a PDB with correct atom naming
        mol=Molecule( xyz )
        Parameterisation._rename_mol(mol) # Canonicalise atom and reside naming
        mol.write( pdb_tmp[1] )

        return {
            "RTF" : rtf_tmp[1],
            "PRM" : prm_tmp[1],
            "PDB" : pdb_tmp[1]
        }


    def plotDihedrals(self, show=True):
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import numpy as np

        fn=[] 
        names=[]
        if not self.completed:
            raise ValueError("Parameterisation is not complete")
        if not show:
            mpl.use('Agg')            
		
        # look for the right files  in the output directory

        dir = os.path.join( self.directory, "999-results" )

        sdl = self.getDihedralIndexes()


        #print(sdl)

 #       fig = plt.figure()

        m=self.getMolecule()
        for i in range(len(sdl)):
            fh = plt.figure()
#            ax1 = fig.add_subplot(
#                 len(sdl), 1, 1+ i ,
#                 ylabel="kcal/mol"
#            )
            ax1=fh.gca()
            ax1.set_xlim( -180, 180)
            ax1.set_xticks( [ -180, -135, -90 , -45, 0, 45, 90, 135, 180 ] )

#            if( i == len(sdl)-1):
            ax1.set_xlabel("phi")
            ax1.set_ylabel("kcal/mol")
            dl=sdl[i]
            title= m.name[ dl[0]-1 ] +" " + m.name[ dl[1]-1 ] + " " + m.name[ dl[2]-1 ] + " " + m.name[dl[3]-1]
            ax1.set_title(title)
            names.append(title)
#            qm =np.loadtxt( os.path.join( dir, "fitting-1d-" + str(i+1) + ".dat" ) )
            qm =np.loadtxt( os.path.join( dir, "1d-qm-mm-" + str(i+1) + ".dat" ) )
            org=np.loadtxt( os.path.join( dir, "org-1d-" + str(i+1) + ".dat") )

            # Remove any values which are > 20kcal/mol from the minimum
            # and so would have been excluded from the fit
            min = np.amin( qm[:,1] )
            qmx=[]
            qmy=[]
            for i in range( qm[:,1].shape[0] ):
               if qm[i,1] <= (min +20):
                 qmx.append( qm[i,0] )
                 qmy.append( qm[i,1] )

            ax1.plot( qmx, qmy, label="QM", color="r", marker="o" )
            ax1.plot( qm[:,0], qm[:,2], label="MM Fitted", color="g", marker="o" )
            ax1.plot( org[:,0], org[:,2], label="MM Original", color="b", marker="o" )
           # if i==0:
            ax1.legend( prop={'size':8} )
#        plt.tight_layout()
            if show:
               plt.show()
            tf = tempfile.mkstemp( suffix=".svg" )
            os.close(tf[0])
            plt.savefig( tf[1], format="svg" )
            fn.append(tf[1])
        return (fn, names)
   
    def run_qm_jobs( self, directory ):
        if self._config.ExecutionMode == "Inline":
           self._run_qm_jobs_inline( directory )
        elif self._config.ExecutionMode == 'PBS':
           self._run_qm_jobs_pbs( directory )
        elif self._config.ExecutionMode == 'LSF':
           self._run_qm_jobs_lsf( directory )
        else:
           raise ValueError( "Unknown exection mode" )

    def _run_qm_jobs_inline( self, directory ):

        cmd = self._env['BIN_G09']

        fni=[]
        fno=[]
        for root, dirs, files in os.walk( directory ):
             for f in files:
                 if(f.endswith(".gjf")):
                      op = f.replace(".gjf", ".out")
                      if not os.path.exists( os.path.join( root, op ) ):
                            fni.append( os.path.join( root, f  ) )
                            fno.append( os.path.join( root, op ) )

        if( len(fni) ):
          bar = ProgressBar( len(fni), description = "Running QM Calculations" )
          for i in range(len(fni)):
             subprocess.check_output( [ cmd, fni[i], fno[i] ] )
             bar.progress()
          bar.stop()


    def _map_requested_torsions( self, ligand, torsionlist, outputfile ):
#config.FileName, config.Torsions, os.path.join( indir, "soft-dih-list.txt" )
          torsions = Parameterisation.listDihedrals( ligand )
          fh=open( outputfile, "w" )
          for t in torsionlist:
            try:
              ii = torsions[1].index(t)
            except:
              raise ValueError( "Torsion [" + str(t) + "] doesn't exist" );
            index = torsions[0][ii] 
            print( "%d %d %d %d" % ( index[0], index[1], index[2], index[3] ), file=fh )
              
          fh.close()

    def _run_qm_jobs_pbs( self, directory ):

        cmd = self._env['BIN_G09']

        fni=[]
        fno=[]
        for root, dirs, files in os.walk( directory ):
             for f in files:
                 if(f.endswith(".gjf")):
                      op = f.replace(".gjf", ".out")
                      if not os.path.exists(op):
                            fni.append( os.path.join( root, f  ) )
                            fno.append( os.path.join( root, op ) )

        if( len(fni) ):
          # Make a PBS script
          for i in range(len(fni)):
             fpbs = fni[i] + ".pbs"
             if not os.path.exists(fpbs):
               f = open(fpbs, "w")
               print("#PBS -lselect=1:ncpus=%d:mem=%dgb" % ( self._config.NCORES, self._config.MEMORY), file=f )
               print("#PBS -lwalltime=24:0:0\n", file=f )
               print("cd \"%s\"" % ( directory ), file=f )
               print("\"%s\" %s %s" % ( cmd, fni[i], fno[i] ), file=f )
               f.close()


             # Look to see if there is already a job submitted
             # If not, qsub it
             fpbsstate= fni[i] + ".jobid"
             if not os.path.exists(fpbsstate):
                # Qsub, saving jobid to file
                subprocess.check_output( "\""+ self._env['BIN_QSUB'] +"\" \"" + fpbs + "\" > \"" + fpbsstate +"\"", shell=True )  

          # Finally monitor progress. Continue until all jobs have produced an output
          # NB TODO FIXME: should also poll qstat to see if job is still live
          bar = ProgressBar( len(fni), description = "Running QM Calculations" )
          complete=False
          lastcount=0
          while not complete:
            count=0
            complete=True
            for i in fno:
               #print(" Checking [" + i +"]" )
               #print( os.path.exists(i) )
               #print( os.access(i, os.R_OK) )
               try:
                 os.stat(i) # Try to flush any cache (for NFS)
               except:
                 pass
               if os.access( i, os.R_OK ):
                 #print("FOUND")
                 count=count+1
               else:
                 complete=False

            #print( str(count) + " completed of " + str(len(fno)) )
            while( lastcount < count ):
               bar.progress()
               lastcount=lastcount+1;
            time.sleep(10);
          bar.stop()
          time.sleep(5) # A bit of time for any outputfile to complete writing
             
    def _run_qm_jobs_lsf( self, directory ):

        cmd = self._env['BIN_G09']

        fni=[]
        fno=[]
        for root, dirs, files in os.walk( directory ):
             for f in files:
                 if(f.endswith(".gjf")):
                      op = f.replace(".gjf", ".out")
                      if not os.path.exists(op):
                            fni.append( os.path.join( root, f  ) )
                            fno.append( os.path.join( root, op ) )

        if( len(fni) ):
          # Make an LSF  script
          for i in range(len(fni)):
             fpbs = fni[i] + ".lsf"
             if not os.path.exists(fpbs):
               f = open(fpbs, "w")

               print("#BSUB -n %d" % ( self._config.NCORES), file=f )
               print("#BSUB -R \"span[ptile=%d]\"" % ( self._config.NCORES), file=f )
               print("#BSUB -W 24:00", file=f )
               print("#BSUB -J gaussian", file=f )
               print("#BSUB -app gaussian", file=f )
               print("#BSUB -o /dev/null", file=f )
               print("#BSUB -M %d000" % (self._config.MEMORY), file=f )
               print("\nmodule load gaussian\n", file=f )

               print("cd \"%s\"" % ( directory ), file=f )
               print("\"%s\" %s %s" % ( cmd, fni[i], fno[i] ), file=f )
               f.close()


             # Look to see if there is already a job submitted
             # If not, qsub it
             fpbsstate= fni[i] + ".jobid"
             if not os.path.exists(fpbsstate):
                # Qsub, saving jobid to file
                subprocess.check_output( "\""+ self._env['BIN_BSUB'] +"\" < \"" + fpbs + "\" > \"" + fpbsstate +"\"", shell=True )  

          # Finally monitor progress. Continue until all jobs have produced an output
          # NB TODO FIXME: should also poll qstat to see if job is still live
          bar = ProgressBar( len(fni), description = "Running QM Calculations" )
          complete=False
          lastcount=0
          while not complete:
            count=0
            complete=True
            for i in fno:
               #print(" Checking [" + i +"]" )
               #print( os.path.exists(i) )
               #print( os.access(i, os.R_OK) )
               try:
                 os.stat(i) # Try to flush any cache (for NFS)
               except:
                 pass
               if os.access( i, os.R_OK ):
                 #print("FOUND")
                 count=count+1
               else:
                 complete=False

            #print( str(count) + " completed of " + str(len(fno)) )
            while( lastcount < count ):
               bar.progress()
               lastcount=lastcount+1;
            time.sleep(10);
          bar.stop()
          time.sleep(5) # A bit of time for any outputfile to complete writing
             




def installExamples():
        import shutil
        import inspect
        d = os.path.dirname(inspect.getfile(Parameterisation))
        d = os.path.join( d, "share" )
        d = os.path.join( d, "examples" )
        if not os.path.exists(d):
            raise ValueError( "Examples not found")
        dest = os.path.join( os.path.expanduser( "~" ), "gaamp_examples" )
        shutil.copytree( d, dest )
        return dest		


	
