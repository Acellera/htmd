from __future__ import print_function

from enum import Enum
import numpy as np
import subprocess
from shutil import which 
import shutil
import os
import time
from math import pi as PI
from math import sqrt
from math import acos
from math import cos
from math import sin
from numpy.random import uniform as rand

from htmd.molecule.vdw import VDW
from htmd.molecule.molecule import Molecule
from htmd.progress.progress import ProgressBar

from htmd.apps.lsf import LSF

class BasisSet(Enum):
 _6_31G_star = 1000
 _cc_pVTZ    = 1001

class Theory(Enum):
 HF = 2000
 DFT= 2001

class Code(Enum):
 PSI4     = 3000
 Gaussian = 3001
 TeraChem = 3002

class Execution(Enum):
 Inline = 4000
 LSF    = 4001

class QMResult:
  completed  = False
  errored    = False
  coords     = None
  esp_points = None
  esp_scalar = None
  energy     = 0.00
  dipole     = None
  quadrupole = None 
  mulliken   = None

class QMCalculation:
  
  def __init__( self, molecule, 
                basis  = BasisSet._6_31G_star,
                theory = Theory.HF,
                charge = 0, 
                multiplicity=1, 
                frozen = None,
                optimize = False,
                esp    = False,
                esp_vdw_radii  = [1.4, 1.6, 1.8, 2.0, 2.2],
                esp_density= 10.,
                code   = None,
                execution = Execution.Inline,
                directory = ".", mem=2, ncpus=-1
              ):

    if not isinstance(basis, BasisSet ): raise ValueError( "basis must of type BasisSet" )
    if not isinstance(theory, Theory ) : raise ValueError( "theory must of type Theory" )
    if code and (not isinstance(code, Code )): raise ValueError( "code must be of type Code" )
    if not isinstance(execution, Execution ) : raise ValueError( "execution must be of type Execution" )

    if not isinstance( mem, int)             : raise ValueError( "mem must be integer" )
    if not isinstance( ncpus, int)           : raise ValueError( "ncpus must be integer" )
    if not isinstance( charge, int)          : raise ValueError( "charge must be integer" )
    if not isinstance( multiplicity, int )   : raise ValueError( "multiplicity must be integer" )
    if not isinstance( optimize, bool )      : raise ValueError( "optimize must be boolean" )
    if not isinstance( esp     , bool ) and not isinstance( esp, np.ndarray )      : raise ValueError( "esp must be boolean or ndarray" )
    if not isinstance( esp_density , float ) : raise ValueError( "esp_density must be float" )

    if( ncpus == -1 ):
      try:
        ncpus = int( os.getenv("NCPUS") )
      except:
        pass
    if( ncpus == -1 ):
      ncpus = os.cpu_count()



    # TODO esp validation, etc

    self.molecule= molecule.copy()
    self.basis   = basis
    self.theory  = theory
    self.charge  = charge
    self.multiplicity = multiplicity
    self.frozen  = None
    self.esp     = esp
    self.esp_vdw_radii = esp_vdw_radii
    self.esp_density = esp_density
    self.optimize= optimize
    self.directory= directory
    self.ncpus   = ncpus
    self.mem     = mem
    self.execution = execution
    self.code    = code
    self._results = []
    self.natoms  = self.molecule.coords.shape[0]

    self.psi4_binary = shutil.which( "psi4", mode = os.X_OK )
    self.gaussian_binary = shutil.which( "g03", mode = os.X_OK )
    if not self.gaussian_binary:
       self.gaussian_binary = shutil.which( "g09", mode = os.X_OK )
    if (not self.gaussian_binary) and ( not self.psi4_binary ):
       raise RuntimeError("Can not find neither Gaussian nor PSI4" )
    if self.code == None:
       if self.gaussian_binary: self.code = Code.Gaussian
       elif self.psi4_binary:   self.code = Code.PSI4
    if self.code == Code.PSI4    and (not self.psi4_binary )     : raise RuntimeError( "PSI4 not found" )
    if self.code == Code.Gaussian and (not self.gaussian_binary ): raise RuntimeError( "Gaussian not found" )

    # Set up point cloud if esp calculation requested
   
    if isinstance( self.esp, np.ndarray ):
       self.points = [ self.esp ]
       if self.points[0].shape[1]!=3: 
          raise ValueError( "ESP point array must be npoints x 3" )
       if self.molecule.coords.shape[2] != 1:
          raise ValueError( "Can only specift ESP point array with a single frame of coords" )
    elif self.esp == True:
       self.points = self._generate_points()
    else:
       self.points = None

    for i in range(self.molecule.coords.shape[2]):
       qmr = QMResult()
       qmr.coords     = self.molecule.coords[:,:,i]
       qmr.esp_points = np.squeeze(np.asarray(self.points))
       self._results.append( qmr )

    # For frozen dihedrals, map atom names to 1-based atom indices    
    if frozen:
       self.frozen=[]
       for i in frozen:
         if type(i[0]) ==  str:
           a1 = self.molecule.name.tolist().index( i[0] ) + 1
           a2 = self.molecule.name.tolist().index( i[1] ) + 1
           a3 = self.molecule.name.tolist().index( i[2] ) + 1
           a4 = self.molecule.name.tolist().index( i[3] ) + 1
           self.frozen.append( [a1,a2,a3,a4] )
         else:
           self.frozen.append( [i[0]+1, i[1]+1, i[2]+1, i[3]+1 ] )
  
    self._prepare()

  def _generate_points( self ):
    points = [] 
#np.zeros((0,3, self.molecule.coords.shape[2]), dtype=np.float32)
    for frame in range(self.molecule.coords.shape[2]):
      pp =  self._points( 
        self.molecule.coords[:,:,frame], 
        self._vdw_radii( self.molecule.element ),  
        self.esp_vdw_radii, self.esp_density 
      )
      points.append( pp )
    return points
  
  def _write_points( self, filename, points ):
    # Save out a list of points to a named file
    f = open(filename, "w" )
    for p in points:
      print( "%f %f %f" %( p[0], p[1], p[2]), file=f )
    f.close()

  def _vdw_radii( self, elements ):
    radii = np.zeros( ( elements.shape[0] ), dtype=np.float32 )
    i=0
    for e in elements:
      radii[i] = VDW.radiusByElement( e ) 
      i=i+1
    return radii

  def _points( self, coords, radii, multipliers, density ):
    points=[]
    np.random.seed(0)
    # Make a set of points in a vdw shell around each atom
    for m in multipliers:
      for i in range(coords.shape[0]):
         p = self._rand_sphere_sample( coords[i,:], radii[i] * m, density )
         # remove any points that are within radii[i]*m of i-th atom
         for pp in p:
           too_close = False
           for j in range(coords.shape[0]):
              if self._dist( coords[j,:], pp ) < radii[j]*m :
                too_close = True
                break
           if not too_close:
              points.append( pp )
             
    return  np.asarray( points, dtype=np.float32 )

  def _dist(self, a, b):
    c = a-b
    return sqrt(c.dot(c))

  def _dist2(self, a, b):
    c = a-b
    return c.dot(c)



  def _rand_sphere_sample( self, centre, r, density ):
    # Produce a set of points on the sphere of radius r centred on centre
    # with ~density points / unit^2


    surface_area = 4./3. * PI * r * r
    n_points     = (int) (density * surface_area)
    area_per_point = 1./density  # surface_area / n_points
    mindist      = sqrt( area_per_point / PI )

    points       = np.zeros( ( n_points, 3 ) )

    i=0;
    mindist2 = mindist * mindist
    pos = np.zeros((3))
    while( i < n_points ):
       z = 2. * rand() -1.
       lon = 2. * PI * rand();
       lat = acos(z)  
       x = cos(lon) * sin(lat)
       y = sin(lon) * sin(lat)

       pos[0] = x * r
       pos[1] = y * r 
       pos[2] = z * r

       # Crudely test to see if it is in range of other points 
       too_close = False
       for j in range(i):
         if self._dist2( points[j, :], pos ) < (mindist2)  : 
           too_close=True
           break
       if not too_close:
         points[i,:] = pos
         i=i+1
    points = points
    points = points + centre  
    return points

  def _prepare( self ):
    dirs = []
    try: 
      os.mkdir( self.directory )
    except:
      pass
    i=0
    for c in range(self.molecule.coords.shape[2]):
      dn = os.path.join( self.directory, "%05d" % ( i ) )
      if not os.path.exists( dn ):
          os.mkdir( dn )
      dirs.append(dn)
      # Write input files for both PSI4 and Gaussian
      self._write_xyz( dn, c )	
      self._write_psi4( dn, c )
      self._write_gaussian( dn, c )
      i=i+1
      self._job_dirs = dirs;
    pass
    self._start( dirs )

  def _start( self, directories ):
     if self.execution == Execution.Inline:
       self._start_inline( directories )
     elif self.execution == Execution.LSF:
       self._start_lsf( directories )
     else:
       raise RuntimeError( "Invalid execution mode" )

  def _start_lsf( self, directories ):
     print( "Running QM Calculations via LSF" )

     to_submit=[]
     for directory  in directories:
       cwd=os.getcwd()
       try: 
           if self.code == Code.Gaussian: 
             if not os.path.exists( os.path.join( directory, "output.gau" ) ):
               to_submit.append( directory ) 
           elif self.code == Code.PSI4: 
             if not os.path.exists( os.path.join( directory, "psi4.out" ) ):
               to_submit.append( directory ) 
       except:
         raise

     if self.code == Code.Gaussian: 
       cmd =  '"' + self.gaussian_binary + '" < input.gjf > output.gau 2>&1'
     elif self.code == Code.PSI4:
       cmd =  '"' + self.psi4_binary + '" -i psi4.in -o psi4.out 2>&1'

     lsf = LSF( ncpus=self.ncpus, executable = cmd, queue = "general", resources = "span[ptile=%d]" % (self.ncpus), app = "gaussian")
     lsf.submit( to_submit )
     time.sleep(5)
     lsf.wait()


  def _start_inline( self, directories ):
     bar = ProgressBar( len(directories), description="Running QM Calculations" )

     for directory  in directories:
       cwd=os.getcwd()
       try: 
           os.chdir( directory )
           if self.code == Code.Gaussian: 
             if not os.path.exists( "output.gau" ):
               subprocess.call( '"' + self.gaussian_binary + '" < input.gjf > output.gau 2>&1', shell=True )
           elif self.code == Code.PSI4: 
             if not os.path.exists( "psi4.out" ):
               subprocess.call( [ self.psi4_binary, "-i", "psi4.in", "-o", "psi4.out" ] )
       except:
         os.chdir(cwd)
         raise
       os.chdir(cwd)
       bar.progress()
     bar.stop()

  def _complete( self ):
    pass

  def results( self ): 
    self._complete()
    i=0
    for dn in self._job_dirs:
       self._results[i].completed=True
       #print(dn)
       if( self.code == Code.PSI4 ):
         ret = self._read_psi4( dn )  
       elif( self.code == Code.Gaussian ):
         ret = self._read_gaussian( dn )  

       if( ret == None ):
         self._results[i].errored = True
       elif( (not ( "energy" in ret )) or (ret['energy'] == 0.) ):    
         self._results[i].errored = True
       else:
         self._results[i].energy = ret['energy'] * 627.509469 # Hartree to kcal
         self._results[i].coords = np.atleast_3d(ret['coords'])
         if "gridesp" in ret:
            self._results[i].esp_scalar = ret['gridesp']
            self._results[i].esp_scalar = np.divide( self._results[i].esp_scalar,  0.529177249  ) # Unit conversion from bohrs to angstoms
         if "dipole" in ret:
            self._results[i].dipole     = ret['dipole']
         if "quadrupole" in ret:
            self._results[i].quadrupole = ret['quadrupole']
         if "mulliken" in ret:
            self._results[i].mulliken   = ret['mulliken']
       i=i+1      

    return self._results

  def _read_gaussian( self, dirname ):
   try:
    f = open( os.path.join( dirname, "output.gau" ), "r" )
    fl=f.readlines()
    data={}

    completed = False

    for i in range(len(fl)):
      if "Dipole moment (field-independent basis, Debye):" in fl[i]:
        s = fl[i+1].split()
        data['dipole'] = [ float(s[1]), float(s[3]), float(s[5]), float(s[7]) ]
      if "Traceless Quadrupole moment (field-independent basis, Debye-Ang):" in fl[i]:
        s1 = fl[i+1].split()
        s2 = fl[i+2].split()
        data['quadrupole'] = [ float(s1[1]), float(s1[3]), float(s1[5]), float(s2[1]), float(s2[3]), float(s2[5]) ]
      if "Mulliken atomic charges:" in fl[i]:
        data['mulliken'] = []
        for j in range(self.natoms):
          data['mulliken'].append( float( fl[i+1+j].split()[2] ) )

    for l in fl:
        if "SCF Done:  E(RHF) = " in l:
            ff = l.split()
            data['energy'] = float(ff[4])
            completed=True
    i=0
    while(i<len(fl)):
        if "Number     Number       Type             X           Y           Z" in fl[i]:
           i=i+2
           data['coords'] = np.zeros((self.natoms, 3))
           for j in range( self.natoms ):
              ff = fl[j+i].split()
              data['coords'][j,0]=float( ff[3] ) 
              data['coords'][j,1]=float( ff[4] ) 
              data['coords'][j,2]=float( ff[5] ) 
        else: i=i+1

    i=0
    while(i<len(fl)):
        if "               Potential          X             Y             Z"  in fl[i]:
           i= i + 2 + self.natoms
           data['gridesp'] = []
           while not ("----" in fl[i]):
              ff = fl[i].split()
              data['gridesp'].append( float( ff[1] ) )
              i=i+1
           data['gridesp'] = np.asarray( data['gridesp'], dtype=np.float )

        else: i=i+1

    f.close()
    return data
   except:
    return None

       
  def _read_psi4( self, dirname ):
   try:
    data={}

    f = open( os.path.join( dirname, "psi4.out" ), "r" )
    fl= f.readlines()
    f.close()
    dipole=None
    quadrupole=None
    mulliken=None

    for i in range( len(fl) ):
      if "Mulliken Charges: (a.u.)" in fl[i]:
        mulliken=[]
        for j in range(self.natoms):
           mulliken.append( fl[i+2+j].split()[5] )
        data['mulliken'] = mulliken

      if "Dipole Moment: (Debye)" in fl[i]:
        dipole = fl[i+1]

      if " Traceless Quadrupole Moment: (Debye Ang)" in fl[i]:
        quadrupole = [ fl[i+1], fl[i+2] ]

    if( dipole ):
      s = dipole.split()
      data['dipole'] = [ float(s[1]), float(s[3]), float(s[5]), float(s[7]) ]
    if( quadrupole ):
      s1 = quadrupole[0].split()
      s2 = quadrupole[1].split()
      data['quadrupole'] = [ float(s1[1]), float(s1[3]), float(s1[5]), float(s2[1]), float(s2[3]), float(s2[5]) ]


    f = open( os.path.join( dirname, "psi4out.xyz" ), "r" )
    fl=f.readlines()

    natoms=int(fl[0].split()[0]);
    energy=float(fl[0].split()[1]);
    atoms=np.zeros(( natoms, 3));
  
    for i in range(natoms):
      ff = ( fl[i+2].split() )
      atoms[i,0] = float(ff[1])
      atoms[i,1] = float(ff[2])
      atoms[i,2] = float(ff[3])
    data['coords']= atoms
    data['energy']= energy

    gp = os.path.join( dirname, "grid_esp.dat" )
    if os.path.exists( gp ):
        dd = np.loadtxt( gp )
        data['gridesp'] = np.zeros((len(dd))) 
        i=0
        for d in dd:
           data['gridesp'][i]=float(d)
           i=i+1


    f.close()
    return data
   except:
    return None

  def _write_psi4( self, dirname, frame ):
    coords = self.molecule.coords[:,:,frame]

    nrealatoms = coords.shape[0] # TODO: compensate for dummy atoms

    f = open( os.path.join( dirname, "psi4.in" ), "w" )
    basis = "unknown"
    # If the charge is < 0, need to use a diffuse basis set
    if   self.basis == BasisSet._6_31G_star: 
      if self.charge < 0:
        basis = "6-31+G*"
      else:
        basis = "6-31G*"
    elif self.basis == BasisSet._cc_pVTZ: 
      if self.charge < 0:
        basis = "aug-cc-pvtz"
      else:  
        basis = "cc-pvtz"

    else: raise ValueError( "Unknown basis set %s" % (self.basis) )
    if self.theory==Theory.HF:
      print( "set {\n\treference rhf\n\tbasis %s\n}\n" % ( basis ), file=f )

    print("\nset_num_threads( %d )" % (self.ncpus ), file=f )
    print("memory %d gb" % (self.mem ), file=f )

    print("\nmolecule MOL {\n%d %d\n" % (self.charge, self.multiplicity), file=f )
    for i in range(coords.shape[0] ):
      print("%s\t %f\t %f\t %f" % ( self.molecule.element[i], coords[i,0], coords[i,1], coords[i,2] ), file=f )
    print( "\n\tsymmetry c1\n}", file=f )

    if self.frozen:
      print( "set optking {\n\tfrozen_dihedral = (\"", file=f )
      for i in range(len(self.frozen)):
         bb=""
         if(i < len(self.frozen)-1): bb=","
         print("\t\t%d %d %d %d%s" % ( self.frozen[i][0], self.frozen[i][1], self.frozen[i][2],  self.frozen[i][3], bb), file=f )
      print("\t\")\n}\n", file=f )

    if self.optimize: 
      print("ee,wfn = optimize('scf', return_wfn=True)", file = f )
    else:
      print("ee,wfn = energy('scf', return_wfn=True)", file = f )

    print("oeprop( wfn, 'DIPOLE', 'QUADRUPOLE', 'MULLIKEN_CHARGES')", file=f )
    if self.points is not None:
      print("oeprop( wfn, 'GRID_ESP')", file = f )
      self._write_points( os.path.join( dirname, "grid.dat" ), self.points[frame] )

    print( "f=open( 'psi4out.xyz', 'w' )", file=f )
    print( "f.write( \"%d  \" )" % ( nrealatoms ), file=f )
    print( "f.write( str(ee) + \"\\n\" )", file=f )
    print( "f.write( MOL.save_string_xyz() )", file=f )
    print( "f.close()", file=f )
    f.close()


  def _write_xyz( self, dirname, frame ):
    coords = self.molecule.coords[:,:,frame]
    f = open( os.path.join( dirname, "input.xyz"  ), "w" )
    print( "%d\n" % ( coords.shape[0] ), file=f )
    for i in range(coords.shape[0] ):
      print("%s\t %f\t %f\t %f" % ( self.molecule.element[i], coords[i,0], coords[i,1], coords[i,2] ), file=f )
    f.close()

    pass
  

  def _write_gaussian( self, dirname, frame ):
    coords = self.molecule.coords[:,:,frame]
    f = open( os.path.join( dirname, "input.gjf"  ), "w" )

    print("%%nprocshared=%d" % (self.ncpus), file=f )
    print("%%mem=%dGB" % (self.mem), file=f )
    theory="unknown"
    basis="unknown"
    if self.theory  == Theory.HF       : theory = "HF"

    if self.basis   == BasisSet._6_31G_star: 
     if self.charge < 0:
      basis  = "6-31+G*"  
     else:
      basis  = "6-31G*"  
    elif self.basis == BasisSet._cc_pVTZ: 
      if self.charge < 0:
       basis = "AUG-cc-pVTZ"
      else:
       basis = "cc-pVTZ"
    else: raise ValueError( "Unknown basis set %s" % (self.basis) )

    opt=""
    if self.optimize : opt="opt=ModRedundant"

    print( "#%s/%s nosymm scf=tight %s" % (theory, basis, opt ), file=f )
    if self.points is not None:
      print( "prop=(read,field)", file=f )
    print("\nMol\n\n%d %d" %(self.charge, self.multiplicity), file=f )
    for i in range(coords.shape[0] ):
      print("%s\t %f\t %f\t %f" % ( self.molecule.element[i], coords[i,0], coords[i,1], coords[i,2] ), file=f )
    print("", file=f )
    if self.frozen:
      for i in range(len(self.frozen)):
         print("%s %s %s %s F" % ( self.frozen[i][0], self.frozen[i][1], self.frozen[i][2],  self.frozen[i][3]), file=f )



    if self.points is not None:
      print("@grid.dat /N", file=f )
      self._write_points( os.path.join( dirname, "grid.dat" ), self.points[frame] )

    f.close()

    pass
  
  def complete( self ):
    pass


#if __name__ == "__main__":
#  from htmd import *
#  m = Molecule('ethanol.mol2')
#  p = QMCalculation( m, directory="test/" , esp=True, code=Code.PSI4 )
#  p = p.results()[0]
#  print(p.energy)
#  print(p.coords)
#  print(p.esp_scalar)
