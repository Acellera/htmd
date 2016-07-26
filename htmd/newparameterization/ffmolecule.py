# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from htmd.molecule.molecule import Molecule
from htmd.molecule.vmdparser import guessAnglesAndDihedrals
from htmd.newparameterization.detectsoftdihedrals import detectSoftDihedrals
from htmd.newparameterization.detectequivalents import detectEquivalents
from htmd.newparameterization.fftype import FFTypeMethod,FFType
from htmd.newparameterization.ff import RTF,PRM
from htmd.newparameterization.ffevaluate import FFEvaluate
from htmd.qm.qmcalculation  import *
from htmd.newparameterization.phi import setPhi, getPhi
from htmd.progress.progress import ProgressBar
import scipy.spatial.distance
import re
import shutil
import math
import os
import sys
import tempfile
import scipy.optimize as optimize
import numpy as np
from copy import deepcopy


class QMFittingSet:
 # phi_coords  = []
  def __init__( self ):
    self.phi= []
    self.qm = []
    self.mm_original = []
    self.mm_zeroed   = []
    self.mm_delta    = []
    self.mm_fitted   = []
    self.coords      = []


class FFMolecule(Molecule):
  def __init__(self, filename=None, name=None, rtf=None, prm=None, netcharge=None, method=FFTypeMethod.CGenFF_2b6 ):
    # filename -- a mol2 format input geometry  
    # rtf, prm -- rtf, prm files
    # method  -- if rtf, prm == None, guess atom types according to this method ( of enum FFTypeMethod )

    if( not ( filename.endswith(".mol2")) ):
      raise ValueError( "Input file must be mol2 format" )

    super().__init__(filename=filename, name=name )
    (a,b) = guessAnglesAndDihedrals( self.bonds )
    self.natoms    = self.serial.shape[0]
    self.angles    = a
    self.dihedrals = b
    ee = detectEquivalents( self )
    self._soft_dihedrals = detectSoftDihedrals( self, ee ) 
    self._equivalent_atom_groups= ee[0] # list of groups of equivalent atoms
    self._equivalent_atoms      = ee[1] # list of equivalent atoms, indexed by atom
    self._equivalent_group_by_atom = ee[2] # mapping from atom index to equivalent atom group
    if netcharge==None:
      self.netcharge = int(round(np.sum( self.charge )))
    else:
      self.netcharge = int(round(charge))

    # Canonicalise the atom naming.
    self._rename_mol()

    if( rtf and prm ):
      # If the user has specified explicit RTF and PRM files go ahead and load those
      self._rtf = RTF( rtf )
      self._prm = PRM( prm )
    else:
      # Otherwise make atom types using the specified method
      # (Right now only MATCH)
      fftype = FFType( self, method=method )
      self._rtf = fftype._rtf
      self._prm = fftype._prm
    if( not self._rtf  or not self._prm ):
      raise ValueError("RTF and PRM not defined")

    self.report()

  def report(self):
    print("Net Charge: %d" % ( self.netcharge ) )
    print("Equivalent atom groups:" )
    for i in self._equivalent_atom_groups:
      for j in i:
        print( " %s" % ( self.name[j] ), end="" )
      print("")

    print("Soft dihedrals:" )
    for i in self._soft_dihedrals:
      for j in i.atoms:
        print( " %s" % ( self.name[j] ), end="" )
      print("")

  def _rename_mol(self):
        # This fixes up the atom naming and reside name to be consistent
        # NB this scheme matches what MATCH does. Don't change it
        # Or the naming will be inconsistent with the RTF
        import re

        hh = dict()

        for i in range(len(self.name)):
            t = re.sub('[1234567890]*', "", self.name[i]).upper()
            idx = 0

            if not t in hh:
                hh[t] = idx

            idx = hh[t] + 1;
            hh[t] = idx

            t = t + str(idx)
            self.name[i] = t
            self.resname[i] = "MOL"

  def minimize(self):
    # Kick off a QM calculation -- unconstrained geometry optimization
    qm = QMCalculation( self, charge=self.netcharge, optimize=True, directory="minimize/" )
    results = qm.results()
    if results[0].errored:
      raise RuntimeError("QM Optimization failed")
    # Replace coordinates with the minimized set
    self.coords = np.atleast_3d( results[0].coords )

  def _fitCharges_map_back_to_charges( self, x ):
    charges=np.zeros((self.natoms))

    qsum = 0.
    for i in range(len(x)):
       charges[ self._equivalent_atom_groups[i] ] = x[i]
       qsum = qsum + x[i] * len(self._equivalent_atom_groups[i] )

  #  diff = self.netcharge - qsum;
  #  diff = diff / len(self._equivalent_atom_groups[ len(x) ])
  #  print( self._equivalent_atom_groups[ len(x) ] )
  #  print( diff )
  #  charges[ self._equivalent_atom_groups[ len(x) ] ] = diff

    return charges 


  def _fitCharges_con( self, x ):
    charges = self._fitCharges_map_back_to_charges( x )
    s=np.sum(charges) -  self.netcharge;
    return s



  def _fitCharges_objective( self, x ):
    # Map the fit variables back to per-atom charges
    chisq = 0.
    charges = self._fitCharges_map_back_to_charges( x )



#    if( range_penalty == 1 ): chisq = 1000.

    for i in range( self._fitCharges_grid.shape[0] ):
      ee =  np.sum(charges *  self._fitCharges_distances[i,:])
      delta_ee = self._fitCharges_esp[ i ] - ee
      chisq = chisq + (delta_ee * delta_ee )

    return chisq

  def _fitDihedral_objective( self, x ):
    inv = math.pi / 180.

    # evalaute the torsion with the input params
    # for each of the phi's poses
    chisq = 0.
    for t in range(self._fitDihedral_results.N):
      e  = .0 # FFEvaluate.evaluateTorsion( self._fitDihedral_results["phi_coords"][t], phi )
      for s in range(len(self._fitDihedral_results.phis[t])):
        for j in range(6):
           e = e +  x[j] * (1. + cos( (j+1) * ( self._fitDihedral_results.phis[t][s]*inv ) - x[6+j] * inv ) )
      
      e = e + x[12] 
      diff  = self._fitDihedral_results.mm_delta[t] - e 
      chisq = chisq + diff * diff

    return chisq

  def _removeCOM(self):
    # Relocate centre of mass to the origin
    for f in range(self.coords.shape[2]):
      com=np.zeros(3)
      mass=0.
      for i in range( self.coords.shape[0]):  
        m =  self._rtf.mass_by_type[ self._rtf.type_by_index[i] ]
        mass=mass+m
        com = com + self.coords[i,:,f] * m
      com = com / mass
      self.coords[:,:,f] = self.coords[:,:,f] - com

  def _try_load_pointfile(self):
    # Load a point file if one exists from a previous job
    pointfile = os.path.join( "esp", "00000", "grid.dat" )
    if( os.path.exists(pointfile)):
      f = open(pointfile,"r" );
      fl = f.readlines()
      f.close()
      ret = np.zeros( ( len(fl), 3 ))
      for i in range(len(fl)):
        s = fl[i].split()
        ret[i,0] = float(s[0]) 
        ret[i,1] = float(s[1]) 
        ret[i,2] = float(s[2]) 

      print("Reusing previously-generated point cloud")
      return ret
    return True
 
  def fitCharges(self):
    # Remove the COM from the coords, or PSI4 does it and then the grid is incorrectly centred
    self._removeCOM()
    # Kick off a QM calculation -- unconstrained single point with grid
    points = self._try_load_pointfile() 
     
    qm = QMCalculation( self, charge=self.netcharge, optimize=False, esp=points, directory="esp/"  )
    results = qm.results()
    if results[0].errored:
      raise RuntimeError("QM Calculation failed")
    esp_grid = results[0].esp_points
    esp      = results[0].esp_scalar
    self.coords = results[0].coords
 
#    print(results[0].dipole )
#    print(results[0].quadrupole )
#    print(results[0].mulliken )
 
    self._fitCharges_grid = esp_grid
    self._fitCharges_esp  = esp

    # set up the restraints to fit
   
    N = len(self._equivalent_atom_groups) #- 1
    lb= np.ones((N)) * -1.25 
    ub= np.ones((N)) * +1.25 
    # If the restraint relates to an H, set the lower bound to 0
    for i in range(N):
      if "H" == self.element[ self._equivalent_atom_groups[i][0] ]:
        lb[i]=0.001

    bounds=[]
    for a in range(len(lb)):
      bounds.append( ( lb[a], ub[a] ) )   

    # Start off by equally distributing the mol's charge
    start = np.zeros(N)

    # Precompute the 1/r distances
    self._fitCharges_distances = np.zeros( ( self._fitCharges_grid.shape[0], self.coords.shape[0] ) )

    for i in range( self._fitCharges_grid.shape[0] ):
      p1 = self._fitCharges_grid[i,:]
      for j in range( self.coords.shape[0] ):
        p2 = self.coords[j,:,0]
        r  = np.linalg.norm( p1 - p2 )
        self._fitCharges_distances[i,j] = 1./r



#    initial_chisq = self._fitCharges_objective( start )

    xopt = optimize.minimize( self._fitCharges_objective, start, method="SLSQP", bounds = bounds, options={"disp":False}, constraints= { 'type':'eq', 'fun':self._fitCharges_con } )
#    xopt = optimize.minimize( self._fitCharges_objective, start, method="L-BFGS-B", bounds = bounds, options={"disp":False} )

    charges = self._fitCharges_map_back_to_charges( xopt.x )
  
    # Calculate the dipole from the fitted charges

    dpx=dpy=dpz=0.
    nc=0.
    for i in range(len(charges)):
      dpx = dpx + charges[i] * self.coords[i,0,0]
      dpy = dpy + charges[i] * self.coords[i,1,0]
      dpz = dpz + charges[i] * self.coords[i,2,0]
      nc=nc + charges[i]
    fac = (2.541766 / 0.529177249)
    dpx = dpx * fac
    dpy = dpy * fac
    dpr = dpz * fac
    dp = math.sqrt( dpx*dpx + dpy*dpy + dpz*dpz )


    fit_chisq = self._fitCharges_objective( xopt.x )

    self.charges = charges
    self._rtf.updateCharges( charges )



    return (fit_chisq, results[0].dipole, [dpx, dpy, dpz, dp] )

  def getSoftDihedrals(self):
    dd=[]
    for d in self._soft_dihedrals:
      dd.append( d.atoms.copy() )
    return dd

  def scanSoftDihedral(self, phi, directory="dihedral", step=10):
    found=False
    phi_to_fit = None
    frozens=[]
    dih_index=0
    i=0
    for d in self._soft_dihedrals:
      if (d.atoms == phi).all():  
         phi_to_fit = d
         dih_index=i
      else:
         frozens.append(d.atoms)
      i=i+1
    if not phi_to_fit: raise ValueError( "specified phi is not a recognised soft dihedral" )

    atoms = phi_to_fit.atoms 
    left  = phi_to_fit.left 
    right = phi_to_fit.right 
    equivs= phi_to_fit.equivalents
 
#    step  = 10 # degrees
    nstep = (int)(360/step)
    cset  = np.zeros( ( self.natoms, 3, nstep ) )

    i=0
    for phi in range( -180, 180, step ):
      cset[:,:,i] = setPhi( self.coords[:,:,0], atoms, left, right, phi )
      i=i+1

    mol        = self.copy()
    mol.coords = cset
    try:
      os.mkdir( directory )
    except:
      pass
    dih_name = "%s-%s-%s-%s" % ( self.name[atoms[0]], self.name[atoms[1]], self.name[atoms[2]], self.name[atoms[3]] )
    qmset   = QMCalculation( mol, charge=self.netcharge, directory="%s/%s" % (directory, dih_name), frozen=frozens )
    r = qmset.results()
    x=0
    ret=[]
    for phi in range( -180, 180, step ):
      r[x].phi = phi
      if r[x].errored == False:
        ret.append(r[x])
      x=x+1
    return ret

  def fitSoftDihedral( self, phi ):
    found=False
    phi_to_fit = None
    frozens=[]
    dih_index=0
    i=0
    for d in self._soft_dihedrals:
      if (d.atoms == phi).all():  
         phi_to_fit = d
         dih_index=i
      else:
         frozens.append(d.atoms)
      i=i+1
    if not phi_to_fit: raise ValueError( "specified phi is not a recognised soft dihedral" )
    self._makeDihedralUnique( phi_to_fit )

    atoms = phi_to_fit.atoms 
    left  = phi_to_fit.left 
    right = phi_to_fit.right 
    equivs= phi_to_fit.equivalents
 
    step  = 10 # degrees
    nstep = (int)(360/step)
    cset  = np.zeros( ( self.natoms, 3, nstep ) )

    i=0
    for phi in range( -180, 180, step ):
      cset[:,:,i] = setPhi( self.coords[:,:,0], atoms, left, right, phi )
      i=i+1

    mol        = self.copy()
    mol.coords = cset
    try:
      os.mkdir("dihedral")
    except:
      pass
    dih_name = "%s-%s-%s-%s" % ( self.name[atoms[0]], self.name[atoms[1]], self.name[atoms[2]], self.name[atoms[3]] )
    qmset   = QMCalculation( mol, charge=self.netcharge, directory="dihedral/%s" % (dih_name), frozen=frozens )

    ret = self._makeDihedralFittingSetFromQMResults( atoms, qmset.results() )

    # Get the initial parameters of the dihedral we are going to fit

    param = self._prm.dihedralParam( self._rtf.type_by_index[ atoms[0] ],
                             self._rtf.type_by_index[ atoms[1] ],
                             self._rtf.type_by_index[ atoms[2] ],
                             self._rtf.type_by_index[ atoms[3] ] )

    # Save these parameters as the best fit (fit to beat)
    best_param=np.zeros((13))
    for t in range(6):
       best_param[t]   = param[t].k0
       best_param[t+6] = param[t].phi0
    best_param[12] = 0.

#    print(param)

    # Evalaute the mm potential with this dihedral zeroed out
    # The objective function will try to fit to the delta between
    # The QM potential and the this modified mm potential

    for t in param:
       t.k0 = t.phi0 = 0.
    self._prm.updateDihedral( param )
 
    ffe = FFEvaluate( self )
  #  print(ffe.evaluate( ret.coords[0] ) )
  #  input
    # Now evaluate the ff without the dihedral being fitted
    for t in range(ret.N):
       mm_zeroed    = ( ffe.evaluate( ret.coords[t] )["total"])
       ret.mm_delta.append( ret.qm[t] - mm_zeroed )
       ret.mm_zeroed.append( mm_zeroed )  

    mmin1 = min( ret.mm_zeroed )
    mmin2 = min( ret.mm_delta )
    for t in range(ret.N): 
      ret.mm_zeroed[t] = ret.mm_zeroed[t] - mmin1
      ret.mm_delta[t]  = ret.mm_delta[t]  - mmin2

    self._fitDihedral_results = ret
    self._fitDihedral_phi     = param

    # Now measure all of the soft dihedrals phis that are mapped to this dihedral 
    ret.phis= []
    for i in range(ret.N):
      ret.phis.append( [ ret.phi[i] ] )
      for e in equivs:
         ret.phis[i].append( getPhi( ret.coords[i], e ) )
#    print ("EQUIVALENT DIHEDRALS FOR THIS DIHEDRAL" )
#    print(equivs)
#    print ("PHI VALUES TO FIT")
#    print (ret.phis)
    #  Set up the NOLOPT fit
    #  There are 13 parameters, k,phi for n=1,2,3,4,5,6 and a shift
    N = 13
    # initial guess,
    st= np.zeros(13)
    # bounds

    best_chisq = self._fitDihedral_objective( best_param )
#    print("CHISQ of initial = %f" % ( best_chisq ) )

    # Now zero out the terms of the dihedral we are going to fit
    bar=ProgressBar(64, description="Fitting")
    for i in range(64):
      
      ( bounds, start ) = self._fitDihedral_make_bounds( i )

      xopt = optimize.minimize( self._fitDihedral_objective, start, method="L-BFGS-B", bounds = bounds  , options={'disp': False } )

      chisq = self._fitDihedral_objective( xopt.x )
#      print( "CHISQ of fit = %f " % (chisq) )
      if( chisq < best_chisq ):
         best_chisq = chisq
         best_param = xopt.x
      bar.progress()
    bar.stop()
#    print("Best ChiSQ = %f" %(best_chisq) )

    # Update the target dihedral with the optimized parameters
    # print(param)
    # print(best_param )
    for i in range(6):
      param[i].k0   = best_param[0+i]
      param[i].phi0 = best_param[6+i]
    
    self._prm.updateDihedral( param )
    # print(param)
    param = self._prm.dihedralParam( self._rtf.type_by_index[ atoms[0] ],
                             self._rtf.type_by_index[ atoms[1] ],
                             self._rtf.type_by_index[ atoms[2] ],
                             self._rtf.type_by_index[ atoms[3] ] )
    # print(param)

    # Finally evaluate the fitted potential
    ffe = FFEvaluate( self )
    for t in range(ret.N):
       ret.mm_fitted.append( ffe.evaluate( ret.coords[t] )["total"] )
    mmin = min(ret.mm_fitted )
    chisq=0.

#    print( "QM energies" )
#    print( ret.qm )

    for t in range(ret.N):
       ret.mm_fitted[t] = ret.mm_fitted[t] - mmin
       delta = ret.mm_fitted[t] - ret.qm[t]
       chisq = chisq + (delta * delta )
    ret.chisq = chisq

# TODO Score it

    return ret

  def _fitDihedral_make_bounds( self, i ): 
    lb= np.zeros(13)
    ub= np.zeros(13)
    start = np.zeros(13)

    bounds = []

    for j in range(6):
       start[j] = 0.;
       bounds.append( ( -20., 20. ) )

    for j in range( 6 ):
      if( i & (2**j) ):
         bounds.append( ( 180., 180. ) )
         start[6+j] = 180.
      else:
         bounds.append( ( 0., 0. ) )
         start[6+j] = 0.

    bounds.append( ( -10., 10. ) )
    return (bounds, start)

   
  def _makeDihedralFittingSetFromQMResults( self, atoms, results ): 
    # Extract the valid QM poses and energies from the QM result set
    # Evaluate the MM on those poses
    ffe = FFEvaluate( self )

    ret=QMFittingSet()
    ret.name        = "%s-%s-%s-%s" % ( self._rtf.names[atoms[0]], self._rtf.names[atoms[1]], self._rtf.names[atoms[2]], self._rtf.names[atoms[3]]  )

    completed = 0 

    qmin = 1.e100
    for q in results:
      if( q.completed and not q.errored ):
        if( q.energy < qmin ) : qmin = q.energy
    
    completed = 0
    for q in results:
      if( q.completed and not q.errored ):
        if ( q.energy - qmin ) < 20. : # Only fit against QM points < 20kcal above the minimum 
          completed = completed + 1
          ret.phi.append( getPhi( q.coords, atoms ) )
          
          ret.qm.append( q.energy - qmin  )
          ret.mm_original.append( ffe.evaluate( q.coords)['total'] )
          ret.coords.append( q.coords )
          if(ret.phi[0] > 175. ) : ret.phi[0] = ret.phi[0] - 360. 
    mmin = min( ret.mm_original )
    # roughly align the qm with the mm
    for q in range(completed):
       ret.mm_original[q] = ret.mm_original[q] - mmin
    ret.N = completed

    if completed < 5:
      raise RuntimeError( "Fewer than 5 valid QM points. Not enough to fit!" )

    return ret
   

  def _makeDihedralUnique( self,  phi_to_fit ):
#    (number_of_uses, uses) = self._countUsesOfDihedral( phi_to_fit.atoms )
#    if( number_of_uses > 1 ):
       # Create a new type for (arbitrarily) a middle atom of the dihedral
       # So that the dihedral we are going to modify is unique
       # TODO -- check symmetry
   #    print( "Dihedral term is not unique. Copying type.." ) # Used %d times, by:" % ( number_of_uses ) )
       # print( uses )

    # Duplicate the dihedrals types so this modified term is unique
    print( "Duplicating types.." )
    for i in range(4):
         if not ( "_" in self._rtf.type_by_index[phi_to_fit.atoms[i]]):
           self._duplicateTypeOfAtom( phi_to_fit.atoms[i] )

    (number_of_uses, uses) = self._countUsesOfDihedral( phi_to_fit.atoms )
    if( number_of_uses > 1 ):
         print( phi_to_fit.atoms )
         print( number_of_uses )
         print( uses )
         raise ValueError( "Dihedral term still not unique after duplication" )





  def _countUsesOfDihedral( self, aidx ):
  
    # Return the number of uses of the dihedral
    # specified by the types of the 4 atom indices in the aidx list
    # 
 
#    print( "countUsesOfDihedral in " )
    t1 = self._rtf.type_by_index[aidx[0]]
    t2 = self._rtf.type_by_index[aidx[1]]
    t3 = self._rtf.type_by_index[aidx[2]]
    t4 = self._rtf.type_by_index[aidx[3]]

    count=0
    uses=[]
    for d in self.dihedrals:
      s1 = self._rtf.type_by_index[d[0]] 
      s2 = self._rtf.type_by_index[d[1]] 
      s3 = self._rtf.type_by_index[d[2]] 
      s4 = self._rtf.type_by_index[d[3]] 
      if  ( s1==t1 and s2==t2 and s3==t3 and s4==t4 ) : 
        count = count + 1
        uses.append( d )
      elif( s1==t4 and s2==t3 and s3==t2 and s4==t1 ) : 
        count = count + 1
        uses.append( d )

#    return(count, uses )
#    print(uses)

    # Now for each of the uses, remove any which are equivalent
    c=1
    unique_uses=[ aidx ]
    g1 = self._equivalent_group_by_atom[ aidx[0] ]
    g2 = self._equivalent_group_by_atom[ aidx[1] ]
    g3 = self._equivalent_group_by_atom[ aidx[2] ]
    g4 = self._equivalent_group_by_atom[ aidx[3] ]
    for u in uses:
      h1 = self._equivalent_group_by_atom[ u[0] ]
      h2 = self._equivalent_group_by_atom[ u[1] ]
      h3 = self._equivalent_group_by_atom[ u[2] ]
      h4 = self._equivalent_group_by_atom[ u[3] ]
      equiv=False
      if g1==h1 and g2==h2 and g3==h3 and g4==h4: equiv=True
      if g1==h4 and g2==h3 and g3==h2 and g4==h1: equiv=True
      if( equiv==False):
         c= c+1
         unique_uses.append(u)
      else:
         print(" Dih %d-%d-%d-%d and %d-%d-%d-%d are equivalent " % ( aidx[0], aidx[1], aidx[2], aidx[3], u[0], u[1], u[2], u[3] ) )
         pass
 
  #  return(count, uses )
#    print( c )
#    print( unique_uses )

    return (c, unique_uses)

     
  def _duplicateTypeOfAtom( self, aidx ):
    # This duplicates the type of the specified atom
    # First get the type
    atype = self._rtf.type_by_index[ aidx ]

    # perhaps the type is already a duplicate? if so
    # remove the duplicated suffix
    atype = re.sub( "_[0123456789]+$", "", atype )
    i=0
    # make the new type name
    while ("%s_%d" % ( atype, i )) in self._rtf.types: i=i+1

    newtype = "%s_%d" % (atype, i)
    print("Creating new type %s from %s for atom %s" %( newtype, atype, self._rtf.names[ aidx ] ) )

    # duplicate the type in the fields RTF -- todo: move to a method in the RTF
    self._rtf.type_by_index[ aidx ] = newtype
    self._rtf.mass_by_type[ newtype ] = self._rtf.mass_by_type[ atype ]
    self._rtf.types.append( newtype )
    self._rtf.type_by_name[ self._rtf.names[ aidx ] ] = newtype
    self._rtf.type_by_index[ aidx ] = newtype
    self._rtf.typeindex_by_type[ newtype ] = self._rtf.typeindex_by_type[ atype ] + 1000
    self._rtf.element_by_type[ newtype ] = self._rtf.element_by_type[ atype ]

#    # Now also reset the type of  any atoms that share equivalency
    for bidx in  self._equivalent_atoms[aidx]:
      if aidx != bidx:
         if "_" in self._rtf.type_by_index[ bidx ]:
            raise RuntimeError("Equivalent atom already has a duplicated type: %d %d" % ( bidx, self._rtf.type_by_index[ bidx ] ) )
         self._rtf.type_by_index[ bidx ] = newtype
         self._rtf.type_by_name[ self._rtf.names[ bidx ] ] = newtype

    # the PRM parameters will be automatically duplicated by forcing an ff evaluation
    ffe = FFEvaluate( self )
    ffe.evaluate( self.coords )

    # 

  def plotDihedralFit( self, fit, show=True, directory="." ):
    import matplotlib as mpl
    import matplotlib.pyplot as plt 

    
    if not show:
      mpl.use('Agg')

   
    fh = plt.figure()
    ax1 = fh.gca()
    ax1.set_xlim( -180., 180. )
    ax1.set_xticks([ -180, -135, -90, -45, 0, 45, 90, 135, 180 ] )
    ax1.set_xlabel("Phi")
    ax1.set_ylabel( "kcal/mol" )
    ax1.set_title( fit.name )
    ax1.plot( fit.phi , fit.qm         , label="QM", color="r", marker="o" )
    ax1.plot( fit.phi , fit.mm_original, label="MM Original", color="b", marker="o" )  
    ax1.plot( fit.phi , fit.mm_zeroed  , label="MM With phi zeroed", color="black", marker="x" )  
    ax1.plot( fit.phi , fit.mm_delta   , label="QM-MM target", color="magenta", marker="x" )  
    ax1.plot( fit.phi , fit.mm_fitted  , label="MM Fitted", color="g", marker="o" )  
    ax1.legend(prop={'size': 8})
    if show:
      print("SHOW") # TODO FIXME doesn't work for some reason; nothing shown?!
      plt.show()
    else:
      try:
        os.mkdir( directory )
      except:
        pass
      tf = os.path.join( directory, fit.name ) + ".svg"
      plt.savefig( tf, format="svg")
      return tf

#if __name__ == "__main__":
#  os.chdir("tests/benzamidine" )
#
#  mol = FFMolecule( filename="benzamidine.mol2", method=FFTypeMethod.CGenFF_2b6 )
#
#  print("Minimizing")
#  mol.minimize()
#
#  print("Charge fitting")
#  (score, qm_dipole, mm_dipole) = mol.fitCharges()
#
#  print( "Chi^2 score : %f" %( score ) )
#  print( "QM Dipole   : %f %f %f ; %f" % ( qm_dipole[0], qm_dipole[1], qm_dipole[2], qm_dipole[3] ) )
#  print( "MM Dipole   : %f %f %f ; %f" % ( mm_dipole[0], mm_dipole[1], mm_dipole[2], mm_dipole[3] ) )
#
#  dihedrals = mol.getSoftDihedrals()
#  for d in dihedrals:
#    print("\nFitting dihedral %s-%s-%s-%s" % ( mol.name[ d[0] ], mol.name[ d[1] ], mol.name[ d[2] ], mol.name[ d[3] ]  ))
#    ret = mol.fitSoftDihedral( d )
#
#    print("Chi^2 score : %f" % ( ret.chisq ) )
#
#    fn  = mol.plotDihedralFit( ret, show=False, directory="plots" )
#    #print(fn)
#
#  try:
#    os.mkdir("parameters")
#  except:
#    pass
#  mol._rtf.write( "parameters/mol.rtf" )
#  mol._prm.write( "parameters/mol.prm" )
#  mol.write( "parameters/mol.psf" )
#  mol.write( "parameters/mol.xyz" )
#  mol.write( "parameters/mol.pdb" )
#  
#
