# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import numpy as np
import scipy as sp
import math

class FFEvaluate:
  def __init__( self, ffmol ):
    self.prm = ffmol._prm
    self.rtf = ffmol._rtf
    self.mol = ffmol
    self.natoms = ffmol.natoms

    # set up two sparse matrices
    # that indicate with atom pairs have 1-2,1-3 exclusions
    # and 1-4 scaling

    self.s14 = sp.sparse.lil_matrix( (ffmol.natoms, ffmol.natoms) )
    self.excl= sp.sparse.lil_matrix( (ffmol.natoms, ffmol.natoms) )

    for d in ffmol.dihedrals:
       self.s14[ d[0], d[3] ] = 1
       self.s14[ d[3], d[0] ] = 1
       self.excl[ d[0], d[3] ] = 1
       self.excl[ d[3], d[0] ] = 1

    for d in ffmol.bonds:
       self.excl[ d[0], d[1] ] = 1
       self.excl[ d[1], d[0] ] = 1

    for d in ffmol.angles:
       self.excl[ d[0], d[2] ] = 1
       self.excl[ d[2], d[0] ] = 1
    
    for d in range(ffmol.natoms):
       self.excl[d,d] = 1

  def _evaluate_elec( self, coords ):
    ee = 0.
    for i in range( 0, self.natoms ):
      qi = self.mol.charge[i]
      for j in range( i+1, self.natoms ):
         if not self.excl[ i, j ]:
             qj = self.mol.charge[j]

             dr = np.linalg.norm( coords[j,:] - coords[i,:] )
             e = qi * qj * 332.0636 / dr
             ee = ee + e

    return ee;

  def _evaluate_vdw( self, coords ):
    ee = 0.
    for i in range( 0, self.natoms ):
      for j in range( i+1, self.natoms ):
         if not self.excl[i,j]:
           (A,B) = self.prm.vdwParam( self.rtf.type_by_index[i], self.rtf.type_by_index[i], self.s14[i,j] )
           dr = np.linalg.norm( coords[j,:] - coords[i,:] )
           e = ( A / math.pow( dr, 12 ) ) - ( B / math.pow( dr, 6 ) )
           ee = ee + e
 
    return ee;
 
     
  def _evaluate_bonds( self, coords ):
    ee = 0.
    for b in self.mol.bonds:
      b1 = b[0]
      b2 = b[1]
      r12len = np.linalg.norm( coords[b1,:] - coords[b2,: ] )
      p      = self.prm.bondParam( self.rtf.type_by_index[b1], self.rtf.type_by_index[b2] )
      dist = r12len - p.r0
      coef = -2.0 * p.k0 * dist / r12len
      e = p.k0 * dist * dist
      ee = ee + e
    return ee

  def _evaluate_angles( self, coords ):
    ee = 0.
    for b in self.mol.angles:
      b1 = b[0]
      b2 = b[1]
      b3 = b[2]
   #   dr = np.linalg.norm( coords[b1,:] - coords[b2,: ] )
      p = self.prm.angleParam( self.rtf.type_by_index[b1], self.rtf.type_by_index[b2], self.rtf.type_by_index[b3] )
      theta0 = p.theta0 * math.pi / 180.
      r23 = coords[b3,:] - coords[b2,:]
      r21 = coords[b1,:] - coords[b2,:]

      r23 = r23.flatten()
      r21 = r21.flatten()

      inv_r23len = 1. / np.linalg.norm(r23)
      inv_r21len = 1. / np.linalg.norm(r21)

      cos_theta = r21.dot(r23) * inv_r21len * inv_r23len
      
      cos_theta = min( cos_theta, 1.0 )
      cos_theta = max( cos_theta,-1.0 )
 
      delta_theta = math.acos( cos_theta) - theta0
   
      sin_theta = math.sqrt( 1. - cos_theta * cos_theta )

      e = p.k0 * delta_theta * delta_theta

      ee = ee + e

      if (p.kUB!=None) and (p.kUB != 0.):
         r13 = r23 - r21
         r12len = np.linalg.norm( r13 )
         dist = r12len / p.rUB
         e = p.kUB * dist / r12len
         ee = ee + e

    return ee

  def _evaluate_dihedrals( self, coords ):
    ee = 0.
    for b in self.mol.dihedrals:
      b1 = b[0]
      b2 = b[1]
      b3 = b[2]
      b4 = b[3]
   #   dr = np.linalg.norm( coords[b1,:] - coords[b2,: ] )
      (phi) = self.prm.dihedralParam( self.rtf.type_by_index[b1], self.rtf.type_by_index[b2], self.rtf.type_by_index[b3], self.rtf.type_by_index[b4] )
      ee = ee + self.evaluateTorsion( coords[b,:], phi )


    return ee

  def _evaluate_impropers( self, coords ):
    ee = 0.
    for b in self.rtf.impropers:

      b1 = b[0]
      b2 = b[1]
      b3 = b[2]
      b4 = b[3]
   #   dr = np.linalg.norm( coords[b1,:] - coords[b2,: ] )
      (phi) = self.prm.improperParam( self.rtf.type_by_index[b1], self.rtf.type_by_index[b2], self.rtf.type_by_index[b3], self.rtf.type_by_index[b4] )

      ee = ee + self.evaluateTorsion( coords[b,:], phi )

    return ee

  @staticmethod
  def evaluateTorsion(  coords, phis ):
    pos1 = coords[0,:]
    pos2 = coords[1,:]
    pos3 = coords[2,:]
    pos4 = coords[3,:]
   
    r12 = pos1 - pos2
    r23 = pos2 - pos3
    r34 = pos3 - pos4

    r12 = r12.flatten()
    r23 = r23.flatten()
    r34 = r34.flatten()

    A = np.cross( r12, r23 )
    B = np.cross( r23, r34 )
    C = np.cross( r23, A )

    rA = 1. / np.linalg.norm(A)
    rB = 1. / np.linalg.norm(B)
    rC = 1. / np.linalg.norm(C)

    B = B.dot(rB)
   
    cos_phi = A.dot(B) * rA
    sin_phi = C.dot(B) * rC

    phi = - math.atan2( sin_phi, cos_phi )
    
    K = 0.
    K1= 0.
    for pp in phis:
      k    = pp.k0
      n    = pp.n
      
      phi0 = pp.phi0 * math.pi / 180.
      if( n ):
        K = K + k * ( 1. + math.cos( n * phi - phi0 ) )
      else: # it's an improper
        diff = phi - phi0
        K = K + k * diff * diff
      #  K1= K1 -n * k * math.sin( n * phi - phi0 )

    return K


  def evaluate( self, coords ):
    ee = dict()
    ee['elec'] = self._evaluate_elec( coords )
    ee['vdw'] = self._evaluate_vdw( coords )
    ee['bond'] = self._evaluate_bonds( coords )
    ee['angle']= self._evaluate_angles( coords )
    ee['dihedral'] = self._evaluate_dihedrals( coords )
    ee['improper'] = self._evaluate_impropers( coords )
    s=0.
    for e in ee:
      s = s + ee[e]
    ee['total'] = s

    return ee


if __name__ == "__main__":
  from htmd.newparameterization.ffmolecule import FFMolecule
  ff    = FFMolecule( "benzamidine.mol2" )
  ff._prm.write( "benzamidine.prm" )
  ff._rtf.write( "benzamidine.rtf" )
  ffeval=FFEvaluate( ff ) 
  ee = ffeval.evaluate( ff.coords[:,:,0] )
  print(ee)
