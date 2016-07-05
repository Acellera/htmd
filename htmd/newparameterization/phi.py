from htmd.newparameterization.detectsoftdihedrals import detectSoftDihedrals

import numpy as np
import math
  
#   Converted from old Matlab AxelRot function 
# (See acellera/parameterise/matlab/AzelRot.m)
#   MJH (2015)

def _mkaff( R, t=None ):
  if t is None:
    if( R.size == 4  or R.size == 9 ):
      nn = R.shape[0]
      t = np.zeros((nn))
    if( R.size == 2 or R.size == 3 ):
      t=R
      nn=t.shape[0]
      R=np.eye(nn)
  else:
    nn = R.shape[0]
 
  M=np.eye(nn+1)
  M[0:nn,0:nn] = R
  M[0:nn,nn]   = t
  return M


def _R3d( deg, u ):
  R = np.eye(3)

  u = u / np.linalg.norm(u)
  x = deg

#  print("R3D")
#  print(deg)
#  print(R)

  for ii in range(3):
    v = R[:,ii]

    
    t = v.dot(math.cos(x)) + np.cross(u,v).dot(math.sin(x)) 
    t = t + (( np.transpose(u).dot(v) )* ((1 - math.cos(x) ))) * ( u )
    R[:,ii] = t
  return R


def _AxelRot( coords, deg, u,  x0 ):
  # To radians
  deg = math.pi * deg/180.
  u = u / np.linalg.norm(u)

  AxisShift =  x0 - ( np.transpose(x0).dot(u) ) * u  
#  print("AxisShift")
#  print(AxisShift)
  Mshift    = _mkaff( np.eye(3), -AxisShift )
#  print("MShift")
#  print(Mshift)
  rr = _R3d(deg,u)
#  print("R3d")
#  print(rr)
  Mroto     = _mkaff( rr )
#  print("Mroto")
#  print(Mroto)
  M         = np.linalg.inv(Mshift).dot( Mroto  ).dot(  Mshift )
#  print("M")
#  print(M)

  R = M[0:3,0:3]
  t = M[0:3,3]
 
  for i in range( coords.shape[0] ):
    coords[i,:]  = R.dot(coords[i,:]) + t
  #print(coords.shape)
  return coords
 
def setPhi( coords, atoms, left, right, phi ):
  coords = np.atleast_3d(coords)
  pos1 = coords[ atoms[0], :, 0 ] 
  pos2 = coords[ atoms[1], :, 0 ] 
  pos3 = coords[ atoms[2], :, 0 ] 
  pos4 = coords[ atoms[3], :, 0 ] 

  line = pos3 - pos2
  line = line / np.linalg.norm(line)

  centre = pos2
  old_phi= getPhi( coords, atoms )
  delta_phi = - (phi - old_phi)

#  print( "ATOMS TO BE ROTATED" )
#  print(left)
#  print("CENTRE LINE %d -> %d\n" % ( atoms[2], atoms[1] ) )
#  print("OLD PHI= %f" %(old_phi ))
#  print("NEW PHI= %f" %(phi ))
#  print("DEL PHI= %f" %(delta_phi) )
  coords_old = coords[ left, :, 0 ]
  coords_new = _AxelRot( coords_old, delta_phi, line, centre )
  coords[ left, :, 0 ] = coords_new

  return coords[:,:,0]

def getPhi( coords, atoms ):
  coords = np.atleast_3d(coords)
  pos1 = coords[ atoms[0], :, 0 ] 
  pos2 = coords[ atoms[1], :, 0 ] 
  pos3 = coords[ atoms[2], :, 0 ] 
  pos4 = coords[ atoms[3], :, 0 ] 

  r12 = pos1 - pos2
  r23 = pos2 - pos3
  r34 = pos3 - pos4

  A = np.cross( r12, r23 )
  B = np.cross( r23, r34 )
  C = np.cross( r23, A )

  rA = 1. / np.linalg.norm( A )
  rB = 1. / np.linalg.norm( B )
  rC = 1. / np.linalg.norm( C )

  B = B * rB

  cos_phi = ( np.dot( A, B ) ) * rA
  sin_phi = ( np.dot( C, B ) ) * rC
  
  phi = - math.atan2( sin_phi, cos_phi )
  phi = phi * 180. / math.pi

  return phi


if __name__ == "__main__":
  from htmd.newparameterization.ffmolecule import FFMolecule
  m = FFMolecule( "benzamidine.mol2" )
  
  d = m.getSoftDihedrals() 
  print(d)
  m.fitSoftDihedral( d[0] )
  

