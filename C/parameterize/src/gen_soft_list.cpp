/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ff.h"
#include <unistd.h>

#define MAX_DIHEDRAL_LIST	(4096)
#define MAX_STEP	(800000)
#define GAP			(200)
//#define SIG_CUT		(30.0)	// for the dihedrals involved in aromatic molecules
#define SIG_CUT		(70.0)	// for the dihedrals involved in aromatic molecules and other ring 
#define MAX_ATOM_LOCAL	(256)
#define MAX_BOND		(6)

#define radianInv	(0.017453292519943295)

char szForceFiled[]="mol.prm";
char szXpsfFile[]="mol.xpsf";
//char szCrdFile[]="mol-opt.crd";

CMol Mol;
CForceField ForceField;
double T_Sim=2000.0;

int nDihedral;

double dih_List[MAX_STEP/GAP][MAX_DIHEDRAL_LIST];
double dih_Mean[MAX_DIHEDRAL_LIST], dih_Sig[MAX_DIHEDRAL_LIST];
int IsRigidDih[MAX_DIHEDRAL_LIST], Dihedral_Selected[MAX_DIHEDRAL_LIST];

int Bond_List[MAX_ATOM_LOCAL][MAX_BOND], Bond_Count[MAX_ATOM_LOCAL];


FILE *fFile_Run_Log;
void Quit_With_Error_Msg(char szMsg[]);
void Cal_Sig(void);
void Cal_Sig_Shift_PI2(void);

int To_Exclude_Methyl_Torsion;
int n_Soft_Dih=0, Dih_Bond[MAX_DIHEDRAL_LIST][4];

void Gen_Soft_Dihedral_List_MD_High_T(void);
void Gen_Soft_Dihedral_List_Excluding_Ring(void);
int Is_In_Soft_Dihedral_List(int b, int c);
int Check_Is_Methyl_Involved(int iDih);
void Setup_Bond_List(void);
int Get_Number_of_H_Bonded(int iAtom);
int Get_Number_of_F_Bonded(int iAtom);
int To_Refine_Soft_Dihedral_List(void);
void Get_Representative_14_Atoms(int iDih);

int main(int argc, char *argv[])
{

  timebomb();

	if(argc == 2)	{
		To_Exclude_Methyl_Torsion = 0;
	}
	else	{
		To_Exclude_Methyl_Torsion = 1;
	}

	memset(IsRigidDih, 0, sizeof(int)*MAX_DIHEDRAL_LIST);

	fFile_Run_Log = fopen("identify-soft.log", "w");

	ForceField.ReadForceField(szForceFiled);
	Mol.ReadPSF(szXpsfFile, 0);

	Mol.AssignForceFieldParameters(&ForceField);

	if( !access( "mol-opt.xyz", R_OK ) ) {
		Mol.ReadXYZ( "mol-opt.xyz" );
	}
	else {
		printf("Could not open .xyz file\n" );
		exit(1);

	}
	Setup_Bond_List();

//	Gen_Soft_Dihedral_List_MD_High_T();
//	printf("\n\n");

	Gen_Soft_Dihedral_List_Excluding_Ring();
	printf("\n\n");

	To_Refine_Soft_Dihedral_List();

/*	
	fData = fopen("dihedral.txt", "w");
	for(Step=0; Step<MAX_STEP/GAP; Step++)	{
		fprintf(fData, "%8d ", Step);
		for(Idx=0; Idx<nDihedral; Idx++)	{
			if(IsRigidDih[Idx] == 0)	{
				fprintf(fData, " %7.2lf", dih_List[Step][Idx]);
			}
		}
		fprintf(fData, "\n");
	}
	fclose(fData);
*/
	fclose(fFile_Run_Log);


	return 0;
}

void Get_Representative_14_Atoms(int iDih)
{
	int ia, ib, ic, id, i, Idx, nAtom_Connected, nAtom_Connected_Max=0;

	ib = Dih_Bond[iDih][1];
	ic = Dih_Bond[iDih][2];

	//start	to find the representative atom for the left side, ia
	nAtom_Connected_Max=-100;
	ia = -1;
	for(i=0; i<Bond_Count[ib]; i++)	{
		Idx = Bond_List[ib][i];
		if(Idx == ic)	{
			continue;
		}
		nAtom_Connected = Mol.Count_All_Atoms_Connected(Idx, ib);
		if(nAtom_Connected > nAtom_Connected_Max)	{
			nAtom_Connected_Max = nAtom_Connected;
			ia = Idx;
		}
		else if( (i!=0) && (nAtom_Connected == nAtom_Connected_Max) && (Mol.mass[Idx] > Mol.mass[ia]) )	{
			ia = Idx;
		}
	}
	if(ia>=0)	{
		Dih_Bond[iDih][0] = ia;
	}
	else	{
		Quit_With_Error_Msg("Fail to find the representative atom (heavy atom), ia, in dihedral.\nQuit\n");
	}
	//end	to find the representative atom for the left side, ia


	//start	to find the representative atom for the left side, ia
	nAtom_Connected_Max=-100;
	id = -1;
	for(i=0; i<Bond_Count[ic]; i++)	{
		Idx = Bond_List[ic][i];
		if(Idx == ib)	{
			continue;
		}
		nAtom_Connected = Mol.Count_All_Atoms_Connected(Idx, ic);
		if(nAtom_Connected > nAtom_Connected_Max)	{
			nAtom_Connected_Max = nAtom_Connected;
			id = Idx;
		}
		else if( (i!=0) && (nAtom_Connected == nAtom_Connected_Max) && (Mol.mass[Idx] > Mol.mass[id]) )	{
			id = Idx;
		}
	}
	if(id>=0)	{
		Dih_Bond[iDih][3] = id;
	}
	else	{
		Quit_With_Error_Msg("Fail to find the representative atom (heavy atom), id, in dihedral.\nQuit\n");
	}
	//end	to find the representative atom for the left side, ia

	return;
}


double Normalizex(double& vect_x, double& vect_y, double& vect_z)
{
	double dist, dist_Inv;

	dist = sqrt(vect_x*vect_x + vect_y*vect_y + vect_z*vect_z);
	if(dist < 1.0E-100)	{
		printf("");
	}
	dist_Inv = 1.0/dist;
	vect_x *= dist_Inv;
	vect_y *= dist_Inv;
	vect_z *= dist_Inv;
	return dist;
}

inline double Dot_Product(double x_a, double y_a, double z_a, double x_b, double y_b, double z_b)
{
	return (x_a*x_b + y_a*y_b + z_a*z_b);
}


int Is_Three_Atoms_Colinear(int ia, int ib, int ic)
{
	double v1_x, v1_y, v1_z;
	double v2_x, v2_y, v2_z;
	double dot_v1_v2;
	double *x_Mol, *y_Mol, *z_Mol;

	x_Mol = Mol.x;
	y_Mol = Mol.y;
	z_Mol = Mol.z;

	v1_x = x_Mol[ib] - x_Mol[ia];
	v1_y = y_Mol[ib] - y_Mol[ia];
	v1_z = z_Mol[ib] - z_Mol[ia];

	v2_x = x_Mol[ic] - x_Mol[ib];
	v2_y = y_Mol[ic] - y_Mol[ib];
	v2_z = z_Mol[ic] - z_Mol[ib];

	Normalizex(v1_x, v1_y, v1_z);
	Normalizex(v2_x, v2_y, v2_z);

	dot_v1_v2 = fabs(Dot_Product(v1_x, v1_y, v1_z, v2_x, v2_y, v2_z));
	if(dot_v1_v2 > cos(15.0*radianInv))	{	// smaller than 10 degree
		return 1;
	}
	else	{
		return 0;
	}
}


int To_Refine_Soft_Dihedral_List(void)
{
	int i, Count=0;
	FILE *fOut;

	fOut = fopen("soft-dih-list.txt", "w");
	for(i=0; i<n_Soft_Dih; i++)	{
		if(Check_Is_Methyl_Involved(i) == 0)	{
			Get_Representative_14_Atoms(i);

			if( (Is_Three_Atoms_Colinear(Dih_Bond[i][0], Dih_Bond[i][1], Dih_Bond[i][2])==0) && (Is_Three_Atoms_Colinear(Dih_Bond[i][1], Dih_Bond[i][2], Dih_Bond[i][3])==0) )	{
				Dihedral_Selected[i] = 1;
				printf("%3d %3d %3d %3d \n", Dih_Bond[i][0]+1,  Dih_Bond[i][1]+1,  Dih_Bond[i][2]+1,  Dih_Bond[i][3]+1);
				fprintf(fOut, "%3d %3d %3d %3d \n", Dih_Bond[i][0]+1,  Dih_Bond[i][1]+1,  Dih_Bond[i][2]+1,  Dih_Bond[i][3]+1);
				Count++;
			}
			else	{
				Dihedral_Selected[i] = 0;
			}
		}
		else	{
			Dihedral_Selected[i] = 0;
		}
	}
	fclose(fOut);
	
	return Count;
}

int Get_Number_of_H_Bonded(int iAtom)
{
	int i, Count=0;
	double mass;

	for(i=0; i<Bond_Count[iAtom]; i++)	{
		mass = Mol.mass[Bond_List[iAtom][i]];
		if( (mass > 0.9) && (mass < 1.1) )	{
			Count++;
		}
	}
	return Count;
}

int Get_Number_of_F_Bonded(int iAtom)
{
	int i, Count=0;
	double mass;

	for(i=0; i<Bond_Count[iAtom]; i++)	{
		mass = Mol.mass[Bond_List[iAtom][i]];
		if( (mass > 18.8) && (mass < 19.2) )	{	// F
			Count++;
		}
	}
	return Count;
}

int Check_Is_Methyl_Involved(int iDih)
{
	int ib, ic, HCount=0, FCount=0;

	ib = Dih_Bond[iDih][1];
	ic = Dih_Bond[iDih][2];

	if( (Mol.mass[ib] > 11.5) && (Mol.mass[ib] < 12.5) )	{	// C atom
		HCount = Get_Number_of_H_Bonded(ib);
		FCount = Get_Number_of_F_Bonded(ib);
		if( ((HCount+FCount) == 3) && (To_Exclude_Methyl_Torsion==1) )	{
			return 1;
		}
	}
	if( (Mol.mass[ic] > 11.5) && (Mol.mass[ic] < 12.5) )	{	// C atom
		HCount = Get_Number_of_H_Bonded(ic);
		FCount = Get_Number_of_F_Bonded(ic);
		if( ((HCount+FCount) == 3) && (To_Exclude_Methyl_Torsion==1))	{
			return 1;
		}
	}
	
	return 0;
}


void Setup_Bond_List(void)
{
	int i, iPos, ia, ib, nAtom, nBond, *pBondList;

	nAtom = Mol.nAtom;
	nBond = Mol.nBond;
	pBondList = Mol.BondList;

	for(i=0; i<nAtom; i++)	{
		Bond_Count[i] = 0;
	}

	for(i=0; i<nBond; i++)	{
		iPos = 2*i;
		ia = pBondList[iPos];
		ib = pBondList[iPos+1];

		Bond_List[ia][Bond_Count[ia]] = ib;
		Bond_Count[ia]++;

		Bond_List[ib][Bond_Count[ib]] = ia;
		Bond_Count[ib]++;
	}

}

void Gen_Soft_Dihedral_List_Excluding_Ring(void)
{
	int i, nDih, iPos, ia, ib, ic, id;

	nDih       = Mol.nDihedral;
	n_Soft_Dih = 0;

	for(i=0; i<nDih; i++)	{
		if(Mol.Is_A_Dihedrals_In_A_Ring(i) == 0)	{
			iPos = i*4;
			ia = Mol.DihedralList[iPos  ];
			ib = Mol.DihedralList[iPos+1];
			ic = Mol.DihedralList[iPos+2];
			id = Mol.DihedralList[iPos+3];
			if(ib <= ic)	{
				if(!Is_In_Soft_Dihedral_List(ib, ic))	{
					Dih_Bond[n_Soft_Dih][0] = ia;
					Dih_Bond[n_Soft_Dih][1] = ib;
					Dih_Bond[n_Soft_Dih][2] = ic;
					Dih_Bond[n_Soft_Dih][3] = id;
					n_Soft_Dih++;
					printf("%3d %3d %3d %3d \n", Mol.DihedralList[iPos]+1, Mol.DihedralList[iPos+1]+1, Mol.DihedralList[iPos+2]+1, Mol.DihedralList[iPos+3]+1);
				}
			}
			else	{
				if(!Is_In_Soft_Dihedral_List(ic, ib))	{
					Dih_Bond[n_Soft_Dih][0] = id;
					Dih_Bond[n_Soft_Dih][1] = ic;
					Dih_Bond[n_Soft_Dih][2] = ib;
					Dih_Bond[n_Soft_Dih][3] = ia;
					n_Soft_Dih++;
					printf("%3d %3d %3d %3d \n", Mol.DihedralList[iPos]+1, Mol.DihedralList[iPos+1]+1, Mol.DihedralList[iPos+2]+1, Mol.DihedralList[iPos+3]+1);
				}
			}
		}
	}
}

void Gen_Soft_Dihedral_List_MD_High_T(void)
{
	int Step, Idx, RecIdx, iPos, ia, ib, ic, id;

	Mol.Init_LangevinDynamics(T_Sim);

	nDihedral = Mol.nDihedral;

	for(Step=1; Step<=MAX_STEP; Step++)	{
		Mol.LangevinDynamics(Step);
		if(Step%GAP==0)	{
			RecIdx = Step/GAP - 1;
			for(Idx=0; Idx<nDihedral; Idx++)	{
				dih_List[RecIdx][Idx] = Mol.dih_Phi_List[Idx];
			}
		}
	}

	Cal_Sig();
	Cal_Sig_Shift_PI2();

	n_Soft_Dih = 0;
	for(Idx=0; Idx<nDihedral; Idx++)	{
		if(IsRigidDih[Idx] == 0)	{
			iPos = Idx*4;
			ia = Mol.DihedralList[iPos  ];
			ib = Mol.DihedralList[iPos+1];
			ic = Mol.DihedralList[iPos+2];
			id = Mol.DihedralList[iPos+3];
			if(ib <= ic)	{
				if(!Is_In_Soft_Dihedral_List(ib, ic))	{
					Dih_Bond[n_Soft_Dih][0] = ia;
					Dih_Bond[n_Soft_Dih][1] = ib;
					Dih_Bond[n_Soft_Dih][2] = ic;
					Dih_Bond[n_Soft_Dih][3] = id;
					n_Soft_Dih++;
					printf("%3d %3d %3d %3d \n", Mol.DihedralList[iPos]+1, Mol.DihedralList[iPos+1]+1, Mol.DihedralList[iPos+2]+1, Mol.DihedralList[iPos+3]+1);
				}
			}
			else	{
				if(!Is_In_Soft_Dihedral_List(ic, ib))	{
					Dih_Bond[n_Soft_Dih][0] = id;
					Dih_Bond[n_Soft_Dih][1] = ic;
					Dih_Bond[n_Soft_Dih][2] = ib;
					Dih_Bond[n_Soft_Dih][3] = ia;
					n_Soft_Dih++;
					printf("%3d %3d %3d %3d \n", Mol.DihedralList[iPos]+1, Mol.DihedralList[iPos+1]+1, Mol.DihedralList[iPos+2]+1, Mol.DihedralList[iPos+3]+1);
				}
			}
		}
	}
}

int Is_In_Soft_Dihedral_List(int b, int c)
{
	int i;

	for(i=0; i<n_Soft_Dih; i++)	{
		if( (b==Dih_Bond[i][1]) && (c==Dih_Bond[i][2]) )	{
			return 1;
		}
	}
	return 0;
}

void Quit_With_Error_Msg(char szMsg[])
{
	fprintf(fFile_Run_Log, "%s", szMsg);
	fflush(fFile_Run_Log);

	FILE *fOut;
	fOut = fopen("error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in gen_soft_list.cpp\n");
	fclose(fOut);

	exit(1);
}

void Cal_Sig(void)
{
	int i, j;
	double d_dih;

	for(i=0; i<nDihedral; i++)	{
		dih_Mean[i] = 0.0;

		for(j=0; j<MAX_STEP/GAP; j++)	{
			dih_Mean[i] += dih_List[j][i];
		}
		dih_Mean[i] /= (MAX_STEP/GAP);
	}

	for(i=0; i<nDihedral; i++)	{
		dih_Sig[i] = 0.0;

		for(j=0; j<MAX_STEP/GAP; j++)	{
			d_dih = dih_List[j][i] - dih_Mean[i];
			dih_Sig[i] += (d_dih*d_dih);
		}
		dih_Sig[i] = sqrt(dih_Sig[i]/(MAX_STEP/GAP));

		if(dih_Sig[i] < SIG_CUT)	{
			IsRigidDih[i] = 1;
		}
	}
}

void Cal_Sig_Shift_PI2(void)
{
	int i, j;
	double d_dih;

	for(i=0; i<nDihedral; i++)	{
		dih_Mean[i] = 0.0;

		for(j=0; j<MAX_STEP/GAP; j++)	{
			if(dih_List[j][i] < 0.0)	{
				dih_List[j][i] += 360.0;	// shift by 2*PI
			}
			dih_Mean[i] += dih_List[j][i];
		}
		dih_Mean[i] /= (MAX_STEP/GAP);
	}

	for(i=0; i<nDihedral; i++)	{
		dih_Sig[i] = 0.0;

		for(j=0; j<MAX_STEP/GAP; j++)	{
			d_dih = dih_List[j][i] - dih_Mean[i];
			dih_Sig[i] += (d_dih*d_dih);
		}
		dih_Sig[i] = sqrt(dih_Sig[i]/(MAX_STEP/GAP));

		if(dih_Sig[i] < SIG_CUT)	{
			IsRigidDih[i] = 1;
		}
	}

	for(i=0; i<nDihedral; i++)	{	// restore the original dihedral
		for(j=0; j<MAX_STEP/GAP; j++)	{
			if(dih_List[j][i] >= 180.0)	{
				dih_List[j][i] -= 360.0;	// shift by -2*PI
			}
		}
	}

}

