/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ff.h"

#define radian	(57.29577951308232088)
#define radianInv	(0.017453292519943295)
#define MAX_ENTRY	(4096)
#define b0_DEV_MAX	(0.05)
#define theta0_DEV_MAX	(8.0*radianInv)
#define SCALE		(0.4)


#define MAX_BOND_ADD	(MAX_ATOM*16)
#define MAX_ANGLE_ADD	(MAX_ATOM*16)
#define MAX_DIHEDRAL_ADD	(MAX_ATOM*144)
#define MAX_IMPROPER_ADD	(MAX_ATOM*96)
#define MAX_LJ_ADD		(MAX_ATOM)

#define MAX_LINE_RTF	(4096)
#define MAX_LINE_PRM	(4096)
#define MAX_LINE_XPSF	(4096)


BOND_REC Bond_Rec_Add[MAX_BOND_ADD];
ANGLE_REC Angle_Rec_Add[MAX_ANGLE_ADD];
DIHEDRAL_REC Dihedral_Rec_Add[MAX_DIHEDRAL_ADD];
IMPRODIHEDRAL_REC ImproDihedral_Rec_Add[MAX_IMPROPER_ADD];

int nBond_Add=0, nAngle_Add=0, nDihedral_Add=0, nImproper_Add=0;
int To_Add_Bond[MAX_BOND_ADD], To_Add_Angle[MAX_ANGLE_ADD], To_Add_Dihedral[MAX_DIHEDRAL_ADD], To_Add_Improper[MAX_IMPROPER_ADD];

char szAddTxt_LJ[MAX_ATOM*128];


CMol Mol;
CForceField ForceField;
FILE *fFile_Run_Log;
double b0_QM[MAX_ENTRY], theta0_QM[MAX_ENTRY];
int New_Atom_Type[MAX_ATOM];
char szOldChemName[MAX_ATOM][12], szNewChemName[MAX_ATOM][12];
int nNewAtomType=0, Atom_To_Update[MAX_ATOM];


void Quit_With_Error_Msg(char szMsg[]);
void Build_List_of_Atom_To_Update(void);
int Is_ChemName_Used(char szChemName[], int Exclude);
void Generate_A_New_ChemName(int iAtom);
void Add_New_Atom_Type_Prm_File(void);
void Update_Rtf_File_with_New_Atom_Type(void);
void Update_Xpsf_File_with_New_Atom_Type(void);
void Update_Water_Xpsf_File_with_New_Atom_Type(void);
void Determine_Bond_Angle_Dihedral_Improper_Add(int ActiveAtom);

int main(void)
{
	FILE *fOut;
	int i, nBond, nAngle, CountChangedBond=0, CountChangedAngle=0;
	int ia, ib, ic, iPos;
	int *pBondList, *pAngleList;
	fFile_Run_Log = fopen("log-check_b0_theta0.txt", "w");

  timebomb();

	memset(New_Atom_Type, 0, sizeof(int)*MAX_ATOM);
	ForceField.ReadForceField("mol.prm");

	Mol.ReadPSF("mol.xpsf", 0);
	Mol.AssignForceFieldParameters(&ForceField);
	pBondList = Mol.BondList;
	pAngleList = Mol.AngleList;

	nBond=Mol.nBond;
	nAngle=Mol.nAngle;

//	Mol.ReadCRD("mol-short.crd");
	Mol.ReadXYZ("mol-opt.xyz");
	Mol.Cal_E(1);
//	Mol.SavePdb("mol-opt.pdb");
	for(i=0; i<nBond; i++)	{
		b0_QM[i] = Mol.b0_MM[i];
	}
	for(i=0; i<nAngle; i++)	{
		theta0_QM[i] = Mol.theta0_MM[i];
	}


	Mol.FullGeometryOptimization_LBFGS();
	Mol.Cal_E(1);


	fOut = fopen("bond.txt", "w");
	for(i=0; i<nBond; i++)	{
		fprintf(fOut, "%3d %8.5lf %8.5lf %8.5lf\n", i+1, Mol.Para_b0[i], Mol.b0_MM[i], b0_QM[i]);
		if(fabs(Mol.b0_MM[i]-b0_QM[i]) > b0_DEV_MAX)	{
			iPos = i * 2;
			ia = pBondList[iPos  ];
			ib = pBondList[iPos+1];
			printf("Bond %8s %8s   %.4lf %.4lf %.4lf\n", Mol.ChemName[ia], Mol.ChemName[ib], Mol.Para_b0[i], Mol.b0_MM[i], b0_QM[i]);

			if(fabs(Mol.Para_b0[i]-b0_QM[i]) > b0_DEV_MAX*SCALE)	{
				Mol.Para_b0[i] = b0_QM[i];
				New_Atom_Type[ia] = 1;
				New_Atom_Type[ib] = 1;
				printf("changed: bond  %d\n", i+1);
				CountChangedBond++;
			}
		}
	}
	fclose(fOut);

	if(CountChangedBond > 0)	{
		printf("\n\n");
		Mol.FullGeometryOptimization_LBFGS();	// Geom. Opt. after adjusting bond length 
		Mol.Cal_E(1);
	}

	printf("\n\n");

	fOut = fopen("angle.txt", "w");
	for(i=0; i<nAngle; i++)	{
		fprintf(fOut, "%3d %10.2lf %10.2lf %10.2lf\n", i+1, Mol.Para_theta0[i]*radian, Mol.theta0_MM[i]*radian, theta0_QM[i]*radian);
		if(fabs(Mol.theta0_MM[i]-theta0_QM[i]) > theta0_DEV_MAX)	{
			iPos = i * 3;
			ia = pAngleList[iPos  ];
			ib = pAngleList[iPos+1];
			ic = pAngleList[iPos+2];
			printf("Angle %8s %8s %8s   %.2lf %.2lf %.2lf\n", Mol.ChemName[ia], Mol.ChemName[ib], Mol.ChemName[ic], Mol.Para_theta0[i]*radian, Mol.theta0_MM[i]*radian, theta0_QM[i]*radian);

			if(fabs(Mol.Para_theta0[i]-theta0_QM[i]) > theta0_DEV_MAX*SCALE)	{
				Mol.Para_theta0[i] = theta0_QM[i];
				New_Atom_Type[ia] = 1;
				New_Atom_Type[ic] = 1;
				printf("changed: angle %d\n", i+1);
				CountChangedAngle++;
			}
		}
	}
	fclose(fOut);

	if( (CountChangedBond==0) && (CountChangedAngle==0) )	{	// nothing is changed
		printf("No bond length or angle parameter is changed.\n");
		return 0;
	}

	printf("\n\nAfter adjusting theta_0:\n\n");

	Mol.FullGeometryOptimization_LBFGS();
	Mol.Cal_E(1);

	fOut = fopen("bond-new.txt", "w");
	for(i=0; i<nBond; i++)	{
		fprintf(fOut, "%3d %8.5lf %8.5lf %8.5lf\n", i+1, Mol.Para_b0[i], Mol.b0_MM[i], b0_QM[i]);
		if(fabs(Mol.b0_MM[i]-b0_QM[i]) > b0_DEV_MAX)	{
			iPos = i * 2;
			ia = pBondList[iPos  ];
			ib = pBondList[iPos+1];
			printf("Bond %8s %8s\n", Mol.ChemName[ia], Mol.ChemName[ib]);
		}
	}
	fclose(fOut);

	fOut = fopen("angle-new.txt", "w");
	for(i=0; i<nAngle; i++)	{
		fprintf(fOut, "%3d %10.2lf %10.2lf %10.2lf\n", i+1, Mol.Para_theta0[i]*radian, Mol.theta0_MM[i]*radian, theta0_QM[i]*radian);
		if(fabs(Mol.theta0_MM[i]-theta0_QM[i]) > theta0_DEV_MAX)	{
			iPos = i * 3;
			ia = pAngleList[iPos  ];
			ib = pAngleList[iPos+1];
			ic = pAngleList[iPos+2];
			printf("Angle %8s %8s %8s\n", Mol.ChemName[ia], Mol.ChemName[ib], Mol.ChemName[ic]);

			Mol.Para_theta0[i] = theta0_QM[i];
		}
	}
	fclose(fOut);

//	Mol.WriteCRDFile("mol-short.crd");
//	Mol.WriteMyCRDFile("mol-test.crd");
//	Mol.ReadCRD("mol-test.crd");

	printf("\n");
	Mol.Cal_E(1);

	Build_List_of_Atom_To_Update();

	for(i=0; i<nNewAtomType; i++)	{
		Determine_Bond_Angle_Dihedral_Improper_Add(Atom_To_Update[i]);
	}

	Add_New_Atom_Type_Prm_File();
	Update_Rtf_File_with_New_Atom_Type();
	Update_Xpsf_File_with_New_Atom_Type();
	Update_Water_Xpsf_File_with_New_Atom_Type();


	printf("\n\nRead in modified prm and xpsf files, \n");
	ForceField.ReadForceField("mol.prm");
	Mol.ReadPSF("mol.xpsf", 0);
	Mol.AssignForceFieldParameters(&ForceField);
	Mol.ReadXYZ("mol-opt.xyz");
	Mol.Cal_E(1);

	fclose(fFile_Run_Log);

	return 0;
}


void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in check_b0_theta0.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}

void Build_List_of_Atom_To_Update(void)
{
	int i, nAtom;
	char szChemPara[2][12], szBuff[256];
	double Para_List[6];

	nNewAtomType=0;
	strcpy(szAddTxt_LJ, "");

	nAtom = Mol.nAtom;
	for(i=0; i<nAtom; i++)	{
		if(New_Atom_Type[i] == 1)	{	// to be updated. A new atom type to be added
			Atom_To_Update[nNewAtomType] = i;
			Generate_A_New_ChemName(i);

			strcpy(szChemPara[0], szOldChemName[i]);
			ForceField.GetPara_LJ(szChemPara, Para_List);

			if( (Para_List[1] == Para_List[4]) && (Para_List[2] == Para_List[5]) )	{
				sprintf(szBuff, "%-8s   0.00 %9.4lf %9.4lf\n", szNewChemName[i], Para_List[1], Para_List[2]);
			}
			else	{
				sprintf(szBuff, "%-8s   0.00 %9.4lf %9.4lf   0.00 %9.4lf %9.4lf\n", szNewChemName[i], Para_List[1], Para_List[2], Para_List[4], Para_List[5]);
			}
			strcat(szAddTxt_LJ, szBuff);
			
			nNewAtomType++;
		}
	}
}

int Is_ChemName_Used(char szChemName[], int Exclude)
{
	int nAtom, i;

	nAtom = Mol.nAtom;

	for(i=0; i<nAtom; i++)	{
		if(i == Exclude)	{
			continue;
		}
		if(strcmp(szChemName, Mol.ChemName[i])==0)	{
			return 1;
		}
	}
	return 0;
}

void Generate_A_New_ChemName(int iAtom)
{
	int i, PoolSize;
	char CharPool[]="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ", CharEnd, NewChemName[256], BaseName[256], szErrorMsg[256];

	PoolSize = strlen(CharPool);
	strcpy(szOldChemName[iAtom], Mol.ChemName[iAtom]);

	strcpy(BaseName, szOldChemName[iAtom]);
	BaseName[4]=0;	// only keep maximum 4 lettters in the original name

	for(i=0; i<PoolSize; i++)	{
		CharEnd = CharPool[i];
		sprintf(NewChemName, "%sx%c", BaseName, CharEnd);
		if(Is_ChemName_Used(NewChemName, iAtom) == 0)	{	// not used. Done.
			strcpy(Mol.ChemName[iAtom], NewChemName);
			strcpy(szNewChemName[iAtom], NewChemName);
			return;
		}
		else	{	// used. Keep trying
		}
	}

	sprintf(szErrorMsg, "Fail to create a new name for atom %d, %s\nQuit\n", iAtom+1, szOldChemName[iAtom]);
	Quit_With_Error_Msg(szErrorMsg);
}

void Determine_Bond_Angle_Dihedral_Improper_Add(int ActiveAtom)
{
	int i, j, iPos, ia, ib, ic, id, LocalCount;
	int nBond, nAngle, nDihedral, nImproper;

	nBond = Mol.nBond;
	nAngle = Mol.nAngle;
	nDihedral = Mol.nDihedral;
	nImproper = Mol.nImpro;

	//start	to generate the list of new bond entries
	for(i=0; i<nBond; i++)	{
		iPos = 2*i;
		ia = Mol.BondList[iPos  ];
		ib = Mol.BondList[iPos+1];

		if( (ActiveAtom == ia) || (ActiveAtom == ib) )	{
			strcpy(Bond_Rec_Add[nBond_Add].Chem[0], Mol.ChemName[ia]);
			strcpy(Bond_Rec_Add[nBond_Add].Chem[1], Mol.ChemName[ib]);

			Bond_Rec_Add[nBond_Add].para[0] = Mol.Para_k_b[i];
			Bond_Rec_Add[nBond_Add].para[1] = Mol.Para_b0[i];
			To_Add_Bond[nBond_Add] = 1;
			nBond_Add++;
		}
	}
	//end	to generate the list of new bond entries

	//start	to re-organize the new bond entries, delete redundant entries
	for(i=0; i<nBond_Add; i++)	{
		for(j=i+1; j<nBond_Add; j++)	{	// to find out those redundant entries
			if(To_Add_Bond[j] == 0)	{	// already deleted
				continue;
			}
			if( ( (strcmp(Bond_Rec_Add[i].Chem[0], Bond_Rec_Add[j].Chem[0])==0) && (strcmp(Bond_Rec_Add[i].Chem[1], Bond_Rec_Add[j].Chem[1])==0) ) || 
				( (strcmp(Bond_Rec_Add[i].Chem[0], Bond_Rec_Add[j].Chem[1])==0) && (strcmp(Bond_Rec_Add[i].Chem[1], Bond_Rec_Add[j].Chem[0])==0) ) )	{
				if( (fabs(Bond_Rec_Add[i].para[0]-Bond_Rec_Add[j].para[0]) < 1.0E-7) && (fabs(Bond_Rec_Add[i].para[1]-Bond_Rec_Add[j].para[1]) < 1.0E-7) )	{
					To_Add_Bond[j] = 0;
				}
			}
		}
	}
	//end	to re-organize the new bond entries, delete redundant entries


	//start	to generate the list of new angle entries
	for(i=0; i<nAngle; i++)	{
		iPos = 3*i;
		ia = Mol.AngleList[iPos  ];
		ib = Mol.AngleList[iPos+1];
		ic = Mol.AngleList[iPos+2];

		if( (ActiveAtom == ia) || (ActiveAtom == ib) || (ActiveAtom == ic) )	{
			strcpy(Angle_Rec_Add[nAngle_Add].Chem[0], Mol.ChemName[ia]);
			strcpy(Angle_Rec_Add[nAngle_Add].Chem[1], Mol.ChemName[ib]);
			strcpy(Angle_Rec_Add[nAngle_Add].Chem[2], Mol.ChemName[ic]);

			Angle_Rec_Add[nAngle_Add].para[0] = Mol.Para_k_a[i];
			Angle_Rec_Add[nAngle_Add].para[1] = Mol.Para_theta0[i]*radian;
			Angle_Rec_Add[nAngle_Add].para[2] = Mol.Para_k_Urey[i];
			Angle_Rec_Add[nAngle_Add].para[3] = Mol.Para_b0_Urey[i];
			To_Add_Angle[nAngle_Add] = 1;
			nAngle_Add++;
		}
	}
	//end	to generate the list of new angle entries

	//start	to re-organize the new angle entries, delete redundant entries
	for(i=0; i<nAngle_Add; i++)	{
		for(j=i+1; j<nAngle_Add; j++)	{	// to find out those redundant entries
			if(To_Add_Angle[j] == 0)	{	// already deleted
				continue;
			}
			if( ( (strcmp(Angle_Rec_Add[i].Chem[0], Angle_Rec_Add[j].Chem[0])==0) && (strcmp(Angle_Rec_Add[i].Chem[1], Angle_Rec_Add[j].Chem[1])==0) && (strcmp(Angle_Rec_Add[i].Chem[2], Angle_Rec_Add[j].Chem[2])==0) ) || 
				( (strcmp(Angle_Rec_Add[i].Chem[0], Angle_Rec_Add[j].Chem[2])==0) && (strcmp(Angle_Rec_Add[i].Chem[1], Angle_Rec_Add[j].Chem[1])==0) && (strcmp(Angle_Rec_Add[i].Chem[2], Angle_Rec_Add[j].Chem[0])==0) ) )	{
//				if( (fabs(Angle_Rec_Add[i].para[0]-Angle_Rec_Add[j].para[0]) < 1.0E-7) && (fabs(Angle_Rec_Add[i].para[1]-Angle_Rec_Add[j].para[1]) < 1.0E-7) )	{
//					To_Add_Angle[j] = 0;
//				}
				if( (fabs(Angle_Rec_Add[i].para[0]-Angle_Rec_Add[j].para[0]) < 1.0E-7) && (fabs(Angle_Rec_Add[i].para[1]-Angle_Rec_Add[j].para[1]) < 3.0) )	{	// assume angles for three atoms with same atom types having same angle
					Angle_Rec_Add[i].para[1] = 0.5*(Angle_Rec_Add[i].para[1] + Angle_Rec_Add[j].para[1]);
					To_Add_Angle[j] = 0;
				}
				else	{
					Quit_With_Error_Msg("Unhandled situation in updating angle parameters in check_b0_theta0.cpp.\nQuit\n\n");
				}
			}
		}
	}
	//end	to re-organize the new Angle entries, delete redundant entries


	//start	to generate the list of new Dihedral entries
	for(i=0; i<nDihedral; i++)	{
		iPos = 4*i;
		ia = Mol.DihedralList[iPos  ];
		ib = Mol.DihedralList[iPos+1];
		ic = Mol.DihedralList[iPos+2];
		id = Mol.DihedralList[iPos+3];
		
		if( (ActiveAtom == ia) || (ActiveAtom == ib) || (ActiveAtom == ic)  || (ActiveAtom == id) )	{
			LocalCount = 0;
			for(j=1; j<=6; j++)	{
				Dihedral_Rec_Add[nDihedral_Add].para[0] = Mol.Para_k_Dih[i][j];
				Dihedral_Rec_Add[nDihedral_Add].para[1] = 1.0*j;
				Dihedral_Rec_Add[nDihedral_Add].para[2] = Mol.Para_phi[i][j]*radian;
				if(Dihedral_Rec_Add[nDihedral_Add].para[0] > 1.0E-6)	{	// a valid entry
					strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[0], Mol.ChemName[ia]);
					strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[1], Mol.ChemName[ib]);
					strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[2], Mol.ChemName[ic]);
					strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[3], Mol.ChemName[id]);

					To_Add_Dihedral[nDihedral_Add] = 1;
					nDihedral_Add++;
					LocalCount++;
				}
			}
			if(LocalCount == 0)	{	// We need to add a dummy entry with k=0. We need to find the entry in rtf file
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[0], Mol.ChemName[ia]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[1], Mol.ChemName[ib]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[2], Mol.ChemName[ic]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[3], Mol.ChemName[id]);
				
				Dihedral_Rec_Add[nDihedral_Add].para[0] = 0.0;
				Dihedral_Rec_Add[nDihedral_Add].para[1] = 2;
				Dihedral_Rec_Add[nDihedral_Add].para[2] = 180.0;
				To_Add_Dihedral[nDihedral_Add] = 1;
				nDihedral_Add++;
//				Quit_With_Error_Msg("LocalCount == 0\nCode is not completed yet in Determine_Bond_Angle_Dihedral_Improper_Add_1().\nQuit\n\n");
			}
		}
	}
	//end	to generate the list of new Dihedral entries

	//start	to re-organize the new dihedral entries, delete redundant entries
	for(i=0; i<nDihedral_Add; i++)	{
		for(j=i+1; j<nDihedral_Add; j++)	{	// to find out those redundant entries
			if(To_Add_Dihedral[j] == 0)	{	// already deleted
				continue;
			}
			if( ( (strcmp(Dihedral_Rec_Add[i].Chem[0], Dihedral_Rec_Add[j].Chem[0])==0) && (strcmp(Dihedral_Rec_Add[i].Chem[1], Dihedral_Rec_Add[j].Chem[1])==0) && (strcmp(Dihedral_Rec_Add[i].Chem[2], Dihedral_Rec_Add[j].Chem[2])==0) && (strcmp(Dihedral_Rec_Add[i].Chem[3], Dihedral_Rec_Add[j].Chem[3])==0) ) || 
				( (strcmp(Dihedral_Rec_Add[i].Chem[0], Dihedral_Rec_Add[j].Chem[3])==0) && (strcmp(Dihedral_Rec_Add[i].Chem[1], Dihedral_Rec_Add[j].Chem[2])==0) && (strcmp(Dihedral_Rec_Add[i].Chem[2], Dihedral_Rec_Add[j].Chem[1])==0) && (strcmp(Dihedral_Rec_Add[i].Chem[3], Dihedral_Rec_Add[j].Chem[0])==0) ) )	{
				if( (fabs(Dihedral_Rec_Add[i].para[0]-Dihedral_Rec_Add[j].para[0]) < 1.0E-3) && (fabs(Dihedral_Rec_Add[i].para[1]-Dihedral_Rec_Add[j].para[1]) < 1.0E-3) && (fabs(Dihedral_Rec_Add[i].para[2]-Dihedral_Rec_Add[j].para[2]) < 1.0E-3) )	{
					To_Add_Dihedral[j] = 0;
				}
			}
		}
	}
	//end	to re-organize the new Dihedral entries, delete redundant entries

	

	//start	to generate the list of new Improper entries
	for(i=0; i<nImproper; i++)	{
		iPos = 4*i;
		ia = Mol.ImprDihedralList[iPos  ];
		ib = Mol.ImprDihedralList[iPos+1];
		ic = Mol.ImprDihedralList[iPos+2];
		id = Mol.ImprDihedralList[iPos+3];
		
		if( (ActiveAtom == ia) || (ActiveAtom == ib) || (ActiveAtom == ic)  || (ActiveAtom == id) )	{
			ImproDihedral_Rec_Add[nImproper_Add].para[0] = Mol.Para_k_ImpDih[i];
			ImproDihedral_Rec_Add[nImproper_Add].para[1] = Mol.Para_Imp_phi[i]*radian;
			ImproDihedral_Rec_Add[nImproper_Add].para[2] = Mol.Para_Type_ImpDih[i];
			if(Mol.Para_k_ImpDih[i] > 1.0E-3)	{	// a valid entry
				strcpy(ImproDihedral_Rec_Add[nImproper_Add].Chem[0], Mol.ChemName[ia]);
				strcpy(ImproDihedral_Rec_Add[nImproper_Add].Chem[1], Mol.ChemName[ib]);
				strcpy(ImproDihedral_Rec_Add[nImproper_Add].Chem[2], Mol.ChemName[ic]);
				strcpy(ImproDihedral_Rec_Add[nImproper_Add].Chem[3], Mol.ChemName[id]);
				
				To_Add_Improper[nImproper_Add] = 1;
				nImproper_Add++;
			}
		}
	}
	//end	to generate the list of new Improper entries

	//start	to re-organize the new improper entries, delete redundant entries
	for(i=0; i<nImproper_Add; i++)	{
		for(j=i+1; j<nImproper_Add; j++)	{	// to find out those redundant entries
			if(To_Add_Improper[j] == 0)	{	// already deleted
				continue;
			}
			if( ( (strcmp(ImproDihedral_Rec_Add[i].Chem[0], ImproDihedral_Rec_Add[j].Chem[0])==0) && (strcmp(ImproDihedral_Rec_Add[i].Chem[1], ImproDihedral_Rec_Add[j].Chem[1])==0) && (strcmp(ImproDihedral_Rec_Add[i].Chem[2], ImproDihedral_Rec_Add[j].Chem[2])==0) && (strcmp(ImproDihedral_Rec_Add[i].Chem[3], ImproDihedral_Rec_Add[j].Chem[3])==0) ) || 
				( (strcmp(ImproDihedral_Rec_Add[i].Chem[0], ImproDihedral_Rec_Add[j].Chem[3])==0) && (strcmp(ImproDihedral_Rec_Add[i].Chem[1], ImproDihedral_Rec_Add[j].Chem[2])==0) && (strcmp(ImproDihedral_Rec_Add[i].Chem[2], ImproDihedral_Rec_Add[j].Chem[1])==0) && (strcmp(ImproDihedral_Rec_Add[i].Chem[3], ImproDihedral_Rec_Add[j].Chem[0])==0) ) )	{
				if( (fabs(ImproDihedral_Rec_Add[i].para[0]-ImproDihedral_Rec_Add[j].para[0]) < 1.0E-3) && (fabs(ImproDihedral_Rec_Add[i].para[1]-ImproDihedral_Rec_Add[j].para[1]) < 1.0E-3) && (fabs(ImproDihedral_Rec_Add[i].para[2]-ImproDihedral_Rec_Add[j].para[2]) < 1.0E-3) )	{
					To_Add_Improper[j] = 0;
				}
			}
		}
	}
	//end	to re-organize the new improper entries, delete redundant entries
}




char szXpsfTxt[MAX_LINE_XPSF][256];
void Update_Xpsf_File_with_New_Atom_Type(void)
{
	FILE *fIn, *fOut;
	int i, nLine, LineToUpdate, ReadItem, nAtom, Index, ResID, Fixed;
	char *ReadLine, szTag[256], MolName[24], ResName[24], AtomName[24], ChemName[24];
	double charge, mass, last_1, last_2;

	if(nNewAtomType == 0)	{
		return;
	}

	fIn = fopen("mol.xpsf", "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open file: mol.xpsf in Update_Xpsf_File_with_New_Atom_Type().\nQuit\n");
	}
	nLine = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szXpsfTxt[nLine], 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		else	{
			nLine++;
			if(nLine >= MAX_LINE_XPSF)	{
				fclose(fIn);
				Quit_With_Error_Msg("nLine >= MAX_LINE_XPSF in Update_Xpsf_File_with_New_Atom_Type().\nQuit\n");
			}
		}
	}
	fclose(fIn);

	LineToUpdate = -1;
	for(i=0; i<nLine; i++)	{
		ReadItem = sscanf(szXpsfTxt[i], "%d%s", &nAtom, szTag);
		if(ReadItem == 2)	{	// the beginning of atoms
			if( (nAtom == Mol.nAtom) || (strcmp(szTag, "!NATOM")==0) )	{
				LineToUpdate = i+1;
				break;
			}
		}
	}

	for(i=0; i<nAtom; i++)	{
		ReadItem = sscanf(szXpsfTxt[LineToUpdate+i], "%d %s %d %s%s%s%lf%lf%d%lf%lf", 
			&Index, MolName, &ResID, ResName, AtomName, ChemName, &charge, &mass, &Fixed, &last_2, &last_1);
//		sprintf(szXpsfTxt[LineToUpdate+i], "%8d %-4s %2d    %-5s %-5s %-5s %10.6lf %9.4lf  %10d%10.5lf    %.6E\n", 
//			Index, MolName, ResID, ResName, AtomName, Mol.ChemName[i], charge, mass, Fixed, last_2, last_1);
		sprintf(szXpsfTxt[LineToUpdate+i], "%10d %-6s %3d        %-8s %-8s %-8s %12.5lf%10.5lf           0%10.5lf %14.5lf\n", 
			Index, MolName, ResID, ResName, AtomName, Mol.ChemName[i], charge, mass, last_2, last_1);
	}

	fOut = fopen("mol.xpsf", "w");
	for(i=0; i<nLine; i++)	{
		fprintf(fOut, "%s", szXpsfTxt[i]);
	}
	fclose(fOut);

}

void Update_Water_Xpsf_File_with_New_Atom_Type(void)
{
	FILE *fIn, *fOut;
	int i, nLine, LineToUpdate, ReadItem, nAtom, Index, ResID, Fixed;
	char *ReadLine, szTag[256], MolName[24], ResName[24], AtomName[24], ChemName[24];
	double charge, mass, last_1, last_2;

	if(nNewAtomType == 0)	{
		return;
	}

	fIn = fopen("mol-wat.xpsf", "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open file: mol-wat.xpsf in Update_Water_Xpsf_File_with_New_Atom_Type().\nQuit\n");
	}
	nLine = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szXpsfTxt[nLine], 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		else	{
			nLine++;
			if(nLine >= MAX_LINE_XPSF)	{
				fclose(fIn);
				Quit_With_Error_Msg("nLine >= MAX_LINE_XPSF in Update_Water_Xpsf_File_with_New_Atom_Type().\nQuit\n");
			}
		}
	}
	fclose(fIn);

	LineToUpdate = -1;
	for(i=0; i<nLine; i++)	{
		ReadItem = sscanf(szXpsfTxt[i], "%d%s", &nAtom, szTag);
		if(ReadItem == 2)	{	// the beginning of atoms
			if( (nAtom == (Mol.nAtom+3)) || (strcmp(szTag, "!NATOM")==0) )	{
				LineToUpdate = i+1;
				break;
			}
		}
	}

	for(i=0; i<Mol.nAtom; i++)	{
		ReadItem = sscanf(szXpsfTxt[LineToUpdate+i], "%d %s %d %s%s%s%lf%lf%d%lf%lf", 
			&Index, MolName, &ResID, ResName, AtomName, ChemName, &charge, &mass, &Fixed, &last_2, &last_1);
//		sprintf(szXpsfTxt[LineToUpdate+i], "%8d %-4s %2d    %-5s %-5s %-5s %10.6lf %9.4lf  %10d%10.5lf    %.6E\n", 
//			Index, MolName, ResID, ResName, AtomName, Mol.ChemName[i], charge, mass, Fixed, last_2, last_1);
		sprintf(szXpsfTxt[LineToUpdate+i], "%10d %-6s %3d        %-8s %-8s %-8s %12.5lf%10.5lf           0%10.5lf %14.5lf\n", 
			Index, MolName, ResID, ResName, AtomName, Mol.ChemName[i], charge, mass, last_2, last_1);
	}

	fOut = fopen("mol-wat.xpsf", "w");
	for(i=0; i<nLine; i++)	{
		fprintf(fOut, "%s", szXpsfTxt[i]);
	}
	fclose(fOut);

}


char szRtfTxt[MAX_LINE_RTF][256];
void Update_Rtf_File_with_New_Atom_Type(void)
{
	FILE *fIn, *fOut;
	char *ReadLine, szLineSave[256];
	int i, j, nLine, LineLastMassMol, AtomeTypeCount, ActiveAtom, AtomCount, Idx_Line_To_Change;

	if(nNewAtomType > 0)	{
		fIn = fopen("mol.rtf", "r");
		if(fIn == NULL)	{
			Quit_With_Error_Msg("Fail to open file: mol.rtf\nQuit\n");
		}
		
		nLine = 0;
		while(1)	{
			if(feof(fIn))	{
				break;
			}
			ReadLine = fgets(szRtfTxt[nLine], 256, fIn);
			if(ReadLine == NULL)	{
				break;
			}
			else	{
				nLine++;
				if(nLine >= MAX_LINE_RTF)	{
					fclose(fIn);
					Quit_With_Error_Msg("nLine >= MAX_LINE_RTF in Update_Rtf_File_with_New_Atom_Type().\nQuit\n");
				}
			}
		}
		fclose(fIn);
		
		//start to find the index of the last line of the entry of MASS (atom type) for MOL
		LineLastMassMol = -1;
		AtomeTypeCount = 0;
		for(i=0; i<nLine; i++)	{
			if(strncmp(szRtfTxt[i], "RESI ", 5)==0)	{
				break;
			}
			if(strncmp(szRtfTxt[i], "MASS", 4)==0)	{
				LineLastMassMol = i;	// the place where the new atom type should be added
				AtomeTypeCount++;
			}
		}
		if(LineLastMassMol < 0)	{
			Quit_With_Error_Msg("Fail to find an entry for MASS in mol.rtf.\nQuit\n");
		}
		//end to find the index of the last line of the entry of MASS (atom type) for MOL
		
		//start	to find the line of the atom needed change
		AtomCount = 0;
		Idx_Line_To_Change = -1;
		for(; i<nLine; i++)	{
			if(strncmp(szRtfTxt[i], "BOND", 4)==0)	{	// the end of the definition of atoms in MOL
				break;
			}
			if(strncmp(szRtfTxt[i], "ATOM", 4)==0)	{
				AtomCount++;
				
				if(New_Atom_Type[AtomCount-1] == 1)	{	// the atom to be changed
					ActiveAtom = AtomCount-1;
					strcpy(szLineSave, szRtfTxt[i]);
					//				sprintf(szRtfTxt[i], "ATOM %-6s%-6s%7.3lf\n", Mol.AtomName[ActiveAtom], Mol.ChemName[ActiveAtom], Mol.CG[ActiveAtom]);
					sprintf(szRtfTxt[i], "ATOM %-6s %-8s%10.6lf    !  %s\n", Mol.AtomName[ActiveAtom], Mol.ChemName[ActiveAtom], Mol.CG[ActiveAtom], szOldChemName[ActiveAtom]);
					printf("Replacing the atom type of atom %d from %8s to %8s\n", ActiveAtom+1, szOldChemName[ActiveAtom], Mol.ChemName[ActiveAtom]);
					printf("%s%s", szLineSave, szRtfTxt[i]);
					Idx_Line_To_Change = i;
				}
			}
		}
		if(Idx_Line_To_Change < 0)	{
			Quit_With_Error_Msg("Fail to indentify the entry of ATOM for changing the atom type in To_Add_New_Atom_Type_Rtf().\nQuit\n");
		}
		//end	to find the line of the atom needed change
		
		fOut = fopen("mol.rtf", "w");
		for(i=0; i<nLine; i++)	{
			fprintf(fOut, "%s", szRtfTxt[i]);
			if(i==LineLastMassMol)	{	// add one line for the new atom type
				for(j=0; j<nNewAtomType; j++)	{
					fprintf(fOut, "MASS %5d %-8s %10.6lf\n", AtomeTypeCount+j+1, Mol.ChemName[Atom_To_Update[j]], Mol.mass[Atom_To_Update[j]]);
				}
			}
		}
		fclose(fOut);
	}
	

	//start	to prepare the rtf file used for drude model
	fIn = fopen("mol.rtf", "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open file: mol.rtf\nQuit\n");
	}

	nLine = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szRtfTxt[nLine], 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		else	{
			if(strncmp(szRtfTxt[nLine], "MASS 101   HTXW", 15) == 0)	{	// the entry for TIP3
//				strcpy(szRtfTxt[nLine], "END\n");
//				nLine++;
				break;
			}
			else	{
				nLine++;
				if(nLine >= MAX_LINE_RTF)	{
					fclose(fIn);
					Quit_With_Error_Msg("nLine >= MAX_LINE_RTF in Update_Rtf_File_with_New_Atom_Type().\nQuit\n");
				}
			}
		}
	}
	fclose(fIn);

	fOut = fopen("org-mol.rtf", "w");
	for(i=0; i<nLine; i++)	{
		fprintf(fOut, "%s", szRtfTxt[i]);
	}
	fclose(fOut);
	//end	to prepare the rtf file used for drude model
}

char szPrmTxt[MAX_LINE_PRM][256];
void Add_New_Atom_Type_Prm_File(void)
{
	FILE *fIn, *fOut;
	char szName[]="mol.prm", *ReadLine, szChemRead[4][16];
	int i, j, nLine, To_Output[MAX_LINE_PRM], ReadItem;
	int Idx_Last_Entry_Bond, Idx_Last_Entry_Angle, Idx_Last_Entry_Dihedral, Idx_Last_Entry_Improper, iDummy_1;
	double fDummy_1, fDummy_2; 

	if(nNewAtomType > 0)	{
		//start	to read all lines
		fIn = fopen("mol.prm", "r");
		if(fIn == NULL)	{
			Quit_With_Error_Msg("Fail to open mol.prm for read in Add_New_Atom_Type_Prm_File().\nQuit\n");
		}
		nLine=0;
		while(1)	{
			if(feof(fIn))	{
				break;
			}
			ReadLine = fgets(szPrmTxt[nLine], 256, fIn);
			if(ReadLine == NULL)	{
				break;
			}
			else	{
				To_Output[nLine] = 1;
				nLine++;
				if(nLine >= MAX_LINE_PRM)	{
					fclose(fIn);
					Quit_With_Error_Msg("nLine >= MAX_LINE_PRM in Add_New_Atom_Type_Prm_File().\nQuit\n");
				}
			}
		}
		fclose(fIn);
		//end	to read all lines
		
		//start	to determine the line index of the last entry in bond
		Idx_Last_Entry_Bond = -1;
		for(i=0; i<nLine; i++)	{
			if(strncmp(szPrmTxt[i], "ANGLES", 6)==0)	{
				Idx_Last_Entry_Angle = i;
				break;
			}
			ReadItem = sscanf(szPrmTxt[i], "%s%s%lf%lf", szChemRead[0], szChemRead[1], &fDummy_1, &fDummy_2);
			if(ReadItem == 4)	{
				Idx_Last_Entry_Bond = i;
			}
		}
		//end	to determine the line index of the last entry in bond
		
		
		//start	to determine the line index of the last entry in Angle
		for(; i<nLine; i++)	{
			if(strncmp(szPrmTxt[i], "DIHEDRALS", 9)==0)	{
				Idx_Last_Entry_Dihedral = i;
				break;
			}
			ReadItem = sscanf(szPrmTxt[i], "%s%s%s%lf%lf", szChemRead[0], szChemRead[1], szChemRead[2], &fDummy_1, &fDummy_2);
			if(ReadItem == 5)	{
				Idx_Last_Entry_Angle = i;
			}
		}
		//end	to determine the line index of the last entry in Angle
		
		
		
		//start	to determine the line index of the last entry in dihedral
		for(; i<nLine; i++)	{
			if(strncmp(szPrmTxt[i], "IMPROPERS", 9)==0)	{
				Idx_Last_Entry_Improper = i;
				break;
			}
			ReadItem = sscanf(szPrmTxt[i], "%s%s%s%s%lf%d%lf", 
				szChemRead[0], szChemRead[1], szChemRead[2], szChemRead[3], &fDummy_1, &iDummy_1, &fDummy_2);
			if(ReadItem == 7)	{
				Idx_Last_Entry_Dihedral = i;
			}
		}
		//end	to determine the line index of the last entry in dihedral
		
		//start	to determine the line index of the last entry in improper
		for(; i<nLine; i++)	{
			if(strncmp(szPrmTxt[i], "NONBONDED", 9)==0)	{
				break;
			}
			ReadItem = sscanf(szPrmTxt[i], "%s%s%s%s%lf%d%lf", 
				szChemRead[0], szChemRead[1], szChemRead[2], szChemRead[3], &fDummy_1, &iDummy_1, &fDummy_2);
			if(ReadItem == 7)	{
				Idx_Last_Entry_Improper = i;
			}
		}
		//end	to determine the line index of the last entry in improper
		
		
		
		fOut = fopen(szName, "w");
		
		for(i=0; i<nLine; i++)	{
			fprintf(fOut, "%s", szPrmTxt[i]);
			
			if(i==Idx_Last_Entry_Bond)	{
				for(j=0; j<nBond_Add; j++)	{
					if(To_Add_Bond[j])	{
						//					fprintf(fOut, "%-5s %-5s  %8.3lf %7.3lf\n", Bond_Rec_Add[j].Chem[0], Bond_Rec_Add[j].Chem[1], Bond_Rec_Add[j].para[0], Bond_Rec_Add[j].para[1]);
						fprintf(fOut, "%-8s %-8s  %8.3lf %.15lf\n", Bond_Rec_Add[j].Chem[0], Bond_Rec_Add[j].Chem[1], Bond_Rec_Add[j].para[0], Bond_Rec_Add[j].para[1]);
					}
				}
			}
			else if(i==Idx_Last_Entry_Angle)	{
				for(j=0; j<nAngle_Add; j++)	{
					if(To_Add_Angle[j])	{
						if(Angle_Rec_Add[j].para[3] == 0.0)	{
							fprintf(fOut, "%-8s %-8s %-8s  %8.3lf %.15lf\n", 
								Angle_Rec_Add[j].Chem[0], Angle_Rec_Add[j].Chem[1], Angle_Rec_Add[j].Chem[2], Angle_Rec_Add[j].para[0], Angle_Rec_Add[j].para[1]);
						}
						else	{
							fprintf(fOut, "%-8s %-8s %-8s  %8.3lf %.15lf %9.4lf %.15lf\n", 
								Angle_Rec_Add[j].Chem[0], Angle_Rec_Add[j].Chem[1], Angle_Rec_Add[j].Chem[2], Angle_Rec_Add[j].para[0], Angle_Rec_Add[j].para[1], 
								Angle_Rec_Add[j].para[2], Angle_Rec_Add[j].para[3]);
						}
					}
				}
			}
			else if(i==Idx_Last_Entry_Dihedral)	{
				for(j=0; j<nDihedral_Add; j++)	{
					if(To_Add_Dihedral[j])	{
						// to control the precision of torsion parameters
						fprintf(fOut, "%-8s %-8s %-8s %-8s  %7.3lf %3.0lf %9.1lf\n", 
							Dihedral_Rec_Add[j].Chem[0], Dihedral_Rec_Add[j].Chem[1], Dihedral_Rec_Add[j].Chem[2], Dihedral_Rec_Add[j].Chem[3], 
							Dihedral_Rec_Add[j].para[0], Dihedral_Rec_Add[j].para[1], Dihedral_Rec_Add[j].para[2]);
						//					fprintf(fOut, "%-5s %-5s %-5s %-5s  %7.3lf %3.0lf %9.1lf\n", 
						//						Dihedral_Rec_Add[j].Chem[0], Dihedral_Rec_Add[j].Chem[1], Dihedral_Rec_Add[j].Chem[2], Dihedral_Rec_Add[j].Chem[3], 
						//						Dihedral_Rec_Add[j].para[0], Dihedral_Rec_Add[j].para[1], Dihedral_Rec_Add[j].para[2]);
					}
				}
			}
			else if(i==Idx_Last_Entry_Improper)	{
				for(j=0; j<nImproper_Add; j++)	{
					if(To_Add_Improper[j])	{
						fprintf(fOut, "%-8s %-8s %-8s %-8s  %7.3lf   %3.0lf %9.1lf\n", 
							ImproDihedral_Rec_Add[j].Chem[0], ImproDihedral_Rec_Add[j].Chem[1], ImproDihedral_Rec_Add[j].Chem[2], ImproDihedral_Rec_Add[j].Chem[3], 
							ImproDihedral_Rec_Add[j].para[0], ImproDihedral_Rec_Add[j].para[2], ImproDihedral_Rec_Add[j].para[1]);
					}
				}
			}
			
		}
		fprintf(fOut, "%s", szAddTxt_LJ);
		
		fclose(fOut);
		
	}



	//start	to prepare the rtf file used for drude model
	fIn = fopen("mol.prm", "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open mol.prm for read in Add_New_Atom_Type_Prm_File().\nQuit\n");
	}
	nLine=0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szPrmTxt[nLine], 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		else if(strncmp(szPrmTxt[nLine],"HTXW   HTXW ", 12)==0)	{
		}
		else if(strncmp(szPrmTxt[nLine],"OTXW   HTXW ", 12)==0)	{
		}
		else if(strncmp(szPrmTxt[nLine],"HTXW   OTXW   HTXW", 18)==0)	{
		}
		else if(strncmp(szPrmTxt[nLine],"HTXW     0.000000", 17)==0)	{
		}
		else if(strncmp(szPrmTxt[nLine],"OTXW     0.000000", 17)==0)	{
		}
		else	{
			To_Output[nLine] = 1;
			nLine++;
			if(nLine >= MAX_LINE_PRM)	{
				fclose(fIn);
				Quit_With_Error_Msg("nLine >= MAX_LINE_PRM in Add_New_Atom_Type_Prm_File().\nQuit\n");
			}
		}
	}
	fclose(fIn);

	fOut = fopen("org-mol.prm", "w");
	for(i=0; i<nLine; i++)	{
		fprintf(fOut, "%s", szPrmTxt[i]);
	}
	fclose(fOut);
	//end	to prepare the rtf file used for drude model

}

