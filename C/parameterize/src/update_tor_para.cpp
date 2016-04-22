/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ff.h"

#define N_MAX_DIH	(20)
#define N_PARA_MAX	(N_MAX_DIH*10+1)
#define MAX_DIH_REC	(2048)
#define MAX_LEN_NAME	(16)
#define radian	(57.29577951308232088)

#define TO_ADD				(100)
#define TO_DELETE_AND_ADD	(101)

#define MAX_BOND_ADD	(8)
#define MAX_ANGLE_ADD	(64)
#define MAX_DIHEDRAL_ADD	(1536)
#define MAX_IMPROPER_ADD	(1536)
#define MAX_LJ_ADD		(2)

#define MAX_ATOM_TYPE	(256)

int nAtomType_Rtf;
char szChemName_Rtf[MAX_ATOM_TYPE][16];


BOND_REC Bond_Rec_Add[MAX_BOND_ADD];
ANGLE_REC Angle_Rec_Add[MAX_ANGLE_ADD];
DIHEDRAL_REC Dihedral_Rec_Add[MAX_DIHEDRAL_ADD];
IMPRODIHEDRAL_REC ImproDihedral_Rec_Add[MAX_IMPROPER_ADD];
char szAddTxt_LJ[256];
int nBond_Add, nAngle_Add, nDihedral_Add, nImproper_Add;
int To_Add_Bond[MAX_BOND_ADD], To_Add_Angle[MAX_ANGLE_ADD], To_Add_Dihedral[MAX_DIHEDRAL_ADD], To_Add_Improper[MAX_IMPROPER_ADD];
int nBond, nAngle, nDihedral, nImproper;

int New_Type_To_Add;	// 0 / 1 / 2
char szNewChem_1[16], szOldChem_1[16];
char szNewChem_4[16], szOldChem_4[16];

int n_Phi=0;
int DihList[N_MAX_DIH][4], IdxDihSelect[N_MAX_DIH];

char szPhiToScan[]="soft-dih-list.txt";
char szDihedralRec[]="dihedral-rec.txt";
char szRtfFile[]="mol-tor.rtf";

int n_Dihedral;
char Dih_ChemName[MAX_DIH_REC][4][MAX_LEN_NAME], Dih_ParaName[MAX_DIH_REC][4][MAX_LEN_NAME];
double Dih_Para_k[MAX_DIH_REC][7], Dih_Para_Phi[MAX_DIH_REC][7];

void Quit_With_Error_Msg(char szMsg[]);
int Read_Soft_DihedralList(void);
void Read_Dihedral_Record(void);
void BackupCurrent_Torsion_Parameters(void);
void Update_rtf_prm_xpsf(int Idx);
void Rewrite_Prm_File(int Flag, char szAddTxt[], char szChemName[][256]);
int Find_Dihedral_Entries_By_Name(char szChemName[][256], int List[]);
void Generate_New_Prm_Text(char szText[], char szChemName[][256], int Phi_Idx);
int IsValid_Updated_Prm(void);
void Determine_Bond_Angle_Dihedral_Improper_Add_1(int ActiveAtom, int Active_Phi);
void Determine_Bond_Angle_Dihedral_Improper_Add_14(int ActiveAtom_1, int ActiveAtom_4, int Active_Phi);
void To_Add_New_Atom_Type_Rtf(int ActiveAtom_1, char szOldChem_1[], int ActiveAtom_4, char szOldChem_4[]);
void Add_New_Atom_Type_Prm_File(void);
void Update_Xpsf_File_with_New_Atom_Type(int ActiveAtom_1, int ActiveAtom_4);
void To_Find_The_Org_Dihedral_Entry(char szChem_1[], char szChem_2[], char szChem_3[], char szChem_4[], double& k_Phi, double& Order, double& Phi0);
extern int Compare_Chem_Type(char Type_1[], char Type_2[]);
int Is_ChemName_Used(char szChemName[], int Exclude);
void Read_Atom_Type_RTF(void);

double rmsfit(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);
void quatfit(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);
void jacobi(int n, int np, double a[5][5], double d[5], double v[5][5], double b[5], double z[5]);
double CalRMSD(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);


CForceField ForceField;
CMol Mol;

double E_Correct;

FILE *fFile_Run_Log, *fDihedral=NULL;

#define MAX_LINE_XPSF	(2048)
void Update_Xpsf_File_with_New_Atom_Type(int ActiveAtom_1, int ActiveAtom_4)
{
	FILE *fIn, *fOut;
	int i, nLine, LineToUpdate, ReadItem, nAtom, Index, ResID, Fixed;
	char *ReadLine, szXpsfTxt[MAX_LINE_XPSF][256], szTag[256], MolName[256], ResName[256], AtomName[256], ChemName[256];
	double charge, mass, last_1, last_2;

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
				LineToUpdate = i+1+ActiveAtom_1;
				break;
			}
		}
	}
	ReadItem = sscanf(szXpsfTxt[LineToUpdate], "%d%s%d%s%s%s%lf%lf%d%lf%lf", 
		&Index, MolName, &ResID, ResName, AtomName, ChemName, &charge, &mass, &Fixed, &last_2, &last_1);
	
//	sprintf(szXpsfTxt[LineToUpdate], "%8d %-4s%2d    %-5s%-5s%-6s%8.3lf %9.4lf  %10d%10.5lf    %.6E\n", 
//		Index, MolName, ResID, ResName, AtomName, Mol.ChemName[ActiveAtom_1], charge, mass, Fixed, last_2, last_1);
	sprintf(szXpsfTxt[LineToUpdate], "%10d %-6s %3d        %-8s %-8s %-8s %12.5lf%10.5lf           0%10.5lf %14.5lf\n", 
		Index, MolName, ResID, ResName, AtomName, Mol.ChemName[ActiveAtom_1], charge, mass, last_2, last_1);

	if(ActiveAtom_4 >= 0)	{
		LineToUpdate = LineToUpdate - ActiveAtom_1 + ActiveAtom_4;
		ReadItem = sscanf(szXpsfTxt[LineToUpdate], "%d %s %d %s%s%s%lf%lf%d%lf%lf", 
			&Index, MolName, &ResID, ResName, AtomName, ChemName, &charge, &mass, &Fixed, &last_2, &last_1);
		
//		sprintf(szXpsfTxt[LineToUpdate], "%8d %-4s%2d    %-5s%-5s%-6s%8.3lf %9.4lf  %10d%10.5lf    %.6E\n", 
//			Index, MolName, ResID, ResName, AtomName, Mol.ChemName[ActiveAtom_4], charge, mass, Fixed, last_2, last_1);
		sprintf(szXpsfTxt[LineToUpdate], "%10d %-6s %3d        %-8s %-8s %-8s %12.5lf%10.5lf           0%10.5lf %14.5lf\n", 
			Index, MolName, ResID, ResName, AtomName, Mol.ChemName[ActiveAtom_4], charge, mass, last_2, last_1);
	}

	fOut = fopen("mol.xpsf", "w");
	for(i=0; i<nLine; i++)	{
		fprintf(fOut, "%s", szXpsfTxt[i]);
	}
	fclose(fOut);

}
#undef MAX_LINE_XPSF

int main(void)
{
  timebomb();

	int i;

	fFile_Run_Log = stdout;
	printf("Writing [%s]\n", szDihedralRec );

	fDihedral = fopen(szDihedralRec, "w");
	

	Read_Atom_Type_RTF();
	ForceField.ReadForceField("mol.prm");	// read the force field. Parameter file. It is fine if the file includes rtf info. 
	
	Mol.ReadPSF("mol.xpsf", 0);	//the first parameter is the xpsf file name; the second parameter is 0 (only one molecule in xpsf) or 1 (two molecules in xpsf). 
	n_Dihedral = Mol.nDihedral;
	if(n_Dihedral > MAX_DIH_REC)	{
		Quit_With_Error_Msg("n_Dihedral > MAX_DIH_REC.\nQuit\n");
	}

	printf( "calling -- %p\n", fDihedral );

	Mol.AssignForceFieldParameters(&ForceField, fDihedral );
	fclose(fDihedral);

	Read_Soft_DihedralList();
	BackupCurrent_Torsion_Parameters();	// to save the desired parameters used for validating new rtf, prm and xpsf files


	printf("There are nPhi=%d \n", n_Phi );

	for(i=0; i<n_Phi; i++)	{	// already read in previous modified prm / xpsf file
		printf("Writing to [%s]\n", szDihedralRec );
		fDihedral = fopen(szDihedralRec, "w");
		ForceField.ReadForceField("mol.prm");	// read the force field. Parameter file. It is fine if the file includes rtf info. 
		Mol.AssignForceFieldParameters(&ForceField, fDihedral );
		fclose(fDihedral);

		Read_Dihedral_Record();
		Update_rtf_prm_xpsf(i);
	}

	IsValid_Updated_Prm();

	return 0;
}

void Generate_A_New_ChemName(int iAtom, char szOldChemName[])
{
	int i, PoolSize;
	char CharPool[]="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ", CharEnd, NewChemName[256], BaseName[256], szErrorMsg[256];

	PoolSize = strlen(CharPool);
	strcpy(szOldChemName, Mol.ChemName[iAtom]);

	strcpy(BaseName, szOldChemName);
	BaseName[4]=0;	// only keep maximum 4 lettters in the original name

	for(i=0; i<PoolSize; i++)	{
		CharEnd = CharPool[i];
		sprintf(NewChemName, "%sx%c", BaseName, CharEnd);
		if(Is_ChemName_Used(NewChemName, iAtom) == 0)	{	// not used. Done.
			strcpy(Mol.ChemName[iAtom], NewChemName);
			return;
		}
		else	{	// used. Keep trying
		}
	}

	sprintf(szErrorMsg, "Fail to create a new name for atom %d, %s\nQuit\n", iAtom+1, szOldChemName);
	Quit_With_Error_Msg(szErrorMsg);
}

void Update_rtf_prm_xpsf(int Idx)
{
	int Idx_Phi, iPos, iPos_Check;
	int i, j, k, ia, ib, ic, id, Idx_Check, n_ChemName_Hit_List[MAX_DIH_REC], I_abcd_Count[4], n_Chem_Hit, ActiveAtom_1, ActiveAtom_4, CountMin;
	char ChemName[4][256], szTxtAdd[1024], szChemPara[2][12], szBuff[256];
	double Para_List[6];

	Idx_Phi = IdxDihSelect[Idx];	// the dihedral with updated parameters
	iPos = Idx_Phi*4;
	ia = Mol.DihedralList[iPos  ];
	ib = Mol.DihedralList[iPos+1];
	ic = Mol.DihedralList[iPos+2];
	id = Mol.DihedralList[iPos+3];

	strcpy(ChemName[0], Mol.ChemName[ia]);
	strcpy(ChemName[1], Mol.ChemName[ib]);
	strcpy(ChemName[2], Mol.ChemName[ic]);
	strcpy(ChemName[3], Mol.ChemName[id]);

	n_Chem_Hit = Find_Dihedral_Entries_By_Name(ChemName, n_ChemName_Hit_List);
	if(n_Chem_Hit == 1)	{	// easy case. 
		Generate_New_Prm_Text(szTxtAdd, ChemName, Idx_Phi);

		if( (strcmp(Dih_ParaName[Idx_Phi][0], "X")==0) || (strcmp(Dih_ParaName[Idx_Phi][3], "X")==0) )	{	// contained wild character? Add the new entries. 
			Rewrite_Prm_File(TO_ADD, szTxtAdd, ChemName);
		}
		else	{	// just delete the old entries and add the new one
			Rewrite_Prm_File(TO_DELETE_AND_ADD, szTxtAdd, ChemName);
		}
	}
	else	{	// multi-entries. New atom type has to be added
		//start	to determine the active atoms
		memset(I_abcd_Count, 0, sizeof(int)*4);
		for(i=0; i<n_Chem_Hit; i++)	{
			Idx_Check = n_ChemName_Hit_List[i];
			iPos_Check = Idx_Check*4;
			for(j=0; j<4; j++)	{
				for(k=0; k<4; k++)	{
					if( Mol.DihedralList[iPos_Check+j] == Mol.DihedralList[iPos+k])	{
						I_abcd_Count[k]++;
					}
				}
			}
		}
		CountMin = I_abcd_Count[0];
		ActiveAtom_1 = 0;
		for(k=1; k<4; k++)	{
			if(I_abcd_Count[k] < CountMin)	{
				CountMin = I_abcd_Count[k];
				ActiveAtom_1 = k;
			}
		}

		if(CountMin > 1)	{	// to change 1-4 two atoms. Mol.DihedralList[iPos+0] and Mol.DihedralList[iPos+3]
			ActiveAtom_1 = Mol.DihedralList[iPos+0];
			ActiveAtom_4 = Mol.DihedralList[iPos+3];
//			strcpy(szOldChem_1, Mol.ChemName[ActiveAtom_1]);
//			strcpy(szOldChem_4, Mol.ChemName[ActiveAtom_4]);

			Generate_A_New_ChemName(ActiveAtom_1, szOldChem_1);
			strcpy(szNewChem_1, Mol.ChemName[ActiveAtom_1]);	// update the chem name for the active atom
//			sprintf(szNewChem_1, "%sT%d", szOldChem_1, Idx);
//			strcpy(Mol.ChemName[ActiveAtom_1], szNewChem_1);	// update the chem name for the active atom

			Generate_A_New_ChemName(ActiveAtom_4, szOldChem_4);
			strcpy(szNewChem_4, Mol.ChemName[ActiveAtom_4]);	// update the chem name for the active atom
//			sprintf(szNewChem_4, "%sT%d", szOldChem_4, Idx);
//			strcpy(Mol.ChemName[ActiveAtom_4], szNewChem_4);	// update the chem name for the active atom
			
			strcpy(szChemPara[0], szOldChem_1);
			ForceField.GetPara_LJ(szChemPara, Para_List);

			if( (fabs(Para_List[4]-Para_List[1]) > 1.0E-4) || (fabs(Para_List[5]-Para_List[2]) > 1.0E-4) )	{	// to output 1-4 parameters
				sprintf(szAddTxt_LJ, "%-8s   0.00 %9.4lf %9.4lf      0.00 %9.4lf %9.4lf\n", 
					szNewChem_1, Para_List[1], Para_List[2], Para_List[4], Para_List[5]);
			}
			else	{	// no 1-4 parameters
				sprintf(szAddTxt_LJ, "%-8s   0.00 %9.4lf %9.4lf\n", 
					szNewChem_1, Para_List[1], Para_List[2]);
			}

			strcpy(szChemPara[0], szOldChem_4);
			ForceField.GetPara_LJ(szChemPara, Para_List);

			if( (fabs(Para_List[4]-Para_List[1]) > 1.0E-4) || (fabs(Para_List[5]-Para_List[2]) > 1.0E-4) )	{	// to output 1-4 parameters
				sprintf(szBuff, "%-8s   0.00 %9.4lf %9.4lf      0.00 %9.4lf %9.4lf\n", 
					szNewChem_4, Para_List[1], Para_List[2], Para_List[4], Para_List[5]);
			}
			else	{	// no 1-4 parameters
				sprintf(szBuff, "%-8s   0.00 %9.4lf %9.4lf\n", 
					szNewChem_4, Para_List[1], Para_List[2]);
			}
			strcat(szAddTxt_LJ, szBuff);
			
			
			nBond_Add = nAngle_Add = nDihedral_Add = nImproper_Add = 0;
			Determine_Bond_Angle_Dihedral_Improper_Add_14(ActiveAtom_1, ActiveAtom_4, Idx_Phi);
			
			Add_New_Atom_Type_Prm_File();
			
			To_Add_New_Atom_Type_Rtf(ActiveAtom_1, szOldChem_1, ActiveAtom_4, szOldChem_4);
			Update_Xpsf_File_with_New_Atom_Type(ActiveAtom_1, ActiveAtom_4);
		}
		else	{
			ActiveAtom_1 = Mol.DihedralList[iPos+ActiveAtom_1];
//			strcpy(szOldChem_1, Mol.ChemName[ActiveAtom_1]);

			Generate_A_New_ChemName(ActiveAtom_1, szOldChem_1);
			strcpy(szNewChem_1, Mol.ChemName[ActiveAtom_1]);	// update the chem name for the active atom
			
			strcpy(szChemPara[0], szOldChem_1);
			ForceField.GetPara_LJ(szChemPara, Para_List);

//			sprintf(szNewChem_1, "%sT%d", szOldChem_1, Idx);
//			strcpy(Mol.ChemName[ActiveAtom_1], szNewChem_1);	// update the chem name for the active atom

			if( (fabs(Para_List[4]-Para_List[1]) > 1.0E-4) || (fabs(Para_List[5]-Para_List[2]) > 1.0E-4) )	{	// to output 1-4 parameters
				sprintf(szAddTxt_LJ, "%-8s   0.00 %9.4lf %9.4lf      0.00 %9.4lf %9.4lf\n", 
					szNewChem_1, Para_List[1], Para_List[2], Para_List[4], Para_List[5]);
			}
			else	{	// no 1-4 parameters
				sprintf(szAddTxt_LJ, "%-8s   0.00 %9.4lf %9.4lf\n", 
					szNewChem_1, Para_List[1], Para_List[2]);
			}
			
			nBond_Add = nAngle_Add = nDihedral_Add = nImproper_Add = 0;
			Determine_Bond_Angle_Dihedral_Improper_Add_1(ActiveAtom_1, Idx_Phi);
			
			Add_New_Atom_Type_Prm_File();
			
			To_Add_New_Atom_Type_Rtf(ActiveAtom_1, szOldChem_1, -1, NULL);
			Update_Xpsf_File_with_New_Atom_Type(ActiveAtom_1, -1);
		}

		//end	to determine the active atoms
	}
}

#define MAX_LINE_PRM	(2048)
void Rewrite_Prm_File(int Flag, char szAddTxt[], char szChemName[][256])
{
	FILE *fIn, *fOut;
	char szPrmTxt[MAX_LINE_PRM][256], szName[]="mol.prm", *ReadLine, szChemRead[4][256];
	int i, nLine, To_Output[MAX_LINE_PRM], Idx_Line_Dihedral, Idx_Last_Entry_Dihedral, Idx_Line_Delete, ReadItem, Order;
	double k_Dih, Phi0; 

	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open mol.prm for read in Rewrite_Prm_File().\nQuit\n");
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
				Quit_With_Error_Msg("nLine >= MAX_LINE_PRM in Rewrite_Prm_File().\nQuit\n");
			}
		}
	}
	fclose(fIn);

	Idx_Line_Dihedral = -1;
	for(i=0; i<nLine; i++)	{
		if(strncmp(szPrmTxt[i], "DIHEDRALS", 9)==0)	{
			Idx_Line_Dihedral = i;
			break;
		}
	}
	if(Idx_Line_Dihedral < 0)	{
		Quit_With_Error_Msg("Fail to find the tag DIHEDRALS in mol.prm.\nQuit\n");
	}

	Idx_Line_Delete = -1;	// to record the last line to be deleted
	Idx_Last_Entry_Dihedral = -1;
	for(; i<nLine; i++)	{
		if(strncmp(szPrmTxt[i], "IMPROPER", 8)==0)	{
			break;
		}
		else	{
			ReadItem = sscanf(szPrmTxt[i], "%s%s%s%s%lf%d%lf", szChemRead[0], szChemRead[1], szChemRead[2], szChemRead[3], &k_Dih, &Order, &Phi0);
			if(ReadItem == 7)	{
				Idx_Last_Entry_Dihedral = i;
				if(Flag == TO_DELETE_AND_ADD)	{	// to delete existing entries
					if( (strcmp(szChemName[0], szChemRead[0])==0) && (strcmp(szChemName[1], szChemRead[1])==0) && (strcmp(szChemName[2], szChemRead[2])==0) && (strcmp(szChemName[3], szChemRead[3])==0) )	{
						To_Output[i] = 0;	// to be deleted
						Idx_Line_Delete = i;
					}
					else if( (strcmp(szChemName[0], szChemRead[3])==0) && (strcmp(szChemName[1], szChemRead[2])==0) && (strcmp(szChemName[2], szChemRead[1])==0) && (strcmp(szChemName[3], szChemRead[0])==0) )	{
						To_Output[i] = 0;	// to be deleted
						Idx_Line_Delete = i;
					}
				}
			}
		}
	}
	if(Idx_Last_Entry_Dihedral < 0)	{
		Quit_With_Error_Msg("Fail to find the last entry of DIHEDRALS in mol.prm.\nQuit\n");
	}


	fOut = fopen(szName, "w");

	for(i=0; i<nLine; i++)	{
		if(To_Output[i])	{
			fprintf(fOut, "%s", szPrmTxt[i]);
		}
		if( (Flag==TO_ADD) && (i==Idx_Last_Entry_Dihedral) )	{
			fprintf(fOut, "%s", szAddTxt);
		}
		if( (Flag==TO_DELETE_AND_ADD) && (i==Idx_Line_Delete))	{
			fprintf(fOut, "%s", szAddTxt);
		}
	}

	fclose(fOut);
}

void To_Find_The_Org_Dihedral_Entry(char szChem_1[], char szChem_2[], char szChem_3[], char szChem_4[], double& k_Phi, double& Order, double& Phi0)
{
	FILE *fIn;
	char szPrmTxt[MAX_LINE_PRM][256], szName[]="mol.prm", *ReadLine, szChemRead[4][256], ErrorMsg[256];
	int i, nLine, Line_Dihedral_Start, Line_Dihedral_End, ReadItem;

	//start	to read all lines
	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open mol.prm for read in Rewrite_Prm_File().\nQuit\n");
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
			nLine++;
			if(nLine >= MAX_LINE_PRM)	{
				fclose(fIn);
				Quit_With_Error_Msg("nLine >= MAX_LINE_PRM in To_Find_The_Org_Dihedral_Entry().\nQuit\n");
			}
		}
	}
	fclose(fIn);
	//end	to read all lines

	Line_Dihedral_Start = -1;
	for(i=0; i<nLine; i++)	{
		if(strncmp(szPrmTxt[i], "DIHEDRALS", 9)==0)	{
			Line_Dihedral_Start = i+1;
			break;
		}
	}
	if(Line_Dihedral_Start < 0)	{
		Quit_With_Error_Msg("Line_Dihedral_Start < 0 in To_Find_The_Org_Dihedral_Entry().\nQuit\n");
	}

	Line_Dihedral_End = -1;
	for(; i<nLine; i++)	{
		if(strncmp(szPrmTxt[i], "IMPROPERS", 9)==0)	{
			Line_Dihedral_End = i-1;
			break;
		}
	}
	if(Line_Dihedral_End < 0)	{
		Quit_With_Error_Msg("Line_Dihedral_End < 0 in To_Find_The_Org_Dihedral_Entry().\nQuit\n");
	}

	for(i=Line_Dihedral_Start; i<=Line_Dihedral_End; i++)	{
		ReadItem = sscanf(szPrmTxt[i], "%s%s%s%s%lf%lf%lf", szChemRead[0], szChemRead[1], szChemRead[2], szChemRead[3], &k_Phi, &Order, &Phi0);
		if(ReadItem == 7)	{
			if( (strcmp(szChemRead[0], szChem_1)==0) && (strcmp(szChemRead[1], szChem_2)==0) && (strcmp(szChemRead[2], szChem_3)==0) && (strcmp(szChemRead[3], szChem_4)==0) )	{
				return;
			}
			if( (strcmp(szChemRead[3], szChem_1)==0) && (strcmp(szChemRead[2], szChem_2)==0) && (strcmp(szChemRead[1], szChem_3)==0) && (strcmp(szChemRead[0], szChem_4)==0) )	{
				return;
			}
		}
	}

	for(i=Line_Dihedral_Start; i<=Line_Dihedral_End; i++)	{
		ReadItem = sscanf(szPrmTxt[i], "%s%s%s%s%lf%lf%lf", szChemRead[0], szChemRead[1], szChemRead[2], szChemRead[3], &k_Phi, &Order, &Phi0);
		if(ReadItem == 7)	{
			if( (Compare_Chem_Type(szChemRead[0], szChem_1)==0) && (Compare_Chem_Type(szChemRead[1], szChem_2)==0) && (Compare_Chem_Type(szChemRead[2], szChem_3)==0) && (Compare_Chem_Type(szChemRead[3], szChem_4)==0) )	{
				return;
			}
			if( (Compare_Chem_Type(szChemRead[3], szChem_1)==0) && (Compare_Chem_Type(szChemRead[2], szChem_2)==0) && (Compare_Chem_Type(szChemRead[1], szChem_3)==0) && (Compare_Chem_Type(szChemRead[0], szChem_4)==0) )	{
				return;
			}
		}
	}

	sprintf(ErrorMsg, "Fail to find the dihedral parameters for, \n%-8s %-8s %-8s %-8s\nQuit\n", szChem_1, szChem_2, szChem_3, szChem_4);
	Quit_With_Error_Msg(ErrorMsg);
	return;
}

void Add_New_Atom_Type_Prm_File(void)
{
	FILE *fIn, *fOut;
	char szPrmTxt[MAX_LINE_PRM][256], szName[]="mol.prm", *ReadLine, szChemRead[4][256];
	int i, j, nLine, To_Output[MAX_LINE_PRM], ReadItem;
	int Idx_Last_Entry_Bond, Idx_Last_Entry_Angle, Idx_Last_Entry_Dihedral, Idx_Last_Entry_Improper, Idx_Last_Entry_LJ=-1, iDummy_1;
	double fDummy_1, fDummy_2; 

	//start	to read all lines
	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open mol.prm for read in Rewrite_Prm_File().\nQuit\n");
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
		if(strncmp(szPrmTxt[i], "IMPROPER", 8)==0)	{
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

/*
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
*/
	
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

	//start	to determine the line index of the last entry in LJ
	for(; i<nLine; i++)	{
		if(strncmp(szPrmTxt[i], "NBFIX", 5)==0)	{
			Idx_Last_Entry_LJ = i - 1;
			break;
		}
		if(strncmp(szPrmTxt[i], "HBOND", 5)==0)	{
			Idx_Last_Entry_LJ = i - 1;
			break;
		}
		if(strncmp(szPrmTxt[i], "END", 3)==0)	{
			Idx_Last_Entry_LJ = i - 1;
			break;
		}
	}
	if(Idx_Last_Entry_LJ < 0)	{
		Idx_Last_Entry_LJ = nLine-1;
	}
	//end	to determine the line index of the last entry in LJ
	


	printf("Writing [%s]\n", szName );
	fOut = fopen(szName, "w");

	for(i=0; i<nLine; i++)	{
		fprintf(fOut, "%s", szPrmTxt[i]);
		
		if(i==Idx_Last_Entry_Bond)	{
			for(j=0; j<nBond_Add; j++)	{
				if(To_Add_Bond[j])	{
					fprintf(fOut, "%-8s %-8s  %8.3lf %7.3lf\n", Bond_Rec_Add[j].Chem[0], Bond_Rec_Add[j].Chem[1], Bond_Rec_Add[j].para[0], Bond_Rec_Add[j].para[1]);
				}
			}
		}
		else if(i==Idx_Last_Entry_Angle)	{
			for(j=0; j<nAngle_Add; j++)	{
				if(To_Add_Angle[j])	{
					if(Angle_Rec_Add[j].para[3] == 0.0)	{
						fprintf(fOut, "%-8s %-8s %-8s  %8.3lf %12.3lf\n", 
							Angle_Rec_Add[j].Chem[0], Angle_Rec_Add[j].Chem[1], Angle_Rec_Add[j].Chem[2], Angle_Rec_Add[j].para[0], Angle_Rec_Add[j].para[1]);
					}
					else	{
						fprintf(fOut, "%-8s %-8s %-8s  %8.3lf %12.3lf %9.3lf %10.5lf\n", 
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
					fprintf(fOut, "%-8s %-8s %-8s %-8s  %.15lf %3.0lf %9.1lf\n", 
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
		else if(i==Idx_Last_Entry_LJ)	{
			fprintf(fOut, "%s", szAddTxt_LJ);
		}

	}

	fclose(fOut);
}


#undef MAX_LINE_PRM

void BackupCurrent_Torsion_Parameters(void)
{
	int i, j, ReadItem, Idx_Phi;
	FILE *fIn;
	double Para_Read[10];

	for(i=0; i<n_Dihedral; i++)	{
		for(j=1; j<=6; j++)	{
			Dih_Para_k[i][j] = Mol.Para_k_Dih[i][j];
			Dih_Para_Phi[i][j] = Mol.Para_phi[i][j];
		}
	}

	fIn = fopen("saved-para.dat", "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Error in opening file: saved-para.dat in BackupCurrent_Torsion_Parameters()\nQuit\n");
	}

	for(i=0; i<n_Phi; i++)	{
		Idx_Phi = IdxDihSelect[i];
		for(j=0; j<10; j++)	{
			ReadItem = fscanf(fIn, "%lf", &(Para_Read[j]));
			if(ReadItem != 1)	{
				fclose(fIn);
				Quit_With_Error_Msg("Error in reading saved-para.dat in BackupCurrent_Torsion_Parameters().\nQuit\n");
			}
		}
		Dih_Para_k[Idx_Phi][1] = Para_Read[0];
		Dih_Para_Phi[Idx_Phi][1] = Para_Read[1];
		Dih_Para_k[Idx_Phi][2] = Para_Read[2];
		Dih_Para_Phi[Idx_Phi][2] = Para_Read[3];
		Dih_Para_k[Idx_Phi][3] = Para_Read[4];
		Dih_Para_Phi[Idx_Phi][3] = Para_Read[5];
		Dih_Para_k[Idx_Phi][4] = Para_Read[6];
		Dih_Para_Phi[Idx_Phi][4] = Para_Read[7];
		Dih_Para_k[Idx_Phi][6] = Para_Read[8];
		Dih_Para_Phi[Idx_Phi][6] = Para_Read[9];

		for(j=1; j<=6; j++)	{	// to assign the torsion parameters
			Mol.Para_k_Dih[Idx_Phi][j] = Dih_Para_k[Idx_Phi][j];
			Mol.Para_phi[Idx_Phi][j] = Dih_Para_Phi[Idx_Phi][j];
		}
	}

	fclose(fIn);

	Mol.ReadXYZ("mol-opt.xyz");
	printf("----------------------------------------\nUsing the original prm/xpsf files\n");
	Mol.Cal_E(1);
	printf("----------------------------------------\n");

	E_Correct = Mol.E_Total;

	//start	to test only
	FILE *fOut;
	fOut = fopen("before.dat", "w");
	for(i=0; i<n_Dihedral; i++)	{
		for(j=1; j<=6; j++)	{
			fprintf(fOut, "%4d %2d %.15lf %.15lf\n", i, j, Mol.Para_k_Dih[i][j], Mol.Para_phi[i][j]);
		}
	}
	fclose(fOut);
	//end	to test only
}

void Read_Dihedral_Record(void)
{
	int i, ReadItem, iTmp;

	fDihedral = fopen(szDihedralRec, "r");
	if(fDihedral == NULL)	{
		Quit_With_Error_Msg("Error in open file dihedral-rec.txt for reading in Read_Dihedral_Record().\nQuit\n");
	}

	for(i=0; i<n_Dihedral; i++)	{
		ReadItem = fscanf(fDihedral, "%d%s%s%s%s%s%s%s%s", 
//			&iTmp, &(Dih_ChemName[i][0]), &(Dih_ChemName[i][1]), &(Dih_ChemName[i][2]), &(Dih_ChemName[i][3]), 
//			&(Dih_ParaName[i][0]), &(Dih_ParaName[i][1]), &(Dih_ParaName[i][2]), &(Dih_ParaName[i][3]));
			&iTmp, Dih_ChemName[i][0], Dih_ChemName[i][1], Dih_ChemName[i][2], Dih_ChemName[i][3], 
			Dih_ParaName[i][0], Dih_ParaName[i][1], Dih_ParaName[i][2], Dih_ParaName[i][3]);
		if(ReadItem != 9)	{
			fclose(fDihedral);
			Quit_With_Error_Msg("Error in reading file: dihedral-rec.txt \nQuit\n");
		}
	}
	fclose(fDihedral);
}

void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in update-torsion-para.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}


int Read_Soft_DihedralList(void)
{
	FILE *fIn;
	int ReadItem;
	char szLine[256], *ReadLine, ErrorMsg[256];

	n_Phi = 0;

	fIn = fopen(szPhiToScan, "r");

	while(1)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Fail to open file: %s\nQuit\n", szPhiToScan);
			Quit_With_Error_Msg(ErrorMsg);
		}
		ReadLine = fgets(szLine, 128, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		ReadItem = sscanf(szLine, "%d %d %d %d", 
			&(DihList[n_Phi][0]), &(DihList[n_Phi][1]), &(DihList[n_Phi][2]), &(DihList[n_Phi][3]));

		if(ReadItem == 4)	{
			DihList[n_Phi][0]--;
			DihList[n_Phi][1]--;
			DihList[n_Phi][2]--;
			DihList[n_Phi][3]--;
			IdxDihSelect[n_Phi] = Mol.Query_Dihedral_Index(DihList[n_Phi][0], DihList[n_Phi][1], DihList[n_Phi][2], DihList[n_Phi][3]);

			if(IdxDihSelect[n_Phi] < 0)	{
				Quit_With_Error_Msg("Fail to identify the index of one soft dihedral.\n");
			}
			n_Phi++;
		}
		else	{
			break;
		}
	}

	fclose(fIn);

	return n_Phi;
}

int IsValid_Updated_Prm(void)
{
	int Idx, i, iPos, ia, ib, ic, id, nAtom;
	double rmsd, x_QM[MAX_ATOM], y_QM[MAX_ATOM], z_QM[MAX_ATOM];

	fDihedral = NULL;
	ForceField.ReadForceField("mol.prm");	// read the force field. Parameter file. It is fine if the file includes rtf info. 
	Mol.ReadPSF("mol.xpsf", 0);	//the first parameter is the xpsf file name; the second parameter is 0 (only one molecule in xpsf) or 1 (two molecules in xpsf). 
	Mol.AssignForceFieldParameters(&ForceField, fDihedral );

	for(Idx=0; Idx<n_Dihedral; Idx++)	{
		iPos = Idx*4;
		ia = Mol.DihedralList[iPos  ];
		ib = Mol.DihedralList[iPos+1];
		ic = Mol.DihedralList[iPos+2];
		id = Mol.DihedralList[iPos+3];

		for(i=1; i<=6; i++)	{
			if( fabs(Dih_Para_k[Idx][i] - Mol.Para_k_Dih[Idx][i]) > 1.0E-3 )	{	// something wrong?
				printf("Please check the torsion parameters for %3d dihedral: %-8s %-8s %-8s %-8s\n", Idx+1, Mol.ChemName[ia], Mol.ChemName[ib], Mol.ChemName[ic], Mol.ChemName[id]);
				Quit_With_Error_Msg("Error in generate new prm/xpsf file in IsValid_Updated_Prm().\nQuit\n");
				return 0;
			}
		}
	}

	printf("----------------------------------------\nUsing the new prm/xpsf files\n");
	Mol.Cal_E(1);
	printf("----------------------------------------\n");

	if(fabs(E_Correct - Mol.E_Total) > 1.0E-4)	{	// something wrong!!!
    printf(" WARNING: Delta between inital and modified ffs is > 1.0e-4\n" );
    printf(" Correct    : %f\n", E_Correct );
    printf(" Calculated : %f\n", Mol.E_Total );

//		Quit_With_Error_Msg("Energies are different. Error in generate new prm/xpsf file in IsValid_Updated_Prm().\nQuit\n");
	}


	//start	to test only
	FILE *fOut;
	int j;
	fOut = fopen("after.dat", "w");
	for(i=0; i<n_Dihedral; i++)	{
		for(j=1; j<=6; j++)	{
			fprintf(fOut, "%4d %2d %.15lf %.15lf\n", i, j, Mol.Para_k_Dih[i][j], Mol.Para_phi[i][j]);
		}
	}
	fclose(fOut);
	//end	to test only

	nAtom = Mol.nAtom;
	memcpy(x_QM, Mol.x, sizeof(double)*nAtom);
	memcpy(y_QM, Mol.y, sizeof(double)*nAtom);
	memcpy(z_QM, Mol.z, sizeof(double)*nAtom);

	Mol.WriteXYZ("mol-qm.xyz");

	Mol.FullGeometryOptimization_SD(3000,0.005);
	Mol.FullGeometryOptimization_SD(3000,0.002);
	Mol.FullGeometryOptimization_SD(3000,0.001);
	Mol.FullGeometryOptimization_SD(3000,0.0005);
	Mol.FullGeometryOptimization_SD(3000,0.0002);
	Mol.FullGeometryOptimization_SD(3000,0.0001);
	Mol.FullGeometryOptimization_SD(3000,0.00005);
	Mol.FullGeometryOptimization_SD(3000,0.00002);
	Mol.FullGeometryOptimization_SD(3000,0.00001);
	Mol.FullGeometryOptimization_SD(3000,0.000005);
	Mol.FullGeometryOptimization_SD(3000,0.000002);
	Mol.FullGeometryOptimization_SD(3000,0.000001);
//	Mol.FullGeometryOptimization_LBFGS(1);	// use small step size first
//	Mol.FullGeometryOptimization_LBFGS(0);	// use larger step size

	rmsd = CalRMSD(x_QM, y_QM, z_QM, Mol.x, Mol.y, Mol.z, nAtom);
	Mol.WriteXYZ("mol-mm.xyz");

	printf("Saved mol-mm.pdb and mol-qm.pdb\n" );
	printf("mol-qm.pdb is mol-opt.xyz in PDB form\n" );
	printf("mol-mm.pdb is a geomopt starting from that conf \n" );

	printf("RMSD between QM and MM: %.3lf\n", rmsd);

	return 1;
}

int Find_Dihedral_Entries_By_Name(char szChemName[][256], int List[])
{
	int nEntries=0, i;
	char ErrorMsg[256];

	for(i=0; i<n_Dihedral; i++)	{
		if( ((strcmp(Dih_ChemName[i][0], szChemName[0])==0) && (strcmp(Dih_ChemName[i][1], szChemName[1])==0) && (strcmp(Dih_ChemName[i][2], szChemName[2])==0) && (strcmp(Dih_ChemName[i][3], szChemName[3])==0)) || 
			((strcmp(Dih_ChemName[i][0], szChemName[3])==0) && (strcmp(Dih_ChemName[i][1], szChemName[2])==0) && (strcmp(Dih_ChemName[i][2], szChemName[1])==0) && (strcmp(Dih_ChemName[i][3], szChemName[0])==0)) )	{
			List[nEntries] = i;
			nEntries++;
		}
	}

	if(nEntries == 0)	{
		sprintf(ErrorMsg, "Fail to find dihedral, %s %s %s %s \nQuit\n", szChemName[0], szChemName[1], szChemName[2], szChemName[3]);
		Quit_With_Error_Msg(ErrorMsg);
	}

	return nEntries;
}

void Generate_New_Prm_Text(char szText[], char szChemName[][256], int Phi_Idx)
{
	int i;
	char szLine[256];

	strcpy(szText, "");
	for(i=1; i<=6; i++)	{
		if(fabs(Dih_Para_k[Phi_Idx][i]) > 1.0E-6)	{	// non-zero
			// to control the precision of torsion parameters
			sprintf(szLine, "%-8s %-8s %-8s %-8s    %.15lf         %d    %6.1lf\n", 
				szChemName[0], szChemName[1], szChemName[2], szChemName[3], Dih_Para_k[Phi_Idx][i], i, Dih_Para_Phi[Phi_Idx][i]*radian);
			strcat(szText, szLine);
		}
	}
}


#define MAX_LINE_RTF	(2048)
void To_Add_New_Atom_Type_Rtf(int ActiveAtom_1, char szOldChem_1[], int ActiveAtom_4, char szOldChem_4[])
{
	FILE *fIn, *fOut;
	char *ReadLine, szLineSave[256];
	char szRtfTxt[MAX_LINE_RTF][256];
	int i, nLine, LineLastMassMol, AtomeTypeCount, IdxAtom_1, IdxAtom_4, AtomCount, Idx_Line_To_Change;

	fIn = fopen(szRtfFile, "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open file: mol-tor.rtf\nQuit\n");
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
				Quit_With_Error_Msg("nLine >= MAX_LINE_RTF in To_Add_New_Atom_Type_Rtf().\nQuit\n");
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
		Quit_With_Error_Msg("Fail to find an entry for MASS in mol-tor.rtf.\nQuit\n");
	}
	//end to find the index of the last line of the entry of MASS (atom type) for MOL

	//start	to find the line of the atom needed change
	IdxAtom_1 = ActiveAtom_1 + 1;
	IdxAtom_4 = ActiveAtom_4 + 1;
	AtomCount = 0;
	Idx_Line_To_Change = -1;
	for(; i<nLine; i++)	{
		if(strncmp(szRtfTxt[i], "BOND", 4)==0)	{	// the end of the definition of atoms in MOL
			break;
		}
    printf("MJH parse line [%s]\n", szRtfTxt[i] );
		if(strncmp(szRtfTxt[i], "ATOM", 4)==0)	{

      char old_atom_name[20];
			old_atom_name[0]='\0';
      strcpy( old_atom_name, "XXX" );

      if(strlen( szRtfTxt[i] ) > 5 ) {
				sscanf( szRtfTxt[i]+5, "%19s", old_atom_name );
      }
			AtomCount++;

			if(AtomCount == IdxAtom_1)	{	// the atom to be changed
				strcpy(szLineSave, szRtfTxt[i]);
			//	sprintf(szRtfTxt[i], "ATOM %-6s %-8s %7.3lf\n", Mol.AtomName[ActiveAtom_1], Mol.ChemName[ActiveAtom_1], Mol.CG[ActiveAtom_1]);
				sprintf(szRtfTxt[i], "ATOM %-6s %-8s %7.3lf\n", old_atom_name, Mol.ChemName[ActiveAtom_1], Mol.CG[ActiveAtom_1]);
				printf("1: Replacing the atom type of atom %d from %6s to %6s\n", ActiveAtom_1+1, szOldChem_1, Mol.ChemName[ActiveAtom_1]);
				printf("[%s] to [%s] old atom name [%s]\n", szLineSave, szRtfTxt[i], old_atom_name );
				Idx_Line_To_Change = i;
//				break;
			}
			else if(AtomCount == IdxAtom_4)	{	// the atom to be changed
				strcpy(szLineSave, szRtfTxt[i]);
			//	sprintf(szRtfTxt[i], "ATOM %-6s %-8s %7.3lf\n", Mol.AtomName[ActiveAtom_4], Mol.ChemName[ActiveAtom_4], Mol.CG[ActiveAtom_4]);
				sprintf(szRtfTxt[i], "ATOM %-6s %-8s %7.3lf\n", old_atom_name, Mol.ChemName[ActiveAtom_4], Mol.CG[ActiveAtom_4]);
				printf("2: Replacing the atom type of atom %d from %6s to %6s\n", ActiveAtom_4+1, szOldChem_4, Mol.ChemName[ActiveAtom_4]);
				printf("[%s] to [%s] old atom name [%s]\n", szLineSave, szRtfTxt[i], old_atom_name );
				Idx_Line_To_Change = i;
//				break;
			}
		}
	}
	if(Idx_Line_To_Change < 0)	{
		Quit_With_Error_Msg("Fail to indentify the entry of ATOM for changing the atom type in To_Add_New_Atom_Type_Rtf().\nQuit\n");
	}
	//end	to find the line of the atom needed change

	fOut = fopen(szRtfFile, "w");
	for(i=0; i<nLine; i++)	{
		fprintf(fOut, "%s", szRtfTxt[i]);
		if(i==LineLastMassMol)	{	// add one line for the new atom type
			fprintf(fOut, "MASS %5d %-8s %10.6lf\n", AtomeTypeCount+1, Mol.ChemName[ActiveAtom_1], Mol.mass[ActiveAtom_1]);
			if(ActiveAtom_4 >= 0)	{
				fprintf(fOut, "MASS %5d %-8s %10.6lf\n", AtomeTypeCount+2, Mol.ChemName[ActiveAtom_4], Mol.mass[ActiveAtom_4]);
			}
		}
	}
	fclose(fOut);
}
#undef MAX_LINE_RTF

void Determine_Bond_Angle_Dihedral_Improper_Add_1(int ActiveAtom, int Active_Phi)
{
	int i, j, iPos, ia, ib, ic, id, LocalCount;
	double k_Phi, Order, Phi0;

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
				if( (fabs(Angle_Rec_Add[i].para[0]-Angle_Rec_Add[j].para[0]) < 1.0E-7) && (fabs(Angle_Rec_Add[i].para[1]-Angle_Rec_Add[j].para[1]) < 1.0E-7) )	{
					To_Add_Angle[j] = 0;
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

		if(i == Active_Phi)	{	// to add all parameters for this soft dihedral
			LocalCount = 0;
			for(j=1; j<=6; j++)	{
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[0], Mol.ChemName[ia]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[1], Mol.ChemName[ib]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[2], Mol.ChemName[ic]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[3], Mol.ChemName[id]);
				
				Dihedral_Rec_Add[nDihedral_Add].para[0] = Dih_Para_k[i][j];
				Dihedral_Rec_Add[nDihedral_Add].para[1] = 1.0*j;
				Dihedral_Rec_Add[nDihedral_Add].para[2] = Dih_Para_Phi[i][j]*radian;
				if(fabs(Dihedral_Rec_Add[nDihedral_Add].para[0]) > 1.0E-6)	{	// a valid entry
					To_Add_Dihedral[nDihedral_Add] = 1;
					nDihedral_Add++;
					LocalCount++;
				}
			}
			if(LocalCount == 0)	{	// chem names are set already
				j = 3;
				Dihedral_Rec_Add[nDihedral_Add].para[0] = 0.0;
				Dihedral_Rec_Add[nDihedral_Add].para[1] = 1.0*j;
				Dihedral_Rec_Add[nDihedral_Add].para[2] = Dih_Para_Phi[i][j]*radian;
				To_Add_Dihedral[nDihedral_Add] = 1;
				nDihedral_Add++;
			}
		}
		else if( (ActiveAtom == ia) || (ActiveAtom == ib) || (ActiveAtom == ic)  || (ActiveAtom == id) )	{
			LocalCount = 0;
			for(j=1; j<=6; j++)	{
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[0], Mol.ChemName[ia]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[1], Mol.ChemName[ib]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[2], Mol.ChemName[ic]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[3], Mol.ChemName[id]);
				
				Dihedral_Rec_Add[nDihedral_Add].para[0] = Mol.Para_k_Dih[i][j];
				Dihedral_Rec_Add[nDihedral_Add].para[1] = 1.0*j;
				Dihedral_Rec_Add[nDihedral_Add].para[2] = Mol.Para_phi[i][j]*radian;
				if(fabs(Dihedral_Rec_Add[nDihedral_Add].para[0]) > 1.0E-6)	{	// a valid entry
					To_Add_Dihedral[nDihedral_Add] = 1;
					nDihedral_Add++;
					LocalCount++;
				}
			}
			if(LocalCount == 0)	{	// We need to add a dummy entry with k=0. We need to find the entry in rtf file
				if(ActiveAtom == ia)	{
					To_Find_The_Org_Dihedral_Entry(szOldChem_1, Mol.ChemName[ib], Mol.ChemName[ic], Mol.ChemName[id], k_Phi, Order, Phi0);
				}
				else if(ActiveAtom == ib)	{
					To_Find_The_Org_Dihedral_Entry(Mol.ChemName[ia], szOldChem_1, Mol.ChemName[ic], Mol.ChemName[id], k_Phi, Order, Phi0);
				}
				else if(ActiveAtom == ic)	{
					To_Find_The_Org_Dihedral_Entry(Mol.ChemName[ia], Mol.ChemName[ib], szOldChem_1, Mol.ChemName[id], k_Phi, Order, Phi0);
				}
				else if(ActiveAtom == id)	{
					To_Find_The_Org_Dihedral_Entry(Mol.ChemName[ia], Mol.ChemName[ib], Mol.ChemName[ic], szOldChem_1, k_Phi, Order, Phi0);
				}

				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[0], Mol.ChemName[ia]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[1], Mol.ChemName[ib]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[2], Mol.ChemName[ic]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[3], Mol.ChemName[id]);
				
				Dihedral_Rec_Add[nDihedral_Add].para[0] = k_Phi;
				Dihedral_Rec_Add[nDihedral_Add].para[1] = Order;
				Dihedral_Rec_Add[nDihedral_Add].para[2] = Phi0;
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
			strcpy(ImproDihedral_Rec_Add[nImproper_Add].Chem[0], Mol.ChemName[ia]);
			strcpy(ImproDihedral_Rec_Add[nImproper_Add].Chem[1], Mol.ChemName[ib]);
			strcpy(ImproDihedral_Rec_Add[nImproper_Add].Chem[2], Mol.ChemName[ic]);
			strcpy(ImproDihedral_Rec_Add[nImproper_Add].Chem[3], Mol.ChemName[id]);
			
			ImproDihedral_Rec_Add[nImproper_Add].para[0] = Mol.Para_k_ImpDih[i];
			ImproDihedral_Rec_Add[nImproper_Add].para[1] = Mol.Para_Imp_phi[i]*radian;
			ImproDihedral_Rec_Add[nImproper_Add].para[2] = Mol.Para_Type_ImpDih[i];
			if(Mol.Para_k_ImpDih[i] > 1.0E-3)	{	// a valid entry
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

void Determine_Bond_Angle_Dihedral_Improper_Add_14(int ActiveAtom_1, int ActiveAtom_4, int Active_Phi)
{
	int i, j, iPos, ia, ib, ic, id, LocalCount;
	double k_Phi, Order, Phi0;
	char ErrorMsg[256];

	nBond = Mol.nBond;
	nAngle = Mol.nAngle;
	nDihedral = Mol.nDihedral;
	nImproper = Mol.nImpro;

	//start	to generate the list of new bond entries
	for(i=0; i<nBond; i++)	{
		iPos = 2*i;
		ia = Mol.BondList[iPos  ];
		ib = Mol.BondList[iPos+1];

		if( (ActiveAtom_1 == ia) || (ActiveAtom_1 == ib) || (ActiveAtom_4 == ia) || (ActiveAtom_4 == ib) )	{
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

		if( (ActiveAtom_1 == ia) || (ActiveAtom_1 == ib) || (ActiveAtom_1 == ic) || (ActiveAtom_4 == ia) || (ActiveAtom_4 == ib) || (ActiveAtom_4 == ic))	{
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
				if( (fabs(Angle_Rec_Add[i].para[0]-Angle_Rec_Add[j].para[0]) < 1.0E-7) && (fabs(Angle_Rec_Add[i].para[1]-Angle_Rec_Add[j].para[1]) < 1.0E-7) )	{
					To_Add_Angle[j] = 0;
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

		if(i == Active_Phi)	{	// to add all parameters for this soft dihedral. We don't delete the original parameter
			LocalCount = 0;
			for(j=1; j<=6; j++)	{
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[0], Mol.ChemName[ia]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[1], Mol.ChemName[ib]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[2], Mol.ChemName[ic]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[3], Mol.ChemName[id]);
				
				Dihedral_Rec_Add[nDihedral_Add].para[0] = Dih_Para_k[i][j];
				Dihedral_Rec_Add[nDihedral_Add].para[1] = 1.0*j;
				Dihedral_Rec_Add[nDihedral_Add].para[2] = Dih_Para_Phi[i][j]*radian;
				if(fabs(Dihedral_Rec_Add[nDihedral_Add].para[0]) > 1.0E-6)	{	// a valid entry
					To_Add_Dihedral[nDihedral_Add] = 1;
					nDihedral_Add++;
					LocalCount++;
				}
			}
			if(LocalCount == 0)	{	// chem names are set already
				j = 3;
				Dihedral_Rec_Add[nDihedral_Add].para[0] = 0.0;
				Dihedral_Rec_Add[nDihedral_Add].para[1] = 1.0*j;
				Dihedral_Rec_Add[nDihedral_Add].para[2] = Dih_Para_Phi[i][j]*radian;
				To_Add_Dihedral[nDihedral_Add] = 1;
				nDihedral_Add++;
			}
		}
		else if( (ActiveAtom_1 == ia) || (ActiveAtom_1 == ib) || (ActiveAtom_1 == ic)  || (ActiveAtom_1 == id) || (ActiveAtom_4 == ia) || (ActiveAtom_4 == ib) || (ActiveAtom_4 == ic)  || (ActiveAtom_4 == id) )	{
			LocalCount = 0;
			for(j=1; j<=6; j++)	{
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[0], Mol.ChemName[ia]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[1], Mol.ChemName[ib]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[2], Mol.ChemName[ic]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[3], Mol.ChemName[id]);
				
				Dihedral_Rec_Add[nDihedral_Add].para[0] = Mol.Para_k_Dih[i][j];
				Dihedral_Rec_Add[nDihedral_Add].para[1] = 1.0*j;
				Dihedral_Rec_Add[nDihedral_Add].para[2] = Mol.Para_phi[i][j]*radian;
				if(fabs(Dihedral_Rec_Add[nDihedral_Add].para[0]) > 1.0E-6)	{	// a valid entry
					To_Add_Dihedral[nDihedral_Add] = 1;
					nDihedral_Add++;
					LocalCount++;
				}
			}
			if(LocalCount == 0)	{	// We need to add a dummy entry with k=0. We need to find the entry in rtf file
				if( ActiveAtom_1 == ia )	{
					To_Find_The_Org_Dihedral_Entry(szOldChem_1, Mol.ChemName[ib], Mol.ChemName[ic], Mol.ChemName[id], k_Phi, Order, Phi0);
				}
				else if( ActiveAtom_1 == ib )	{
					To_Find_The_Org_Dihedral_Entry(Mol.ChemName[ia], szOldChem_1, Mol.ChemName[ic], Mol.ChemName[id], k_Phi, Order, Phi0);
				}
				else if( ActiveAtom_1 == ic )	{
					To_Find_The_Org_Dihedral_Entry(Mol.ChemName[ia], Mol.ChemName[ib], szOldChem_1, Mol.ChemName[id], k_Phi, Order, Phi0);
				}
				else if( ActiveAtom_1 == id )	{
					To_Find_The_Org_Dihedral_Entry(Mol.ChemName[ia], Mol.ChemName[ib], Mol.ChemName[ic], szOldChem_1, k_Phi, Order, Phi0);
				}
				else if(ActiveAtom_4 == ia)	{
					To_Find_The_Org_Dihedral_Entry(szOldChem_4, Mol.ChemName[ib], Mol.ChemName[ic], Mol.ChemName[id], k_Phi, Order, Phi0);
				}
				else if(ActiveAtom_4 == ib)	{
					To_Find_The_Org_Dihedral_Entry(Mol.ChemName[ia], szOldChem_4, Mol.ChemName[ic], Mol.ChemName[id], k_Phi, Order, Phi0);
				}
				else if(ActiveAtom_4 == ic)	{
					To_Find_The_Org_Dihedral_Entry(Mol.ChemName[ia], Mol.ChemName[ib], szOldChem_4, Mol.ChemName[id], k_Phi, Order, Phi0);
				}
				else if(ActiveAtom_4 == id)	{
					To_Find_The_Org_Dihedral_Entry(Mol.ChemName[ia], Mol.ChemName[ib], Mol.ChemName[ic], szOldChem_4, k_Phi, Order, Phi0);
				}
				else	{
					sprintf(ErrorMsg, "Unexpected situation in Determine_Bond_Angle_Dihedral_Improper_Add_14().\nQuit\n");
					Quit_With_Error_Msg(ErrorMsg);
				}

				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[0], Mol.ChemName[ia]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[1], Mol.ChemName[ib]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[2], Mol.ChemName[ic]);
				strcpy(Dihedral_Rec_Add[nDihedral_Add].Chem[3], Mol.ChemName[id]);
				
				Dihedral_Rec_Add[nDihedral_Add].para[0] = k_Phi;
				Dihedral_Rec_Add[nDihedral_Add].para[1] = Order;
				Dihedral_Rec_Add[nDihedral_Add].para[2] = Phi0;
				To_Add_Dihedral[nDihedral_Add] = 1;
				nDihedral_Add++;
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
		
		if( (ActiveAtom_1 == ia) || (ActiveAtom_1 == ib) || (ActiveAtom_1 == ic)  || (ActiveAtom_1 == id) || (ActiveAtom_4 == ia) || (ActiveAtom_4 == ib) || (ActiveAtom_4 == ic)  || (ActiveAtom_4 == id) )	{
			strcpy(ImproDihedral_Rec_Add[nImproper_Add].Chem[0], Mol.ChemName[ia]);
			strcpy(ImproDihedral_Rec_Add[nImproper_Add].Chem[1], Mol.ChemName[ib]);
			strcpy(ImproDihedral_Rec_Add[nImproper_Add].Chem[2], Mol.ChemName[ic]);
			strcpy(ImproDihedral_Rec_Add[nImproper_Add].Chem[3], Mol.ChemName[id]);
			
			ImproDihedral_Rec_Add[nImproper_Add].para[0] = Mol.Para_k_ImpDih[i];
			ImproDihedral_Rec_Add[nImproper_Add].para[1] = Mol.Para_Imp_phi[i]*radian;
			ImproDihedral_Rec_Add[nImproper_Add].para[2] = Mol.Para_Type_ImpDih[i];
			if(Mol.Para_k_ImpDih[i] > 1.0E-3)	{	// a valid entry
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

	for(i=0; i<nAtomType_Rtf; i++)	{
		if(strcmp(szChemName, szChemName_Rtf[i])==0)	{
			return 1;
		}
	}

	return 0;
}


void Read_Atom_Type_RTF(void)
{
	FILE *fIn;
	char szLine[256], *ReadLine, szTmp[256], szChemName[256];
	int ChemID, ReadItem;
	double mass=0.0;

	fIn = fopen(szRtfFile, "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open file: mol-tor.rtf\nQuit\n");
	}

	nAtomType_Rtf = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		if( (strncmp(szLine,"MASS",4)==0) || (strncmp(szLine,"mass",4)==0) )	{
			ReadItem = sscanf(szLine, "%s%d%s%lf", szTmp, &ChemID, szChemName, &mass);
			if(ReadItem == 4)	{
				strcpy(szChemName_Rtf[nAtomType_Rtf], szChemName);
				nAtomType_Rtf++;
			}
		}
	}

	fclose(fIn);
}



double CalRMSD(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum)
{
	double RMSD, midxa, midya, midza, midxb, midyb, midzb;
	double xa_tmp[MAX_ATOM], ya_tmp[MAX_ATOM], za_tmp[MAX_ATOM];
	double xb_tmp[MAX_ATOM], yb_tmp[MAX_ATOM], zb_tmp[MAX_ATOM];
	int i;

	memcpy(xa_tmp+1, xa, sizeof(double)*AtomNum);
	memcpy(ya_tmp+1, ya, sizeof(double)*AtomNum);
	memcpy(za_tmp+1, za, sizeof(double)*AtomNum);

	memcpy(xb_tmp+1, xb, sizeof(double)*AtomNum);
	memcpy(yb_tmp+1, yb, sizeof(double)*AtomNum);
	memcpy(zb_tmp+1, zb, sizeof(double)*AtomNum);

	midxa=midya=midza=midxb=midyb=midzb=0.0;
	for(i=1; i<=AtomNum; i++)	{
		midxa += (xa_tmp[i]);
		midya += (ya_tmp[i]);
		midza += (za_tmp[i]);

		midxb += (xb_tmp[i]);
		midyb += (yb_tmp[i]);
		midzb += (zb_tmp[i]);
	}
	midxa/=AtomNum;
	midya/=AtomNum;
	midza/=AtomNum;
	midxb/=AtomNum;
	midyb/=AtomNum;
	midzb/=AtomNum;

	for(i=1; i<=AtomNum; i++)	{
		xa_tmp[i]-=midxa;
		ya_tmp[i]-=midya;
		za_tmp[i]-=midza;
		xb_tmp[i]-=midxb;
		yb_tmp[i]-=midyb;
		zb_tmp[i]-=midzb;
	}

	quatfit(xa_tmp, ya_tmp, za_tmp, xb_tmp, yb_tmp, zb_tmp, AtomNum);
	RMSD=rmsfit(xa_tmp, ya_tmp, za_tmp, xb_tmp, yb_tmp, zb_tmp, AtomNum);

	for(i=1; i<=AtomNum; i++)	{	// generate the rotated coodinates
		xb[i-1]=midxb+xb_tmp[i];
		yb[i-1]=midyb+yb_tmp[i];
		zb[i-1]=midzb+zb_tmp[i];
	}


	return RMSD;
}

double rmsfit (double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n)
{
	int ia;
	double rmsfit, dist2, dx, dy, dz;

	rmsfit=0.0;
	for(ia=1; ia<=n; ia++)	{	//store xi-x0 into xStd[]
		dx = x1[ia] - x2[ia];
		dy = y1[ia] - y2[ia];
		dz = z1[ia] - z2[ia];
		dist2 = dx*dx + dy*dy + dz*dz;
		rmsfit += dist2;
	}
	rmsfit = sqrt(rmsfit/n);

	return rmsfit;
}

void quatfit (double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n)
{
	double xxyx,xxyy,xxyz, xyyx,xyyy,xyyz, xzyx,xzyy,xzyz;
	double xrot,yrot,zrot, work1[5],work2[5], q[5],d[5], rot[4][4], c[5][5], v[5][5];
	int ia;

	xxyx = 0.0;
	xxyy = 0.0;
	xxyz = 0.0;
	xyyx = 0.0;
	xyyy = 0.0;
	xyyz = 0.0;
	xzyx = 0.0;
	xzyy = 0.0;
	xzyz = 0.0;

	for(ia=1; ia<=n; ia++)	{
		xxyx = xxyx + x1[ia]*x2[ia];
		xxyy = xxyy + y1[ia]*x2[ia];
		xxyz = xxyz + z1[ia]*x2[ia];
		xyyx = xyyx + x1[ia]*y2[ia];
		xyyy = xyyy + y1[ia]*y2[ia];
		xyyz = xyyz + z1[ia]*y2[ia];
		xzyx = xzyx + x1[ia]*z2[ia];
		xzyy = xzyy + y1[ia]*z2[ia];
		xzyz = xzyz + z1[ia]*z2[ia];
	}

	c[1][1] = xxyx + xyyy + xzyz;
	c[1][2] = xzyy - xyyz;
	c[2][2] = xxyx - xyyy - xzyz;
	c[1][3] = xxyz - xzyx;
	c[2][3] = xxyy + xyyx;
	c[3][3] = xyyy - xzyz - xxyx;
	c[1][4] = xyyx - xxyy;
	c[2][4] = xzyx + xxyz;
	c[3][4] = xyyz + xzyy;
	c[4][4] = xzyz - xxyx - xyyy;
	
	jacobi(4,4,c,d,v,work1,work2);
	q[1] = v[1][4];
	q[2] = v[2][4];
	q[3] = v[3][4];
	q[4] = v[4][4];

	rot[1][1] = q[1]*q[1] + q[2]*q[2] - q[3]*q[3] - q[4]*q[4];
	rot[2][1] = 2.0 * (q[2] * q[3] - q[1] * q[4]);
	rot[3][1] = 2.0 * (q[2] * q[4] + q[1] * q[3]);
	rot[1][2] = 2.0 * (q[3] * q[2] + q[1] * q[4]);
	rot[2][2] = q[1]*q[1] - q[2]*q[2] + q[3]*q[3] - q[4]*q[4];
	rot[3][2] = 2.0 * (q[3] * q[4] - q[1] * q[2]);
	rot[1][3] = 2.0 * (q[4] * q[2] - q[1] * q[3]);
	rot[2][3] = 2.0 * (q[4] * q[3] + q[1] * q[2]);
	rot[3][3] = q[1]*q[1] - q[2]*q[2] - q[3]*q[3] + q[4]*q[4];

	for(ia=1; ia<=n; ia++)	{
		xrot = x2[ia]*rot[1][1] + y2[ia]*rot[1][2] + z2[ia]*rot[1][3];
		yrot = x2[ia]*rot[2][1] + y2[ia]*rot[2][2] + z2[ia]*rot[2][3];
		zrot = x2[ia]*rot[3][1] + y2[ia]*rot[3][2] + z2[ia]*rot[3][3];
		x2[ia] = xrot;
		y2[ia] = yrot;
		z2[ia] = zrot;
	}

	return;
}


#define ZERO			(1.0E-100)
void jacobi (int n, int np, double a[5][5], double d[5], double v[5][5], double b[5], double z[5])
{
	int i,j,k;
	int ip,iq;
	int nrot,maxrot;
	double sm,tresh,s,c,t;
	double theta,tau,h,g,p;
	
	maxrot = 100;
	nrot = 0;
	
	for(ip=1; ip<=n; ip++)	{
		for(iq=1; iq<=n; iq++)	{
			v[ip][iq] = 0.0;
		}
		v[ip][ip] = 1.0;
	}
	for(ip=1; ip<=n; ip++)	{
		b[ip] = a[ip][ip];
		d[ip] = b[ip];
		z[ip] = 0.0;
	}
	
	
	for(i=1; i<=maxrot; i++)	{
		sm=0.0;
		for(ip=1; ip<=(n-1); ip++)	{
			for(iq=ip+1; iq<=n; iq++)	{
				sm = sm + fabs(a[ip][iq]);
			}
		}
		
		if(sm < ZERO)	{	//converge
			break;
		}
		if(i < 4)	{
			tresh = 0.2*sm / (n*n);
		}
		else	{
			tresh=0.0;
		}
		
		for(ip=1; ip<=(n-1); ip++)	{
			for(iq=ip+1; iq<=n; iq++)	{
				g = 100.0 * fabs(a[ip][iq]);
				if((i > 4) && (fabs(g) < ZERO))	{
					a[ip][iq] = 0.0;
				}
				else if(fabs(a[ip][iq]) > tresh)	{
					h = d[iq] - d[ip];
					if(fabs(g) < ZERO)	{
						t = a[ip][iq] / h;
					}
					else	{
						theta = 0.5*h / a[ip][iq];
						t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0)  t = -t;
					}
					
					c = 1.0 / sqrt(1.0+t*t);
					s = t * c;
					tau = s / (1.0+c);
					h = t * a[ip][iq];
					z[ip] = z[ip] - h;
					z[iq] = z[iq] + h;
					d[ip] = d[ip] - h;
					d[iq] = d[iq] + h;
					a[ip][iq] = 0.0;
					
					for(j=1; j<=(ip-1); j++)	{
						g = a[j][ip];
						h = a[j][iq];
						a[j][ip] = g - s*(h+g*tau);
						a[j][iq] = h + s*(g-h*tau);
					}
					for(j=ip+1; j<=(iq-1); j++)	{
						g = a[ip][j];
						h = a[j][iq];
						a[ip][j] = g - s*(h+g*tau);
						a[j][iq] = h + s*(g-h*tau);
					}
					for(j=iq+1; j<=n; j++)	{
						g = a[ip][j];
						h = a[iq][j];
						a[ip][j] = g - s*(h+g*tau);
						a[iq][j] = h + s*(g-h*tau);
					}
					for(j=1; j<=n; j++)	{
						g = v[j][ip];
						h = v[j][iq];
						v[j][ip] = g - s*(h+g*tau);
						v[j][iq] = h + s*(g-h*tau);
					}
					nrot = nrot + 1;
				}
			}
		}
		for(ip=1; ip<=n; ip++)	{
            b[ip] = b[ip] + z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
		}
	}
	if(nrot >= maxrot)	{
		printf("JACOBI  --  Matrix Diagonalization not Converged.\n");
	}

	for(i=1; i<= (n-1); i++)	{
         k = i;
         p = d[i];
		 for(j=i+1; j<=n; j++)	{
			 if(d[j] < p)	{
               k = j;
               p = d[j];
			 }
		 }
		 if(k != i)	{
            d[k] = d[i];
            d[i] = p;
			for(j=1; j<=n; j++)	{
               p = v[j][i];
               v[j][i] = v[j][k];
               v[j][k] = p;
			}
		 }
	}
}
#undef	ZERO

