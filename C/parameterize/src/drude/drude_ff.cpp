/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ff.h"

#define MAX_LINE	(2048)
#define MAX_LEN		(256)
#define MAX_N_ATOM	(256)
#define THOLE_0		(1.3)	// the target value for thole

char szTxtPrm[MAX_LINE][MAX_LEN];
char szTxtRtf[MAX_LINE][MAX_LEN];
char szTxtRtf_Bond[MAX_LINE*MAX_LEN/2]="";
char szLonePairInfo[65536]="";
char szLonePair_Bond_Para[MAX_N_ATOM][16];
char szRtfFile[MAX_LEN];
char szPrmFile[MAX_LEN];

void Add_Tip3_Rtf_File(void);
void Add_Tip3_Prm_File(void);
void To_Upper_Case(char szBuff[]);
int FindString(char szBuff[], const char szTag[]);
void Reorganize_Lonepair_Bond_Para(void);

int n_LonePair_Bond_Para=0, LonePair_Bond_Para_Output[MAX_N_ATOM];
int Line_Insert_Bond, Line_Insert_Angle, Line_Insert_Impropers;
int Line_Insert_Atom_RTF=-1, LastLine_Mass=-1;
int Found_Improper;
char szXpsfName[]="mol.xpsf";
char szFileElemList[]="elem-list.txt";

#define N_ELEM	(109)
char Name_Elem[N_ELEM][16]={"H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "K", "AR", "CA", 
		"SC", "TI", "V", "CR", "MN", "FE", "NI", "CO", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", 
		"RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "I", "TE", "XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", 
		"HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", 
		"PA", "TH", "NP", "U", "AM", "PU", "BK", "CM", "CF", "ES", "FM", "MD", "NO", "RF", "LR", "DB", "BH", "SG", "MT", "HS"};
double Mass_List[N_ELEM]={1.0079, 4.0026, 6.941, 9.0122, 10.811, 12.0107, 14.0067, 15.9994, 18.9984, 20.1797, 22.9897, 24.305, 26.9815, 28.0855, 
		30.9738, 32.065, 35.453, 39.0983, 39.948, 40.078, 44.9559, 47.867, 50.9415, 51.9961, 54.938, 55.845, 58.6934, 58.9332, 63.546, 65.39, 
		69.723, 72.64, 74.9216, 78.96, 79.904, 83.8, 85.4678, 87.62, 88.9059, 91.224, 92.9064, 95.94, 98, 101.07, 102.9055, 106.42, 107.8682, 
		112.411, 114.818, 118.71, 121.76, 126.9045, 127.6, 131.293, 132.9055, 137.327, 138.9055, 140.116, 140.9077, 144.24, 145, 150.36, 151.964, 
		157.25, 158.9253, 162.5, 164.9303, 167.259, 168.9342, 173.04, 174.967, 178.49, 180.9479, 183.84, 186.207, 190.23, 192.217, 195.078, 196.9665, 
		200.59, 204.3833, 207.2, 208.9804, 209, 210, 222, 223, 226, 227, 231.0359, 232.0381, 237, 238.0289, 243, 244, 247, 247, 251, 252, 257, 258, 
		259, 261, 262, 262, 264, 266, 268, 277};

void GetAtomTypeFromMass(void);
int DetermineAtomType(char szName[], double mass);
void Count_Bond_and_H(void);
int Count_Bonded_SP2_Atoms(int Idx, int IsRestricted);
void Assign_Alpha_Miller(void);
int QueryAtomIndex(char szAtmName[]);
void Add_Lonepair_Anisotropy(FILE *fOut);
int CountHeavyAtoms(int Idx);
int Get_Next_Neighbor(int Idx);

CMol Mol;
int AtomType[MAX_N_ATOM];
int H_Count[MAX_N_ATOM], SP2_O_Count[MAX_N_ATOM], BondCount[MAX_N_ATOM], BondList[MAX_N_ATOM][6], SP2_Flag[MAX_N_ATOM];
double Alpha_Miller[MAX_N_ATOM];

double E14FAC=1.0;

FILE *fFile_Run_Log;
void Quit_With_Error_Msg(char szMsg[]);


int main(int argc, char *argv[])
{

  timebomb();

	if(argc < 4)	{
		printf("Usage: add-tip3 your-rtf your-prm your-psf E14FAC\nQuit\n");
		return 1;
	}

	strcpy(szRtfFile, argv[1]);
	strcpy(szPrmFile, argv[2]);
	strcpy(szXpsfName, argv[3]);

	if(argc == 5)	{
		E14FAC = atof(argv[4]);
	}

	fFile_Run_Log = fopen("log-modify-rtf.txt", "w");
	Mol.ReadPSF(szXpsfName, 0);
	GetAtomTypeFromMass();
	Count_Bond_and_H();
	Assign_Alpha_Miller();

	Add_Tip3_Rtf_File();

	Add_Tip3_Prm_File();

	fclose(fFile_Run_Log);

	return 0;
}

void Reorganize_Lonepair_Bond_Para(void)
{
	int i, j;

	for(i=0; i<n_LonePair_Bond_Para; i++)	{
		LonePair_Bond_Para_Output[i] = 1;
	}

	for(i=0; i<n_LonePair_Bond_Para; i++)	{
		if(LonePair_Bond_Para_Output[i] != 1)	{
			continue;
		}
		for(j=i+1; j<n_LonePair_Bond_Para; j++)	{
			if(strcmp(szLonePair_Bond_Para[j], szLonePair_Bond_Para[i])==0)	{	// delete the entry
				LonePair_Bond_Para_Output[j] = 0;
			}
		}
	}
}

int CountHeavyAtoms(int Idx)	// for the atoms only connected with two other atoms
{
	int Count = 0;

	if(Mol.mass[BondList[Idx][0]] > 1.5)	{
		Count++;
	}
	if(Mol.mass[BondList[Idx][1]] > 1.5)	{
		Count++;
	}
	return Count;
}

void Count_Bond_and_H(void)
{
	int i, nBond, Atom_1, Atom_2, nAtom;
	char *AtomName;

	memset(BondCount, 0, sizeof(int)*MAX_N_ATOM);
	memset(H_Count, 0, sizeof(int)*MAX_N_ATOM);
	memset(SP2_O_Count, 0, sizeof(int)*MAX_N_ATOM);
	memset(BondList, 0, sizeof(int)*MAX_N_ATOM*4);

	nBond = Mol.nBond;
	nAtom = Mol.nAtom;

	for(i=0; i<nBond; i++)	{
		Atom_1 = Mol.BondList[2*i];
		Atom_2 = Mol.BondList[2*i+1];

		BondList[Atom_1][BondCount[Atom_1]] = Atom_2;
		BondCount[Atom_1]++;
		BondList[Atom_2][BondCount[Atom_2]] = Atom_1;
		BondCount[Atom_2]++;

//		if( (BondCount[Atom_1] > 4) || (BondCount[Atom_2] > 4) )	{
//			Quit_With_Error_Msg("The number of bonded atoms is larger than 4.\nQuit\n");
//		}

		if(Mol.mass[Atom_1] < 1.2)	{	// H atom
			H_Count[Atom_2]++;
		}
		if(Mol.mass[Atom_2] < 1.2)	{	// H atom
			H_Count[Atom_1]++;
		}
	}

	memset(SP2_Flag, 0, sizeof(int)*MAX_N_ATOM);

	for(i=0; i<nAtom; i++)	{
		AtomName = Name_Elem[AtomType[i]];

		if(strcmp(AtomName, "C")==0)	{	// Carbon
			if(BondCount[i] == 3)	{
				SP2_Flag[i] = 1;
			}
		}
		else if(strcmp(AtomName, "N")==0)	{
			if(BondCount[i] == 2)	{
				SP2_Flag[i] = 1;
			}
		}
		else if(strcmp(AtomName, "O")==0)	{
			if(BondCount[i] == 1)	{
				SP2_Flag[i] = 1;
				SP2_O_Count[BondList[i][0]] ++;
			}
		}
		else if(strcmp(AtomName, "S")==0)	{
			if(BondCount[i] == 1)	{
				SP2_Flag[i] = 1;
			}
		}
	}

	for(i=0; i<nAtom; i++)	{
		AtomName = Name_Elem[AtomType[i]];

		if(strcmp(AtomName, "N")==0)	{
			if(BondCount[i] == 3)	{
				if(Count_Bonded_SP2_Atoms(i, 1) >= 1)	{	// directly connected with at least one sp2 atom
					SP2_Flag[i] = 2;	// promoted sp2 atom
				}
			}
		}
		else if(strcmp(AtomName, "O")==0)	{
			if(BondCount[i] == 2)	{
				if(Count_Bonded_SP2_Atoms(i, 1) == 2)	{	// directly connected with two sp2 atoms
					SP2_Flag[i] = 2;	// promoted sp2 atom
				}
			}
		}
		else if(strcmp(AtomName, "S")==0)	{
			if(BondCount[i] == 2)	{
				if(Count_Bonded_SP2_Atoms(i, 1) == 2)	{	// directly connected with two sp2 atoms
					SP2_Flag[i] = 2;	// promoted sp2 atom
				}
			}
		}
	}

}

int Get_Next_Neighbor(int Idx)
{
	int Neighbor, i;

	Neighbor = BondList[Idx][0];

	for(i=0; i<BondCount[Neighbor]; i++)	{
		if(BondList[Neighbor][i] != Idx)	{
			return BondList[Neighbor][i];
		}
	}

	Quit_With_Error_Msg("Fail to find a valid 2nd neighbor atom.\nQuit\n");
	return -1;
}

#define CG_LP	(-0.15)
void Add_Lonepair_Anisotropy(FILE *fOut)
{
	int i, nAtom, nBondedHeavyAtom, Atom_L, Atom_R, Next_Neighbor, PrintOut;
	char *AtomName, szLPName_1[256], szLPName_2[256], szTxtAdd[256];

	strcpy(szLonePairInfo, "");

	nAtom = Mol.nAtom;
	for(i=0; i<nAtom; i++)	{
		AtomName = Name_Elem[AtomType[i]];
		PrintOut = 0;

		if(strcmp(AtomName, "N")==0)	{
			if(BondCount[i] == 2)	{	// only sp2 N and connected with two sp2 atoms
				if(Count_Bonded_SP2_Atoms(i, 1) == 2)	{	// the lonepair will be added. "bisector"
					Atom_L = BondList[i][0];
					Atom_R = BondList[i][1];
					sprintf(szLPName_1, "LP%s", Mol.AtomName[i]);

					sprintf(szTxtAdd, "LONEPAIR bisector %-7s%-7s%-7s%-7s distance 0.35 angle 179.99 dihe 179.99\n", 
						szLPName_1, Mol.AtomName[i], Mol.AtomName[Atom_L], Mol.AtomName[Atom_R]);
					strcat(szLonePairInfo, szTxtAdd);
					sprintf(szTxtAdd, "ANISOTROPY        %-7s%-7s%-7s%-7s A11 1.1611   A22 0.6778\n",
						Mol.AtomName[i], szLPName_1, Mol.AtomName[Atom_L], Mol.AtomName[Atom_R]);
					strcat(szLonePairInfo, szTxtAdd);

					sprintf(szTxtAdd, "BOND  %7s %7s\n", Mol.AtomName[i], szLPName_1);
					strcat(szTxtRtf_Bond, szTxtAdd);

					fprintf(fOut, "ATOM %-7s%-7s%7.3lf  ALPHA%7.3lf  THOLE%6.3lf\n", Mol.AtomName[i], Mol.ChemName[i], Mol.CG[i]-CG_LP, -Alpha_Miller[i], THOLE_0);
					fprintf(fOut, "ATOM %-7s%-7s%7.3lf\n", szLPName_1, "LP", CG_LP);

					strcpy(szLonePair_Bond_Para[n_LonePair_Bond_Para], Mol.ChemName[i]);
					n_LonePair_Bond_Para++;

					PrintOut = 1;
				}
			}
		}
		else if(strcmp(AtomName, "S")==0)	{
			if(BondCount[i] == 2)	{	// only sp3 S
				Atom_L = BondList[i][0];
				Atom_R = BondList[i][1];
				sprintf(szLPName_1, "LPA%s", Mol.AtomName[i]);
				sprintf(szLPName_2, "LPB%s", Mol.AtomName[i]);
				nBondedHeavyAtom = CountHeavyAtoms(i);
				
				if(nBondedHeavyAtom==1)	{	// "relative"
					if(Mol.mass[Atom_L] < 1.5)	{	// Hydrogen
						Atom_R = BondList[i][0];
						Atom_L = BondList[i][1];
					}
					sprintf(szTxtAdd, "LONEPAIR relative %-7s%-7s%-7s%-7s distance 0.75 angle  95.00 dihe 100.00\n", 
						szLPName_1, Mol.AtomName[i], Mol.AtomName[Atom_L], Mol.AtomName[Atom_R]);
					strcat(szLonePairInfo, szTxtAdd);
					sprintf(szTxtAdd, "LONEPAIR relative %-7s%-7s%-7s%-7s distance 0.75 angle  95.00 dihe 260.00\n", 
						szLPName_2, Mol.AtomName[i], Mol.AtomName[Atom_L], Mol.AtomName[Atom_R]);
					strcat(szLonePairInfo, szTxtAdd);
					sprintf(szTxtAdd, "ANISOTROPY        %-7s%-7s%-7s%-7s A11 0.8800   A22 1.3200\n",
						Mol.AtomName[i], Mol.AtomName[Atom_L], szLPName_1, szLPName_2);
					strcat(szLonePairInfo, szTxtAdd);
				}
				else	{	// "bisector"
					sprintf(szTxtAdd, "LONEPAIR bisector %-7s%-7s%-7s%-7s distance 0.75 angle  95.00 dihe 100.00\n", 
						szLPName_1, Mol.AtomName[i], Mol.AtomName[Atom_L], Mol.AtomName[Atom_R]);
					strcat(szLonePairInfo, szTxtAdd);
					sprintf(szTxtAdd, "LONEPAIR bisector %-7s%-7s%-7s%-7s distance 0.75 angle  95.00 dihe 100.00\n", 
						szLPName_2, Mol.AtomName[i], Mol.AtomName[Atom_R], Mol.AtomName[Atom_L]);
					strcat(szLonePairInfo, szTxtAdd);
					sprintf(szTxtAdd, "ANISOTROPY        %-7s%-7s%-7s%-7s A11 0.8800   A22 1.3200\n",
						Mol.AtomName[i], Mol.AtomName[Atom_L], szLPName_1, szLPName_2);
					strcat(szLonePairInfo, szTxtAdd);
				}
				sprintf(szTxtAdd, "BOND  %7s %7s  %7s %7s\n", Mol.AtomName[i], szLPName_1, Mol.AtomName[i], szLPName_2);
				strcat(szTxtRtf_Bond, szTxtAdd);
				
				fprintf(fOut, "ATOM %-7s%-7s%7.3lf  ALPHA%7.3lf  THOLE%6.3lf\n", Mol.AtomName[i], Mol.ChemName[i], Mol.CG[i]-CG_LP*2.0, -Alpha_Miller[i], THOLE_0);
				fprintf(fOut, "ATOM %-7s%-7s%7.3lf\n", szLPName_1, "LP", CG_LP);
				fprintf(fOut, "ATOM %-7s%-7s%7.3lf\n", szLPName_2, "LP", CG_LP);

				strcpy(szLonePair_Bond_Para[n_LonePair_Bond_Para], Mol.ChemName[i]);
				n_LonePair_Bond_Para++;
				PrintOut = 1;
			}
		}
		else if(strcmp(AtomName, "O")==0)	{
			if(BondCount[i] == 1)	{	// sp2 O
				Atom_L = BondList[i][0];
				Next_Neighbor = Get_Next_Neighbor(i);
				sprintf(szLPName_1, "LPA%s", Mol.AtomName[i]);
				sprintf(szLPName_2, "LPB%s", Mol.AtomName[i]);

				sprintf(szTxtAdd, "LONEPAIR relative %-7s%-7s%-7s%-7s distance 0.35 angle 110.00 dihe   0.00\n", 
					szLPName_1, Mol.AtomName[i], Mol.AtomName[Atom_L], Mol.AtomName[Next_Neighbor]);
				strcat(szLonePairInfo, szTxtAdd);
				sprintf(szTxtAdd, "LONEPAIR relative %-7s%-7s%-7s%-7s distance 0.35 angle 110.00 dihe 180.00\n", 
					szLPName_2, Mol.AtomName[i], Mol.AtomName[Atom_L], Mol.AtomName[Next_Neighbor]);
				strcat(szLonePairInfo, szTxtAdd);
				sprintf(szTxtAdd, "ANISOTROPY        %-7s%-7s%-7s%-7s A11 0.8800   A22 1.3200\n",
					Mol.AtomName[i], Mol.AtomName[Atom_L], szLPName_1, szLPName_2);
				strcat(szLonePairInfo, szTxtAdd);
			}
			else if(BondCount[i] == 2)	{	// only sp3 O
				Atom_L = BondList[i][0];
				Atom_R = BondList[i][1];
				sprintf(szLPName_1, "LPA%s", Mol.AtomName[i]);
				sprintf(szLPName_2, "LPB%s", Mol.AtomName[i]);
				nBondedHeavyAtom = CountHeavyAtoms(i);
				
				if(nBondedHeavyAtom==1)	{	// "relative"
					if(Mol.mass[Atom_L] < 1.5)	{	// Hydrogen
						Atom_R = BondList[i][0];
						Atom_L = BondList[i][1];
					}
					sprintf(szTxtAdd, "LONEPAIR relative %-7s%-7s%-7s%-7s distance 0.35 angle 110.00 dihe  90.00\n", 
						szLPName_1, Mol.AtomName[i], Mol.AtomName[Atom_L], Mol.AtomName[Atom_R]);
					strcat(szLonePairInfo, szTxtAdd);
					sprintf(szTxtAdd, "LONEPAIR relative %-7s%-7s%-7s%-7s distance 0.35 angle 110.00 dihe 270.00\n", 
						szLPName_2, Mol.AtomName[i], Mol.AtomName[Atom_L], Mol.AtomName[Atom_R]);
					strcat(szLonePairInfo, szTxtAdd);
					sprintf(szTxtAdd, "ANISOTROPY        %-7s%-7s%-7s%-7s A11 0.8800   A22 1.3200\n",
						Mol.AtomName[i], Mol.AtomName[Atom_L], szLPName_1, szLPName_2);
					strcat(szLonePairInfo, szTxtAdd);
				}
				else	{	// "bisector"
					sprintf(szTxtAdd, "LONEPAIR bisector %-7s%-7s%-7s%-7s distance 0.35 angle 110.00 dihe  90.00\n", 
						szLPName_1, Mol.AtomName[i], Mol.AtomName[Atom_L], Mol.AtomName[Atom_R]);
					strcat(szLonePairInfo, szTxtAdd);
					sprintf(szTxtAdd, "LONEPAIR bisector %-7s%-7s%-7s%-7s distance 0.35 angle 110.00 dihe 270.00\n", 
						szLPName_2, Mol.AtomName[i], Mol.AtomName[Atom_L], Mol.AtomName[Atom_R]);
					strcat(szLonePairInfo, szTxtAdd);
					sprintf(szTxtAdd, "ANISOTROPY        %-7s%-7s%-7s%-7s A11 0.8800   A22 1.3200\n",
						Mol.AtomName[i], Mol.AtomName[Atom_L], szLPName_1, szLPName_2);
					strcat(szLonePairInfo, szTxtAdd);
				}
//				sprintf(szTxtAdd, "BOND  %7s %7s   %7s %7s\n", Mol.AtomName[i], szLPName_1, Mol.AtomName[i], szLPName_2);
//				strcat(szTxtRtf_Bond, szTxtAdd);
			}
			else	{
				Quit_With_Error_Msg("Add_Lonepair_Anisotropy() > Unknown Oxygen atom type.\nQuit\n");
			}
			sprintf(szTxtAdd, "BOND  %7s %7s   %7s %7s\n", Mol.AtomName[i], szLPName_1, Mol.AtomName[i], szLPName_2);
			strcat(szTxtRtf_Bond, szTxtAdd);
			
			fprintf(fOut, "ATOM %-7s%-7s%7.3lf  ALPHA%7.3lf  THOLE%6.3lf\n", Mol.AtomName[i], Mol.ChemName[i], Mol.CG[i]-CG_LP*2.0, -Alpha_Miller[i], THOLE_0);
			fprintf(fOut, "ATOM %-7s%-7s%7.3lf\n", szLPName_1, "LP", CG_LP);
			fprintf(fOut, "ATOM %-7s%-7s%7.3lf\n", szLPName_2, "LP", CG_LP);
			
			strcpy(szLonePair_Bond_Para[n_LonePair_Bond_Para], Mol.ChemName[i]);
			n_LonePair_Bond_Para++;
			PrintOut = 1;
		}

		if(PrintOut == 0)	{
			if(Alpha_Miller[i] > 1.0E-10)	{
				fprintf(fOut, "ATOM %-7s%-7s%7.3lf  ALPHA%7.3lf  THOLE%6.3lf\n", Mol.AtomName[i], Mol.ChemName[i], Mol.CG[i], -Alpha_Miller[i], THOLE_0);
			}
			else	{
				fprintf(fOut, "ATOM %-7s%-7s%7.3lf\n", Mol.AtomName[i], Mol.ChemName[i], Mol.CG[i]);
			}
		}
	}

	fprintf(fOut, "\n%s\n%s", szTxtRtf_Bond, szLonePairInfo);

	Reorganize_Lonepair_Bond_Para();
}

int Count_Bonded_SP2_Atoms(int Idx, int IsRestricted)
{
	int i, Count=0;

	if(IsRestricted)	{
		for(i=0; i<BondCount[Idx]; i++)	{
			if( SP2_Flag[BondList[Idx][i]] == 1)	{	// real SP2 atom
				Count++;
			}
		}
	}
	else	{
		for(i=0; i<BondCount[Idx]; i++)	{
			if( SP2_Flag[BondList[Idx][i]] >= 1)	{	// real SP2 atom + promoted SP2 atom
				Count++;
			}
		}
	}

	return Count;
}

#define ALPHA_H	(0.387)
void Assign_Alpha_Miller(void)
{
	int i, nAtom;
	char *AtomName;

	nAtom = Mol.nAtom;
	for(i=0; i<nAtom; i++)	{
		AtomName = Name_Elem[AtomType[i]];

		if(strcmp(AtomName, "H")==0)	{
			Alpha_Miller[i] = 0.0;
		}
		else if(strcmp(AtomName, "F")==0)	{
			Alpha_Miller[i] = 0.296;
		}
		else if(strcmp(AtomName, "CL")==0)	{
			Alpha_Miller[i] = 2.315;
		}
		else if(strcmp(AtomName, "BR")==0)	{
			Alpha_Miller[i] = 3.013;
		}
		else if(strcmp(AtomName, "I")==0)	{
			Alpha_Miller[i] = 5.415;
		}
		else if(strcmp(AtomName, "P")==0)	{	// Benoit's value
			Alpha_Miller[i] = 2.063;
		}
		else if(strcmp(AtomName, "C")==0)	{
			if(BondCount[i] == 4)	{	// sp3
				Alpha_Miller[i] = 1.061 + ALPHA_H*H_Count[i];
			}
			else if(BondCount[i] == 3)	{	// sp2
				if(Count_Bonded_SP2_Atoms(i, 0)==3)	{
					Alpha_Miller[i] = 1.896;
				}
				else	{
					Alpha_Miller[i] = 1.352 + ALPHA_H*H_Count[i];
				}
			}
			else if(BondCount[i] == 2)	{	// sp
				Alpha_Miller[i] = 1.283 + ALPHA_H*H_Count[i];
			}
			else	{
				fprintf(fFile_Run_Log, "Unsupported type for C atom. Index = %d %s\nQuit\n", i+1, AtomName);
				fflush(fFile_Run_Log);
				exit(1);
			}
		}
		else if(strcmp(AtomName, "N")==0)	{
			if(BondCount[i] == 4)	{
				Alpha_Miller[i] = 1.0;	// Added by Lei
			}
			else if(BondCount[i] == 3)	{
				if(SP2_Flag[i] > 0)	{
					Alpha_Miller[i] = 1.090 + ALPHA_H*H_Count[i];
				}
				else	{
					Alpha_Miller[i] = 0.964 + ALPHA_H*H_Count[i];
				}
			}
			else if(BondCount[i] == 2)	{	// sp2
				Alpha_Miller[i] = 1.030 + ALPHA_H*H_Count[i];
			}
			else if(BondCount[i] == 1)	{	// sp
				Alpha_Miller[i] = 0.956 + ALPHA_H*H_Count[i];
			}
			else	{
				fprintf(fFile_Run_Log, "Unsupported type for N atom. Index = %d %s\nQuit\n", i+1, AtomName);
				fflush(fFile_Run_Log);
				exit(1);
			}
		}
		else if(strcmp(AtomName, "O")==0)	{
			if(BondCount[i] == 2)	{
				if(SP2_Flag[i] > 0)	{
					Alpha_Miller[i] = 0.274;
				}
				else	{
					Alpha_Miller[i] = 0.637 + ALPHA_H*H_Count[i];
				}
			}
			else if(BondCount[i] == 1)	{	// sp2
				if(SP2_O_Count[BondList[i][0]] == 2)	{	// anionic oxygen
					Alpha_Miller[i] = 0.858;
				}
				else	{	// regular C=O 
					Alpha_Miller[i] = 0.569;
				}
			}
			else	{
				Alpha_Miller[i] = 0.57 + 0.17*H_Count[i];	// others
//				fprintf(fFile_Run_Log, "Unsupported type for O atom. Index = %d %s\nQuit\n", i+1, AtomName);
//				fflush(fFile_Run_Log);
//				exit(1);
			}
		}
		else if(strcmp(AtomName, "S")==0)	{
			if(BondCount[i] == 2)	{
				if(SP2_Flag[i] > 0)	{
					Alpha_Miller[i] = 2.700;
				}
				else	{
					Alpha_Miller[i] = 3.000 + ALPHA_H*H_Count[i];
				}
			}
			else if(BondCount[i] == 1)	{	// sp2
				Alpha_Miller[i] = 3.729;
			}
			else	{
				Alpha_Miller[i] = 2.99 + 0.17*H_Count[i];	// others
//				fprintf(fFile_Run_Log, "Unsupported type for S atom. Index = %d %s\nQuit\n", i+1, AtomName);
//				fflush(fFile_Run_Log);
//				exit(1);
			}
		}
		else	{
			fprintf(fFile_Run_Log, "Unsupported atom type: atom %s, index = %d \nQuit\n", Mol.AtomName[i], i+1);
			fflush(fFile_Run_Log);
			exit(1);
		}

		if(Alpha_Miller[i] > 1.0E-10)	{
			printf("Alpha(%4s) = %7.3lf\n", Mol.AtomName[i], Alpha_Miller[i]);
		}
	}
}
#undef ALPHA_H

void Add_Tip3_Prm_File(void)
{
	FILE *fOut;
	char *ReadLine, szNewPrm[]="drude-mol.prm";
	int i, j, n_Line;

	Found_Improper = 0;
	fOut = fopen(szPrmFile, "r");
	if(fOut == NULL)	{
		printf("Error in open file: %s\nQuit\n", szPrmFile);
		exit(1);
	}
	else	{
		n_Line = 0;

		while(1)	{
			if(feof(fOut))	{
				break;
			}
			ReadLine = fgets(szTxtPrm[n_Line], MAX_LEN, fOut);
			To_Upper_Case(szTxtPrm[n_Line]);
			if(ReadLine == NULL)	{
				break;
			}
			else	{
				if(strncmp(szTxtPrm[n_Line], "BOND", 4)==0)	{
					sprintf(szTxtPrm[n_Line], "BONDS\n");			// rename "BOND" to "BONDS"
				}
				else if(strncmp(szTxtPrm[n_Line], "ANGLE", 5)==0)	{
					sprintf(szTxtPrm[n_Line], "ANGLES\n");			// rename "ANGLE" to "ANGLES"
					Line_Insert_Bond = n_Line;
				}
				else if(strncmp(szTxtPrm[n_Line], "DIHEDRAL", 8)==0)	{
					sprintf(szTxtPrm[n_Line], "DIHEDRALS\n");			// rename "DIHEDRAL" to "DIHEDRALS"
					Line_Insert_Angle = n_Line;
				}
				else if(strncmp(szTxtPrm[n_Line], "IMPROPER", 8)==0)	{
					Found_Improper = 1;
					sprintf(szTxtPrm[n_Line], "IMPROPERS\n");			// rename "IMPHI" to "IMPROPERS"
				}
				else if(strncmp(szTxtPrm[n_Line], "IMPHI", 5)==0)	{
					Found_Improper = 1;
					sprintf(szTxtPrm[n_Line], "IMPROPERS\n");			// rename "IMPHI" to "IMPROPERS"
				}
				else if(strncmp(szTxtPrm[n_Line], "NONBONDED", 9)==0)	{
					Line_Insert_Impropers = n_Line;
				}
				else if(strncmp(szTxtPrm[n_Line], "END", 3)==0)	{	// NOT saved
					n_Line --;
				}


				n_Line++;
				if(n_Line >= MAX_LINE)	{
					printf("n_Line >= MAX_LINE\nQuit\n");
					fclose(fOut);
					exit(1);
				}
			}
		}
		fclose(fOut);
	}

	fOut = fopen(szNewPrm, "w");
	for(i=0; i<n_Line; i++)	{
		if(i == Line_Insert_Bond)	{
			for(j=0; j<n_LonePair_Bond_Para; j++)	{
				if(LonePair_Bond_Para_Output[j])	{	// to add an entry for lone pair bond parameter
					fprintf(fOut, "%-9s%-9s  0.00      0.000\n", szLonePair_Bond_Para[j], "LP");
				}
			}
			
//			fprintf(fOut, "X        LP*        0.00      0.000\n");
			fprintf(fOut, "X        DRUD     500.00      0.000\n");
			fprintf(fOut, "ODW      HDW      450.00      0.9572     ! SWM4, SWM4-NDP water, Guillaume 2005\n");
			fprintf(fOut, "ODW      LP         0.00      0.24034492 ! SWM4, SWM4-NDP water, Guillaume 2005\n");
			fprintf(fOut, "ODW      DOH2     500.00      0.000      ! SWM4, SWM4-NDP water, Guillaume 2005\n");
			fprintf(fOut, "HDW      HDW        0.00      1.5139     ! SWM4, SWM4-NDP water, Guillaume 2005\n\n");
		}
		else if(i == Line_Insert_Angle)	{
			fprintf(fOut, "HDW      ODW      HDW       55.000   104.52 ! SWM4-NDP water Guillaume 2005\n\n");
		}
		else if( (Found_Improper==0) && (i==Line_Insert_Impropers) )	{	// add the entry for IMPROPERS
			fprintf(fOut, "IMPROPERS\n\n");
		}

		if(FindString(szTxtPrm[i], "NONBONDED") >= 0)	{
			fprintf(fOut, "NONBONDED  E14FAC %9.6lf\n", E14FAC);
			continue;
		}
		else if(FindString(szTxtPrm[i], "NBXMOD") >= 0)	{
			continue;
		}
		else if(FindString(szTxtPrm[i], "CUTNB") >= 0)	{
			continue;
		}
		else if(FindString(szTxtPrm[i], "E14FAC") >= 0)	{
			continue;
		}
		else if(FindString(szTxtPrm[i], "WMIN") >= 0)	{
			continue;
		}

		fprintf(fOut, "%s", szTxtPrm[i]);
	}
	fprintf(fOut, "\nHDW      0.0   -0.0000     0.0000 ! SWM4, SWM4-NDP water, GL, 2005\n");
	fprintf(fOut, "ODW      0.0   -0.21094325 1.78692899 ! SWM4, SWM4-NDP water, GL, 2005\n\n");
	fprintf(fOut, "D*       0.0   -0.0000    0.0000 ! Wildcard for Drudes and dummy atom\n");
	fprintf(fOut, "DRUD     0.0   -0.0000    0.0000 ! Wildcard for Drudes and dummy atom\n");
	fprintf(fOut, "LP       0.0   -0.0000    0.0000\n");
	fprintf(fOut, "CALD     0.0   -0.2100000 1.2708552 ! Ca\n\n");	// used for the perturbation charge



	fclose(fOut);
}


void Add_Tip3_Rtf_File(void)
{
	FILE *fOut;
	int nLine=0, i, ReadItem;
	char *ReadLine, szNewRtf[]="drude-mol.rtf", szBuff[256];
	double net_cg_local=0.0;

	fOut = fopen(szRtfFile, "r");
	if(fOut == NULL)	{
		printf("Error in open file: %s\nQuit\n", szRtfFile);
		exit(1);
	}
	else	{
		nLine = 0;

		while(1)	{
			if(feof(fOut))	{
				break;
			}
			ReadLine = fgets(szTxtRtf[nLine], MAX_LEN, fOut);
			To_Upper_Case(szTxtRtf[nLine]);
			if(ReadLine == NULL)	{
				break;
			}
			else	{
				if(strncmp(szTxtRtf[nLine], "GROUP", 5)==0)	{
					Line_Insert_Atom_RTF=nLine+1;
				}
				else if(strncmp(szTxtRtf[nLine], "RESI", 4)==0)	{
					ReadItem = sscanf(szTxtRtf[nLine], "%s%s%lf", szBuff, szBuff, &net_cg_local);
					if(ReadItem == 3)	{
						sprintf(szTxtRtf[nLine], "RESI  MOL  %6.3lf\n", net_cg_local);
					}
					else	{
						ReadItem = sscanf(szTxtRtf[nLine], "%s%lf", szBuff, &net_cg_local);
						if(ReadItem == 2)	{
							sprintf(szTxtRtf[nLine], "RESI  MOL  %6.3lf\n", net_cg_local);
						}
					}
					
					Line_Insert_Atom_RTF=nLine+1;
				}
				else if(strncmp(szTxtRtf[nLine], "MASS", 4)==0)	{
					LastLine_Mass = nLine+1;
				}
				else if(strncmp(szTxtRtf[nLine], "ATOM", 4)==0)	{
					nLine--;	// not saved in szTxtRtf
				}
				else if(strncmp(szTxtRtf[nLine], "BOND", 4)==0)	{
					strcat(szTxtRtf_Bond, szTxtRtf[nLine]);
					nLine--;	// not saved in szTxtRtf
				}
				else if(strncmp(szTxtRtf[nLine], "ANGL", 4)==0)	{
					nLine--;	// not saved in szTxtRtf
				}
				else if(strncmp(szTxtRtf[nLine], "DIHE", 4)==0)	{
					nLine--;	// not saved in szTxtRtf
				}
				else if(strncmp(szTxtRtf[nLine], "END", 3)==0)	{
					nLine--;	// not saved in szTxtRtf
				}
				nLine++;
				if(nLine >= MAX_LINE)	{
					printf("n_Line >= MAX_LINE\nQuit\n");
					fclose(fOut);
					exit(1);
				}
			}
		}
		fclose(fOut);
	}

	fOut = fopen(szNewRtf, "w");

	for(i=0; i<nLine; i++)	{
		if(i==Line_Insert_Atom_RTF)	{	// to insert the entries for ATOM and BOND
			Add_Lonepair_Anisotropy(fOut);
			if(strlen(szTxtRtf[i])<2)	{	// skip one black line
				i++;
			}
		}
		if(i==LastLine_Mass)	{
			fprintf(fOut, "MASS   113 LP        0.00000 H  ! general lone pair\n");
			fprintf(fOut, "MASS   114 DRUD      0.00000 H  ! drude particle\n\n");
//			fprintf(fOut, "AUTOGENERATE ANGLES DIHEDRALS DRUDE  !note use of DRUDE\n\n");
			fprintf(fOut, "AUTOGENERATE DRUDE  !note use of DRUDE\n\n");
		}
		if(i<nLine) fprintf(fOut, "%s", szTxtRtf[i]);
	}

	fprintf(fOut, "\n\n!! SWM4 TIP4P polarizable water model\n");
	fprintf(fOut, "MASS   151 ODW      15.99940 O  ! water oxygen\n");
	fprintf(fOut, "MASS   152 HDW       1.00800 H  ! water hydrogen\n");
	fprintf(fOut, "MASS   153 DOH2      0.00000 H  ! water Drude\n\n");

	fprintf(fOut, "! SWM4-NDP water\n");
	fprintf(fOut, "RESI SWM4          0.000\n");
	fprintf(fOut, "GROUP\n");
	fprintf(fOut, "ATOM OH2  ODW      0.00000  TYPE DOH2    ALPHA -0.97825258 THOLE 1.3   \n");
	fprintf(fOut, "ATOM OM   LP      -1.11466\n");
	fprintf(fOut, "ATOM H1   HDW      0.55733\n");
	fprintf(fOut, "ATOM H2   HDW      0.55733\n");
	fprintf(fOut, "BOND OH2 H1\n");
	fprintf(fOut, "BOND OH2 H2\n");
	fprintf(fOut, "BOND OH2 OM\n");
	fprintf(fOut, "BOND H1  H2 ! for SHAKE\n");
	fprintf(fOut, "ANGLE H1 OH2 H2\n");
	fprintf(fOut, "LONEPAIR bisector OM OH2 H1 H2 distance 0.24034492 angle 0.0 dihe 0.0\n");
	fprintf(fOut, "IC   H1   OH2   H2    H1    0.9572 104.52     0.00   37.74   1.5139\n");
	fprintf(fOut, "IC   H1   OM   *OH2   H2    0.9572  52.26   180.00   52.26   0.9572\n");
	fprintf(fOut, "IC   H2   H1    OH2   OM    1.5139  37.74     0.01   52.26   0.24034492\n");
	fprintf(fOut, "PATCH FIRST NONE LAST NONE\n");

	fprintf(fOut, "\n\nEND\n");

	fclose(fOut);
}

void To_Upper_Case(char szBuff[])
{
	char gap;
	int i=0;

	gap = 'A'-'a';

	while(1)	{
		if(szBuff[i]==0)	{
			break;
		}
		if( (szBuff[i]>='a') && (szBuff[i]<='z') )	{
			szBuff[i] += gap;
		}
		i++;
	}
}

int FindString(char szBuff[], const char szTag[])
{
	int Len_Buff, Len_Tag, iPos, iPos_Max;

	Len_Tag = strlen(szTag);
	Len_Buff = strlen(szBuff);

	iPos_Max = Len_Buff-Len_Tag;

	for(iPos=0; iPos<iPos_Max; iPos++)	{
		if(strncmp(szBuff+iPos, szTag, Len_Tag) == 0)	{
			return iPos;
		}
	}
	return (-1);
}

void GetAtomTypeFromMass(void)
{
	FILE *fOut;
	int i, nAtom;
	double *mass;
	
	nAtom = Mol.nAtom;
	mass = Mol.mass;

	fOut = fopen(szFileElemList, "w");
	for(i=0; i<nAtom; i++)	{
		AtomType[i] = DetermineAtomType(Mol.AtomName[i], mass[i]);

		fprintf(fOut, "%s\n", Name_Elem[AtomType[i]]);
	}

	fclose(fOut);
}

int DetermineAtomType(char szName[], double mass)
{
	int i;

	for(i=0; i<N_ELEM; i++)	{
		if( (Name_Elem[i][0] == szName[0]) && (fabs(mass - Mass_List[i]) < 0.2) )	{
			return i;
		}
	}

	printf("Fail to identify atom: %s with mass = %8.3lf\n", szName, mass);

	exit(1);

	return -1;
}

void Quit_With_Error_Msg(char szMsg[])
{
	fprintf(fFile_Run_Log, "%s", szMsg);
	fflush(fFile_Run_Log);
	exit(1);
}

int QueryAtomIndex(char szAtmName[])
{
	int Idx, nAtom;

	nAtom = Mol.nAtom;
	for(Idx=0; Idx<nAtom; Idx++)	{
		if(strcmp(szAtmName, Mol.AtomName[Idx])==0)	{
			return Idx;
		}
	}
	fprintf(fFile_Run_Log, "Fail to identify the atom name from %s\nQuit\n", szAtmName);
	fflush(fFile_Run_Log);
	exit(1);

	return -1;
}


