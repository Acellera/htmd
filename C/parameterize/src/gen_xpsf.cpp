/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ff.h"

#define MAX_ATOM	(512)
#define MAX_N_LP		(100)
#define MAX_N_ANISO		(100)
#define MAX_RES_REC		(1024)
#define MAX_MASS_REC	(1024)
#define MAX_LINE_RTF	(65536)
#define MAX_SIZE_TXT	(65536)
#define MASS_DRUDE		(0.4)

#define MAX_ANGLE		(65536)
#define MAX_DIHEDRAL	(65536)

#define MD_COULOMB  (332.0716)	//(in CHARMM)
#define K_DRUDE		(500.0)

#define MAX_BOND_PER_LINE	(32)
#define MAX_DIHEDRAL_PER_LINE	(32)
#define MAX_IMPROPER_PER_LINE	(32)
#define MAX_ANGLE_PER_LINE	(32)
#define LARGE_d		(100)

#ifndef MINAB
#define MINAB 
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif

FILE *fFile_Run_Log;     // will be shared by other source code just for compiling
void Quit_With_Error_Msg(char szMsg[]);


int FindString(char szBuff[], const char szTag[]);
int FindTags_As_First_String(char szBuff[], char szTag[]);
void To_Upper_Case(char szBuff[], char szNew[]);
int QueryResidue(char ResName[]);
int QueryChemName(char AtomName[]);
int Get_Alpha_Thole(char szBuff[], double& alpha, double& thole);
int Get_Global_Atom_Index(char szName[], int iBegin, int iEnd);
int SplitString(char szBuff[], char ItemList[][16]);
int Extract_Entry_Per_Line(char szBuff[], int Entry_List[], int Atm_First, int Atm_Last);
void Para_Aniso_A_K(double& K11, double& K22, double& K33, double A11, double A22);
int Query_Angle(int a, int b, int c);
int Query_Dihedral(int a, int b, int c, int d);

char ItemList[MAX_BOND_PER_LINE][16];
int Bond_in_Line[MAX_BOND_PER_LINE], BondCount[MAX_ATOM], BondList[MAX_ATOM][6];
int Dihedral_in_Line[MAX_DIHEDRAL_PER_LINE], Improper_in_Line[MAX_IMPROPER_PER_LINE], Angle_in_Line[MAX_ANGLE_PER_LINE];


typedef struct	{
	int IsHeavyAtom;
	int Type;
	char Name[16];
	char ChemName[16];
	double mass, CG, alpha, thole;
}MYATOM, *PMYATOM;

typedef struct	{
	char ResName[16];
	int nAtom;
	double CG;
	int Line_First, Line_Last, Line_Last_Bond;
}RESIDUE_REC, *PRESIDUE_REC;

typedef struct	{
	int id;
	double mass;
	char Name[16];
	char Elem[4];
}MASS_REC, *PMASS_REC;

typedef struct{
	char szName[4][16];
	int Atom[4];
	double r, theta, phi;
}LONEPAIR;

typedef struct{
	char szName[4][16];
	int Atom[4];
	double A11, A22;
	double k11, k22, k33;
}ANISOTROPY;


char szTxt_Atom[MAX_SIZE_TXT], szTxt_Bond[MAX_SIZE_TXT], szTxt_Angle[MAX_SIZE_TXT], szTxt_Dihedral[MAX_SIZE_TXT], szTxt_Improper[MAX_SIZE_TXT];
char szTxt_LP[MAX_SIZE_TXT], szTxt_Aniso[MAX_SIZE_TXT];

char All_Atom[MAX_ATOM][16];
int IsHeavyAtom[MAX_ATOM];

int nRes, ResList[MAX_RES_REC];
int nAtom, nBond, nAngle, nDihedral, nImproper, nLP, nAniso;	// C-Map will be added later !!!
int nRes_Rec=0, nMass_Rec=0;
int IsLP[MAX_ATOM], IsDrude[MAX_ATOM];
MYATOM AtomList[MAX_ATOM];
RESIDUE_REC Res_Rec[MAX_RES_REC];
MASS_REC Mass_Rec[MAX_MASS_REC];
LONEPAIR LP_List[MAX_N_LP];
ANISOTROPY Aniso_List[MAX_N_ANISO];

int Angle_List[MAX_ANGLE][4];
int Dihedral_List[MAX_DIHEDRAL][4];

int ActiveLine[MAX_LINE_RTF];
char Dist_Mat[MAX_ATOM][MAX_ATOM];


int main(int argc, char *argv[])
{

  timebomb();

	FILE *fIn, *fOut;
	char szLine[256], szLine_Upper[256], *ReadLine, szTmp[256], szAddLine[256], szLPType[256];
	int i, j, k, l, ReadItem, LineCount=0, ActivePara, IdxRes, LineNum, IdxLine;
	int Atom_i, Atom_j, Atom_k, Atom_l;
	int nAtomLocal, To_Add_Test_Charge=0, Atom_Begin, Atom_End, Bond_Line, Angle_Line, Improper_Line;
	double fTmp=0.0, alpha, thole, CG_Host, CG_Drude;

	if(argc < 4)	{
		printf("Usage: genxpsf your-rtf output.xpsf res_1 res_2 res_3 ...\n");
		exit(1);
	}

	szTxt_Atom[0] = szTxt_Bond[0] = szTxt_Angle[0] = szTxt_Dihedral[0] = szTxt_Improper[0] = szTxt_LP[0] = szTxt_Aniso[0] = 0;


	fIn = fopen(argv[1], "r");

	if(fIn == NULL)	{
		printf("Fail to open the force field (or rtf) file: %s\nQuit\n", argv[1]);
		exit(1);
	}

	for(i=0; i<MAX_ATOM; i++)	{
		for(j=0; j<MAX_ATOM; j++)	{
			Dist_Mat[i][j] = LARGE_d;
		}
		Dist_Mat[i][i] = 0;
	}

	memset(ActiveLine, 0, sizeof(int)*MAX_LINE_RTF);
	memset(IsHeavyAtom, 0, sizeof(int)*MAX_ATOM);
	memset(IsLP, 0, sizeof(int)*MAX_ATOM);
	memset(IsDrude, 0, sizeof(int)*MAX_ATOM);
	memset(BondCount, 0, sizeof(int)*MAX_ATOM);
	memset(BondList, 0, sizeof(int)*MAX_ATOM*4);
	
	nRes_Rec = -1;
	while(1)	{
		if(feof(fIn))	{
			if(nRes_Rec == 0)	{
				Res_Rec[nRes_Rec].Line_Last = LineCount - 1;
				nRes_Rec++;
			}
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			if(nRes_Rec == 0)	{
				Res_Rec[nRes_Rec].Line_Last = LineCount - 1;
				nRes_Rec++;
			}
			break;
		}
		else	{
			To_Upper_Case(szLine, szLine_Upper);
			if(FindTags_As_First_String(szLine_Upper, "MASS"))	{
//				ReadItem = sscanf(szLine, "%s %d %s %lf %s", szTmp, &(Mass_Rec[nMass_Rec].id), Mass_Rec[nMass_Rec].Name, &(Mass_Rec[nMass_Rec].mass), Mass_Rec[nMass_Rec].Elem);	// charmm format
				ReadItem = sscanf(szLine_Upper, "%s %d %s %lf", szTmp, &(Mass_Rec[nMass_Rec].id), Mass_Rec[nMass_Rec].Name, &(Mass_Rec[nMass_Rec].mass));	// the RTF from Antechamber
//				if(ReadItem == 5)	{
				if(ReadItem == 4)	{
					nMass_Rec++;
				}
			}
			else if(FindTags_As_First_String(szLine_Upper, "RESI"))	{
				nRes_Rec++;
				Res_Rec[nRes_Rec].Line_First = LineCount;

				ReadItem = sscanf(szLine_Upper, "%s %s %lf", szTmp, Res_Rec[nRes_Rec].ResName, &(Res_Rec[nRes_Rec].CG));
				if(ReadItem != 3)	{
					printf("Error in reading an entry of RESI.\n%s\nQuit\n", szLine);
					exit(1);
				}

				if(nRes_Rec > 0)	{
					Res_Rec[nRes_Rec-1].Line_Last = LineCount - 1;
				}
			}
			else if(FindTags_As_First_String(szLine_Upper, "ATOM"))	{
				Res_Rec[nRes_Rec].nAtom ++;
			}
			else if(FindTags_As_First_String(szLine_Upper, "BOND") && (FindTags_As_First_String(szLine_Upper, "BONDS")==0))	{	// "BOND "
				Res_Rec[nRes_Rec].Line_Last_Bond = LineCount+1;
			}
			else if(FindTags_As_First_String(szLine_Upper, "END"))	{	// the beginning of parameters or the end of RTF file
				Res_Rec[nRes_Rec].Line_Last = LineCount - 1;
				nRes_Rec++;
				LineCount++;
				break;
			}
			LineCount++;
		}
	}
	LineNum = LineCount;


	nRes = 0;
	ActivePara = 3;
	while(ActivePara < argc)	{
		To_Upper_Case(argv[ActivePara], szLine_Upper);
		
		if(strcmp(szLine_Upper, "TEST")==0)	{
			To_Add_Test_Charge = 1;
			ActivePara++;
		}
		else	{
			IdxRes = QueryResidue(argv[ActivePara]);
			ActivePara++;
			ResList[nRes] = IdxRes;
			nRes++;
		}
	}


	nAtom = nBond = nAngle = nDihedral = nImproper = nLP = nAniso = 0;

	for(IdxRes=0; IdxRes<nRes; IdxRes++)	{
		fseek(fIn, 0, SEEK_SET);	// starting from the beginning of the rtf file
		nAtomLocal = 0;
		Atom_Begin = nAtom;
		for(IdxLine=0; IdxLine<LineNum; IdxLine++)	{
			fgets(szLine, 256, fIn);
			if( (IdxLine>=Res_Rec[ResList[IdxRes]].Line_First) && (IdxLine<=Res_Rec[ResList[IdxRes]].Line_Last) )	{
				To_Upper_Case(szLine, szLine_Upper);

				//start	to parse the entry for atoms
				if(FindTags_As_First_String(szLine_Upper, "ATOM"))	{
					ReadItem = sscanf(szLine_Upper, "%s %s %s %lf", szTmp, AtomList[nAtomLocal].Name, AtomList[nAtomLocal].ChemName, &(AtomList[nAtomLocal].CG));
					if(ReadItem != 4)	{
						printf("Error in read the entry of atom: %s\nQuit\n", szLine);
						exit(1);
					}
					AtomList[nAtomLocal].Type = QueryChemName(AtomList[nAtomLocal].ChemName);
					AtomList[nAtomLocal].IsHeavyAtom = Get_Alpha_Thole(szLine_Upper, alpha, thole);
					AtomList[nAtomLocal].mass = Mass_Rec[AtomList[nAtomLocal].Type].mass;

					if(AtomList[nAtomLocal].IsHeavyAtom)	{	// the host and drude
						CG_Drude = -sqrt(2.0*K_DRUDE*fabs(alpha)/MD_COULOMB);
						CG_Host = AtomList[nAtomLocal].CG - CG_Drude;
						sprintf(szAddLine, "%10d %-6s %3d        %-8s %-8s %-8s %12.5lf%10.5lf           0%10.5lf %14.5lf\n", 
							nAtom+1, Res_Rec[ResList[IdxRes]].ResName, IdxRes+1, Res_Rec[ResList[IdxRes]].ResName, 
							AtomList[nAtomLocal].Name, AtomList[nAtomLocal].ChemName, CG_Host, AtomList[nAtomLocal].mass-MASS_DRUDE, 
							alpha, thole);
						strcat(szTxt_Atom, szAddLine);
						IsHeavyAtom[nAtom] = 1;
						strcpy(All_Atom[nAtom], AtomList[nAtomLocal].Name);
						sprintf(szAddLine, "%10d %-6s %3d        %-8s D%-7s %-8s %12.5lf%10.5lf           0%10.5lf %14.5lf\n", 
							nAtom+2, Res_Rec[ResList[IdxRes]].ResName, IdxRes+1, Res_Rec[ResList[IdxRes]].ResName, 
							AtomList[nAtomLocal].Name, "DRUD", CG_Drude, MASS_DRUDE, 
							0.0, 0.0);
						strcat(szTxt_Atom, szAddLine);
						nAtom++;	// add one more for the drude
						strcpy(All_Atom[nAtom], AtomList[nAtomLocal].Name);
						IsDrude[nAtom] = 1;	// this is a drude
					}
					else	{
						if(AtomList[nAtomLocal].mass < 1.0E-10)	{	// lone pair
							sprintf(szAddLine, "%10d %-6s %3d        %-8s %-8s %-8s %12.5lf%10.5lf          -1%10.5lf %14.5lf\n", 
								nAtom+1, Res_Rec[ResList[IdxRes]].ResName, IdxRes+1, Res_Rec[ResList[IdxRes]].ResName, 
								AtomList[nAtomLocal].Name, AtomList[nAtomLocal].ChemName, AtomList[nAtomLocal].CG, AtomList[nAtomLocal].mass, 
								0.0, 0.0);
							strcat(szTxt_Atom, szAddLine);
							IsLP[nAtom] = 1;
						}
						else	{
							sprintf(szAddLine, "%10d %-6s %3d        %-8s %-8s %-8s %12.5lf%10.5lf           0%10.5lf %14.5lf\n", 
								nAtom+1, Res_Rec[ResList[IdxRes]].ResName, IdxRes+1, Res_Rec[ResList[IdxRes]].ResName, 
								AtomList[nAtomLocal].Name, AtomList[nAtomLocal].ChemName, AtomList[nAtomLocal].CG, AtomList[nAtomLocal].mass, 
								0.0, 0.0);
							strcat(szTxt_Atom, szAddLine);
						}
						strcpy(All_Atom[nAtom], AtomList[nAtomLocal].Name);
					}
					nAtom++;
					Atom_End = nAtom;
				}
				//end	to parse the entry for atoms

				//start	to parse the entry for bonds
				if(FindTags_As_First_String(szLine_Upper, "BOND"))	{
					Bond_Line = Extract_Entry_Per_Line(szLine_Upper, Bond_in_Line, Atom_Begin, Atom_End);
					for(i=0; i<Bond_Line; i+=2)	{
						BondList[Bond_in_Line[i]][BondCount[Bond_in_Line[i]]] = Bond_in_Line[i+1];
						BondCount[Bond_in_Line[i]] ++;
						BondList[Bond_in_Line[i+1]][BondCount[Bond_in_Line[i+1]]] = Bond_in_Line[i];
						BondCount[Bond_in_Line[i+1]] ++;

						sprintf(szAddLine, "%10d%10d", 
							min(Bond_in_Line[i], Bond_in_Line[i+1])+1, max(Bond_in_Line[i], Bond_in_Line[i+1])+1);
						strcat(szTxt_Bond, szAddLine);
						nBond++;
						if(nBond % 4 == 0)	{
							strcat(szTxt_Bond, "\n");
						}
					}
				}
				//end	to parse the entry for bonds


				//start	to parse the entry for angles
				if(strncmp(szLine_Upper, "ANGL", 4)==0)	{
					if( (strcmp(Res_Rec[ResList[IdxRes]].ResName, "TIP3")==0) || (strcmp(Res_Rec[ResList[IdxRes]].ResName, "SWM4")==0) )	{	// ONLY record the angle in water !!
						Angle_Line = Extract_Entry_Per_Line(szLine_Upper, Angle_in_Line, Atom_Begin, Atom_End);;
						for(i=0; i<Angle_Line; i+=3)	{
							sprintf(szAddLine, "%10d%10d%10d", 
								Angle_in_Line[i]+1, Angle_in_Line[i+1]+1, Angle_in_Line[i+2]+1);
							strcat(szTxt_Angle, szAddLine);
							Angle_List[nAngle][0] = Angle_in_Line[i];	Angle_List[nAngle][1] = Angle_in_Line[i+1];	Angle_List[nAngle][2] = Angle_in_Line[i+2];	
							nAngle++;
							if(nAngle % 3 == 0)	{
								strcat(szTxt_Angle, "\n");
							}
						}
					}
				}
				//end	to parse the entry for angles
			

				if(IdxLine == Res_Rec[ResList[IdxRes]].Line_Last_Bond)	{
					//start	to parse the entry for angles, auto angle
					//start	to constrcut distance matrix
					for(i=Atom_Begin; i<Atom_End; i++)	{	//start	enumeration
						Atom_i = i;
						
						for(j=0; j<BondCount[Atom_i]; j++)	{	//i->j
							Atom_j = BondList[Atom_i][j];
							
							if(Dist_Mat[Atom_i][Atom_j] > 1)	{
								Dist_Mat[Atom_i][Atom_j] = 1;
								Dist_Mat[Atom_j][Atom_i] = 1;
							}
							
							for(k=0; k<BondCount[Atom_j]; k++)	{	//i->j->k
								Atom_k = BondList[Atom_j][k];
								
								if( Atom_k != Atom_i )	{
									if(Dist_Mat[Atom_i][Atom_k] >= 2)	{
										Dist_Mat[Atom_i][Atom_k] = 2;
										Dist_Mat[Atom_k][Atom_i] = 2;
									}

									if( (strcmp(Res_Rec[ResList[IdxRes]].ResName, "TIP3")!=0) && (strcmp(Res_Rec[ResList[IdxRes]].ResName, "SWM4")!=0) )	{	// skip the angle in water, TIP3 and SWM4 model
										if( (IsLP[Atom_i]==0) && (IsLP[Atom_j]==0) && (IsLP[Atom_k]==0))	{	// skip those angles with LP involved
											if( Query_Angle(Atom_i, Atom_j, Atom_k) < 0 )	{
												sprintf(szAddLine, "%10d%10d%10d", Atom_i+1, Atom_j+1, Atom_k+1);
												strcat(szTxt_Angle, szAddLine);
												Angle_List[nAngle][0] = Atom_i;	Angle_List[nAngle][1] = Atom_j;	Angle_List[nAngle][2] = Atom_k;
												nAngle++;
												if(nAngle % 3 == 0)	{
													strcat(szTxt_Angle, "\n");
												}
											}
											
											
											for(l=0; l<BondCount[Atom_k]; l++)	{	//i->j->k->l
												Atom_l = BondList[Atom_k][l];
												
												if( (Atom_l != Atom_j) && (Atom_l != Atom_i) )	{
													if(Dist_Mat[Atom_i][Atom_l] >= 3)	{
														Dist_Mat[Atom_i][Atom_l] = 3;
														Dist_Mat[Atom_l][Atom_i] = 3;
													}
													
													if( (IsLP[Atom_i]==0) && (IsLP[Atom_j]==0) && (IsLP[Atom_k]==0) && (IsLP[Atom_l]==0))	{	// skip those angles with LP involved
														if( Query_Dihedral(Atom_i, Atom_j, Atom_k, Atom_l) < 0 )	{
															sprintf(szAddLine, "%10d%10d%10d%10d", Atom_i+1, Atom_j+1, Atom_k+1, Atom_l+1);
															strcat(szTxt_Dihedral, szAddLine);
															Dihedral_List[nDihedral][0] = Atom_i;	Dihedral_List[nDihedral][1] = Atom_j;	Dihedral_List[nDihedral][2] = Atom_k;	Dihedral_List[nDihedral][3] = Atom_l;	
															nDihedral++;
															if(nDihedral % 2 == 0)	{
																strcat(szTxt_Dihedral, "\n");
															}
														}
													}
												}
											}

										}
									}
								}
							}
						}
					}
					//end	to constrcut distance matrix
					//end	to parse the entry for angles, auto angle
				}


/*
				//start	to parse the entry for dihedral
				if(FindTags_As_First_String(szLine_Upper, "DIHE"))	{
					Dihedral_Line = Extract_Entry_Per_Line(szLine_Upper, Dihedral_in_Line, Atom_Begin, Atom_End);;
					for(i=0; i<Dihedral_Line; i+=4)	{
						sprintf(szAddLine, "%10d%10d%10d%10d", 
							Dihedral_in_Line[i]+1, Dihedral_in_Line[i+1]+1, Dihedral_in_Line[i+2]+1, Dihedral_in_Line[i+3]+1);
						strcat(szTxt_Dihedral, szAddLine);
						nDihedral++;
						if(nDihedral % 2 == 0)	{
							strcat(szTxt_Dihedral, "\n");
						}
					}
				}
				//end	to parse the entry for dihedral
*/

				//start	to parse the entry for improper
				if(FindTags_As_First_String(szLine_Upper, "IMPR") || FindTags_As_First_String(szLine_Upper, "IMPH"))	{
					Improper_Line = Extract_Entry_Per_Line(szLine_Upper, Improper_in_Line, Atom_Begin, Atom_End);;
					for(i=0; i<Improper_Line; i+=4)	{
						sprintf(szAddLine, "%10d%10d%10d%10d", 
							Improper_in_Line[i]+1, Improper_in_Line[i+1]+1, Improper_in_Line[i+2]+1, Improper_in_Line[i+3]+1);
						strcat(szTxt_Improper, szAddLine);
						nImproper++;
						if(nImproper % 2 == 0)	{
							strcat(szTxt_Improper, "\n");
						}
					}
				}
				//end	to parse the entry for improper

				//start	to parse the entry for lone pair
				if(FindTags_As_First_String(szLine_Upper, "LONEPAIR"))	{
					ReadItem = sscanf(szLine_Upper, "%s %s %s %s %s %s %s %lf %s %lf %s %lf", 
						szTmp, szLPType, LP_List[nLP].szName[0], LP_List[nLP].szName[1], LP_List[nLP].szName[2], LP_List[nLP].szName[3], 
						szTmp, &(LP_List[nLP].r), szTmp, &(LP_List[nLP].theta), szTmp, &(LP_List[nLP].phi));
					if(ReadItem == 12)	{
						LP_List[nLP].Atom[0] = Get_Global_Atom_Index(LP_List[nLP].szName[0], Atom_Begin, Atom_End);
						LP_List[nLP].Atom[1] = Get_Global_Atom_Index(LP_List[nLP].szName[1], Atom_Begin, Atom_End);
						LP_List[nLP].Atom[2] = Get_Global_Atom_Index(LP_List[nLP].szName[2], Atom_Begin, Atom_End);
						LP_List[nLP].Atom[3] = Get_Global_Atom_Index(LP_List[nLP].szName[3], Atom_Begin, Atom_End);
						if(strcmp(szLPType, "BISECTOR")==0)	{
							LP_List[nLP].r = -LP_List[nLP].r;
						}
						nLP++;
					}
					else	{
						printf("Unrecogonized definition for lonepair, \n%s\nQuit\n", szLine);
						exit(1);
					}
				}
				//end	to parse the entry for lone pair

				//start	to parse the entry for lone pair
				if(FindTags_As_First_String(szLine_Upper, "ANISOTROPY"))	{
					ReadItem = sscanf(szLine_Upper, "%s %s %s %s %s %s %lf %s %lf", 
						szTmp, Aniso_List[nAniso].szName[0], Aniso_List[nAniso].szName[1], Aniso_List[nAniso].szName[2], Aniso_List[nAniso].szName[3], 
						szTmp, &(Aniso_List[nAniso].A11), szTmp, &(Aniso_List[nAniso].A22));
					if(ReadItem == 9)	{
						Aniso_List[nAniso].Atom[0] = Get_Global_Atom_Index(Aniso_List[nAniso].szName[0], Atom_Begin, Atom_End);
						Aniso_List[nAniso].Atom[1] = Get_Global_Atom_Index(Aniso_List[nAniso].szName[1], Atom_Begin, Atom_End);
						Aniso_List[nAniso].Atom[2] = Get_Global_Atom_Index(Aniso_List[nAniso].szName[2], Atom_Begin, Atom_End);
						Aniso_List[nAniso].Atom[3] = Get_Global_Atom_Index(Aniso_List[nAniso].szName[3], Atom_Begin, Atom_End);
						nAniso++;
					}
					else	{
						printf("Unrecogonized definition for ANISOTROPY, \n%s\nQuit\n", szLine);
						exit(1);
					}
				}
				//end	to parse the entry for lone pair

			}
		}

		for(i=Atom_Begin; i<Atom_End; i++)	{
			if(IsHeavyAtom[i])	{	// to add the bond between the heavy atom and its drude
				sprintf(szAddLine, "%10d%10d", i+1, i+2);
				strcat(szTxt_Bond, szAddLine);
				nBond++;
				if(nBond % 4 == 0)	{
					strcat(szTxt_Bond, "\n");
				}
			}
		}
	}
	fclose(fIn);

	if(To_Add_Test_Charge)	{
		nAtom++;
		sprintf(szAddLine, "%10d CAL      2        CAL      CAL      CALD          0.50000  40.08000           0   0.00000        0.00000\n", nAtom, nRes+1);
		strcat(szTxt_Atom, szAddLine);
	}


	for(i=0; i<nLP; i++)	{
		sprintf(szAddLine, "%10d%10d   F%10.6lf%14.3lf%14.3lf\n", 3, 4*i+1, LP_List[i].r, LP_List[i].theta, LP_List[i].phi);
		strcat(szTxt_LP, szAddLine);
	}
	for(i=0; i<nLP; i++)	{
		for(j=0; j<4; j++)	{
			sprintf(szAddLine, "%10d", LP_List[i].Atom[j]+1);
			strcat(szTxt_LP, szAddLine);
		}
		if((i+1)%2 == 0)	{
			strcat(szTxt_LP, "\n");
		}
	}

	for(i=0; i<nAniso; i++)	{
		Para_Aniso_A_K(Aniso_List[i].k11, Aniso_List[i].k22, Aniso_List[i].k33, Aniso_List[i].A11, Aniso_List[i].A22);
		sprintf(szAddLine, "      %14.4lf%14.4lf%14.4lf\n", Aniso_List[i].k11, Aniso_List[i].k22, Aniso_List[i].k33);
		strcat(szTxt_Aniso, szAddLine);
	}
	for(i=0; i<nAniso; i++)	{
		for(j=0; j<4; j++)	{
			sprintf(szAddLine, "%10d", Aniso_List[i].Atom[j]+1);
			strcat(szTxt_Aniso, szAddLine);
		}
		if((i+1)%2 == 0)	{
			strcat(szTxt_Aniso, "\n");
		}
	}


	fOut = fopen(argv[2], "w");
	fprintf(fOut, "PSF EXT CMAP CHEQ\n\n         1 !NTITLE\n* X-plor file generated by GAAMP.\n\n");
	fprintf(fOut, "%10d !NATOM\n", nAtom);
	fprintf(fOut, "%s\n", szTxt_Atom);

	fprintf(fOut, "%10d !NBOND: bonds\n", nBond);
	fprintf(fOut, "%s\n\n", szTxt_Bond);

	fprintf(fOut, "%10d !NTHETA: angles\n", nAngle);
	fprintf(fOut, "%s\n\n", szTxt_Angle);

	fprintf(fOut, "%10d !NPHI: dihedrals\n", nDihedral);
	fprintf(fOut, "%s\n\n", szTxt_Dihedral);

	fprintf(fOut, "%10d !NIMPHI: impropers\n", nImproper);
	fprintf(fOut, "%s\n\n", szTxt_Improper);

	fprintf(fOut, "         0 !NDON: donors\n\n\n");
	fprintf(fOut, "         0 !NACC: acceptors\n\n\n");
	fprintf(fOut, "         0 !NNB\n\n\n");
	fprintf(fOut, "         0         0 !NGRP NST2\n\n\n");
	fprintf(fOut, "         0         0 !NGRP NST2\n\n\n");
	fprintf(fOut, "         0 !MOLNT\n\n\n");

	fprintf(fOut, "%10d%10d !NUMLP NUMLPH\n%s\n\n", nLP, 4*nLP, szTxt_LP);
	fprintf(fOut, "%10d !NUMANISO\n%s\n\n", nAniso, szTxt_Aniso);
	fprintf(fOut, "         0 !NCRTERM: cross-terms\n\n");

	fclose(fOut);




	return 0;
}



int FindTags_As_First_String(char szBuff[], char szTag[])
{
	int ReadItem;
	char szFirstStr[256];

	ReadItem = sscanf(szBuff, "%s", szFirstStr);;
	if(ReadItem == 1)	{
		if(strcmp(szFirstStr, szTag) == 0)	{
			return 1; 
		}
	}
	return 0;
}

void To_Upper_Case(char szBuff[], char szNew[])
{
	char gap;
	int i=0;

	gap = 'A'-'a';

	while(1)	{
		if(szBuff[i]==0)	{
			break;
		}
		if( (szBuff[i]>='a') && (szBuff[i]<='z') )	{
			szNew[i] = szBuff[i] + gap;
		}
		else	{
			szNew[i] = szBuff[i];
		}
		i++;
	}
	szNew[i]=0;
}

int QueryResidue(char ResName[])
{
	int i;
	char szName_upper[256];

	To_Upper_Case(ResName, szName_upper);
	for(i=0; i<nRes_Rec; i++)	{
		if(strcmp(Res_Rec[i].ResName, szName_upper)==0)	{
			return i;
		}
	}
	printf("Fail to find the residue you inputed: %s\nQuit\n", ResName);
	exit(1);

	return -1;
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

int QueryChemName(char AtomName[])
{
	int i;

	for(i=0; i<nMass_Rec; i++)	{
		if(strcmp(Mass_Rec[i].Name, AtomName)==0)	{
			return i;
		}
	}
	printf("Fail to find the chemical name of atom %s\nQuit\n", AtomName);
	exit(1);
	return -1;
}

int Get_Alpha_Thole(char szBuff[], double& alpha, double& thole)
{
	int Len, iMax, i, ReadItem;
	char szThole[256];

	Len = strlen(szBuff);
	iMax = Len - 5;

	for(i=0; i<iMax; i++)	{
		if(strncmp(szBuff+i, "ALPHA", 5)==0)	{
			break;
		}
	}
	if(i>=iMax)	{
		return 0;
	}

	ReadItem = sscanf(szBuff+i+6, "%lf %s %lf", &alpha, szThole, &thole);
	if(ReadItem == 1)	{
		thole = 1.3;	// the default thole
		return 1;
	}
	else if(ReadItem == 3)	{
		return 1;
	}
	else	{
		printf("Find ALPHA, but fail to get its value. \n%s\nQuit\n", szBuff);
		exit(1);

		return 0;
	}
}

int SplitString(char szBuff[], char ItemList[][16])
{
	int n_Item, ReadItem, Unfinished=1;
	int iPos, nLen;
	char szItem[256];
	
	n_Item = 0;
	iPos = 0;
	while(Unfinished)	{
		while( (szBuff[iPos] == ' ') || (szBuff[iPos] == '\t') )	{	//to find the first non-blank character
			if(szBuff[iPos] == 0)	{
				Unfinished = 0;
				break;
			}
			iPos++;
		}
		if(Unfinished == 0)	{
			break;
		}
		
		ReadItem = sscanf(szBuff+iPos, "%s", szItem);
		if(ReadItem != 1)	{
			break;
		}
		strcpy(ItemList[n_Item], szItem);
		nLen = strlen(szItem);
		iPos += nLen;
		n_Item++;
	}
	
	return n_Item;
}



int Get_Global_Atom_Index(char szName[], int iBegin, int iEnd)
{
	int i;
	char szName_Upper[16];

	To_Upper_Case(szName, szName_Upper);
	for(i=iBegin; i<iEnd; i++)	{
		if(strcmp(All_Atom[i], szName_Upper)==0)	{
			return i;
		}
	}

	printf("Fail to find atom %s within [%d, %d)\nQuit\n", szName, iBegin, iEnd);
	exit(1);

	return -1;
}

int Extract_Entry_Per_Line(char szBuff[], int Entry_List[], int Atm_First, int Atm_Last)
{
	int nItem, i;
	
	nItem = SplitString(szBuff, ItemList);

	for(i=0; i<nItem; i++)	{	// stop when find !
		if(ItemList[i][0]=='!')	{
			nItem = i;
			break;
		}
	}

	for(i=1; i<nItem; i++)	{
		Entry_List[i-1] = Get_Global_Atom_Index(ItemList[i], Atm_First, Atm_Last);
	}
	return (nItem-1);
}

void Para_Aniso_A_K(double& K11, double& K22, double& K33, double A11, double A22)
{
	K11 = 1.0/A11;
	K22 = 1.0/A22;
	K33 = 1.0/(3.0-A11-A22);

	K11 *= K_DRUDE;
	K22 *= K_DRUDE;
	K33 *= K_DRUDE;

	K33 -= K_DRUDE;
	K11 -= (K_DRUDE+K33);
	K22 -= (K_DRUDE+K33);
}

int Query_Angle(int a, int b, int c)
{
	int i;

	for(i=0; i<nAngle; i++)	{
		if( (a==Angle_List[i][0]) && (b==Angle_List[i][1]) && (c==Angle_List[i][2]) )	{
			return i;
		}
		if( (c==Angle_List[i][0]) && (b==Angle_List[i][1]) && (a==Angle_List[i][2]) )	{
			return i;
		}
	}
	return (-1);	// a not existing angle
}

int Query_Dihedral(int a, int b, int c, int d)
{
	int i;

	for(i=0; i<nDihedral; i++)	{
		if( (a==Dihedral_List[i][0]) && (b==Dihedral_List[i][1]) && (c==Dihedral_List[i][2]) && (d==Dihedral_List[i][3]) )	{
			return i;
		}
		if( (d==Dihedral_List[i][0]) && (c==Dihedral_List[i][1]) && (b==Dihedral_List[i][2]) && (a==Dihedral_List[i][3]) )	{
			return i;
		}
	}
	return (-1);	// a not existing dihedral
}


void Quit_With_Error_Msg(char szMsg[])
{
        FILE *fOut;
        fOut = fopen("error.txt", "a+");
        fseek(fOut, 0, SEEK_END);
        fprintf(fOut, "Error check_lj.cpp.\nQuit\n");
        fprintf(fOut, "%s\n", szMsg);
        fclose(fOut);

        exit(1);
}

