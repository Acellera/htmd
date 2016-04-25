/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include "ff.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAX_LINE	(2048)
#define MAX_LEN		(256)

char szTxtPrm[MAX_LINE][MAX_LEN];
char szRtfFile[MAX_LEN];
char szPrmFile[MAX_LEN];

void Add_Tip3_Rtf_File(void);
void Add_Tip3_Prm_File(void);
void To_Upper_Case(char szBuff[]);
int FindString(char szBuff[], const char szTag[]);
void Quit_With_Error_Msg(char szMsg[]);

int n_Line;
int Line_Insert_Bond, Line_Insert_Angle, Line_Insert_Impropers;
int Found_Improper;
double E14FAC=1.0;

int net_charge;

FILE *fFile_Run_Log;     // will be shared by other source code
///void Quit_With_Error_Msg(char szMsg[]);


int main(int argc, char *argv[])
{
	if(argc < 4)	{
		printf("Usage: add-tip3 your-rtf your-prm net_charge E14FAC\nQuit\n");
		return 1;
	}

  timebomb();

	strcpy(szRtfFile, argv[1]);
	strcpy(szPrmFile, argv[2]);
	net_charge = atoi(argv[3]);
	if(argc == 5)	{
		E14FAC = atof(argv[4]);
	}

	Add_Tip3_Rtf_File();

	Add_Tip3_Prm_File();


	return 0;
}

void Add_Tip3_Prm_File(void)
{
	FILE *fOut;
	char *ReadLine, ErrorMsg[256];
	int i;

	Found_Improper = 0;
	fOut = fopen(szPrmFile, "r");
	if(fOut == NULL)	{
		sprintf(ErrorMsg, "Error in open file: %s\nQuit\n", szPrmFile);
		Quit_With_Error_Msg(ErrorMsg);
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
				else if(strncmp(szTxtPrm[n_Line], "IMPHI", 5)==0)	{
					Found_Improper = 1;
					sprintf(szTxtPrm[n_Line], "IMPROPERS\n");			// rename "IMPHI" to "IMPROPERS"
				}
				else if(strncmp(szTxtPrm[n_Line], "IMPROPER", 8)==0)	{
					Found_Improper = 1;
					sprintf(szTxtPrm[n_Line], "IMPROPERS\n");			// rename "IMPHI" to "IMPROPERS"
				}
				else if(strncmp(szTxtPrm[n_Line], "NONBONDED", 9)==0)	{
					Line_Insert_Impropers = n_Line;
				}

				if(strncmp(szTxtPrm[n_Line], "END", 3)==0)	{
					continue;
				}

				n_Line++;
				if(n_Line >= MAX_LINE)	{
					sprintf(ErrorMsg, "n_Line >= MAX_LINE\nQuit\n");
					fclose(fOut);
					Quit_With_Error_Msg(ErrorMsg);
				}
			}
		}
		fclose(fOut);
	}

	fOut = fopen(szPrmFile, "w");
	for(i=0; i<n_Line; i++)	{
		if(i == Line_Insert_Bond)	{
			fprintf(fOut, "HTXW   HTXW      0.000     1.5139 ! FROM TIPS3P GEOMETRY (FOR SHAKE/W PARAM)\n");
			fprintf(fOut, "OTXW   HTXW    450.000     0.9572 ! FROM TIPS3P GEOM\n\n");
		}
		else if(i == Line_Insert_Angle)	{
			fprintf(fOut, "HTXW   OTXW   HTXW     55.000   104.5200 ! TIP3P GEOMETRY\n\n");
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
	fprintf(fOut, "\nHTXW     0.000000  -0.046000     0.224500 ! TIP3P HYDROGEN PARAMETERS, adm jr., NBFIX obsolete\n");
	fprintf(fOut, "OTXW     0.000000  -0.152100     1.768200 ! TIP3P OXYGEN PARAMETERS, adm jr., NBFIX obsolete\n\n");
//	fprintf(fOut, "\n\nEND\n");

	fclose(fOut);
}

void Extract_Comment(char szBuff[], char szComment[])
{
	int i=0;

	while(szBuff[i])	{
		if(szBuff[i] == '!')	{
			break;
		}
		i++;
	}

	if(szBuff[i] == 0)	{
		strcpy(szComment, "\n");
	}
	else	{
		strcpy(szComment, szBuff+i);
	}
}

char szTxt_Rtf[MAX_LINE*MAX_LEN];
void Add_Tip3_Rtf_File(void)
{
	FILE *fOut, *fIn;
	char ErrorMsg[256], szLine[256], *ReadLine, szBuff[256], szComment[256];
	int ReadItem;
	double net_cg_local=0.0;

	szTxt_Rtf[0] = 0;
	fIn = fopen(szRtfFile, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Error in open file: %s\nQuit\n", szRtfFile);
		Quit_With_Error_Msg(ErrorMsg);
	}
	else	{
		while(1)	{
			if(feof(fIn))	{
				break;
			}
			ReadLine = fgets(szLine, 256, fIn);
			if(ReadLine == NULL)	{
				break;
			}
			else	{
				if(strncmp(szLine,"RESI", 4)==0)	{
					ReadItem = sscanf(szLine, "%s%s%lf", szBuff, szBuff, &net_cg_local);
					if(ReadItem == 3)	{
						if(fabs(net_cg_local-net_charge) > 1.0E-1)	{
							Quit_With_Error_Msg("The charge in rtf is NOT consistent with the charge provided by user.\nQuit\n");
						}
						Extract_Comment(szLine, szComment);
						sprintf(szLine, "RESI  MOL  %6.3lf %s", net_cg_local, szComment);
					}
					else	{
						ReadItem = sscanf(szLine, "%s%lf", szBuff, &net_cg_local);
						if(ReadItem == 2)	{
							if(fabs(net_cg_local-net_charge) > 1.0E-1)	{
								Quit_With_Error_Msg("The charge in rtf is NOT consistent with the charge provided by user.\nQuit\n");
							}
							Extract_Comment(szLine, szComment);
							sprintf(szLine, "RESI  MOL  %6.3lf %s", net_cg_local, szComment);
						}
					}
				}
				if(strncmp(szLine, "END", 3)==0)	{
					continue;
				}
				else if(strncmp(szLine, "ANGL", 4)==0)	{
					continue;
				}
				else if(strncmp(szLine, "DIHE", 4)==0)	{
					continue;
				}
				strcat(szTxt_Rtf, szLine);
			}
		}
		fclose(fIn);
	}

	fOut = fopen(szRtfFile, "w");
	fprintf(fOut, "%s", szTxt_Rtf);

	fprintf(fOut, "\n\nMASS 101   HTXW    1.008000 H ! TIPS3P WATER HYDROGEN\n");
	fprintf(fOut, "MASS 102   OTXW   15.999400 O ! TIPS3P WATER OXYGEN\n\n");

	fprintf(fOut, "RESI TIP3         0.000 ! tip3p water model, generate using noangle nodihedral\n");
	fprintf(fOut, "GROUP\nATOM OH2  OTXW     -0.834\nATOM H1   HTXW      0.417\nATOM H2   HTXW      0.417\nBOND OH2 H1 OH2 H2 H1 H2    ! the last bond is needed for shake\nANGLE H1 OH2 H2             ! required\n");

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

void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in add_tip3.cpp.\nQuit\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}

