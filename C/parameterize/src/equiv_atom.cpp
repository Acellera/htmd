/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include"ff.h"

#define MAXPATHATOMNUM 5000
#define MAX_RESP_ATOM 1000
#define MAX_LEN	(256)
#define MAX_LINE_RTF	(512)
#define MAX_ATOM	(512)
#define MAX_BOND	(1024)

typedef struct {
	int atomicnum;
	int con[6];
} ATOM;

FILE *fFile_Run_Log;     // will be shared by other source code just for compiling


void Quit_With_Error_Msg(char szMsg[]);
int FindTags_As_First_String(char szBuff[], char szTag[]);
void Find_String_in_File(FILE *fRead, char szTag[]);
void Read_Input(char szName_Input[]);
int Extract_Atom_Info();
int Extract_Bond_Info();
int Split_Str_One_Line(char szBuff[]);
int Query_Atom_Name(char szAtom[]);
void equatom(void);
void Setup_Atom_Bond_Info(void);
void Setup_Atom_Element_Info(void);	// icnum
int DetermineAtomType_from_Mass(double mass);
int From_Chem_Name_to_Elem_Index(char szName[]);
void Reorganized_Equivalency(void);

//start	data for parsing one line
int nItem_Line;
char szItemLine[64][128];
//end	data for parsing one line

int selectindex[MAX_ATOM];
int equatomno[MAX_ATOM];
int pathnum[MAX_ATOM];
int selectelement[MAX_ATOM];
int pathatomnum[MAX_ATOM];
int maxlength = -1;
double *pathscore[MAX_RESP_ATOM];


ATOM atom[MAX_ATOM];
int pathnumindex = 0;
int atomindex = 0;
int atomnum = 0;
int selectnum = 0;

char szRtf[MAX_LINE_RTF][MAX_LEN];
int nAtom, nLineRec, nBond;

int nChemType=0;
int BondList[MAX_BOND][2];
int AtomBond[MAX_ATOM], Bond_Array[MAX_ATOM][6];
char szAtomName[MAX_ATOM][16];
char szChemName[MAX_ATOM][16];
char szChemName_Lib[MAX_ATOM][16];
int ElemIndex_List[MAX_ATOM];

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

int main(int argc, char *argv[])
{

  timebomb();

	int i;
	
	if(argc != 2)	{
		printf("Usage: equiv_atom mol.rtf\n");
		exit(1);
	}

	Read_Input(argv[1]);
	Extract_Atom_Info();
	Extract_Bond_Info();

	atomnum = nAtom;
	Setup_Atom_Bond_Info();

	for (i = 0; i < atomnum; i++) {
		pathatomnum[i] = MAXPATHATOMNUM;
		pathscore[i] = (double *) calloc(pathatomnum[i], sizeof(double));
		if (pathscore == NULL) {
			fprintf(stderr, "memory allocation error for *pathscore[%d]\n", i+1);
			exit(1);
		}
	}


	equatom();

	for (i = 0; i < atomnum; i++) {
		free(pathscore[i]);
	}

	Reorganized_Equivalency();

	return 0;
}

void Reorganized_Equivalency(void)
{
	FILE *fIn, *fOut;
	int i, j, Representative[MAX_ATOM], Output[MAX_ATOM], nAtom, Index, ReadItem;


	fIn = fopen("equiv.txt", "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open file equiv.txt.\nQuit\n");
	}

	memset(Output, 0, sizeof(int)*MAX_ATOM);

	fOut = fopen("equiv-org.txt", "w");

	nAtom = 0;
	while(1)	{
		ReadItem = fscanf(fIn, "%d %d", &Index, &(Representative[nAtom]));
		if(ReadItem == 2)	{
			if(Representative[nAtom] > 0)	{	// a valid new entry for atom equivalency
				Representative[nAtom]--;
				Output[Representative[nAtom]] = 1;
			}
			else	{
				Representative[nAtom] = -1;
			}
			nAtom++;
		}
		else	{
			break;
		}
	}
	fclose(fIn);

	for(i=0; i<nAtom; i++)	{
		if(Output[i] == 1)	{
			fprintf(fOut, "equivalent  %3d", i+1);
			for(j=0; j<nAtom; j++)	{
				if(j==i)	{
					continue;
				}
				if(Representative[j] == i)	{
					fprintf(fOut, " %3d", j+1);
				}
			}
			fprintf(fOut, "\n");
		}
	}



	fclose(fOut);
}


void Setup_Atom_Bond_Info(void)
{
	int i, j, nBond;

	for(i=0; i<nAtom; i++)	{
		nBond = AtomBond[i];
		for(j=0; j<nBond; j++)	{
			atom[i].con[j] = Bond_Array[i][j];
		}
		for(; j<6; j++)	{
			atom[i].con[j] = -1;
		}

		atom[i].atomicnum = From_Chem_Name_to_Elem_Index(szChemName[i]);
	}
}

void Read_Input(char szName_Input[])
{
	FILE *fIn;
	char szError[256], szLine[256], *ReadLine, szResName[24]="RESI MOL              ", szBuff[256], szTmp[256];
	int ToSave, ReadItem, ChemID, ElemIdx;
	double mass=0.0;

	fIn = fopen(szName_Input, "r");
	if(fIn == NULL)	{
		sprintf(szError, "Fail to open file %s\nQuit\n", szName_Input);
		Quit_With_Error_Msg(szError);
	}

	nChemType = 0;
	// to find "RESI "
	while(1)	{
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			fclose(fIn);
			Quit_With_Error_Msg("Fail to find the residue information.\nQuit\n");
		}
		else	{
			if( (FindTags_As_First_String(szLine, "mass")==1) || (FindTags_As_First_String(szLine, "MASS")==1))	{
				ReadItem = sscanf(szLine, "%s%d%s%lf", szTmp, &ChemID, szBuff, &mass);
				if(ReadItem == 4)	{	// a valid chem type
					ElemIdx = DetermineAtomType_from_Mass(mass) + 1;
					ElemIndex_List[nChemType] = ElemIdx;
					strcpy(szChemName_Lib[nChemType], szBuff);
					nChemType++;
				}
			}

			if(FindTags_As_First_String(szLine, "RESI")==1)	{
				break;
			}
		}
	}

	nLineRec = 1;
	memcpy(szLine, szResName, 19);
	strcpy(szRtf[0], szLine);
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		ToSave = 0;

		if(FindTags_As_First_String(szLine, "ATOM")==1)	{	// to 
			ToSave = 1;
		}
		else if(FindTags_As_First_String(szLine, "BOND")==1)	{
			ToSave = 1;
		}
		else if(FindTags_As_First_String(szLine, "DOUBLE")==1)	{
			ToSave = 1;
		}
		else if(FindTags_As_First_String(szLine, "IMPR")==1)	{
			ToSave = 1;
		}

		// auto angle and dihedral

		if(ToSave)	{
			strcpy(szRtf[nLineRec], szLine);
			nLineRec++;
		}

		if(FindTags_As_First_String(szLine, "END")==1)	{
			strcpy(szRtf[nLineRec], "\nEND\n");
			nLineRec++;
			break;
		}

	}


	fclose(fIn);
}

void scorepath(ATOM atm[], int selectnum, int startnum)
{
	int i, j, k;
	int start;
	double score;
	int resetindex;
	start = -1;
	resetindex = -1;
	selectindex[startnum] = selectnum;
	selectelement[selectnum++] = startnum;
	if(maxlength != -1 && selectnum > maxlength)
		return;
	for (i = 0; i < 6; i++) {
		if (atm[startnum].con[i] == -1) {
			score = 0.0;
			for (j = 0; j < selectnum; j++) {
				/* printf("%5d", selectelement[j]); */
				score +=
					(j + 1) * 0.11 +
					atom[selectelement[j]].atomicnum * 0.08;
			}
			pathscore[atomindex][pathnumindex++] = score;
			if (pathnumindex >= pathatomnum[atomindex]) {
				pathatomnum[atomindex] += MAXPATHATOMNUM;
				fprintf(stderr, "\nInfo: the number of the path atoms exceeds MAXPATHATOMNUM(%d) for atom[%d],extend the size and reallocate the memory automatically", pathatomnum[atomindex], atomindex);
				pathscore[atomindex] = (double *) realloc(pathscore[atomindex], pathatomnum[atomindex] * sizeof(double));
				if (pathscore[atomindex] == NULL) {
					fprintf(stderr, " reallocate memory for pathscore[%d] failed\n", atomindex);
					exit(1);
				}
			}
			return;
		}
		start = atm[startnum].con[i];
		for (k = 0; k < selectnum; k++)
			if (start == selectelement[k]) {
				resetindex = 1;
				break;
			}
		if (resetindex == 1) {
			resetindex = -1;
			continue; 
		}
		if (start == -1)
			return;
		scorepath(atm, selectnum, start);
		/* we have already visited this atom */
	}
}

void sort(double array[], int elemnum)
{
	int i, j;
	double tmp;
	for (i = 0; i < elemnum; i++)
		for (j = i + 1; j < elemnum; j++) {
			if (array[j] < array[i]) {
				tmp = array[i];
				array[i] = array[j];
				array[j] = tmp;
			}
		}
}

void equatom(void)
{
	FILE *fOut;
	int i, j, k;
	int equindex;
	double sum;
	for (i = 0; i < atomnum; i++)
		equatomno[i] = -1;
	for (i = 0; i < atomnum; i++) {
		selectnum = 0;
		pathnumindex = 0;
		atomindex = i;
		for (j = 0; j < atomnum; j++) {
			selectindex[j] = -1;
			selectelement[i] = -1;
		}
		scorepath(atom, 0, i);
		pathnum[i] = pathnumindex;
	}
	selectnum = 0;
	for (i = 0; i < atomnum; i++) {
		sum = 0.0;
		for (j = 0; j < pathnum[i]; j++) {
			sum += pathscore[i][j];
		}
	}
	for (i = 0; i < atomnum; i++) {
		sort(pathscore[i], pathnum[i]);
	}
	for (i = 0; i < atomnum; i++) {
		for (j = i + 1; j < atomnum; j++) {
			if (equatomno[j] != -1)
				continue;
			equindex = 1;
			if (pathnum[i] != pathnum[j])
				continue;
			for (k = 0; k < pathnum[i]; k++)
				if (pathscore[i][k] != pathscore[j][k]) {
					equindex = -1;
					break;
				}
			if (equindex == 1)
				equatomno[j] = i;
		}
	}
	
	fOut = fopen("equiv.txt", "w");
  for(i=0;i<atomnum;i++) fprintf(fOut, "%5d%5d\n", i+1, equatomno[i]+1);
  fclose(fOut);
}

void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in extract_mini_ff.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
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

void Find_String_in_File(FILE *fRead, char szTag[])
{
	char szLine[256], *ReadLine, szError[256];
	int nLen;

	nLen = strlen(szTag);
	while(1)	{
		if(feof(fRead))	{
			sprintf(szError, "Fail to find tag: %s\nQuit\n", szTag);
			fclose(fRead);
			Quit_With_Error_Msg(szError);
		}
		ReadLine = fgets(szLine, 256, fRead);
		if(ReadLine == NULL)	{
			sprintf(szError, "Fail to find tag: %s\nQuit\n", szTag);
			fclose(fRead);
			Quit_With_Error_Msg(szError);
		}
		else	{
			if(strncmp(szLine, szTag, nLen)==0)	{
				return;
			}
		}
	}
}


int Extract_Bond_Info()
{
	int i, j, nItem, iAtom_1, iAtom_2;

	nBond = 0;

	memset(AtomBond, 0, sizeof(int)*nAtom);

	for(i=0; i<nLineRec; i++)	{
		if( (FindTags_As_First_String(szRtf[i], "BOND")==1) || (FindTags_As_First_String(szRtf[i], "DOUBLE")==1) )	{
			nItem = Split_Str_One_Line(szRtf[i]);
			
			for(j=1; j<nItem; j+=2)	{
				iAtom_1 = Query_Atom_Name(szItemLine[j]);
				iAtom_2 = Query_Atom_Name(szItemLine[j+1]);

				BondList[nBond][0] = iAtom_1;
				BondList[nBond][1] = iAtom_2;

				Bond_Array[iAtom_1][AtomBond[iAtom_1]] = iAtom_2;
				AtomBond[iAtom_1] ++;

				Bond_Array[iAtom_2][AtomBond[iAtom_2]] = iAtom_1;
				AtomBond[iAtom_2] ++;

				nBond++;
			}
		}
	}

	return nBond;
}

int Extract_Atom_Info()
{
	int i, ReadItem;
	double chg=0.0;
	char szTmp[256];

	nAtom = 0;
	for(i=0; i<nLineRec; i++)	{
		if(FindTags_As_First_String(szRtf[i], "ATOM")==1)	{
			ReadItem = sscanf(szRtf[i], "%s%s%s%lf", szTmp, szAtomName[nAtom], szChemName[nAtom], &chg);
			if(ReadItem == 4)	{
				nAtom++;
			}
			else	{
				printf("Error in line: %s\nQuit\n", szRtf[i]);
				exit(1);
			}
		}
	}

	return nAtom;
}


int Split_Str_One_Line(char szBuff[])
{
	int iLen, iPos, iPos_End, CountChange=0, iLen_Atom_Name;
	char c, szBak[256];

	nItem_Line = 0;

	iLen = strlen(szBuff);
	strcpy(szBak, szBuff);

	iPos = 0;
	while(1)	{
		while(1)	{	// to find the beginning of a string
			c = szBuff[iPos];
			if(c == '!')	{
				break;
			}
			if( (c != 0xA) && (c != 0xD) && (c != 0x20)  && (c != '\t') && (c != 0) )	{
				break;
			}
			iPos++;
			if(iPos >= iLen)	{
				break;
			}
		}
		if(iPos >= iLen)	{
			break;
		}
		if(c == '!')	{
			break;
		}

		iPos_End = iPos;
		while(1)	{	// to find the end of a string
			c = szBuff[iPos_End];
			if( (c == 0xA) || (c == 0xD) || (c == 0x20)  || (c == '\t') || (c == 0) || (c == '!') )	{
				break;
			}
			iPos_End++;
			if(iPos >= iLen)	{
				break;
			}
		}

		iLen_Atom_Name = iPos_End-iPos;
		memcpy(szItemLine[nItem_Line], szBuff+iPos, iLen_Atom_Name);
		szItemLine[nItem_Line][iLen_Atom_Name] = 0;	// end
		nItem_Line++;


		iPos = iPos_End;
	}

	return nItem_Line;
}


int Query_Atom_Name(char szAtom[])
{
	int i;

	for(i=0; i<nAtom; i++)	{
		if(strcmp(szAtom, szAtomName[i])==0)	{
			return i;
		}
	}

	printf("Fail to find atom: %s\nQuit\n", szAtom);
	exit(1);

	return -1;
}

int From_Chem_Name_to_Elem_Index(char szName[])
{
	int i;
	char ErrorMsg[256];

	for(i=0; i<nChemType; i++)	{
		if(strcmp(szName, szChemName_Lib[i])==0)	{
			return ElemIndex_List[i];
		}
	}

	sprintf(ErrorMsg, "Fail to identify atom: %s.\n", szName);

	Quit_With_Error_Msg(ErrorMsg);

	return -1;
}

int DetermineAtomType_from_Mass(double mass)
{
	int i;
	char ErrorMsg[256];

	for(i=0; i<N_ELEM; i++)	{
		if( fabs(mass - Mass_List[i]) < 0.1 )	{
			return i;
		}
	}

	sprintf(ErrorMsg, "Fail to identify atom with mass = %8.3lf\n", mass);

	Quit_With_Error_Msg(ErrorMsg);

	return -1;
}

