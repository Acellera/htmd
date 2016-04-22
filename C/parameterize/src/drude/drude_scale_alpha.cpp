/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ff.h"

#define MAX_N	(8192)
#define MAX_LEN	(256)

#define MD_COULOMB  (332.0716)	//(in CHARMM)
#define K_DRUDE		(500.0)

#ifndef MINAB
#define MINAB 
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif

char szRec[MAX_N][MAX_LEN];
double cg[MAX_N], alpha[MAX_N];

int nAtom=0, nLine=0, Line_List[MAX_N];

void Read_xpsf(char szName[]);
void Scale_alpha(double scale);
void Update_charge(int Idx, double charge);
void Update_alpha(int Idx, double polarizability);
void Save_xpsf(char szName[]);

FILE *fFile_Run_Log;     // will be shared by other source code just for compiling
void Quit_With_Error_Msg(char szMsg[]);


int main(int argc, char *argv[])
{

  timebomb();

	if(argc != 3)	{
		printf("Usage: scale-alpha drude-mol.xpsf scale-factor\n");
		exit(1);
	}

	Read_xpsf(argv[1]);

	Scale_alpha(atof(argv[2]));

	Save_xpsf(argv[1]);

	return 0;
}

void Save_xpsf(char szName[])
{
	FILE *fOut;
	int i;

	fOut = fopen(szName, "w");

	for(i=0; i<nLine; i++)	{
		fprintf(fOut, "%s", szRec[i]);
	}
	
	fclose(fOut);
}

void Update_alpha(int Idx, double polarizability)
{
	char szBuff[256];
	
	sprintf(szBuff, "%10.5lf", polarizability);
	memcpy(szRec[Idx]+90, szBuff, 10);
}

void Update_charge(int Idx, double charge)
{
	char szBuff[256];
	
	sprintf(szBuff, "%10.5lf", charge);
	memcpy(szRec[Idx]+58, szBuff, 10);
}

void Scale_alpha(double scale)
{
	int i, j, Idx;
	double cg_sum, polarizability, alpha_new;

	for(i=0; i<nAtom; i++)	{
		if(fabs(alpha[i]) > 0.01)	{	// a heavy atom
			j = i+1;
			polarizability = -MD_COULOMB*cg[j]*cg[j]/(2.0*K_DRUDE);
//			printf("%d  %lf\n", i+1, polarizability-alpha[i]);
			if( fabs(polarizability-alpha[i]) > 0.02 )	{	// something wrong?
				printf("Error in alpha: %d  %lf\n", i+1, fabs(polarizability-alpha[i]));
				exit(1);
			}

			cg_sum = cg[i] + cg[j];
			alpha_new = polarizability * scale;
			alpha_new = min(2.40, alpha_new);

			cg[j] = -sqrt(2.0*K_DRUDE*fabs(alpha_new)/MD_COULOMB);
			cg[i] = cg_sum - cg[j];

			Idx = Line_List[i];	// the line of the heavy atom
			Update_alpha(Idx, alpha_new);
			Update_charge(Idx, cg[i]);
			Update_charge(Idx+1, cg[j]);
		}
	}
}

void Read_xpsf(char szName[])
{
	int ReadItem, iTmp, nDrude=0;
	FILE *fIn;
	char szLine[256], *ReadLine, szResName[256], szAtom[256], szChem[256], szTmp[256];
	double chg=0.0, polarizability, mass, thole;

	nAtom = 0;
	nLine = 0;

	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		printf("Fail to open file: %s\nQuit.\n", szName);
		exit(1);
	}

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		ReadItem = sscanf(szLine, "%d%s%d%s%s%s%lf%lf%d%lf%lf", 
			&iTmp, &szTmp, &iTmp, szResName, szAtom, szChem, &chg, &mass, &iTmp, &polarizability, &thole);

		if(ReadItem == 11)	{
			cg[nAtom] = chg;
			alpha[nAtom] = polarizability;
			Line_List[nAtom] = nLine;
			
			nAtom++;
		}

		strcpy(szRec[nLine], szLine);
		nLine++;
	}

	fclose(fIn);

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

