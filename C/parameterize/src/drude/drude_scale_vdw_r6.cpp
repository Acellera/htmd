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

#ifndef MINAB
#define MINAB 
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif

char szRec[MAX_N][MAX_LEN];
double k_Emin, k_Rmin;

int nLine=0;

void Read_Prm(char szName[]);
void Save_Prm(char szName[]);

FILE *fFile_Run_Log;     // will be shared by other source code just for compiling
void Quit_With_Error_Msg(char szMsg[]);


int main(int argc, char *argv[])
{

  timebomb();

	double Coeff_R6=1.0;

	if(argc != 3)	{
		printf("Usage: scale-vdw-r6 ff.str scale-factor\n");
		exit(1);
	}

	Coeff_R6 = atof(argv[2]);
	k_Emin = Coeff_R6 * Coeff_R6;
	k_Rmin = pow(Coeff_R6, -1.0/6.0);


	Read_Prm(argv[1]);

	Save_Prm(argv[1]);

	return 0;
}

void Save_Prm(char szName[])
{
	FILE *fOut;
	int i;

	fOut = fopen(szName, "w");

	for(i=0; i<nLine; i++)	{
		fprintf(fOut, "%s", szRec[i]);
	}
	
	fclose(fOut);
}

void Read_Prm(char szName[])
{
	int ReadItem, nDrude=0, FoundVDW=0;
	FILE *fIn;
	char szLine[256], *ReadLine, szChem[256];
	double Emin, Rmin, fTmp, fTmp_14, Emin_14, Rmin_14;

	nLine = 0;

	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		printf("Fail to open file: %s\nQuit.\n", szName);
		exit(1);
	}

	while(1)	{	// to find the entry of VDW parameters
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		if(strncmp(szLine, "NONBONDED", 9) == 0)	{
			FoundVDW = 1;
		}

		strcpy(szRec[nLine], szLine);
		nLine++;

		if(FoundVDW)	{
			break;
		}
	}

	if(FoundVDW == 0)	{
		printf("Fail to find the entry of VDW parameters in %s\nQuit\n", szName);
		fclose(fIn);
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

		ReadItem = sscanf(szLine, "%s%lf%lf%lf%lf%lf%lf", 
			szChem, &fTmp, &Emin, &Rmin, &fTmp_14, &Emin_14, &Rmin_14);

		if( (strcmp(szChem, "HDW")==0) || (strcmp(szChem, "ODW")==0) || (strcmp(szChem, "D*")==0) || (strcmp(szChem, "DRUD")==0) || (strcmp(szChem, "LP")==0) || (strcmp(szChem, "CALD")==0) )	{
		}
		else	{
			if(ReadItem == 7)	{
				if(fabs(fTmp) <1.0E-10)	{
					Emin *= k_Emin;
					Rmin *= k_Rmin;
					sprintf(szLine, "%-8s0.00%11.5lf %11.6lf       0.00%11.5lf %11.6lf\n", 
						szChem, Emin, Rmin, Emin, Rmin);
				}
			}
			else if(ReadItem == 4)	{
				if(fabs(fTmp) <1.0E-10)	{
					Emin *= k_Emin;
					Rmin *= k_Rmin;
					sprintf(szLine, "%-8s0.00%11.5lf %11.6lf\n", 
						szChem, Emin, Rmin);
				}
			}
		}


		strcpy(szRec[nLine], szLine);	// 
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

