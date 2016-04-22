/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ff.h"


void Quit_With_Error_Msg(char szMsg[]);
void Update_CG_xpsf(void);


FILE *fFile_Run_Log;     // will be shared by other source code just for compiling


int main(void)
{
  timebomb();

	Update_CG_xpsf();


	return 0;
}

void Update_CG_xpsf(void)
{
	char szXpsfFile[]="mol.xpsf", szNewXpsfFile[]="new-mol.xpsf", szCGList[]="cg-list.txt", ErrorMsg[256];
	char szLine[256], *ReadLine, szTag_1[256], szTag_2[256];
	double CG[MAX_ATOM], fTmp=0.0;
	FILE *fIn, *fOut;
	int nAtom, ReadItem, nLine, iTmp;
	int Index, ResID, Fixed;
	char MolName[256], ResName[256], AtomName[256], ChemName[256];
	double charge, mass, last_1, last_2;

	fIn = fopen(szCGList, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file: %s\nQuit\n", szCGList);
		Quit_With_Error_Msg(ErrorMsg);
	}

	nAtom = 0;
	while(1)	{
		ReadItem = fscanf(fIn, "%lf", &(CG[nAtom]));
		if(ReadItem==1)	{
			nAtom++;
		}
		else	{
			break;
		}
	}
	fclose(fIn);
	
	fIn = fopen(szXpsfFile, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file: %s\nQuit\n", szXpsfFile);
		Quit_With_Error_Msg(ErrorMsg);
	}
	fOut = fopen(szNewXpsfFile, "w");

	nLine = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		else	{
			ReadItem = sscanf(szLine, "%d %s %d %s", &iTmp, szTag_1, &iTmp, szTag_2);
			if( (ReadItem==4) && (strcmp(szTag_1, "MOL")==0) && (strcmp(szTag_2, "MOL")==0) )	{
				ReadItem = sscanf(szLine, "%d %s %d %s%s%s%lf%lf%d%lf%lf", 
					&Index, MolName, &ResID, ResName, AtomName, ChemName, &charge, &mass, &Fixed, &last_2, &last_1);
				if(ReadItem == 11)	{
				}
				else	{
					fclose(fIn);
					fclose(fOut);
					sprintf(ErrorMsg, "Error in read file %s\n%s\nQuit\n", szXpsfFile, szLine);
					Quit_With_Error_Msg(ErrorMsg);
				}
//				fprintf(fOut, "%8d %-4s%2d    %-5s%-5s%-5s%10.6lf     %9.4lf  %10d%10.5lf %.6E\n", 
//					Index, MolName, ResID, ResName, AtomName, ChemName, CG[nLine], mass, Fixed, last_2, last_1);
				fprintf(fOut, "%10d %-6s %3d        %-8s %-8s %-8s %12.5lf%10.5lf           0%10.5lf %14.5lf\n", 
					Index, MolName, ResID, ResName, AtomName, ChemName, CG[nLine], mass, last_2, last_1);
				nLine++;
			}
			else	{
				fprintf(fOut, "%s", szLine);
			}
		}
	}

	fclose(fIn);
	fclose(fOut);

	if(nLine != nAtom)	{
		sprintf(ErrorMsg, "Fatal error: nLine != nAtom\nQuit\n");
		Quit_With_Error_Msg(ErrorMsg);
	}
}

void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in update-xpsf.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}
