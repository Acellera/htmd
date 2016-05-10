/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ff.h"

#define MAX_LEN	(512)
#define MAXLINE	(2048)

char szTxtPrm[MAXLINE][MAX_LEN];

void To_Upper_Case(char szBuff[]);

FILE *fFile_Run_Log;     // will be shared by other source code just for compiling
void Quit_With_Error_Msg(char szMsg[]);


int main(int argc, char *argv[])
{

  timebomb();

	FILE *fOut;
	char *ReadLine, szAtomType[256];
	int i, n_Line, Idx_NBond=100000, ReadItem;
	double LJ_Para[3];

	fOut = fopen(argv[1], "r");
	if(fOut == NULL)	{
		printf("Error in open file: %s\nQuit\n", argv[1]);
		exit(1);
	}
	else	{
		n_Line = 0;

		while(1)	{
			if(feof(fOut))	{
				break;
			}
			ReadLine = fgets(szTxtPrm[n_Line], MAX_LEN, fOut);
			if(ReadLine == NULL)	{
				break;
			}
			else	{
				To_Upper_Case(szTxtPrm[n_Line]);

				if(strncmp(szTxtPrm[n_Line], "NONBONDED", 9)==0)	{
					Idx_NBond = n_Line;
				}

				if(n_Line > Idx_NBond)	{	// the entry for LJ parameters
					ReadItem = sscanf(szTxtPrm[n_Line], "%s %lf %lf %lf", szAtomType, &(LJ_Para[0]), &(LJ_Para[1]), &(LJ_Para[2]));
					if(ReadItem == 4)	{
//						if( (strncmp(szAtomType, "HO", 2)==0) && (LJ_Para[0]==0.0) && (LJ_Para[1]==0.0) && (LJ_Para[2]==0.0) ) 	{	//
						if( (strncmp(szAtomType, "HO", 2)==0) && (LJ_Para[0]==0.0) && (LJ_Para[2]==0.0) ) 	{	//
//							sprintf(szTxtPrm[n_Line], "%-8s0.00   -0.0100    0.4000\n", szAtomType);
							sprintf(szTxtPrm[n_Line], "%-8s 0.00   -0.0100    0.5500\n", szAtomType);
						}
						else if( (strncmp(szAtomType, "HO", 2)==0) && (LJ_Para[0]==0.0) && (LJ_Para[2]==0.4) )	{
							sprintf(szTxtPrm[n_Line], "%-8s 0.00   -0.0100    0.5500\n", szAtomType);
						}

//						if( (strncmp(szAtomType, "HN", 2)==0) && ( fabs(LJ_Para[2]-0.2245)<0.02) ) 	{	//
//							sprintf(szTxtPrm[n_Line], "%-8s 0.00   -0.0157    0.6000\n", szAtomType);
//						}

					}
				}

				n_Line++;

				if(n_Line >= MAXLINE)	{
					printf("n_Line >= MAX_LINE\nQuit\n");
					fclose(fOut);
					exit(1);
				}
			}
		}
		fclose(fOut);
	}

	fOut = fopen(argv[1], "w");
	for(i=0; i<n_Line; i++)	{
		fprintf(fOut, "%s", szTxtPrm[i]);
	}
	fclose(fOut);

	return 0;
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

