/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ff.h"

#define MAX_N_ATOM	(256)
#define MAX_ESP_GRID	(5000)

double x_Mol[MAX_N_ATOM], y_Mol[MAX_N_ATOM], z_Mol[MAX_N_ATOM];
char szTxtCrdFile[MAX_N_ATOM][256], szAtomInclude[256];

double Grid_x[MAX_ESP_GRID], Grid_y[MAX_ESP_GRID], Grid_z[MAX_ESP_GRID], Grid_ESP[MAX_ESP_GRID];
int netcharge=0;

char szMyPath[]="../mypath.txt";
char szQMParaFile[]="../QM-para.txt";
char szGaussianOpt[256]="mol-opt.gjf";
char szGaussianOptOut[256]="mol-opt.out";
char szCrdName[256]="mol.crd";
char szOptCrdName[256]="mol-opt.crd";
char szInpName[256];
char szXpsfName[256];
char szGjf_CGrid[256]="mol-cgrid.gjf";
char szConf_CGrid[]="ConnollyGrid.cfg";
char szExe_G09[256];
char szMol_ESP_pot0[256]="mol-esp-pot0.dat";
char szMol_ESP_pert[256]="mol-esp-pert.dat";
char szFileElemList[]="elem-list.txt";

int AtomType[MAX_N_ATOM];

int nAtom;

FILE *fFile_Run_Log;     // will be shared by other source code just for compiling
void Quit_With_Error_Msg(char szMsg[]);


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

char szQM_Level_Opt_Cmd[256]="%mem=200MB\n%nproc=4\n#P HF/6-31G* SCF=Tight opt\n\n";	// preset
char szQM_Level_ESP_Cmd[256]="%mem=200MB\n%nproc=4\n#B3LYP/aug-cc-pVDZ nosymm scf=tight \nprop=(read,field)\n\nMol\n\n0 1\n";


void GenerateCRD(void);
int DetermineAtomType(char szName[]);
void GetAtomTypeFromXpsf(void);
int FindString(char szBuff[], const char szTag[]);
void Output_Gaussian_Opt_Input(void);
void GetOptimizedStructure(void);
void Get_EXE_Path(char szExeFile[], char szExePath[]);
int To_Find_Tag(FILE *fIn, char szFileName[], const char szTag[], char szLine[], int Quit);
int Skip_N_Line(FILE *fIn, char szFileName[], int n);
void Generate_CGrid_Input(void);
int Count_Line_Num_File(char szName[]);
void Get_ESP_Result(char szESPOutput[], int& nGrid, int& nIon);

void Replace_NewLine(char szBuff[]);
int Split_Into_Two_String(char szBuff[], char str1[], char str2[]);
void Setup_QM_Level(void);



int main(int argc, char *argv[])
{

  timebomb();

	if(argc != 4)	{
		printf("Usage: gen-crd your.inp your.xpsf netcharge\nQuit\n");
		return 1;
	}

	strcpy(szInpName, argv[1]);
	strcpy(szXpsfName, argv[2]);
	netcharge = atoi(argv[3]);

	Get_EXE_Path("G09_EXE_PATH", szExe_G09);
	Setup_QM_Level();

	GenerateCRD();
	GetAtomTypeFromXpsf();

	Output_Gaussian_Opt_Input();
	GetOptimizedStructure();
	Generate_CGrid_Input();


	return 0;
}

void GenerateCRD(void)
{
	FILE *fIn, *fOut;
	char szLine[256], *ReadLine, ResName[256], AtomName[256];
	int ReadItem, i, iTmp;
	double fTmp=0.0;

	fIn = fopen(szInpName, "r");
	if(fIn == NULL)	{
		printf("Error in open file: %s\nQuit\n", szInpName);
		exit(1);
	}

	fOut = fopen(szCrdName, "w");
	while(1)	{
		if(feof(fIn))	{
			printf("Error in reading file: %s\nThe coordinate data are incomplete!\nQuit\n", szInpName);
			fclose(fIn);
			fclose(fOut);
			exit(1);
		}

		ReadLine = fgets(szLine, 256, fIn);

		if(ReadLine == NULL)	{
			printf("Error in reading file: %s\nFail to find the head of crd data!\nQuit\n", szInpName);
			fclose(fIn);
			fclose(fOut);
			exit(1);
		}
		else	{
			if(strncmp(szLine, "* Residues coordinate", 21)==0)	{
				ReadLine = fgets(szLine, 256, fIn);

				ReadLine = fgets(szLine, 256, fIn);
				ReadItem = sscanf(szLine, "%d", &nAtom);
				if(ReadItem !=1)	{
					printf("Error in reading file: %s\nFail to find the crd data!\nQuit\n", szInpName);
					fclose(fIn);
					fclose(fOut);
					exit(1);
				}
				fprintf(fOut, "%s", szLine);
				for(i=0; i<nAtom; i++)	{
					ReadLine = fgets(szLine, 256, fIn);
					ReadItem = sscanf(szLine, "%d%d%s%s%lf%lf%lf", &iTmp, &iTmp, ResName, AtomName, &(x_Mol[i]), &(y_Mol[i]), &(z_Mol[i]));
					if(ReadItem != 7)	{
						printf("Error in reading file: %s\nFail to parse the crd data!\n%s\nQuit\n", szInpName, szLine);
						fclose(fIn);
						fclose(fOut);
						exit(1);
					}
					else	{
						strcpy(szTxtCrdFile[i], szLine);	// save the crd file in to array
						fprintf(fOut, "%s", szLine);
					}
				}
				break;
			}
		}
	}

	fclose(fIn);
	fclose(fOut);
}

void GetAtomTypeFromXpsf(void)
{
	FILE *fIn;
	char AtomName[256], szAtomAdd[256];
	int ReadItem, i;
	
	fIn = fopen(szFileElemList, "r");
	if(fIn == NULL)	{
		printf("Fail to open file: elem-list.txt\nQuit\n");
		exit(1);
	}
	strcpy(szAtomInclude, "");
	for(i=0; i<nAtom; i++)	{
		ReadItem = fscanf(fIn, "%s", AtomName);
		if(ReadItem != 1)	{
			printf("Error in reading file: elem-list.txt\nQuit\n");
		}

		AtomType[i] = DetermineAtomType(AtomName);
		sprintf(szAtomAdd, " %d", i+1);

		if(strcmp(AtomName, "N")==0)	{
			strcat(szAtomInclude, szAtomAdd);
		}
		else if(strcmp(AtomName, "O")==0)	{
			strcat(szAtomInclude, szAtomAdd);
		}
		else if(strcmp(AtomName, "F")==0)	{
			strcat(szAtomInclude, szAtomAdd);
		}
		else if(strcmp(AtomName, "P")==0)	{
			strcat(szAtomInclude, szAtomAdd);
		}
		else if(strcmp(AtomName, "S")==0)	{
			strcat(szAtomInclude, szAtomAdd);
		}
		else if(strcmp(AtomName, "CL")==0)	{
			strcat(szAtomInclude, szAtomAdd);
		}
		else if(strcmp(AtomName, "BR")==0)	{
			strcat(szAtomInclude, szAtomAdd);
		}
		else if(strcmp(AtomName, "I")==0)	{
			strcat(szAtomInclude, szAtomAdd);
		}
	}

	fclose(fIn);
}

int DetermineAtomType(char szName[])
{
	int i;

	for(i=0; i<N_ELEM; i++)	{
		if( strcmp(Name_Elem[i], szName) == 0)	{
			return i;
		}
	}

	printf("Fail to identify atom: %s\n", szName);

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

void Output_Gaussian_Opt_Input(void)
{
	FILE *fOut;
	int i;

	fOut = fopen(szGaussianOpt, "w");


	fprintf(fOut, "%s", szQM_Level_Opt_Cmd);
	fprintf(fOut, "Mol \n\n%d,1 \n", netcharge);

	for(i=0; i<nAtom; i++)	{
		fprintf(fOut, "%-2s  %10.5lf%10.5lf%10.5lf\n", Name_Elem[AtomType[i]], x_Mol[i], y_Mol[i], z_Mol[i]);
	}
	fprintf(fOut, "\n\n");
	fclose(fOut);
}

void GetOptimizedStructure(void)
{
	FILE *fIn, *fOut;
	char szCmd[256], szLine[256], *ReadLine;
	int i, ReadItem;

	fIn = fopen(szGaussianOptOut, "r");
	if(fIn == NULL)	{
		sprintf(szCmd, "%s < %s > %s", szExe_G09, szGaussianOpt, szGaussianOptOut);
		system(szCmd);
		printf("Gen-Crd_file> Finished geometry optimization.\n");
	}
	else	{
		printf("Gen-Crd_file> Read previous geometry optimization.\n");
		fclose(fIn);
	}




	fIn = fopen(szGaussianOptOut, "r");
	if(fIn == NULL)	{
		printf("Fail to open file: %s\nQuit\n", szGaussianOptOut);
		exit(1);
	}

	To_Find_Tag(fIn, szGaussianOptOut, " Optimization completed", szLine, 1);
	To_Find_Tag(fIn, szGaussianOptOut, "                         Standard orientation:", szLine, 1);

	Skip_N_Line(fIn, szGaussianOptOut, 4);

	for(i=0; i<nAtom; i++)	{
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			printf("Error in reading the optimized coordinate from Gaussian output file: %s\nQuit\n", szGaussianOptOut);
			exit(1);
		}
		ReadItem = sscanf(szLine+34, "%lf %lf %lf", &(x_Mol[i]), &(y_Mol[i]), &(z_Mol[i]));
		if(ReadItem != 3)	{
			printf("ReadItem != 3\nError in reading the optimized coordinate from Gaussian output file: %s\nQuit\n", szGaussianOptOut);
			exit(1);
		}

		sprintf(szTxtCrdFile[i]+20, "%10.5lf%10.5lf%10.5lf", x_Mol[i], y_Mol[i], z_Mol[i]);	// update the coordinate with the optimized geometry
	}
	fclose(fIn);


	fOut = fopen(szOptCrdName, "w");
	fprintf(fOut, "%5d\n", nAtom);
	for(i=0; i<nAtom; i++)	{
		fprintf(fOut, "%s\n", szTxtCrdFile[i]);
	}
	fclose(fOut);
}


void Get_EXE_Path(char szExeFile[], char szExePath[])
{
	FILE *fIn;
	char szLine[256], szBuff[256], *ReadLine;
	int nLen;

	fIn = fopen(szMyPath, "r");
	if(fIn == NULL)	{
		printf("Fail to open file: %s\nQuit\n", szMyPath);
		exit(1);
	}

	nLen = strlen(szExeFile);
	while(1)	{
		if(feof(fIn))	{
			printf("Fail to find the path of %s in %s\nQuit\n", szExeFile, szMyPath);
			fclose(fIn);
			exit(1);
		}

		ReadLine = fgets(szLine, 256, fIn);

		if(ReadLine == NULL)	{
			printf("Fail to find the path of %s in %s\nQuit\n", szExeFile, szMyPath);
			fclose(fIn);
			exit(1);
		}
		else	{
			if(strncmp(szLine, szExeFile, nLen)==0)	{
				sscanf(szLine, "%s %s", szBuff, szExePath);
				break;
			}
		}
	}
	fclose(fIn);

	fIn = fopen(szExePath, "rb");
	if(fIn == NULL)	{
		printf("The file %s for %s doesn't exist.\nQuit\n", szExePath, szExeFile);
		exit(1);
	}
	else	{
		fclose(fIn);
	}
}

int To_Find_Tag(FILE *fIn, char szFileName[], const char szTag[], char szLine[], int Quit)
{
	char *ReadLine;

	while(1)	{
		if(feof(fIn))	{
			if(Quit)	{
				printf("Fail to find the tag: %s in file %s\nQuit\n", szTag, szFileName);
				fclose(fIn);
				exit(1);
			}
			else	{
				printf("Fail to find the tag: %s in file %s. Reaching the end of this file.\n", szTag, szFileName);
				break;
			}
		}

		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			if(Quit)	{
				printf("Fail to find the tag: %s in file %s\nQuit\n", szTag, szFileName);
				fclose(fIn);
				exit(1);
			}
			else	{
				printf("Fail to find the tag: %s in file %s\n", szTag, szFileName);
				break;
			}
		}
		else	{
			if(FindString(szLine, szTag) >= 0)	{
				return 1;
			}
		}
	}

	return 0;
}

int Skip_N_Line(FILE *fIn, char szFileName[], int n)
{
	int i;
	char szLine[256], *ReadLine;

	for(i=0; i<n; i++)	{
		if(feof(fIn))	{
			printf("Fail in Skip_N_Line(%s, %d).\nQuit\n", szFileName, n);
			fclose(fIn);
			exit(1);
		}

		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			printf("Fail in Skip_N_Line(%s, %d).\nQuit\n", szFileName, n);
			fclose(fIn);
			exit(1);
		}
	}
	return 1;
}

#define N_LAYER	(5)
#define N_ION_MIN	(10)
#define N_ION_MAX	(400)
#define N_GRID_MIN	(800)
#define N_GRID_MAX	(4000)

void Generate_CGrid_Input(void)
{
	FILE *fOut, *fIn;
	int i, nGrid, nGrid_ESP, n_Ion, n_Ion_Read, Ready=0, nIteration=0;
	double Rho_List_Ion[N_LAYER]={5.0, 6.0, 3.0, 0.2, 0.05};
	double Rho_List_ESP[N_LAYER]={7.0, 1.6, 1.2, 0.7, 0.4};

	char szExeCGrid[256], szCmd[256], szCGridOutput[]="pgrid.xyz", szIonPosOutput[]="qpos.all";
	char szGjf_ESP[256], szESP_Gaussian_Output[256]="cal-esp.out";

	fOut = fopen(szGjf_CGrid, "w");
	fprintf(fOut, "%s", szQM_Level_ESP_Cmd);
//	fprintf(fOut, "#B3LYP/aug-cc-pVDZ nosymm scf=tight scrf=(solvent=water) \nprop=(read,field)\n\nMol\n\n0 1\n");
	for(i=0; i<nAtom; i++)	{
		fprintf(fOut, "%-2s  %10.5lf%10.5lf%10.5lf\n", Name_Elem[AtomType[i]], x_Mol[i], y_Mol[i], z_Mol[i]);
	}
	fclose(fOut);

	Get_EXE_Path("CGRID_DRUDE_EXE_PATH", szExeCGrid);
	sprintf(szCmd, "%s %s output > cgrid.log", szExeCGrid, szGjf_CGrid);

	//start	to generate ESP grid points
	while(1)	{
		nIteration++;
		Ready = 0;
		fOut = fopen(szConf_CGrid, "w");
		fprintf(fOut, "NSURF=%d\n", 10);

		fprintf(fOut, "1.3  5.00  1.5 1.6  BI  %s\n", szAtomInclude);
		fprintf(fOut, "1.3  6.00  1.5 2.0  CBI %s\n", szAtomInclude);
		fprintf(fOut, "2.2 %6.3lf  1.5 2.0  AB  \n", Rho_List_Ion[2]);
		fprintf(fOut, "2.2 %6.3lf  1.5 2.0  BC  \n", Rho_List_Ion[3]);
		fprintf(fOut, "4.0 %6.3lf  2.5 2.0  B   \n", Rho_List_Ion[4]);


		fprintf(fOut, "1.3  %6.3lf  0.25 2.0  CI %s\n", Rho_List_ESP[0], szAtomInclude);
		fprintf(fOut, "2.2  %6.3lf  0.7  2.0  C \n", Rho_List_ESP[1]);
		fprintf(fOut, "3.0  %6.3lf  0.7  2.0  C \n", Rho_List_ESP[2]);
		fprintf(fOut, "4.0  %6.3lf  0.8  2.0  C \n", Rho_List_ESP[3]);
		fprintf(fOut, "5.0  %6.3lf  0.9  2.0  C \n", Rho_List_ESP[4]);

		fclose(fOut);

		system(szCmd);

		nGrid = Count_Line_Num_File(szCGridOutput);
		n_Ion = (int)(Count_Line_Num_File(szIonPosOutput)*1.0/(nAtom+1) + 0.1);

		if(nGrid < N_GRID_MIN)	{	// to increase the density
			for(i=0; i<N_LAYER; i++)	{
				Rho_List_ESP[i] *= 1.1;
			}
		}
		else if(nGrid > N_GRID_MAX)	{
			for(i=0; i<N_LAYER; i++)	{
				Rho_List_ESP[i] /= 1.1;
			}
		}
		else	{
			Ready++;
		}

		if(n_Ion < N_ION_MIN)	{	// to increase the density
			for(i=0; i<N_LAYER; i++)	{
				Rho_List_Ion[i] *= 1.1;
			}
		}
		else if(n_Ion > N_ION_MAX)	{
			for(i=0; i<N_LAYER; i++)	{
				Rho_List_Ion[i] /= 1.1;
			}
		}
		else	{
			Ready++;
		}

		if(Ready == 2)	{
			break;
		}
		else if(nIteration > 30)	{	// too many failed iterations
			break;
		}
	}
	//end	to generate ESP grid points

	sprintf(szGjf_ESP, "%s.all", szGjf_CGrid);

	fIn = fopen(szESP_Gaussian_Output, "r");
	if(fIn == NULL)	{
		sprintf(szCmd, "%s < %s > %s ", szExe_G09, szGjf_ESP, szESP_Gaussian_Output);
		system(szCmd);
	}
	else	{
		fclose(fIn);
	}


	Get_ESP_Result(szESP_Gaussian_Output, nGrid_ESP, n_Ion_Read);
	if(nGrid_ESP != nGrid)	{
		printf("Generate_CGrid_Input> nGrid_ESP != nGrid\nQuit\n");
		exit(1);
	}
	if(n_Ion_Read != n_Ion)	{
		printf("Generate_CGrid_Input> n_Ion_Read != n_Ion\nQuit\n");
		exit(1);
	}
}


void Get_ESP_Result(char szESPOutput[], int& nGrid, int& nIon)
{
	FILE *fIn, *fOut;
	char szLine[256], *ReadLine;
	int ReadItem, Idx_Grid, Find_Entry;
	double ESP;

	fIn = fopen(szESPOutput, "r");
	if(fIn == NULL)	{
		printf("Get_ESP_Result> Fail to open Gaussian output file, %s.\nQuit\n", szESPOutput);
		exit(1);
	}

	To_Find_Tag(fIn, szESPOutput, "            Electrostatic Properties Using The SCF Density", szLine, 1);
	Skip_N_Line(fIn, szESPOutput, nAtom+3);

	fOut = fopen(szMol_ESP_pot0, "w");

	//start	to read the unperturbed ESP and grid coordinates
	nGrid = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		else	{
			ReadItem=sscanf(szLine+31, "%lf %lf %lf", &(Grid_x[nGrid]), &(Grid_y[nGrid]), &(Grid_z[nGrid]));
			if(ReadItem == 3)	{
				nGrid++;
			}
			else	{
				break;
			}
		}
	}

	To_Find_Tag(fIn, szESPOutput, "    Center     Electric         -------- Electric Field --------", szLine, 1);
	Skip_N_Line(fIn, szESPOutput, nAtom+2);

	Idx_Grid = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		else	{
			ReadItem=sscanf(szLine+10, "%lf", &(Grid_ESP[Idx_Grid]));
			if(ReadItem == 1)	{
				Grid_ESP[Idx_Grid] /= 0.529177249;	// from Atom Unit to Angstrom
				fprintf(fOut, "%10.5lf %10.5lf %10.5lf %10.6lf\n", Grid_x[Idx_Grid], Grid_y[Idx_Grid], Grid_z[Idx_Grid], Grid_ESP[Idx_Grid]);
				Idx_Grid++;
			}
			else	{
				break;
			}
		}
	}
	fclose(fOut);
	if(Idx_Grid != nGrid)	{
		printf("Get_ESP_Result> Idx_Grid != nGrid\nQuit\n");
		exit(1);
	}
	//end	to read the unperturbed ESP and grid coordinates

	nIon = 0;
	fOut = fopen(szMol_ESP_pert, "w");
	while(1)	{
		if(feof(fIn))	{
			break;
		}

		Find_Entry = To_Find_Tag(fIn, szESPOutput, "    Center     Electric         -------- Electric Field --------", szLine, 0);
		if(Find_Entry == 0)	{
			break;
		}
		Skip_N_Line(fIn, szESPOutput, nAtom+2);

		for(Idx_Grid=0; Idx_Grid<nGrid; Idx_Grid++)	{
			ReadLine = fgets(szLine, 256, fIn);
			if(ReadLine == NULL)	{
				printf("Get_ESP_Result() > Incomplete ESP data.\nQuit\n");
				exit(1);
			}
			else	{
				ReadItem=sscanf(szLine+10, "%lf", &ESP);
				if(ReadItem == 1)	{
					ESP /= 0.529177249;	// from Atom Unit to Angstrom
					fprintf(fOut, "%10.5lf\n", ESP);
				}
				else	{
					break;
				}
			}
		}

		nIon++;
	}
	fclose(fOut);


	fclose(fIn);



	return;
}


int Count_Line_Num_File(char szName[])
{
	FILE *fIn;
	char szLine[256], *ReadLine;
	int Count=0;

	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		printf("Count_Line_Num_File> Fail to open file %s\nQuit\n", szName);
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
		else	{
			Count++;
		}
	}
	fclose(fIn);

	return Count;
}



void Setup_QM_Level(void)
{
	FILE *fIn;
	char *ReadLine, szLine[256], szSubStr_1[256], szSubStr_2[256];
	char szKey_QM_MEM[256], szKey_QM_NPROC[256], szQM_Level_Opt[256], szQM_Level_ESP_Drude[256], szQM_Level_E_Dimer[256], szQM_Level_E_Monomer[256];
	int nStr, KeyCount;

	fIn = fopen(szQMParaFile, "r");
	if(fIn == NULL)	{
		printf("Fail to open file %s\nQuit\n", szQMParaFile);
		exit(1);
	}

	KeyCount = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		Replace_NewLine(szLine);
		if(ReadLine == NULL)	{
			break;
		}
		else	{
			nStr = Split_Into_Two_String(szLine, szSubStr_1, szSubStr_2);
			if(nStr == 2)	{
				if(strcmp(szSubStr_1,"QM_MEM")==0)	{
					strcpy(szKey_QM_MEM, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_NPROC")==0)	{
					strcpy(szKey_QM_NPROC, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_OPT")==0)	{
					strcpy(szQM_Level_Opt, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_ESP_DRUDE")==0)	{
					strcpy(szQM_Level_ESP_Drude, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_E_DIMER")==0)	{
					strcpy(szQM_Level_E_Dimer, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_E_MONOMER")==0)	{
					strcpy(szQM_Level_E_Monomer, szSubStr_2);
					KeyCount++;
				}
			}
		}
	}

	fclose(fIn);


	if(KeyCount < 6)	{	// !!!
		printf("Setup_QM_Level> Error: incomplete entries in %s\nQuit\n", szQMParaFile);
		exit(1);
	}


  // MJH start
  if( getenv("NCORES" ) ) {
      printf(" Overriding QM NCPUS with [%s]\n", getenv("NCORES") );
      strcpy(szKey_QM_NPROC, getenv("NCORES") );
      //nCore_Per_Node = atoi(szKey_QM_NPROC);
  }
  if( getenv( "MEMORY" ) ) {
      printf(" Overriding QM MEMORY with [%s]\n", getenv("MEMORY") );
      strcpy(szKey_QM_MEM, getenv("MEMORY" ) );
      strcat( szKey_QM_MEM, "GB" );
  }
  // MJH stop
  //


	sprintf(szQM_Level_Opt_Cmd, "%%mem=%s\n%%nproc=%s\n%s",
		szKey_QM_MEM, szKey_QM_NPROC, szQM_Level_Opt);
	sprintf(szQM_Level_ESP_Cmd, "%%mem=%s\n%%nproc=%s\n%sMol\n\n%d 1\n",
		szKey_QM_MEM, szKey_QM_NPROC, szQM_Level_ESP_Drude, netcharge);
}

int Split_Into_Two_String(char szBuff[], char str1[], char str2[])
{
	int nLen, i, iBegin_First, Count, iBegin_Second, nLen_1, nLen_2;

	nLen = strlen(szBuff);
	str1[0] = 0;
	str2[0] = 0;

	for(i=0; i<nLen; i++)	{
		if( (szBuff[i] != ' ') && (szBuff[i] != '\t') ) 	{	// To find the first character of the first string
			break;
		}
	}

	iBegin_First = i;

	Count = 0;
	for(i=iBegin_First; i<nLen; i++)	{
		if( (szBuff[i] == ' ') || (szBuff[i] == '\t') ) 	{	// To find the last character of the first string
			break;
		}
		else	{
			str1[Count] = szBuff[i];
			Count++;
		}
	}
	str1[Count] = 0;
	nLen_1 = Count;

	for(; i<nLen; i++)	{
		if( (szBuff[i] != ' ') && (szBuff[i] != '\t')  && (szBuff[i] != 0x22) ) 	{	// To find the first character of the second string
			break;
		}
	}
	iBegin_Second = i;

	Count = 0;
	for(i=iBegin_Second; i<nLen; i++)	{
//		if( (szBuff[i] == 0x0) || (szBuff[i] == 0x0D) || (szBuff[i] == 0x0A) || (szBuff[i] == 0x22) ) 	{	// To find the last character of the second string
		if( (szBuff[i] == 0x0) || (szBuff[i] == 0x22) ) 	{	// To find the last character of the second string
			break;
		}
		else	{
			str2[Count] = szBuff[i];
			Count++;
		}
	}
	str2[Count] = 0;
	nLen_2 = Count;

	if(nLen_2 > 0)	{
		return 2;
	}
	else if(nLen_1 > 0)	{
		return 1;
	}
	else	{
		return 0;
	}
}

void Replace_NewLine(char szBuff[])
{
	int nLen, i;

	nLen =strlen(szBuff);

	for(i=1; i<nLen; i++)	{
		if( (szBuff[i-1]==0x5C) && (szBuff[i]==0x6E) )	{
			szBuff[i-1] = 0x20;
			szBuff[i] = 0x0A;
		}
	}
}

void Quit_With_Error_Msg(char szMsg[])
{
        FILE *fOut;
        fOut = fopen("error.txt", "a+");
        fseek(fOut, 0, SEEK_END);
        fprintf(fOut, "Error in generating drude ESP.\nQuit\n");
        fprintf(fOut, "%s\n", szMsg);
        fclose(fOut);

        exit(1);
}

