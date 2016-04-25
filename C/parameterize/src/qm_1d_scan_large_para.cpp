/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
//#include <mpi.h>
#include "ff.h"

#ifndef WIN32
#include <unistd.h>
#else
#include <windows.h>
#endif

char szForceFiled[]="mol.prm";
char szXpsfFile[]="mol.xpsf";
char szCrdFile[]="mol-opt.crd";
char szPhiToScan[]="soft-dih-list.txt";
char szExe_G09[PATH_MAX];
char szQMParaFile[]="../QM-para.txt";
char szMyPath[]="../mypath.txt";
char szFileElem[]="elem-list.txt";

//const char szGAUSS_SCRDIR[]="GAUSS_SCRDIR";
//char szGAUSS_SCRDIR_Base[PATH_MAX];

char szQM_Level_1D_Scan_Cmd[512]="%mem=1600MB\n%nproc=4\n# HF/6-31G* SCF=Tight nosymm opt=ModRedundant\n\nMol dihedral scan\n\n0 1\n";	// kbt

#define MAX_SCAN	(96)
#define N_MAX_DIH	(32)
#define MAX_N_STATE		(6)
#define MAX_BIN_PHI	(128)
#define BIN_SIZE	(10)
#define BIN_SIZE_QM	(10)
#define N_MAX_QM_PART ( ( 1 + ( 360/BIN_SIZE)) * N_MAX_DIH )


#define QM_LEVEL_HF		(1)
#define QM_LEVEL_MP2	(2)

#define HARTREE_KCAL (627.509)


CMol Mol;
CForceField ForceField;
int nDihedral, QM_Level;

FILE *fFile_Run_Log;
extern int FindString(char szBuff[], char szTag[]);
void Quit_With_Error_Msg(char szMsg[]);
void Replace_D_E(char szStr[]);
int To_Find_Tag(FILE *fIn, char szFileName[], char szTag[], char szLine[]);
int Skip_N_Line(FILE *fIn, char szFileName[], int n);
double Get_Energy(char szLine[]);
void ReadElementList(void);
void Determine_QM_Level(char szQMLevel[]);

int n_Phi=0, Counter=0;
int DihList[N_MAX_DIH][4], IdxDihSelect[N_MAX_DIH], n_State_List[N_MAX_DIH], n_State_List_Save[N_MAX_DIH];
double Phi_Set[N_MAX_DIH][MAX_N_STATE], Phi_Set_Save[N_MAX_DIH][MAX_N_STATE], Phi_To_Set[N_MAX_DIH], Phi_To_Set_Scan[N_MAX_DIH];
double x_Save[MAX_ATOM], y_Save[MAX_ATOM], z_Save[MAX_ATOM];
double E_Phi[MAX_BIN_PHI], E_Scan[MAX_BIN_PHI];

int n_QM_Part=0;
char szQM_Scan[N_MAX_QM_PART][PATH_MAX];
double Phi_QM_Start[N_MAX_QM_PART];
char szElemName[MAX_ATOM][8];
char szOptPdbSave[PATH_MAX];

int	ProgID=0, nProc=1;



int Read_Soft_DihedralList(void);
int Get_Dihedral_Index(int ia, int ib, int ic, int id);
void Enumerate_Dihedrals(int IdxDih, int iState);
void BackupCoordinates(void);
void RestoreCoordinates(void);
double Get_Barrier(double E[], int n, double& E_Min, double& Phi);
void MM_Fixed_1D_Scan(int ReadInPdb);
void SaveOptPdb(char szName[]);
void Output_Gaussian_File(void);
void Get_Netcharge_From_Xpsf(void);
void Determine_QM_Scan_Part(void);
void Extract_Coord_E(char szName[], int Idx_Phi, int Idx_Part);
void Reorganize_QM_1D_Scan_Data(void);
double Cal_Phi_Dist(double Phi_1, double Phi_2);
void Get_EXE_Path(char szExeFile[], char szExePath[]);
void Setup_QM_Level(void);
void Replace_NewLine(char szBuff[]);
int Split_Into_Two_String(char szBuff[], char str1[], char str2[]);

void Read_MM_Rotamers(void);
int Is_This_A_Double_Bond(int ib, int ic);

int netcharge_mol;
int Active_Phi=0, n_Bins;
double E_Min_Org, E_Barrier_Lowest=1.0E100, E_Min, E_Min_Save, Phi_Save;



//start	data and functions related with job scheduler
#define MAX_JOB	(1024)

enum { IDLE, BUSY };
enum { Job_ToRun, Job_Running, Job_Done };

int nJob=0, nJobDone=0, nWorker=0, nCore_Per_Node=1, n_QM_Part_List[N_MAX_DIH];
int JobStatus[MAX_JOB];
char *szInputList[MAX_JOB], *szOutputList[MAX_JOB];

void RunAllJobs(void);
int Get_First_Available_Worker(void);
void Exatract_All_Info_QM_1D_Scan(void);
//end	data and functions related with job scheduler



#define E_RANGE	(5.0)
void Enumerate_Dihedrals(int IdxDih)	// depth
{
	int iState;
//	char szName[256], szIniName[256];
	double E_Barrier, E_Min, Phi_Save;

	if(IdxDih < n_Phi)	{
		for(iState=0; iState<n_State_List[IdxDih]; iState++)	{
			Phi_To_Set[IdxDih] = Phi_Set[IdxDih][iState];
			Enumerate_Dihedrals(IdxDih+1);
		}
	}
	else	{
		RestoreCoordinates();
		MM_Fixed_1D_Scan(0);
		E_Barrier = Get_Barrier(E_Phi, n_Bins, E_Min, Phi_Save);
		if( (E_Barrier + 1.0 < E_Barrier_Lowest) && (E_Min <= E_Min_Org) )	{
			memcpy(E_Scan, E_Phi, sizeof(double)*n_Bins);
			memcpy(Phi_To_Set_Scan, Phi_To_Set, sizeof(double)*n_Phi);
			Phi_To_Set_Scan[Active_Phi] = Phi_Save;
			E_Barrier_Lowest = E_Barrier;
			E_Min_Save = E_Min;
		}


//		sprintf(szIniName, "ini-%d.pdb", Counter);
//		Mol.SavePdb(szIniName);
//		Mol.SavePdb(szName);

		return;
	}
}

void Read_MM_Rotamers(void)
{
	int i, IdxID, MaxID=24, Count;
	FILE *fIn;
	char szName[PATH_MAX];
	double E_Barrier, Phi_Save, E_Min;

	for(IdxID=0; IdxID<MaxID; IdxID++)	{
		Count = 1;

		while(1)	{
			sprintf(szName, "../11-mm-pes/opt-pdb/opt-id=%d-%d.pdb", IdxID, Count);
			fIn = fopen(szName, "r");
			if(fIn == NULL)	{
				break;
			}
			else	{
				fclose(fIn);
			}

			Mol.ReadXYZ(szName);

			for(i=0; i<n_Phi; i++)	{
				Phi_To_Set[i] = Mol.QueryDihedral(IdxDihSelect[i]);
			}

			MM_Fixed_1D_Scan(1);
			E_Barrier = Get_Barrier(E_Phi, n_Bins, E_Min, Phi_Save);
			if( (E_Barrier + 1.0 < E_Barrier_Lowest) && (E_Min <= E_Min_Org) )	{
				memcpy(E_Scan, E_Phi, sizeof(double)*n_Bins);
				memcpy(Phi_To_Set_Scan, Phi_To_Set, sizeof(double)*n_Phi);
				Phi_To_Set_Scan[Active_Phi] = Phi_Save;
				E_Barrier_Lowest = E_Barrier;
				E_Min_Save = E_Min;
				strcpy(szOptPdbSave, szName);
			}

			Count++;
		}
	}
}

int main(int argc, char **argv)
{
  timebomb();

	int i;
	FILE *fOut;
	char szName[PATH_MAX], *szEnv;


	fFile_Run_Log = fopen("qm-1d-scan.log", "w");

//	szEnv = getenv(szGAUSS_SCRDIR);
//	if (szEnv == NULL) {
//		Quit_With_Error_Msg("Environment variable $GAUSS_SCRDIR is NOT set.\nQuit\n");
//	}
//	strcpy(szGAUSS_SCRDIR_Base, szEnv);
	

//	MPI_Init(&argc, &argv);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProgID);
//	MPI_Comm_size(MPI_COMM_WORLD, &nProc);

	Get_EXE_Path("BIN_G09", szExe_G09);

	n_Bins = 360/BIN_SIZE+1;


	ForceField.ReadForceField(szForceFiled);
	Mol.ReadPSF(szXpsfFile, 0);
	Get_Netcharge_From_Xpsf();
	Setup_QM_Level();
	ReadElementList();

	Mol.AssignForceFieldParameters(&ForceField);
	Mol.ReadXYZ(szCrdFile);
	BackupCoordinates();

	Read_Soft_DihedralList();

	nJob = nJobDone = 0;

	for(Active_Phi=0; Active_Phi<n_Phi; Active_Phi++)	{
		szOptPdbSave[0] = 0;

		RestoreCoordinates();
		for(i=0; i<n_Phi; i++)	{
			Phi_To_Set[i] = Mol.QueryDihedral(IdxDihSelect[i]);
			Mol.Edit_Dihedral(IdxDihSelect[i], Phi_To_Set[i]);
		}

		MM_Fixed_1D_Scan(0);
		E_Barrier_Lowest = Get_Barrier(E_Phi, n_Bins, E_Min, Phi_Save);
		E_Min_Save = E_Min;
		E_Min_Org = E_Min + E_RANGE;
		memcpy(E_Scan, E_Phi, sizeof(double)*n_Bins);
		memcpy(Phi_To_Set_Scan, Phi_To_Set, sizeof(double)*n_Phi);
		Phi_To_Set_Scan[Active_Phi] = Phi_Save;

		sprintf(szName, "ini-%d.dat", Active_Phi+1);
		fOut = fopen(szName, "w");
		for(i=0; i<n_Bins; i++)	{
			fprintf(fOut, "%5d %8.5E\n", -180+i*BIN_SIZE, E_Scan[i]-E_Min);
		}
		fclose(fOut);

		sprintf(szName, "ini-%d.xyz", Active_Phi+1);
		SaveOptPdb(szName);

//		n_State_List[Active_Phi] = 1;
//		Enumerate_Dihedrals(0);
		Read_MM_Rotamers();
//		n_State_List[Active_Phi] = n_State_List_Save[Active_Phi];

		if(szOptPdbSave[0])	{
			Mol.ReadXYZ(szOptPdbSave);
		}

		sprintf(szName, "end-%d.dat", Active_Phi+1);
		fOut = fopen(szName, "w");
		for(i=0; i<n_Bins; i++)	{
			fprintf(fOut, "%5d %8.5E\n", -180+i*BIN_SIZE, E_Scan[i]-E_Min_Save);
		}
		fclose(fOut);

		sprintf(szName, "end-%d.xyz", Active_Phi+1);
		if(szOptPdbSave[0])	{
			Mol.WriteXYZ(szName);
		}
		else	{
			SaveOptPdb(szName);
		}

	
		Output_Gaussian_File();
	}


	if( argc>=2 && !strcmp(argv[1], "--prepare" ) ) {
		exit(0);
	}
//  RunAllJobsSynchronously( nJob, szExe_G09, (char**)szInputList, (char**)szOutputList, NULL );

	Exatract_All_Info_QM_1D_Scan();
	Reorganize_QM_1D_Scan_Data();

	fclose(fFile_Run_Log);

//	MPI_Barrier(MPI_COMM_WORLD);
//	MPI_Finalize();

	return 0;
}

void Output_Gaussian_File(void)
{
	FILE *fOut;
	int i, nAtom, IdxScan;
	char szNameGjf[PATH_MAX], szName_Output[PATH_MAX];
  char dirname[PATH_MAX];

	Determine_QM_Scan_Part();

	nAtom = Mol.nAtom;
	for(IdxScan=0; IdxScan<n_QM_Part; IdxScan++)	{
		Mol.QueryDihedral(IdxDihSelect[Active_Phi]);
		Mol.Edit_Dihedral(IdxDihSelect[Active_Phi], Phi_QM_Start[IdxScan]);

        
		sprintf(dirname, "qm-1d-phi-%02d-%02d", Active_Phi+1, IdxScan+1);
    make_directory( dirname, 0700 );
		sprintf(szNameGjf, "%s/qm.gjf", dirname );
		fOut = fopen(szNameGjf, "w");
		fprintf(fOut, "%s", szQM_Level_1D_Scan_Cmd);
		for(i=0; i<nAtom; i++)	{
//			fprintf(fOut, "%c  %9.5lf %9.5lf %9.5lf\n", Mol.AtomName[i][0], Mol.x[i], Mol.y[i], Mol.z[i]);
			fprintf(fOut, "%s  %9.5lf %9.5lf %9.5lf\n", szElemName[i], Mol.x[i], Mol.y[i], Mol.z[i]);
		}
		fprintf(fOut, "\n");
		
		for(i=0; i<n_Phi; i++)	{
			if(i == Active_Phi)	{
				fprintf(fOut, "%s", szQM_Scan[IdxScan]);
			}
			else	{
				fprintf(fOut, "%d %d %d %d F\n", DihList[i][0]+1, DihList[i][1]+1, DihList[i][2]+1, DihList[i][3]+1);
			}
		}
		fprintf(fOut, "\n");
		fclose(fOut);

		sprintf(szName_Output, "%s/qm.out", dirname );
		szInputList[nJob] = strdup(szNameGjf);
		szOutputList[nJob]= strdup(szName_Output);
		nJob++;
	}

	n_QM_Part_List[Active_Phi] = n_QM_Part;
}


void Exatract_All_Info_QM_1D_Scan(void)
{
	int IdxScan;
	char szName_Output[PATH_MAX];
	
	for(Active_Phi=0; Active_Phi<n_Phi; Active_Phi++)	{
		for(IdxScan=0; IdxScan<n_QM_Part_List[Active_Phi]; IdxScan++)	{
			sprintf(szName_Output, "qm-1d-phi-%02d-%02d/qm.out", Active_Phi+1, IdxScan+1);
			Extract_Coord_E(szName_Output, Active_Phi, IdxScan);
		}
	}
}

void Extract_Coord_E(char szName[], int Idx_Phi, int Idx_Part)
{
	FILE *fIn, *fOut;
	int i, nAtom, ToRead=1, ReadItem, n_Rec=0;
	char szOutput[PATH_MAX], szErrorMsg[PATH_MAX], szLine[PATH_MAX], *ReadLine, szTag[PATH_MAX], ErrorMsg[PATH_MAX];
	double E_QM, E_MM, Phi, Phi_Set, x_Save[MAX_ATOM], y_Save[MAX_ATOM], z_Save[MAX_ATOM];

	nAtom = Mol.nAtom;
	memcpy(x_Save, Mol.x, sizeof(double)*nAtom);
	memcpy(y_Save, Mol.y, sizeof(double)*nAtom);
	memcpy(z_Save, Mol.z, sizeof(double)*nAtom);


	sprintf(szOutput, "tor-1D-idx-%d.dat", Idx_Phi+1);
	fOut = fopen(szOutput, "a+");
	fseek(fOut, 0, SEEK_END);

	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		sprintf(szErrorMsg, "Fail to open %s\nQuit\n", szName);
		Quit_With_Error_Msg(szErrorMsg);
	}


	int coordsset=0;
	int energyset=0;
	int optcomplete=0;

	while(ToRead)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, PATH_MAX, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		else	{
//			if(FindString(szLine, " Center     Atomic      Atomic")>=0)	{	// to extract the coordinate
//				Skip_N_Line(fIn, szName, 2);
//			if(FindString(szLine, "     Input orientation:")>=0)	{	// to extract the coordinate
			if(FindString(szLine, " orientation:")>=0)	{	// to extract the coordinate
				Skip_N_Line(fIn, szName, 4);
				for(i=0; i<nAtom; i++)	{
					ReadLine = fgets(szLine, PATH_MAX, fIn);
					if(ReadLine == NULL)	{
						break;
					}
					ReadItem = sscanf(szLine+31, "%lf %lf %lf", &(Mol.x[i]), &(Mol.y[i]), &(Mol.z[i]));
					if(ReadItem != 3)	{
						ToRead = 0;
						break;
					}
				}
				coordsset=1;
			}
			else if( (FindString(szLine, " SCF Done: ")>=0) && (QM_Level == QM_LEVEL_HF) )	{	// HF
				E_QM = Get_Energy(szLine);
				energyset=1;
			}
			else if( (FindString(szLine, "EUMP2 =")>=0) && (QM_Level == QM_LEVEL_MP2) )	{
				E_QM = Get_Energy(szLine+27);
				energyset=1;
			}
			else if(FindString(szLine, " Optimization completed")>=0)	{
				optcomplete=1;
			}

			if( energyset && coordsset && optcomplete ) {
#if 0
				sprintf(szTag, "  D(%d,%d,%d,%d)", DihList[Idx_Phi][0]+1, DihList[Idx_Phi][1]+1, DihList[Idx_Phi][2]+1, DihList[Idx_Phi][3]+1);
				To_Find_Tag(fIn, szName, szTag, szLine);
				ReadItem = sscanf(szLine+28, "%lf", &Phi_Set);	// previous, modredundant
				if(ReadItem != 1)	{
					sprintf(ErrorMsg, "Error in extracting the dihedral.\n%s\nQuit\n", szLine);
					Quit_With_Error_Msg(ErrorMsg);
				}
#endif
				Phi = Mol.QueryDihedral(IdxDihSelect[Idx_Phi]);
				E_MM = Mol.Cal_E(0);

				fprintf(fOut, "E_Scan %.13E %.13E  Phi  %.1lf\n", E_QM, E_MM, Phi);
				fprintf(fOut, "Coordinate\n");
				for(i=0; i<nAtom; i++)	{
					fprintf(fOut, "%12.6lf %12.6lf %12.6lf\n", Mol.x[i], Mol.y[i], Mol.z[i]);
				}

				coordsset=0;
				energyset=0;
				optcomplete=0;
				n_Rec++;
			}
		}
	}




	fclose(fOut);
	fclose(fIn);

	memcpy(Mol.x, x_Save, sizeof(double)*nAtom);
	memcpy(Mol.y, y_Save, sizeof(double)*nAtom);
	memcpy(Mol.z, z_Save, sizeof(double)*nAtom);
}

double Cal_Phi_Dist(double Phi_1, double Phi_2)
{
	double dist_abs;

	dist_abs = fabs(Phi_1-Phi_2);
	if(dist_abs > 180.0)	{
		dist_abs = fabs(360.0 - dist_abs);
	}
	return dist_abs;
}

#define E_BARRIER_DOUBLE_BOND	(12.0)
void Reorganize_QM_1D_Scan_Data(void)
{
	FILE *fIn, *fOut, *fQMStates;
	int i, j, i_Left, i_Right, Idx_Phi, nScan, nScan_Mapped, ReadItem, nAtom, To_Output[MAX_SCAN];
	char szName_Input[PATH_MAX], szName_Output[PATH_MAX], ErrorMsg[PATH_MAX], szLine[PATH_MAX], *ReadLine, szTmp[PATH_MAX];
	double x_List[MAX_SCAN][MAX_ATOM], y_List[MAX_SCAN][MAX_ATOM], z_List[MAX_SCAN][MAX_ATOM];
	double E_QM[MAX_SCAN], E_MM[MAX_SCAN], Phi[MAX_SCAN], Phi_Mapped[MAX_SCAN], E_QM_Mapped[MAX_SCAN], E_Min, E_Max, Phi_Ini, fTmp=0.0;
	int N_Rotamer_Total=1, nState, DoubleBond;

	RestoreCoordinates();

	nAtom = Mol.nAtom;
	fQMStates = fopen("qm-1d-states.dat", "w");
	for(Idx_Phi=0; Idx_Phi<n_Phi; Idx_Phi++)	{
		sprintf(szName_Input, "tor-1D-idx-%d.dat", Idx_Phi+1);
		sprintf(szName_Output, "new-tor-1D-idx-%d.dat", Idx_Phi+1);

		E_Min = 1.0E100;
		E_Max = -1.0E100;
		Phi_Ini = Mol.QueryDihedral(IdxDihSelect[Idx_Phi]);
		DoubleBond = Is_This_A_Double_Bond(DihList[Idx_Phi][1], DihList[Idx_Phi][2]);

		nScan = 0;
		fIn = fopen(szName_Input, "r");
		if(fIn == NULL)	{
			sprintf(ErrorMsg, "Fail to open file for read: %s\nQuit\n", szName_Input);
			Quit_With_Error_Msg(ErrorMsg);
		}
		while(1)	{
			if(feof(fIn))	{
				break;
			}
			ReadLine = fgets(szLine, PATH_MAX, fIn);
			if(ReadLine == NULL)	{
				break;
			}
			if(strncmp(szLine, "E_Scan", 6) == 0)	{
				ReadItem = sscanf(szLine, "%s %lf %lf %s %lf", szTmp, &(E_QM[nScan]), &(E_MM[nScan]), szTmp, &(Phi[nScan]));
				if(ReadItem == 5)	{
					fgets(szLine, PATH_MAX, fIn);
					for(i=0; i<nAtom; i++)	{
						fgets(szLine, PATH_MAX, fIn);
						ReadItem = sscanf(szLine, "%lf %lf %lf", &(x_List[nScan][i]), &(y_List[nScan][i]), &(z_List[nScan][i]));
						if(ReadItem != 3)	{
							fclose(fIn);
							sprintf(ErrorMsg, "Error in readubg file: %s\nQuit\n", szName_Input);
							Quit_With_Error_Msg(ErrorMsg);
						}
					}
					if(E_QM[nScan] < E_Min)	{
						E_Min = E_QM[nScan];
					}
					if(E_QM[nScan] > E_Max)	{
						E_Max = E_QM[nScan];
					}

					nScan++;
				}
				else	{
					break;
				}
			}
		}
		fclose(fIn);


		for(i=0; i<nScan; i++)	{
			To_Output[i] = 1;
		}

		if( (Phi[0] > 0.0) && (Cal_Phi_Dist(Phi[0], -180.0) < 1.0) )	{	// +180.0 at the first point
			Phi[0] -= 360.0;
		}

		//start	to delete those redundant entries
		for(i=0; i<nScan; i++)	{
			for(j=i+1; j<nScan; j++)	{
				if(Cal_Phi_Dist(Phi[i], Phi[j]) < 1.0)	{	// the same dihedral
					if(Cal_Phi_Dist(Phi[j], -180.0) > 1.0)	{	// only -180, 180 can exist, others will be removed
						To_Output[j] = 0;
					}
				}
			}
		}
		//end	to delete those redundant entries

		//start	to set the correct sign of +180 if it exists as -180
		if(Cal_Phi_Dist(Phi[0], -180.0) < 1.0)	{	// starting with -180
			for(i=nScan-1; i>0; i--)	{
				if(To_Output[i] == 0)	{
					continue;
				}
				
				if( (Phi[i] < 0.0) && (Cal_Phi_Dist(Phi[i], -180.0) < 1.0) )	{
					Phi[i] += 360.0;
//					break;
				}
			}
		}
		//end	to set the correct sign of +180 if it exists as -180

//		fOut = fopen(szName_Output, "w");
		fOut = fopen(szName_Input, "w");	// to overwrite 
		for(j=0; j<nScan; j++)	{
			if(To_Output[j]==0)	{
				continue;
			}
			fprintf(fOut, "E_Scan %.13E %.13E  Phi  %.1lf\n", E_QM[j], E_MM[j], Phi[j]);
			fprintf(fOut, "Coordinate\n");
			for(i=0; i<nAtom; i++)	{
				fprintf(fOut, "%12.6lf %12.6lf %12.6lf\n", x_List[j][i], y_List[j][i], z_List[j][i]);
			}
		}
		
		fclose(fOut);

		//start	to identify QM local minima of each dihedral
		nScan_Mapped = 0;
		for(i=0; i<nScan; i++)	{
			if(To_Output[i])	{
				if(fabs(Phi[i]-180.0) > 3.0)	{	// +180 is excluded
					Phi_Mapped[nScan_Mapped] = Phi[i];
					E_QM_Mapped[nScan_Mapped] = E_QM[i];
					nScan_Mapped++;
				}
			}
		}
		if(nScan_Mapped < 3)	{
			sprintf(ErrorMsg, "nScan_Mapped < 3 for dihedral %d \nQuit\n", Idx_Phi+1);
			Quit_With_Error_Msg(ErrorMsg);
		}
		fprintf(fQMStates, "%3d %3d %3d %3d ", DihList[Idx_Phi][0]+1, DihList[Idx_Phi][1]+1, DihList[Idx_Phi][2]+1, DihList[Idx_Phi][3]+1);
		nState = 0;
		for(i=0; i<nScan_Mapped; i++)	{
			i_Left = (i+nScan_Mapped-1)%nScan_Mapped;
			i_Right = (i+nScan_Mapped+1)%nScan_Mapped;

			if( (DoubleBond) && ( ((E_Max-E_Min)*HARTREE_KCAL) > E_BARRIER_DOUBLE_BOND ) )	{
				fprintf(fQMStates, "  %.0lf", Phi_Ini);
				nState = 1;
				break;
			}

			if( (E_QM_Mapped[i] <= E_QM_Mapped[i_Left]) && (E_QM_Mapped[i] <= E_QM_Mapped[i_Right]) )	{	// a local minimum along QM 1D profile
				fprintf(fQMStates, "  %.0lf", Phi_Mapped[i]);
				nState++;
			}
		}
		fprintf(fQMStates, "\n");
		N_Rotamer_Total *= nState;
		//end	to identify QM local minima of each dihedral
	
	}

	fclose(fQMStates);

	fOut = fopen("n_rotamer.txt", "w");
	fprintf(fOut, "%d", N_Rotamer_Total);
	fclose(fOut);
}
#undef E_BARRIER_DOUBLE_BOND

void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in QM_1D_scan.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}

int Read_Soft_DihedralList(void)
{
	FILE *fIn;
	int ReadItem;
	char szLine[PATH_MAX], *ReadLine;

	n_Phi = 0;

	fIn = fopen(szPhiToScan, "r");

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 128, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		ReadItem = sscanf(szLine, "%d %d %d %d %lf %lf %lf %lf %lf %lf", 
			&(DihList[n_Phi][0]), &(DihList[n_Phi][1]), &(DihList[n_Phi][2]), &(DihList[n_Phi][3]), 
			&(Phi_Set[n_Phi][0]), &(Phi_Set[n_Phi][1]), &(Phi_Set[n_Phi][2]), &(Phi_Set[n_Phi][3]), &(Phi_Set[n_Phi][4]), &(Phi_Set[n_Phi][5]));

		if(ReadItem > 4)	{
			DihList[n_Phi][0]--;
			DihList[n_Phi][1]--;
			DihList[n_Phi][2]--;
			DihList[n_Phi][3]--;
			IdxDihSelect[n_Phi] = Mol.Query_Dihedral_Index(DihList[n_Phi][0], DihList[n_Phi][1], DihList[n_Phi][2], DihList[n_Phi][3]);

			if(IdxDihSelect[n_Phi] < 0)	{
				Quit_With_Error_Msg("Fail to identify the index of one soft dihedral.\n");
			}
			Mol.BuildSegmentList_Dihedrals(IdxDihSelect[n_Phi]);

			n_State_List[n_Phi] = ReadItem - 4;
			n_Phi++;
		}
		else	{
			break;
		}
	}

	fclose(fIn);

	memcpy(Phi_Set_Save, Phi_Set, sizeof(double)*N_MAX_DIH*MAX_N_STATE);
	memcpy(n_State_List_Save, n_State_List, sizeof(int)*N_MAX_DIH);

	return n_Phi;
}

void BackupCoordinates(void)
{
	int nAtom;
	nAtom = Mol.nAtom;

	memcpy(x_Save, Mol.x, sizeof(double)*nAtom);
	memcpy(y_Save, Mol.y, sizeof(double)*nAtom);
	memcpy(z_Save, Mol.z, sizeof(double)*nAtom);
}

void RestoreCoordinates(void)
{
	int nAtom;
	nAtom = Mol.nAtom;

	memcpy(Mol.x, x_Save, sizeof(double)*nAtom);
	memcpy(Mol.y, y_Save, sizeof(double)*nAtom);
	memcpy(Mol.z, z_Save, sizeof(double)*nAtom);
}

double Get_Barrier(double E[], int n, double& E_Min, double& Phi)
{
	int i, iMin;
	double E_Max=-1.0E100;

	E_Min=1.0E100;

	for(i=0; i<n; i++)	{
		if(E[i] < E_Min)	{
			E_Min = E[i];
			iMin = i;
		}
		if(E[i] > E_Max)	{
			E_Max = E[i];
		}
	}
	Phi = 1.0*(-180 + iMin*BIN_SIZE);

	return (E_Max-E_Min);
}

void MM_Fixed_1D_Scan(int ReadInPdb)
{
	int i, Phi, Count;

//	RestoreCoordinates();
	if(ReadInPdb == 0)	{
		for(i=0; i<n_Phi; i++)	{
			Mol.QueryDihedral(IdxDihSelect[i]);
			Mol.Edit_Dihedral(IdxDihSelect[i], Phi_To_Set[i]);
		}
	}
	
	Count = 0;
	for(Phi=-180; Phi<=180; Phi+=BIN_SIZE)	{	// 1D scan for this rotamer
		Phi_To_Set[Active_Phi] = Phi*1.0;
		Mol.QueryDihedral(IdxDihSelect[Active_Phi]);
		Mol.Edit_Dihedral(IdxDihSelect[Active_Phi], Phi_To_Set[Active_Phi]);
		E_Phi[Count] = Mol.Cal_E(0);
		Count++;
	}

	return;
}

void SaveOptPdb(char szName[])
{
	int i;
	for(i=0; i<n_Phi; i++)	{
		Mol.QueryDihedral(IdxDihSelect[i]);
		Mol.Edit_Dihedral(IdxDihSelect[i], Phi_To_Set_Scan[i]);
	}
	Mol.WriteXYZ(szName);
}


void Get_Netcharge_From_Xpsf(void)
{
	int i, nAtom;
	double CG_Sum=0.0;
	char szCharge[16];

	nAtom = Mol.nAtom;
	for(i=0; i<nAtom; i++)	{
		CG_Sum += Mol.CG[i];
	}
	sprintf(szCharge, "%.0lf", CG_Sum);
	netcharge_mol = atoi(szCharge);
}


int Is_This_A_Double_Bond(int ib, int ic)
{
	if( (fabs(Mol.mass[ib]-12.0) < 0.2) && (Mol.BondCount[ib] == 3) )	{	// sp2 C
		if( (fabs(Mol.mass[ic]-12.0) < 0.2) && (Mol.BondCount[ic] == 3) )	{	// sp2 C
			return 1;
		}
		if( (fabs(Mol.mass[ic]-14.0) < 0.2) && ( (Mol.BondCount[ic] == 3) || ((Mol.BondCount[ic] == 2)) ) )	{	// sp2 N
			return 1;
		}
	}

	if( (fabs(Mol.mass[ic]-12.0) < 0.2) && (Mol.BondCount[ic] == 3) )	{	// sp2 C
		if( (fabs(Mol.mass[ib]-12.0) < 0.2) && (Mol.BondCount[ib] == 3) )	{	// sp2 C
			return 1;
		}
		if( (fabs(Mol.mass[ib]-14.0) < 0.2) && ( (Mol.BondCount[ib] == 3) || ((Mol.BondCount[ib] == 2)) ) )	{	// sp2 N
			return 1;
		}
	}

	if( (fabs(Mol.mass[ib]-14.0) < 0.2) && (fabs(Mol.mass[ic]-14.0) < 0.2) && (Mol.BondCount[ib] == 2) && (Mol.BondCount[ic] == 2) )	{	// two sp2 N
		return 1;
	}
	

	return 0;
}



void To_Split_1D_Scan(int Phi_Begin, int nScan)
{
	int i, Gap, nJob, nSubScan;
	int Scan_Size;

	if(nScan < 10)	{
		return;
	}

	n_QM_Part--;

	Scan_Size = 8;

	while( nScan > 0)	{
		if(nScan == (Scan_Size+1) )	{
			nSubScan = nScan;
		}
		else if( nScan > (Scan_Size+1) )		{
			nSubScan = Scan_Size;
		}
		else	{
			nSubScan = nScan;
		}

		Phi_QM_Start[n_QM_Part] = 1.0*Phi_Begin;
		
//		sprintf(szQM_Scan[n_QM_Part], "%d %d %d %d %.1lf S %d %.1lf\n", 
//			DihList[Active_Phi][0]+1, DihList[Active_Phi][1]+1, DihList[Active_Phi][2]+1, DihList[Active_Phi][3]+1, 
//			1.0*Phi_Begin, nSubScan, 1.0*BIN_SIZE_QM);
      sprintf(szQM_Scan[n_QM_Part], "%d %d %d %d %.1lf B\n%d %d %d %d S %d %.1lf\n",
            DihList[Active_Phi][0]+1, DihList[Active_Phi][1]+1, DihList[Active_Phi][2]+1, DihList[Active_Phi][3]+1, 1.0*Phi_Begin, 
            DihList[Active_Phi][0]+1, DihList[Active_Phi][1]+1, DihList[Active_Phi][2]+1, DihList[Active_Phi][3]+1, 
            nSubScan, 1.0*BIN_SIZE_QM);

		n_QM_Part++;
		
		nScan -= (nSubScan+1);
		Phi_Begin += (BIN_SIZE_QM*(nSubScan+1));
	}

}

#define MAX_E_To_SCAN	(300.0)
//#define MAX_E_To_SCAN	(120.0)
void Determine_QM_Scan_Part(void)
{
	int i, *Flag_To_Scan;
	int Phi_Begin, Phi_End, nScan, Is_Double_Bond, Phi;

	Is_Double_Bond = Is_This_A_Double_Bond(DihList[Active_Phi][1], DihList[Active_Phi][2]);
	E_Scan[n_Bins] = 1.0E100;	// A stop sign

	Flag_To_Scan = new int[n_Bins+1];
	for(i=0; i<=n_Bins; i++)	{
		if((E_Scan[i]-E_Min_Save) <= MAX_E_To_SCAN)	{
			Flag_To_Scan[i] = 1;
		}
		else	{
			Flag_To_Scan[i] = 0;
		}

		if(Is_Double_Bond)	{
			Phi = -180 + i*BIN_SIZE;
			if( fabs(Phi + 90.0) < 8.0)	{	// not stable
				Flag_To_Scan[i] = 0;
			}
			if( fabs(Phi - 90.0) < 8.0)	{	// not stable
				Flag_To_Scan[i] = 0;
			}
		}
	}


	n_QM_Part = 0;
	
	for(i=0; i<n_Bins; i++)	{
		if(Flag_To_Scan[i])	{
			break;
		}
	}
	Phi_Begin = -180 + i*BIN_SIZE;
	Phi_QM_Start[n_QM_Part] = 1.0*Phi_Begin;

#if 0
	for(; i<=n_Bins; i++)	{
		if( Flag_To_Scan[i] == 0 )	{
			Phi_End = -180 + (i-1)*BIN_SIZE;
			nScan = (Phi_End-Phi_Begin)/BIN_SIZE_QM;

//			sprintf(szQM_Scan[n_QM_Part], "%d %d %d %d %.1lf S %d %.1lf\n", 
//				DihList[Active_Phi][0]+1, DihList[Active_Phi][1]+1, DihList[Active_Phi][2]+1, DihList[Active_Phi][3]+1, 
//				1.0*Phi_Begin, nScan, 1.0*BIN_SIZE_QM);
      sprintf(szQM_Scan[n_QM_Part], "%d %d %d %d %.1lf B\n%d %d %d %d S %d %.1lf\n",
            DihList[Active_Phi][0]+1, DihList[Active_Phi][1]+1, DihList[Active_Phi][2]+1, DihList[Active_Phi][3]+1, 1.0*Phi_Begin, 
            DihList[Active_Phi][0]+1, DihList[Active_Phi][1]+1, DihList[Active_Phi][2]+1, DihList[Active_Phi][3]+1, 
            nScan, 1.0*BIN_SIZE_QM);

			n_QM_Part++;

			To_Split_1D_Scan(Phi_Begin, nScan);

			for(; i<n_Bins; i++)	{	// to find the new beginning
				if(Flag_To_Scan[i])	{
					break;
				}
			}
			Phi_Begin = -180 + i*BIN_SIZE;
			Phi_QM_Start[n_QM_Part] = 1.0*Phi_Begin;
		}
	}
#endif

	while( i<=n_Bins ) {
		if( Flag_To_Scan[i] ) {
			Phi_QM_Start[n_QM_Part] = -180 + i * BIN_SIZE;
      sprintf(szQM_Scan[n_QM_Part], "%d %d %d %d %.1lf B\n%d %d %d %d F\n",
            DihList[Active_Phi][0]+1, DihList[Active_Phi][1]+1, DihList[Active_Phi][2]+1, DihList[Active_Phi][3]+1, Phi_QM_Start[n_QM_Part],
            DihList[Active_Phi][0]+1, DihList[Active_Phi][1]+1, DihList[Active_Phi][2]+1, DihList[Active_Phi][3]+1  );
			n_QM_Part++;
		}
		i++;
		if( n_QM_Part >= N_MAX_QM_PART ) {
			Quit_With_Error_Msg( "Exceeded compile-time constant N_MAX_QM_PART" );
		}

	}


	delete []Flag_To_Scan;
}
#undef MAX_E_To_SCAN



void Replace_D_E(char szStr[])
{
	int nLen, i;

	nLen = strlen(szStr);
	for(i=0; i<nLen; i++)	{
		if(szStr[i] == 'D')	{
			szStr[i] = 'E';
			return;
		}
	}
	return;
}

double Get_Energy(char szLine[])
{
	char szEnergy[PATH_MAX], ErrorMsg[PATH_MAX];
	int i=0, Count=0;

	while(1)	{
		if(szLine[i]=='=')	{
			i++;
			break;
		}
		else if(szLine[i]==0x0)	{
			sprintf(ErrorMsg, "Fail to extract the energy from, \n%s\nQuit\n", szLine);
			Quit_With_Error_Msg(ErrorMsg);
		}
		else	{
			i++;
		}
	}

	while(1)	{
		if(szLine[i] != ' ')	{
			break;
		}
		else if(szLine[i]==0x0)	{
			sprintf(ErrorMsg, "Fail to extract the energy from, \n%s\nQuit\n", szLine);
			Quit_With_Error_Msg(ErrorMsg);
		}
		else	{
			i++;
		}
	}

	while(1)	{
		if(szLine[i] == ' ')	{
			break;
		}
		else if(szLine[i]==0x0)	{
			if(Count == 0)	{
				sprintf(ErrorMsg, "Error in extracting the energy from, \n%s\nQuit\n", szLine);
				Quit_With_Error_Msg(ErrorMsg);
			}
			else	{
				break;
			}
		}
		else	{
			szEnergy[Count] = szLine[i];
			Count++;
			i++;
		}
	}
	szEnergy[Count] = 0;

	Replace_D_E(szEnergy);

	return atof(szEnergy);
}

int To_Find_Tag(FILE *fIn, char szFileName[], char szTag[], char szLine[])
{
	char *ReadLine, ErrorMsg[PATH_MAX];

	while(1)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Fail to find the tag: %s in file %s\nQuit\n", szTag, szFileName);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}

		ReadLine = fgets(szLine, PATH_MAX, fIn);
		if(ReadLine == NULL)	{
			sprintf(ErrorMsg, "Fail to find the tag: %s in file %s\nQuit\n", szTag, szFileName);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
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
	char szLine[PATH_MAX], *ReadLine, ErrorMsg[PATH_MAX];

	for(i=0; i<n; i++)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Fail in Skip_N_Line(%s, %d).\nQuit\n", szFileName, n);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}

		ReadLine = fgets(szLine, PATH_MAX, fIn);
		if(ReadLine == NULL)	{
			sprintf(ErrorMsg, "Fail in Skip_N_Line(%s, %d).\nQuit\n", szFileName, n);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}
	}
	return 1;
}


void Get_EXE_Path(char szExeFile[], char szExePath[])
{

  if(getenv( szExeFile ) ) {
    strcpy( szExePath, getenv( szExeFile ) );
  }
  else {
    fprintf(stderr, "Envvar [%s] not set\n", szExeFile );
    exit(1);
  }

}


void Setup_QM_Level(void)
{
	FILE *fIn;
	char *ReadLine, szLine[PATH_MAX], szSubStr_1[PATH_MAX], szSubStr_2[PATH_MAX], ErrorMsg[PATH_MAX];
	char szKey_QM_MEM[PATH_MAX], szKey_QM_NPROC[PATH_MAX], szQM_Level_Opt[PATH_MAX], szQM_Level_Dimer_Opt[PATH_MAX], szQM_Level_ESP[PATH_MAX], szQM_Level_E_Monomer[PATH_MAX], szQM_Level_1D_Scan[PATH_MAX];
	int nStr, KeyCount;

	fIn = fopen(szQMParaFile, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file %s\nQuit\n", szQMParaFile);
		Quit_With_Error_Msg(ErrorMsg);
	}

	KeyCount = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, PATH_MAX, fIn);
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
					nCore_Per_Node = atoi(szKey_QM_NPROC);
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_OPT")==0)	{
					strcpy(szQM_Level_Opt, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_DIMER_OPT")==0)	{
					strcpy(szQM_Level_Dimer_Opt, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_ESP")==0)	{
					strcpy(szQM_Level_ESP, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_E_MONOMER")==0)	{
					strcpy(szQM_Level_E_Monomer, szSubStr_2);
					KeyCount++;
				}
				else if(strcmp(szSubStr_1,"QM_LEVEL_1D_SCAN")==0)	{
					strcpy(szQM_Level_1D_Scan, szSubStr_2);
					KeyCount++;
				}
			}
		}
	}

	fclose(fIn);


	Determine_QM_Level(szQM_Level_1D_Scan);


	if(KeyCount != 7)	{	// !!!
		sprintf(ErrorMsg, "Setup_QM_Level> Error: incomplete entries in %s\nQuit\n", szQMParaFile);
		Quit_With_Error_Msg(ErrorMsg);
	}


  // MJH start
  if( getenv("NCORES" ) ) {
      printf(" Overriding QM NCPUS with [%s]\n", getenv("NCORES") );
      strcpy(szKey_QM_NPROC, getenv("NCORES") );
      nCore_Per_Node = atoi(szKey_QM_NPROC);
  }
  if( getenv( "MEMORY" ) ) {
      printf(" Overriding QM MEMORY with [%s]\n", getenv("MEMORY") );
      strcpy(szKey_QM_MEM, getenv("MEMORY" ) );
      strcat( szKey_QM_MEM, "GB" );
  }
  // MJH stop
  //


	sprintf(szQM_Level_1D_Scan_Cmd, "%%mem=%s\n%%nproc=%s\n%sMol dihedral scan\n\n%d,1\n",
		szKey_QM_MEM, szKey_QM_NPROC, szQM_Level_1D_Scan, netcharge_mol);
//	printf("%s", szQM_Level_1D_Scan_Cmd);
}

void Determine_QM_Level(char szQMLevel[])
{
	char ErrorMsg[PATH_MAX];

	if(FindString(szQMLevel, "HF") >= 0)	{
		QM_Level = QM_LEVEL_HF;
	}
	else if(FindString(szQMLevel, "MP2") >= 0)	{
		QM_Level = QM_LEVEL_MP2;
	}
	else	{
		sprintf(ErrorMsg, "Fail to determine the QM level.\n%s\nQuit\n", szQMLevel);
		Quit_With_Error_Msg(ErrorMsg);
	}
	return;
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
		if( (szBuff[i] == 0x0) || (szBuff[i] == 0x22) ) 	{	// To find the last character of the second string
//		if( (szBuff[i] == 0x0) || (szBuff[i] == 0x0D) || (szBuff[i] == 0x0A) || (szBuff[i] == 0x22) ) 	{	// To find the last character of the second string
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


void ReadElementList(void)
{
	int n_Atom_In_Mol, i;
	FILE *fIn;
	char ErrorMsg[PATH_MAX];

	n_Atom_In_Mol = Mol.nAtom;

	fIn = fopen(szFileElem, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file %s\nQuit\n", szFileElem);
		Quit_With_Error_Msg(ErrorMsg);
	}

	for(i=0; i<n_Atom_In_Mol; i++)	{
		fscanf(fIn, "%s", szElemName[i]);
	}
	fclose(fIn);

}




