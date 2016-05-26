/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

//#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include "nlopt.h"

#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#else
#include <unistd.h>	//for linux
#endif

#include "ff.h"

char szForceFiled[]="mol.prm";
char szXpsfFile[]="mol.xpsf";
char szCrdFile[]="mol-opt.xyz";
char szPhiToScan[]="qm-1d-states.dat";

char szExe_G09[256];
char szQMParaFile[]="../QM-para.txt";
char szMyPath[]="../mypath.txt";
char szFileElem[]="elem-list.txt";

char szQM_Level_Rotamer_Cmd[512]="%mem=1600MB\n%nproc=4\n# HF/6-31G* nosymm opt\n\nMol rotamer\n\n%d 1\n";	// kbt

const char szGAUSS_SCRDIR[]="GAUSS_SCRDIR";
char szGAUSS_SCRDIR_Base[256];

#ifndef MINAB
#define MINAB 
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif

#define N_MAX_DIH	(16)
#define MAX_N_STATE		(6)
#define MAX_BIN_PHI	(128)
#define BIN_SIZE	(5)
#define BIN_SIZE_QM	(10)
#define N_MAX_QM_PART	(6)
#define MAX_CONF	(32768)
#define MAX_ROTAMER	(16384)	// an approximate maximum
#define N_ROTAMER_PICK	(200)	// the total number of rotamers for QM optimizations
#define E_CUTOFF		(30.0)
//#define RMSD_CUTOFF	(0.5)
#define N_TOR_PARA	(11)

double RMSD_CUTOFF=0.5;

#define QM_LEVEL_HF		(1)
#define QM_LEVEL_MP2	(2)


CMol Mol;
CForceField ForceField;
int nDihedral, QM_Level;
double Ratio_Pick=1.0;

FILE *fFile_Run_Log;
void Quit_With_Error_Msg(char szMsg[]);
double Para_Tor_List[N_MAX_DIH][N_TOR_PARA];

int n_Phi=0, Counter=0;
int DihList[N_MAX_DIH][4], IdxDihSelect[N_MAX_DIH], n_State_List[N_MAX_DIH], n_State_List_Save[N_MAX_DIH];
double Phi_Set[N_MAX_DIH][MAX_N_STATE], Phi_Set_Save[N_MAX_DIH][MAX_N_STATE], Phi_To_Set[N_MAX_DIH], Phi_To_Set_Scan[N_MAX_DIH];
double x_Save[MAX_ATOM], y_Save[MAX_ATOM], z_Save[MAX_ATOM];
double E_Phi[MAX_BIN_PHI], E_Scan[MAX_BIN_PHI];
char szElemName[MAX_ATOM][8];

int nAcceptor, Acceptor[MAX_ATOM], nHDonor, Donor[MAX_ATOM], IsH[MAX_ATOM];
int Bond_Count[MAX_ATOM], BondList[MAX_ATOM][4], dist[MAX_ATOM][MAX_ATOM];

int n_QM_Part=0, Phi_QM_n_Conf[N_MAX_QM_PART];
char szQM_Scan[N_MAX_QM_PART][256];
double Phi_QM_Start[N_MAX_QM_PART];
double E_Conf[MAX_CONF], Phi_Conf[MAX_CONF][N_MAX_DIH];
double x_List[MAX_CONF][MAX_ATOM], y_List[MAX_CONF][MAX_ATOM], z_List[MAX_CONF][MAX_ATOM];
int Rotamer_QM[MAX_CONF], Idx_Sorted[MAX_CONF], ClusterList[MAX_CONF], nCluster;

int	ProgID=0, nProc=1;
int iseed = 7643;

extern int FindString(char szBuff[], char szTag[]);
void Replace_D_E(char szStr[]);
int To_Find_Tag(FILE *fIn, char szFileName[], char szTag[], char szLine[], int Quit);
int Skip_N_Line(FILE *fIn, char szFileName[], int n);
double Get_Energy(char szLine[]);
void ReadElementList(void);

int Read_Soft_DihedralList(void);
int Get_Dihedral_Index(int ia, int ib, int ic, int id);
void Enumerate_Dihedrals(int IdxDih, int iState);
void BackupCoordinates(void);
void RestoreCoordinates(void);
void SaveOptPdb(char szName[]);
void Output_Gaussian_File(void);
void Get_Netcharge_From_Xpsf(void);
double Extract_Rotamer_Coord_E(char szName[]);

void Get_EXE_Path(char szExeFile[], char szExePath[]);
void Setup_QM_Level(void);
int Split_Into_Two_String(char szBuff[], char str1[], char str2[]);
void Replace_NewLine(char szBuff[]);
double rand(int &iseed);

int Iter_Geo_Opt;
double func_Geo_Optimization(unsigned n, const double *x, double *grad, void *my_func_data);
double Do_Geometry_Optimization(CMol *pMol);
int Read_All_MM_Rotamers(int& n_Rotamer);

void Make_List_H_Acceptor(void);
void GetBondList(int Idx);
void Count_Intramolecular_HBond(double x[], double y[], double z[], int nAtom, int& n_HBond, int& n_Far_HBond);
void Exclude_Conf_Longdist_HBond(void);
void SortByEnergy(void);
void ClusteringRotamers(void);
void Read_Tor_Para_1D_Fitting(void);
void Assign_Torsion_Parameters(void);
void WriteIniConf_For_QM(void);

double rmsfit(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);
void quatfit(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);
void jacobi(int n, int np, double a[5][5], double d[5], double v[5][5], double b[5], double z[5]);
double CalRMSD(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);
void Gen_QM_Rotamer_Data(void);
void Determine_QM_Level(char szQMLevel[]);

int netcharge_mol;
int n_Conf=0;
double E_Min=1.0E100;


//start	data and functions related with job scheduler
#define MAX_JOB	(2048)

int nJob=0, nJobDone=0, nWorker=0, nCore_Per_Node=1, n_QM_Part_List[N_MAX_DIH];
int JobStatus[MAX_JOB];
char *szDirList[MAX_JOB], *szInputList[MAX_JOB], *szOutputList[MAX_JOB];

int Count_Finished_Job();
void RunAllJobs(void);
int Get_First_Available_Worker(void);
void Exatract_All_Info_QM_1D_Scan(void);
//end	data and functions related with job scheduler


void Enumerate_Dihedrals(int IdxDih)	// depth
{
	int iState, i;
	FILE *fIn;
	double r;
	char szName[256];

	if(IdxDih < n_Phi)	{
		for(iState=0; iState<n_State_List[IdxDih]; iState++)	{
			Phi_To_Set[IdxDih] = Phi_Set[IdxDih][iState];
			Enumerate_Dihedrals(IdxDih+1);
		}
	}
	else	{
		r = rand(iseed);
		if(r > Ratio_Pick)	{
			return;
		}

		Counter++;
		if(Counter%nProc != ProgID)	{
			return;
		}
		sprintf(szName, "opt-%d.xyz", Counter);
		fIn = fopen(szName, "r");
		if(fIn != NULL)	{	// already done, then skip
			fclose(fIn);
			return;
		}

		RestoreCoordinates();
		for(i=0; i<n_Phi; i++)	{
			Mol.QueryDihedral(IdxDihSelect[i]);
			Mol.Edit_Dihedral(IdxDihSelect[i], Phi_To_Set[i]);
		}
//		sprintf(szNameIni, "ini-%d.pdb", Counter);
//		Mol.WriteXYZ(szNameIni);
		Do_Geometry_Optimization(&Mol);
		Mol.WriteXYZ(szName);	// high precision

		return;
	}
}

int Read_All_MM_Rotamers(int& n_Rotamer)
{
	int i, nAtom;
	char szName[256];
	double E_min_Local=1.0E100;
	
	nAtom = Mol.nAtom;

	for(i=0; i<n_Rotamer; i++)	{
		sprintf(szName, "opt-%d.xyz", i+1);
		printf("Reading %s\n", szName );
		Mol.ReadXYZ(szName);
//		Mol.FullGeometryOptimization_LBFGS();	// should be very fast !
		E_Conf[i] = Mol.Cal_E(0);

		if( E_Conf[i] < E_min_Local)	{
			E_min_Local = E_Conf[i];
		}

		memcpy(x_List[i], Mol.x, sizeof(double)*nAtom);
		memcpy(y_List[i], Mol.y, sizeof(double)*nAtom);
		memcpy(z_List[i], Mol.z, sizeof(double)*nAtom);
		Rotamer_QM[i] = 1;
	}


	RestoreCoordinates();	// append the ini configuration
	memcpy(x_List[i], Mol.x, sizeof(double)*nAtom);
	memcpy(y_List[i], Mol.y, sizeof(double)*nAtom);
	memcpy(z_List[i], Mol.z, sizeof(double)*nAtom);
	E_Conf[i] = E_min_Local - 0.1;	// set as the configuration with lowest energy
	Rotamer_QM[i] = 1;
	n_Rotamer++;


	SortByEnergy();

	return i;
}

char szCurDir[512], szNewDir[512], szGjfFile[]="mol-qm-rotamer.gjf";
int main(int argc, char **argv)
{
  timebomb();

	int i, j, nAtom;
	FILE *fIn, *fOut;
	double E_Rotamer, E_Rotamer_Min=1.0E100;
	char szOutput[]="mol-qm-rotamer.out", *szEnv;

	if(argc < 2)	{
		printf("Usage: qm-rotamer-scan-large rmsd_cutoff [--prepare]\n");
		exit(1);
	}

	RMSD_CUTOFF = atof(argv[1]);


	fFile_Run_Log = fopen("qm-rotamer-scan.log", "w");

	szEnv = getenv(szGAUSS_SCRDIR);
	if (szEnv == NULL) {
		Quit_With_Error_Msg("Environment variable $GAUSS_SCRDIR is NOT set.\nQuit\n");
	}
	strcpy(szGAUSS_SCRDIR_Base, szEnv);

	Get_EXE_Path("BIN_G09", szExe_G09);

	ForceField.ReadForceField(szForceFiled);
	Mol.ReadPSF(szXpsfFile, 0);
	Make_List_H_Acceptor();


	for(i=0; i<Mol.nAtom; i++)	{
		Mol.Dijkstra(i, &(dist[i][0]));
	}
	

	Get_Netcharge_From_Xpsf();

	Setup_QM_Level();
	ReadElementList();

	Mol.AssignForceFieldParameters(&ForceField);
	Mol.ReadXYZ(szCrdFile);
	BackupCoordinates();
	nAtom = Mol.nAtom;

	Read_Soft_DihedralList();


//  The enumeration is deterministic, so can safely repro it in --complete
//  Which is needed to get nJob and the output filename (yuck, should probably just regen only the filenames)
//	fIn = fopen("qm-ini-1.xyz", "r");
//	if(fIn == NULL)	{
		Enumerate_Dihedrals(0);	// do MM optimization
		
		Read_All_MM_Rotamers(Counter);
		system("/bin/tar -zcvf opt-mm.tgz opt-*.xyz > /dev/null ; /bin/rm opt-*.xyz");	// tar and delete pdb files
		Exclude_Conf_Longdist_HBond();	// exclude those configurations with long-distance intra-molecular H-Bond
		ClusteringRotamers();
		WriteIniConf_For_QM();
//	}
//	else	{
//		fclose(fIn);
//	}
	
	getcwd(szCurDir, 512);
	
	Gen_QM_Rotamer_Data();

//	RunAllJobs();
	if( argc>=3 && !strcmp( "--prepare", argv[2] ) ) {
		exit(0);
	}
  //  RunAllJobsSynchronously( nJob, szExe_G09, (char**)szInputList, (char**)szOutputList, (char**)szDirList );
	printf("Completing... nCluster=%d nJob=%d\n", nCluster, nJob );

	chdir(szCurDir);
		
	fOut = fopen("all-rotamer.dat", "w");

	int Finished;
	char szLine[256];
	char path[PATH_MAX];

	for(i=0; i<nJob; i++)	{
    strcpy( path, szDirList[i] );
    strcat( path, "/" );
    strcat( path, szOutputList[i] );

		printf(  "Reading %s\n", szOutputList[i] );
		fIn = fopen( path, "r" ); // szOutputList[i], "r");
		Finished = To_Find_Tag(fIn, szOutput, " Optimization completed", szLine, 0);
		if( fIn) fclose(fIn);
		if(Finished < 0)	{	// un-successful optimization
			continue;
		}

		E_Rotamer = Extract_Rotamer_Coord_E( path ); // szOutputList[i]);
		if(E_Rotamer != 0.0)	{	// a valid energy
			fprintf(fOut, "E_Rotamer %.13E Conf %4d\nCoordinate\n", E_Rotamer, i+1);
			for(j=0; j<nAtom; j++)	{
				fprintf(fOut, "%14.6lf%14.6lf%14.6lf\n", Mol.x[j], Mol.y[j], Mol.z[j]);
			}
		}

		if(E_Rotamer < E_Rotamer_Min)	{
			Mol.WriteXYZ("QM-min.xyz");
			E_Rotamer_Min = E_Rotamer;
		}
	}

	fclose(fOut);

	system("/bin/tar -zcvf qm-ini.tgz qm-ini-*.xyz > /dev/null ; /bin/rm qm-ini-*.xyz");	// tar and delete pdb files

	fclose(fFile_Run_Log);

	return 0;
}

void Gen_QM_Rotamer_Data(void)
{
	FILE *fIn;
	int i, nAtom;
	char szName[256], szCmd[256], szOutput[]="mol-qm-rotamer.out";

	nAtom = Mol.nAtom;

	if(nCluster > MAX_JOB)	{
		Quit_With_Error_Msg("nCluster > MAX_JOB\nQuit\n");
	}
	
	nJob = 0;
	for(i=0; i<nCluster; i++)	{
		chdir(szCurDir);
		sprintf(szName, "qm-ini-%d.xyz", i+1);

		fIn = fopen(szName, "r");
		if(fIn == NULL)	{
			printf("Fail to read the configuration of %d rotamers.\n", i);
			break;
		}
		else	{
			fclose(fIn);
		}

		Mol.ReadXYZ(szName);
		
		sprintf(szNewDir, "%s/rotamer-%d", szCurDir, i+1);
		sprintf(szCmd, "rotamer-%d", i+1);
		make_directory(szCmd, 0700);
		chdir(szNewDir);
		Output_Gaussian_File();

		szDirList[i]   = strdup( szNewDir ); //nnstrcpy(szDirList[i], szNewDir);
		szInputList[i] = strdup( szGjfFile ); // strcpy(szInputList[i], szGjfFile);
		szOutputList[i]= strdup( szOutput ); // sprintf(szOutputList[i], "%s/%s", szNewDir, szOutput);
		nJob++;
	}
	
}

double Extract_Rotamer_Coord_E(char szName[])
{
	FILE *fIn;
	int Finished, j, iTmp, nAtom, ReadItem;
	char szLine[256], ErrorMsg[256], szStr_Energy[256], *ReadLine;
	double E_Rotamer;

	nAtom = Mol.nAtom;
	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Error to open file: %s\nQuit\n", szName);
		Quit_With_Error_Msg(ErrorMsg);
	}
	
	Finished = To_Find_Tag(fIn, szName, " Optimization completed", szLine, 0);
	if(Finished < 0)	{	// skip unfinished jobs
		printf(ErrorMsg, "Warning> Error to reading file: %s\n", szName);
		fclose(fIn);
		return 0.0;	// an invalid energy
	}
	
	To_Find_Tag(fIn, szName, " orientation:", szLine, 1);
	Skip_N_Line(fIn, szName, 4);	// including the coordinates of X1 and X2.
	
	for(j=0; j<nAtom; j++)	{
		if(feof(fIn))	{
			break;
		}
		fgets(szLine, 256, fIn);
		ReadItem = sscanf(szLine, "%d %d %d %lf %lf %lf", &iTmp, &iTmp, &iTmp, &(Mol.x[j]), &(Mol.y[j]), &(Mol.z[j]));
		if(ReadItem!=6)	{
			break;
		}
	}
	
	if(j==nAtom)	{	// correct number of atoms as expected
		fseek(fIn, 0, SEEK_SET);
		while(1)	{
			if(feof(fIn))	{
				break;
			}
			ReadLine = fgets(szLine, 256, fIn);
			if(ReadLine == NULL)	{
				break;
			}
			else	{
				if( (FindString(szLine, "SCF Done:  E(RHF)") >= 0) && (QM_Level == QM_LEVEL_HF) )	{
					sscanf(szLine+20, "%s", szStr_Energy);
					Replace_D_E(szStr_Energy);
					ReadItem = sscanf(szStr_Energy, "%lf", &E_Rotamer);
					if(ReadItem != 1)	{
						sprintf(ErrorMsg, "Fail to extract energy from : %s\nQuit\n", szStr_Energy);
						Quit_With_Error_Msg(ErrorMsg);
					}
				}
				else if( (FindString(szLine, "EUMP2 =") >= 0) && (QM_Level == QM_LEVEL_MP2) )	{
					sscanf(szLine+34, "%s", szStr_Energy);
					Replace_D_E(szStr_Energy);
					ReadItem = sscanf(szStr_Energy, "%lf", &E_Rotamer);
					if(ReadItem != 1)	{
						sprintf(ErrorMsg, "Fail to extract energy from : %s\nQuit\n", szStr_Energy);
						Quit_With_Error_Msg(ErrorMsg);
					}
				}
			}
			
		}
	}
	
	fclose(fIn);

	return E_Rotamer;
}

void Output_Gaussian_File(void)
{
	FILE *fOut;
	int i, nAtom;

	nAtom = Mol.nAtom;
	
	fOut = fopen(szGjfFile, "w");
	
	fprintf(fOut, "%s", szQM_Level_Rotamer_Cmd);
	for(i=0; i<nAtom; i++)	{
//		fprintf(fOut, "%c  %9.5lf %9.5lf %9.5lf\n", Mol.AtomName[i][0], Mol.x[i], Mol.y[i], Mol.z[i]);
		fprintf(fOut, "%s  %9.5lf %9.5lf %9.5lf\n", szElemName[i], Mol.x[i], Mol.y[i], Mol.z[i]);
	}
	fprintf(fOut, "\n");
	
	fclose(fOut);

//	char szName[256];
//	sprintf(szName, "ini-%d.pdb", n_Rotamer);
//	Mol.WriteXYZ(szName);

}

void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in QM_rotamer_scan.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}

int Read_Soft_DihedralList(void)
{
	FILE *fIn;
	int ReadItem, n_Conf=1;
	char szLine[256], *ReadLine;

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
			n_Conf *= n_State_List[n_Phi];
			n_Phi++;
		}
		else	{
			break;
		}
	}

	fclose(fIn);

	memcpy(Phi_Set_Save, Phi_Set, sizeof(double)*N_MAX_DIH*MAX_N_STATE);
	memcpy(n_State_List_Save, n_State_List, sizeof(int)*N_MAX_DIH);

	if(n_Conf >= MAX_ROTAMER)	{	// too many possible combinations, just randomly pick a part of them
		Ratio_Pick = 1.0 * MAX_ROTAMER/n_Conf;
	}
	else	{
		Ratio_Pick = 1.0;
	}

	Read_Tor_Para_1D_Fitting();

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
	double cg_sum=0.0;
	char szCharge[16];

	nAtom = Mol.nAtom;
	for(i=0; i<nAtom; i++)	{
		cg_sum += Mol.CG[i];
	}
	sprintf(szCharge, "%.0lf", cg_sum);

	netcharge_mol = atoi(szCharge);	// assuming the netcharge is an integer !!!
}


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
	char szEnergy[256], ErrorMsg[256];
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

int To_Find_Tag(FILE *fIn, char szFileName[], char szTag[], char szLine[], int Quit)
{
	char *ReadLine, ErrorMsg[256];

	if( !fIn ) return -1;

	while(1)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Fail to find the tag: %s in file %s\nQuit\n", szTag, szFileName);
			if(Quit)	{
				fclose(fIn);
				Quit_With_Error_Msg(ErrorMsg);
			}
			else	{
				break;
			}
		}

		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			sprintf(ErrorMsg, "Fail to find the tag: %s in file %s\nQuit\n", szTag, szFileName);
			if(Quit)	{
				fclose(fIn);
				Quit_With_Error_Msg(ErrorMsg);
			}
			else	{
				break;
			}
		}
		else	{
			if(FindString(szLine, szTag) >= 0)	{
				return 1;
			}
		}
	}

	return -1;
}

int Skip_N_Line(FILE *fIn, char szFileName[], int n)
{
	int i;
	char szLine[256], *ReadLine, ErrorMsg[256];

	for(i=0; i<n; i++)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Fail in Skip_N_Line(%s, %d).\nQuit\n", szFileName, n);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}

		ReadLine = fgets(szLine, 256, fIn);
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

#if 0
void Get_EXE_Path(char szExeFile[], char szExePath[])
{
	FILE *fIn;
	char szLine[256], szBuff[256], *ReadLine, ErrorMsg[256];
	int nLen;

	fIn = fopen(szMyPath, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file: %s\nQuit\n", szMyPath);
		Quit_With_Error_Msg(ErrorMsg);
	}

	nLen = strlen(szExeFile);
	while(1)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Fail to find the path of %s in %s\nQuit\n", szExeFile, szMyPath);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}

		ReadLine = fgets(szLine, 256, fIn);

		if(ReadLine == NULL)	{
			sprintf(ErrorMsg, "Fail to find the path of %s in %s\nQuit\n", szExeFile, szMyPath);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
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
		sprintf(ErrorMsg, "The file %s for %s doesn't exist.\nQuit\n", szExePath, szExeFile);
		Quit_With_Error_Msg(ErrorMsg);
	}
	else	{
		fclose(fIn);
	}
}
#endif


void Setup_QM_Level(void)
{
	FILE *fIn;
	char *ReadLine, szLine[256], szSubStr_1[256], szSubStr_2[256], ErrorMsg[256];
	char szKey_QM_MEM[256], szKey_QM_NPROC[256], szQM_Level_Opt[256], szQM_Level_Dimer_Opt[256], szQM_Level_ESP[256], szQM_Level_E_Monomer[256], szQM_Level_1D_Scan[256], szQM_Level_Rotamer[256];
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
				else if(strcmp(szSubStr_1,"QM_LEVEL_ROTAMER")==0)	{
					strcpy(szQM_Level_Rotamer, szSubStr_2);
					KeyCount++;
				}
			}
		}
	}

	fclose(fIn);

	Determine_QM_Level(szQM_Level_1D_Scan);

	if(KeyCount != 8)	{	// !!!
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


	sprintf(szQM_Level_Rotamer_Cmd, "%%mem=%s\n%%nproc=%s\n%sMol rotamer scan\n\n%d,1\n",
		szKey_QM_MEM, szKey_QM_NPROC, szQM_Level_Rotamer, netcharge_mol);
//	printf("%s", szQM_Level_Rotamer_Cmd);
}


void Determine_QM_Level(char szQMLevel[])
{
	char ErrorMsg[256];

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


#define IADD   453806245
#define IMUL   314159269
#define MASK   2147483647
#define SCALE  0.4656612873e-9

double rand(int &iseed)
{
	iseed = (iseed * IMUL + IADD) & MASK;
	return (iseed * SCALE);
}


double x_Best[MAX_ATOM*3];
nlopt_opt opt_Geo;
double Do_Geometry_Optimization(CMol *pMol)
{
	int i, nAtom, nDim, iPos, Return_Opt;
	double x[MAX_ATOM*3], x_Ini[MAX_ATOM], y_Ini[MAX_ATOM], z_Ini[MAX_ATOM], E_min=1.0E100;

	nAtom = pMol->nAtom;
	nDim = 3*nAtom;

	if(n_Phi <= 5)	{	// this step is not so expensive
		pMol->FullGeometryOptimization_LBFGS();
		E_min = pMol->Cal_E(0);
		return E_min;
	}

	// for a large molecule with many soft dihedrals, we prefer the optimizer in nlopt since it is around 8 times faster

	memcpy(x_Ini, pMol->x, sizeof(double)*nAtom);	// backup the structure
	memcpy(y_Ini, pMol->y, sizeof(double)*nAtom);
	memcpy(z_Ini, pMol->z, sizeof(double)*nAtom);

	opt_Geo = nlopt_create(NLOPT_LD_LBFGS, nDim);	// fast for non-polarizable model. Not applicable for drude model
	nlopt_set_min_objective(opt_Geo, func_Geo_Optimization, pMol);
	nlopt_set_ftol_rel(opt_Geo, 5E-12);

	for(i=0; i<nAtom; i++)	{	// to set up the initial parameters
		iPos = 3*i;
		x[iPos  ] = pMol->x[i];
		x[iPos+1] = pMol->y[i];
		x[iPos+2] = pMol->z[i];
	}
	Iter_Geo_Opt = 0;
	Return_Opt = nlopt_optimize(opt_Geo, x, &E_min);
	nlopt_destroy(opt_Geo);

	if( (Return_Opt == NLOPT_FTOL_REACHED) || (Return_Opt == NLOPT_SUCCESS) )	{	// success in optimization
		return E_min;
	}
	else	{	// try FullGeometryOptimization_LBFGS()
		memcpy(pMol->x, x_Ini, sizeof(double)*nAtom);	// restore the structure
		memcpy(pMol->y, y_Ini, sizeof(double)*nAtom);
		memcpy(pMol->z, z_Ini, sizeof(double)*nAtom);

		pMol->FullGeometryOptimization_LBFGS();
		E_min = pMol->Cal_E(0);
		return E_min;
	}
}

double func_Geo_Optimization(unsigned n, const double *x, double *grad, void *my_func_data)
{
	int i, iPos, nAtom;
	double E_Total, *gx, *gy, *gz, *xx, *yy, *zz;
	CMol* pMol;

	pMol = (CMol*)my_func_data;
	nAtom = pMol->nAtom;
	gx = pMol->grad_x;
	gy = pMol->grad_y;
	gz = pMol->grad_z;
	xx = pMol->x;
	yy = pMol->y;
	zz = pMol->z;

	if(grad)	{
		Iter_Geo_Opt ++;
		for(i=0; i<nAtom; i++)	{	// update the coordinate
			iPos = 3*i;
			xx[i] = x[iPos  ];
			yy[i] = x[iPos+1];
			zz[i] = x[iPos+2];
		}

		E_Total = pMol->Cal_E(0);

		for(i=0; i<nAtom; i++)	{	// get the gradient
			iPos = 3*i;
			grad[iPos  ] = gx[i];
			grad[iPos+1] = gy[i];
			grad[iPos+2] = gz[i];
		}
	}
	else	{
		E_Total = pMol->Cal_E(0);
	}

	return E_Total;
}

void Make_List_H_Acceptor(void)
{
	int i, j, nAtom;
	double *pMass, Mass;

	nAtom = Mol.nAtom;
	pMass = Mol.mass;

	nHDonor = nAcceptor = 0;
	memset(IsH, 0, sizeof(int)*nAtom);

	for(i=0; i<nAtom; i++)	{
		Mass = pMass[i];

		if( (Mass > 0.99) && (Mass < 1.01) )	{
			IsH[i] = 1;
		}
		else if( (Mass > 13.50) && (Mass < 14.50) )	{	// N
			Acceptor[nAcceptor] = i;
			nAcceptor++;
			GetBondList(i);
		}
		else if( (Mass > 15.50) && (Mass < 16.50) )	{	// O
			Acceptor[nAcceptor] = i;
			nAcceptor++;
			GetBondList(i);
		}
		else if( (Mass > 18.50) && (Mass < 19.50) )	{	// F
			Acceptor[nAcceptor] = i;
			nAcceptor++;
			GetBondList(i);
		}
	}

	for(i=0; i<nAcceptor; i++)	{
		for(j=0; j<Bond_Count[Acceptor[i]]; j++)	{
			if(IsH[ BondList[ Acceptor[i] ][j] ])	{
				Donor[nHDonor] = BondList[Acceptor[i]][j];
				nHDonor++;
			}
		}
	}
}

void GetBondList(int Idx)
{
	int i, nBond, iPos, *pBondList;

	nBond = Mol.nBond;
	pBondList = Mol.BondList;

	Bond_Count[Idx] = 0;
	for(i=0; i<nBond; i++)	{
		iPos = 2*i;
		if(pBondList[iPos] == Idx)	{
			BondList[Idx][Bond_Count[Idx]] = pBondList[iPos+1];
			Bond_Count[Idx]++;
		}
		if(pBondList[iPos+1] == Idx)	{
			BondList[Idx][Bond_Count[Idx]] = pBondList[iPos];
			Bond_Count[Idx]++;
		}
	}
	return;
}

void Count_Intramolecular_HBond(double x[], double y[], double z[], int nAtom, int& n_HBond, int& n_Far_HBond)
{
	int i, j, Atom_1, Atom_2, Dist_Cut=12;
	double dx, dy, dz, r, r_Cut_max=2.6, r_Cut_min=2.1;

	n_HBond = n_Far_HBond = 0;
	for(i=0; i<nAcceptor; i++)	{
		Atom_1 = Acceptor[i];
		for(j=0; j<nHDonor; j++)	{
			Atom_2 = Donor[j];
			dx = x[Atom_2] - x[Atom_1];
			dy = y[Atom_2] - y[Atom_1];
			dz = z[Atom_2] - z[Atom_1];
			r = sqrt(dx*dx + dy*dy + dz*dz);
printf("%f %f %f\n",  r , r_Cut_min, r_Cut_max );
			if( (r > r_Cut_min ) && (r < r_Cut_max) )	{
				printf("Count_Intramolecular_HBond(): r(%d, %d) = %.2lf dist=%.2lf \n", Atom_1+1, Atom_2+1, r, dist[Atom_1][Atom_2]);
				if(dist[Atom_1][Atom_2] > Dist_Cut)	{
					n_Far_HBond++;
				}
				n_HBond++;
			}
		}
	}
//	printf("\n");
}

void Exclude_Conf_Longdist_HBond(void)
{
	int i, Idx, n_HBond, n_Far_HBond;

	E_Min = 1.0E100;
	for(i=0; i<Counter; i++)	{
		Idx = Idx_Sorted[i];
		Count_Intramolecular_HBond(x_List[Idx], y_List[Idx], z_List[Idx], Mol.nAtom, n_HBond, n_Far_HBond);
		if(n_Far_HBond > 0)	{	// excluded
			printf("Excluding rotamer %d (case 1)\n", i );
			Rotamer_QM[i] = 0;
		}
		else	{
			if(E_Min > E_Conf[Idx])	{	// to find the global minimum
				E_Min = E_Conf[Idx];
			}
		}
	}

	for(i=0; i<Counter; i++)	{
		if(Rotamer_QM[i])	{
			Idx = Idx_Sorted[i];
			if(E_Conf[Idx] > (E_Min+E_CUTOFF))	{	// delete 
			  printf("Excluding rotamer %d (case 2)\n", i );
				Rotamer_QM[i] = 0;
			}
		}
	}
}



void SortByEnergy(void)
{
	int i, j, i_Tmp;

	for(i=0; i<Counter; i++)	{
		Idx_Sorted[i] = i;
	}

	for(i=0; i<Counter; i++)	{
		for(j=i+1; j<Counter; j++)	{
			if(E_Conf[Idx_Sorted[i]] > E_Conf[Idx_Sorted[j]])	{	// to exchange i <-> j
				i_Tmp = Idx_Sorted[i];	Idx_Sorted[i] = Idx_Sorted[j];	Idx_Sorted[j] = i_Tmp;
			}
		}
	}
}

void ClusteringRotamers(void)
{
	int i, j, k, Idx_1, Idx_2, nAtom;
	double rmsd_1, rmsd_2, rmsd, x_Tmp[MAX_ATOM];

	nCluster = 0;
	nAtom = Mol.nAtom;
	printf("ClusteringRotamers: counter=%d\n", Counter );

	for(i=0; i<Counter; i++)	{
		if(Rotamer_QM[i]==0)	{
			continue;
		}
		Idx_1 = Idx_Sorted[i];
		ClusterList[nCluster] = Idx_1;
//		printf("%3d  %d\n", i+1, ClusterList[nCluster]);
		nCluster++;

#pragma omp parallel for private(Idx_2,rmsd_1,rmsd_2,rmsd,x_Tmp,k)
		for(j=i+1; j<Counter; j++)	{
			if(Rotamer_QM[j]==0)	{
				continue;
			}
			Idx_2 = Idx_Sorted[j];

			rmsd_1 = CalRMSD(x_List[Idx_1], y_List[Idx_1], z_List[Idx_1], x_List[Idx_2], y_List[Idx_2], z_List[Idx_2], nAtom);
			for(k=0; k<nAtom; k++)	{
				x_Tmp[k] = -x_List[Idx_1][k];	// mirror
			}
			rmsd_2 = CalRMSD(x_Tmp, y_List[Idx_1], z_List[Idx_1], x_List[Idx_2], y_List[Idx_2], z_List[Idx_2], nAtom);
			rmsd = min(rmsd_1, rmsd_2);

			if(rmsd < RMSD_CUTOFF)	{
			printf("Excluding rotamer %d (case 3)\n", j );
				Rotamer_QM[j] = 0;
			}
			if( (fabs(E_Conf[Idx_2]-E_Conf[Idx_1]) < 1.0E-5) && (rmsd < 1.5) )	{	// symmetry configurations ?
				printf("Excluding rotamer %d (case 4)\n", j );
				Rotamer_QM[j] = 0;
			}
//			printf( "<T:%d> - %d\n", omp_get_thread_num(), j );
		}
	}

	nCluster = min(nCluster, N_ROTAMER_PICK);
}

void WriteIniConf_For_QM(void)
{
	char szName[256];
	int i, Idx, nAtom;
	
	nAtom = Mol.nAtom;

	printf("Writing %d clusters\n", nCluster );
	for(i=0; i<nCluster; i++)	{
		Idx = ClusterList[i];
		memcpy(Mol.x, x_List[Idx], sizeof(double)*nAtom);
		memcpy(Mol.y, y_List[Idx], sizeof(double)*nAtom);
		memcpy(Mol.z, z_List[Idx], sizeof(double)*nAtom);
		sprintf(szName, "qm-ini-%d.xyz", i+1);
		printf(" Saving [%s]\n", szName );
		Mol.WriteXYZ(szName);
	}
}

void Read_Tor_Para_1D_Fitting(void)
{
	FILE *fIn;
	char szName[256], ErrorMsg[256];
	int j, Idx_Phi, ReadItem;

	for(Idx_Phi=0; Idx_Phi<n_Phi; Idx_Phi++)	{
		sprintf(szName, "torsion-para-%d.dat", Idx_Phi+1);
		fIn = fopen(szName, "r");
		if(fIn == NULL)	{
			sprintf(ErrorMsg, "Fail to open file: %s.\nQuit\n", szName);
			Quit_With_Error_Msg(ErrorMsg);
		}

		for(j=0; j<(N_TOR_PARA-1); j++)	{
			ReadItem = fscanf(fIn, "%lf", &(Para_Tor_List[Idx_Phi][j]));
			if(ReadItem != 1)	{
				break;
			}
		}
		if(j!=(N_TOR_PARA-1))	{	// not a valid record
			fclose(fIn);
			sprintf(ErrorMsg, "Error in reading file: %s.\nQuit\n", szName);
			Quit_With_Error_Msg(ErrorMsg);
		}

		fclose(fIn);
	}

	Assign_Torsion_Parameters();
}

void Assign_Torsion_Parameters(void)
{
	int i_Phi, Idx;
	double *pPara_k_Dih, *pPara_phi;

	for(i_Phi=0; i_Phi<n_Phi; i_Phi++)	{
		Idx = IdxDihSelect[i_Phi];
		pPara_k_Dih = Mol.Para_k_Dih[Idx];
		pPara_phi = Mol.Para_phi[Idx];

		pPara_k_Dih[1] = Para_Tor_List[i_Phi][0];	// k1
		pPara_phi[1] = Para_Tor_List[i_Phi][1];		// sigma_1
		pPara_k_Dih[2] = Para_Tor_List[i_Phi][2];	// k2
		pPara_phi[2] = Para_Tor_List[i_Phi][3];		// sigma_2
		pPara_k_Dih[3] = Para_Tor_List[i_Phi][4];	// k3
		pPara_phi[3] = Para_Tor_List[i_Phi][5];		// sigma_3
		pPara_k_Dih[4] = Para_Tor_List[i_Phi][6];	// k4
		pPara_phi[4] = Para_Tor_List[i_Phi][7];		// sigma_4
		pPara_k_Dih[5] = 0.0;						// k5
		pPara_phi[5] = 0.0;							// sigma_5

		pPara_k_Dih[6] = Para_Tor_List[i_Phi][8];	// k6
		pPara_phi[6] = Para_Tor_List[i_Phi][9];		// sigma_6
	}
}


double CalRMSD(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum)
{
	double RMSD, midxa, midya, midza, midxb, midyb, midzb;
	double xa_tmp[MAX_ATOM], ya_tmp[MAX_ATOM], za_tmp[MAX_ATOM];
	double xb_tmp[MAX_ATOM], yb_tmp[MAX_ATOM], zb_tmp[MAX_ATOM];
	int i;

	memcpy(xa_tmp+1, xa, sizeof(double)*AtomNum);
	memcpy(ya_tmp+1, ya, sizeof(double)*AtomNum);
	memcpy(za_tmp+1, za, sizeof(double)*AtomNum);

	memcpy(xb_tmp+1, xb, sizeof(double)*AtomNum);
	memcpy(yb_tmp+1, yb, sizeof(double)*AtomNum);
	memcpy(zb_tmp+1, zb, sizeof(double)*AtomNum);

	midxa=midya=midza=midxb=midyb=midzb=0.0;
	for(i=1; i<=AtomNum; i++)	{
		midxa += (xa_tmp[i]);
		midya += (ya_tmp[i]);
		midza += (za_tmp[i]);

		midxb += (xb_tmp[i]);
		midyb += (yb_tmp[i]);
		midzb += (zb_tmp[i]);
	}
	midxa/=AtomNum;
	midya/=AtomNum;
	midza/=AtomNum;
	midxb/=AtomNum;
	midyb/=AtomNum;
	midzb/=AtomNum;

	for(i=1; i<=AtomNum; i++)	{
		xa_tmp[i]-=midxa;
		ya_tmp[i]-=midya;
		za_tmp[i]-=midza;
		xb_tmp[i]-=midxb;
		yb_tmp[i]-=midyb;
		zb_tmp[i]-=midzb;
	}

	quatfit(xa_tmp, ya_tmp, za_tmp, xb_tmp, yb_tmp, zb_tmp, AtomNum);
	RMSD=rmsfit(xa_tmp, ya_tmp, za_tmp, xb_tmp, yb_tmp, zb_tmp, AtomNum);

	for(i=1; i<=AtomNum; i++)	{	// generate the rotated coodinates
		xb[i-1]=midxb+xb_tmp[i];
		yb[i-1]=midyb+yb_tmp[i];
		zb[i-1]=midzb+zb_tmp[i];
	}


	return RMSD;
}

double rmsfit (double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n)
{
	int ia;
	double rmsfit, dist2, dx, dy, dz;

	rmsfit=0.0;
	for(ia=1; ia<=n; ia++)	{	//store xi-x0 into xStd[]
		dx = x1[ia] - x2[ia];
		dy = y1[ia] - y2[ia];
		dz = z1[ia] - z2[ia];
		dist2 = dx*dx + dy*dy + dz*dz;
		rmsfit += dist2;
	}
	rmsfit = sqrt(rmsfit/n);

	return rmsfit;
}

void quatfit (double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n)
{
	double xxyx,xxyy,xxyz, xyyx,xyyy,xyyz, xzyx,xzyy,xzyz;
	double xrot,yrot,zrot, work1[5],work2[5], q[5],d[5], rot[4][4], c[5][5], v[5][5];
	int ia;

	xxyx = 0.0;
	xxyy = 0.0;
	xxyz = 0.0;
	xyyx = 0.0;
	xyyy = 0.0;
	xyyz = 0.0;
	xzyx = 0.0;
	xzyy = 0.0;
	xzyz = 0.0;

	for(ia=1; ia<=n; ia++)	{
		xxyx = xxyx + x1[ia]*x2[ia];
		xxyy = xxyy + y1[ia]*x2[ia];
		xxyz = xxyz + z1[ia]*x2[ia];
		xyyx = xyyx + x1[ia]*y2[ia];
		xyyy = xyyy + y1[ia]*y2[ia];
		xyyz = xyyz + z1[ia]*y2[ia];
		xzyx = xzyx + x1[ia]*z2[ia];
		xzyy = xzyy + y1[ia]*z2[ia];
		xzyz = xzyz + z1[ia]*z2[ia];
	}

	c[1][1] = xxyx + xyyy + xzyz;
	c[1][2] = xzyy - xyyz;
	c[2][2] = xxyx - xyyy - xzyz;
	c[1][3] = xxyz - xzyx;
	c[2][3] = xxyy + xyyx;
	c[3][3] = xyyy - xzyz - xxyx;
	c[1][4] = xyyx - xxyy;
	c[2][4] = xzyx + xxyz;
	c[3][4] = xyyz + xzyy;
	c[4][4] = xzyz - xxyx - xyyy;

	jacobi(4,4,c,d,v,work1,work2);
	q[1] = v[1][4];
	q[2] = v[2][4];
	q[3] = v[3][4];
	q[4] = v[4][4];

	rot[1][1] = q[1]*q[1] + q[2]*q[2] - q[3]*q[3] - q[4]*q[4];
	rot[2][1] = 2.0 * (q[2] * q[3] - q[1] * q[4]);
	rot[3][1] = 2.0 * (q[2] * q[4] + q[1] * q[3]);
	rot[1][2] = 2.0 * (q[3] * q[2] + q[1] * q[4]);
	rot[2][2] = q[1]*q[1] - q[2]*q[2] + q[3]*q[3] - q[4]*q[4];
	rot[3][2] = 2.0 * (q[3] * q[4] - q[1] * q[2]);
	rot[1][3] = 2.0 * (q[4] * q[2] - q[1] * q[3]);
	rot[2][3] = 2.0 * (q[4] * q[3] + q[1] * q[2]);
	rot[3][3] = q[1]*q[1] - q[2]*q[2] - q[3]*q[3] + q[4]*q[4];

	for(ia=1; ia<=n; ia++)	{
		xrot = x2[ia]*rot[1][1] + y2[ia]*rot[1][2] + z2[ia]*rot[1][3];
		yrot = x2[ia]*rot[2][1] + y2[ia]*rot[2][2] + z2[ia]*rot[2][3];
		zrot = x2[ia]*rot[3][1] + y2[ia]*rot[3][2] + z2[ia]*rot[3][3];
		x2[ia] = xrot;
		y2[ia] = yrot;
		z2[ia] = zrot;
	}

	return;
}


#define ZERO			(1.0E-100)
void jacobi (int n, int np, double a[5][5], double d[5], double v[5][5], double b[5], double z[5])
{
	int i,j,k;
	int ip,iq;
	int nrot,maxrot;
	double sm,tresh,s,c,t;
	double theta,tau,h,g,p;

	maxrot = 100;
	nrot = 0;

	for(ip=1; ip<=n; ip++)	{
		for(iq=1; iq<=n; iq++)	{
			v[ip][iq] = 0.0;
		}
		v[ip][ip] = 1.0;
	}
	for(ip=1; ip<=n; ip++)	{
		b[ip] = a[ip][ip];
		d[ip] = b[ip];
		z[ip] = 0.0;
	}


	for(i=1; i<=maxrot; i++)	{
		sm=0.0;
		for(ip=1; ip<=(n-1); ip++)	{
			for(iq=ip+1; iq<=n; iq++)	{
				sm = sm + fabs(a[ip][iq]);
			}
		}

		if(sm < ZERO)	{	//converge
			break;
		}
		if(i < 4)	{
			tresh = 0.2*sm / (n*n);
		}
		else	{
			tresh=0.0;
		}

		for(ip=1; ip<=(n-1); ip++)	{
			for(iq=ip+1; iq<=n; iq++)	{
				g = 100.0 * fabs(a[ip][iq]);
				if((i > 4) && (fabs(g) < ZERO))	{
					a[ip][iq] = 0.0;
				}
				else if(fabs(a[ip][iq]) > tresh)	{
					h = d[iq] - d[ip];
					if(fabs(g) < ZERO)	{
						t = a[ip][iq] / h;
					}
					else	{
						theta = 0.5*h / a[ip][iq];
						t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0)  t = -t;
					}

					c = 1.0 / sqrt(1.0+t*t);
					s = t * c;
					tau = s / (1.0+c);
					h = t * a[ip][iq];
					z[ip] = z[ip] - h;
					z[iq] = z[iq] + h;
					d[ip] = d[ip] - h;
					d[iq] = d[iq] + h;
					a[ip][iq] = 0.0;

					for(j=1; j<=(ip-1); j++)	{
						g = a[j][ip];
						h = a[j][iq];
						a[j][ip] = g - s*(h+g*tau);
						a[j][iq] = h + s*(g-h*tau);
					}
					for(j=ip+1; j<=(iq-1); j++)	{
						g = a[ip][j];
						h = a[j][iq];
						a[ip][j] = g - s*(h+g*tau);
						a[j][iq] = h + s*(g-h*tau);
					}
					for(j=iq+1; j<=n; j++)	{
						g = a[ip][j];
						h = a[iq][j];
						a[ip][j] = g - s*(h+g*tau);
						a[iq][j] = h + s*(g-h*tau);
					}
					for(j=1; j<=n; j++)	{
						g = v[j][ip];
						h = v[j][iq];
						v[j][ip] = g - s*(h+g*tau);
						v[j][iq] = h + s*(g-h*tau);
					}
					nrot = nrot + 1;
				}
			}
		}
		for(ip=1; ip<=n; ip++)	{
			b[ip] = b[ip] + z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	if(nrot >= maxrot)	{
		printf("JACOBI  --  Matrix Diagonalization not Converged.\n");
	}

	for(i=1; i<= (n-1); i++)	{
		k = i;
		p = d[i];
		for(j=i+1; j<=n; j++)	{
			if(d[j] < p)	{
				k = j;
				p = d[j];
			}
		}
		if(k != i)	{
			d[k] = d[i];
			d[i] = p;
			for(j=1; j<=n; j++)	{
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}
#undef	ZERO

void ReadElementList(void)
{
	int n_Atom_In_Mol, i;
	FILE *fIn;
	char ErrorMsg[256];

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


