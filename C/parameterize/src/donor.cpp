/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#else
#include <unistd.h>	//for linux
#endif

#include "nlopt.h"
#include "ff.h"

#define MAX_ATOM_MOL	(500)
//#define N_MAX_WAT	(100)
#define MAX_N_CONF	(100)

#define b0_H_O_EXP	(0.9572)	// tip3 and exp
#define theta0_H_O_H_EXP	(104.52)	// tip3, exp?

//#define b0_H_O_MP2	(0.958825)
//#define theta0_H_O_H_MP2	(104.2242)

#define b0_H_O_HF	(0.9474)		// optimized by HF/6-31G*
#define theta0_H_O_H_HF	(105.5395)

#define HARTREE_To_KCAL	(627.509469)
#define radian	(57.29577951308232088)
#define radianInv	(0.017453292519943295)
#define PI	(3.14159265358979323846)
#define PI2	(6.28318530717958647692)

#ifndef MINAB
#define MINAB 
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif



//start	data and functions related with job scheduler
#define MAX_JOB	(1024)

enum { IDLE, BUSY };
enum { Job_ToRun, Job_Running, Job_Done };

int nJob=0, nJobDone=0, nWorker=0, nCore_Per_Node=1;
int JobStatus[MAX_JOB], DonorH_List[MAX_JOB], HeavyAtom_List[MAX_JOB], DonorIdx_List[MAX_JOB], Theta_List[MAX_JOB];

char *szDirList[MAX_JOB], *szInputList[MAX_JOB], *szOutputList[MAX_JOB];

double Para_MM_Save[MAX_JOB][2];
const char szGAUSS_SCRDIR[]="GAUSS_SCRDIR";
char szGAUSS_SCRDIR_Base[PATH_MAX];

void RunAllJobs(void);
int Get_First_Available_Worker(void);
void Exatract_All_Info_QM_1D_Scan(void);
//end	data and functions related with job scheduler


FILE *fFile_Run_Log;
char szLogName[256]="log-gen-HBond-A.txt";
char szQMParaFile[]="../QM-para.txt";
char szMMForceField[]="mol.prm";

char szQM_Level_Int_E_Dimer_Cmd[512]="%mem=200MB\n%nproc=2\n#P HF/6-31G* nosymm SCF=Tight counterpoise=2\n\nMol-Water Dimer\n\n0,1 0,1 0,1 \n";	// kbt - small queue
char szQM_Level_Int_E_Dimer_Opt_Cmd[512]="%mem=200MB\n%nproc=2\n#P HF/6-31G* nosymm SCF=Tight counterpoise=2 opt=ModRedundant\n\nMol-Water Dimer\n\n0,1 0,1 0,1 \n";	// kbt - small queue
char szQM_Level_Int_E_Single_Mol_Cmd[512]="%mem=200MB\n%nproc=2\n#P HF/6-31G* nosymm SCF=Tight\n\nMol/Water\n\n0,1\n";	// kbt - small queue
char szQM_Level_Int_E_Single_Wat_Cmd[512]="%mem=200MB\n%nproc=2\n#P HF/6-31G* nosymm SCF=Tight\n\nMol/Water\n\n0,1\n";	// kbt - small queue

char szExe_G09[256];
char szExe_Cale_MM[256];
char szMyPath[]="../mypath.txt";
char szMolWat_E_Int[]="E-mol-wat.txt";
char szCGList[]="cg-list.txt";
char szCurDir[256];
char szWorkDir[256];

double x_Wat_MM[3], y_Wat_MM[3], z_Wat_MM[3];	// to store the optimized MM water configurations
double x_Wat_QM[3], y_Wat_QM[3], z_Wat_QM[3];	// to store the optimized MM water configurations

int Idx_Donor, ActiveAtom, ActiveHAtom, Neighbor;


#define N_PARA	(2)
double Para_List[N_PARA], Grad_List[N_PARA], Para_Best[N_PARA], Low_Bound[N_PARA], Up_Bound[N_PARA];
char szElemName[MAX_ATOM_MOL][8];

int n_Atom_In_Mol;
int ResIDMol[MAX_ATOM_MOL];
int Idx_O_Wat, Idx_H1_Wat, Idx_H2_Wat;

double x_ReadIn[MAX_ATOM_MOL], y_ReadIn[MAX_ATOM_MOL], z_ReadIn[MAX_ATOM_MOL];

double x_Opt_QM[MAX_ATOM_MOL], y_Opt_QM[MAX_ATOM_MOL], z_Opt_QM[MAX_ATOM_MOL];
double r_Emin_MM=0.0, r_Emin_QM=0.0, E_Min_MM=0.0, E_Min_QM=0.0;

char szNameMolCRD[256];
char szAtomNameMol[MAX_ATOM_MOL][8], szResNameMol[MAX_ATOM_MOL][8];
double *x_Mol, *y_Mol, *z_Mol;
double E_Min=1.0E20;

double H1_Wat_x, H1_Wat_y, H1_Wat_z;
double O_Wat_x, O_Wat_y, O_Wat_z;
double H2_Wat_x, H2_Wat_y, H2_Wat_z;

double E_Dimer_Int, E_Dimer, E_Dimer_Min, E_Dimer_Far;

#define R_DUMMY_1	(1.5)
#define R_DUMMY_2	(1.5)
#define R_DUMMY_3	(5.0)
double x_Dummy_1, y_Dummy_1, z_Dummy_1, x_Dummy_1_Save, y_Dummy_1_Save, z_Dummy_1_Save, x_Dummy_1_Org, y_Dummy_1_Org, z_Dummy_1_Org;
double x_Dummy_2, y_Dummy_2, z_Dummy_2, x_Dummy_2_Save, y_Dummy_2_Save, z_Dummy_2_Save, x_Dummy_2_Org, y_Dummy_2_Org, z_Dummy_2_Org;
double x_Dummy_3, y_Dummy_3, z_Dummy_3, x_Dummy_3_Save, y_Dummy_3_Save, z_Dummy_3_Save, x_Dummy_3_Org, y_Dummy_3_Org, z_Dummy_3_Org;
double v_x_X1_H, v_y_X1_H, v_z_X1_H, v_x_X1_H_Save, v_y_X1_H_Save, v_z_X1_H_Save;
double Para_Save_for_QM[N_PARA];
int theta, theta_Save;

CForceField ForceField;
CMol Mol, Dimer;

int netcharge_mol;
void Get_Netcharge_From_Xpsf(void);

char szName_Gaussian[256]="qm.gjf";
char szFileElem[]="elem-list.txt";
char szFilePsf[]="mol.xpsf";
char szFileDimerPsf[]="mol-wat.xpsf";
char szFileCrd[]="mol-opt.xyz";

int Is_N[MAX_ATOM_MOL], Is_O[MAX_ATOM_MOL], Is_F[MAX_ATOM_MOL], Is_S[MAX_ATOM_MOL];
int Bond_Count[MAX_ATOM_MOL], BondList[MAX_ATOM_MOL][4], BondWithH[MAX_ATOM_MOL][4], nCount_H[MAX_ATOM_MOL];


void ReadMolCRD(void);
void GenerateWater(void);
double Normalizex(double& vect_x, double& vect_y, double& vect_z);
void OutputCRD(void);
void Generate_Gaussian_Input(void);
void Generate_Ini_Configuration(void);
double Optimize_Water_Pose(void);
void Output_Optimize_Water_Pose_QM(void);
double Optimize_Water_Pose_1D(double x_H1, double y_H1, double z_H1);	// optimize the water position along a line

void Quit_With_Error_Msg(char szMsg[]);
void Make_List_H_Acceptor(void);
void GetBondList(int Idx);
int Is_N_Planar(int Idx, double& x, double& y, double& z);

void Generate_Dummy_Atoms(int Idx);
void Rotate_Dummy_Atoms(int Idx, double theta);
void SaveDummyAtoms(void);
void RestoreDummyAtoms(void);

double Extract_E_Gaussian_Output(char szFileName[], char szTag[]);
double Get_E_Mol_Water_Far_QM();
void Get_EXE_Path(char szExeFile[], char szExePath[]);
void Cal_r_Dependent_E(int MM_Orientation);
double Cal_Dist_Two_Atoms(int ia, int ib);

// to generate xyz with r, theta, phi based on three atoms, A, B, C
void Gen_xyz(double x_A, double y_A, double z_A, 
			 double x_B, double y_B, double z_B, 
			 double x_C, double y_C, double z_C, 
			 double r, double theta, double phi, 
			 double& x_D, double& y_D, double& z_D);


void Replace_NewLine(char szBuff[]);
int Split_Into_Two_String(char szBuff[], char str1[], char str2[]);
void Setup_QM_Level(void);

int To_Find_Tag(FILE *fIn, char szFileName[], char szTag[], char szLine[]);
int Try_To_Find_Tag(FILE *fIn, char szFileName[], char szTag[], char szLine[]);
int Skip_N_Line(FILE *fIn, char szFileName[], int n);
void Replace_D_E(char szStr[]);


int Iteration;
static int FailCount=0;
double Callback_Eval_Gradient(unsigned n, const double *x, double *grad, void *my_func_data);
double CalObjectiveFunction(void);
double CalGradient(void);

extern int FindString(char szBuff[], char szTag[]);

void ReadPreviousCG(void);

double rmsfit(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);
void quatfit(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);
void jacobi(int n, int np, double a[5][5], double d[5], double v[5][5], double b[5], double z[5]);
double CalRMSD(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);
double Dot_Product(double x_a, double y_a, double z_a, double x_b, double y_b, double z_b);
void Cross_Product(double x_a, double y_a, double z_a, double x_b, double y_b, double z_b, double& x_Product, double& y_Product, double& z_Product);

void Extract_All_E_Dimer(void);


int main(int argc, char *argv[])
{
  int complete=0;
	int i, j, Phi=0, DonorSelect = -1, Count;
	int ToGenerate_Water;
	FILE *fOut;
	char szCmd[256], *szEnv;

  if( argc==2 && !strcmp( argv[1], "--complete" ) ) {
		printf("Phase: Completion\n" ); complete = 1;
	}
	else {
		printf("Phase: Preparation\n" ); complete = 0;
	}

  timebomb();

	getcwd(szCurDir, 256);

	fFile_Run_Log = fopen(szLogName, "a+");


	fseek(fFile_Run_Log, SEEK_END, 0);
	Get_EXE_Path( "BIN_G09", szExe_G09);
	
	ForceField.ReadForceField(szMMForceField);
	Dimer.ReadPSF(szFileDimerPsf, 1);
	Dimer.AssignForceFieldParameters(&ForceField);
	
	Mol.ReadPSF(szFilePsf, 0);
	Mol.AssignForceFieldParameters(&ForceField);
	Mol.ReadXYZ(szFileCrd);

	for(i=0; i<Mol.nAtom; i++)	{	// Readin QM optimized structure
		x_ReadIn[i] = Mol.x[i];
		y_ReadIn[i] = Mol.y[i];
		z_ReadIn[i] = Mol.z[i];
	}
	
	n_Atom_In_Mol = Mol.nAtom;
	ReadPreviousCG();
	Get_Netcharge_From_Xpsf();

	for(i=0; i<Mol.nAtom; i++)	{	// !!! the number of atom
		Dimer.x[i] = Mol.x[i];
		Dimer.y[i] = Mol.y[i];
		Dimer.z[i] = Mol.z[i];

		x_Opt_QM[i] = Mol.x[i];
		y_Opt_QM[i] = Mol.y[i];
		z_Opt_QM[i] = Mol.z[i];
	}
	x_Mol = Mol.x;
	y_Mol = Mol.y;
	z_Mol = Mol.z;

	Idx_O_Wat = Mol.nAtom;
	Idx_H1_Wat=Idx_O_Wat+1;
	Idx_H2_Wat=Idx_O_Wat+2;

	Setup_QM_Level();

	Make_List_H_Acceptor();


//	if(argc == 3)	{	// only to generate the list of donors
//		DonorSelect = atoi(argv[2]);
//	}
//	else	{
//		DonorSelect = -1;	// a null donor
//	}

	Idx_Donor = 0;
	for(i=0; i<n_Atom_In_Mol; i++)	{
		ToGenerate_Water = 0;
		Neighbor = 0xFFFFFFF;	// an invalid number
		
		if(Is_N[i] && (nCount_H[i]>=0) )	{
			ToGenerate_Water = 1;
		}
		else if(Is_O[i] && (BondWithH[i]>=0) )	{
			ToGenerate_Water = 1;
		}
		else if(Is_F[i] && (BondWithH[i]>=0) )	{
			ToGenerate_Water = 1;
		}
//		else if(Is_S[i] && (Bond_Count[i]<=2) )	{
//			ToGenerate_Water = 1;
//		}
		
		for(Count=0; Count<nCount_H[i]; Count++)	{
			if(ToGenerate_Water)	{
//				if( (Idx_Donor+1) != DonorSelect )	{	// only run for one selected donor
//					Idx_Donor++;
//					continue;
//				}
				
				ActiveAtom = i;
				ActiveHAtom = BondWithH[i][Count];
				
				Generate_Dummy_Atoms(i);
				SaveDummyAtoms();
				x_Dummy_1_Org = x_Dummy_1;	y_Dummy_1_Org = y_Dummy_1;	z_Dummy_1_Org = z_Dummy_1;
				x_Dummy_2_Org = x_Dummy_2;	y_Dummy_2_Org = y_Dummy_2;	z_Dummy_2_Org = z_Dummy_2;
				x_Dummy_3_Org = x_Dummy_3;	y_Dummy_3_Org = y_Dummy_3;	z_Dummy_3_Org = z_Dummy_3;
				
				//				sprintf(szWorkDir, "%s/dat-donor-%d-%d", szCurDir, Idx_Donor+1, ActiveAtom+1);
				sprintf(szWorkDir, "%s/dat-donor-%d-%d", szCurDir, Idx_Donor+1, ActiveHAtom+1);
				mkdir( szWorkDir, 0700 );
				chdir(szWorkDir);
				
				E_Min_MM = 1.0E100;
				r_Emin_MM = 1.0E100;
				
				theta = 0;
				for(Phi=0; Phi<360; Phi+=60)	{
					Para_List[0] = 1.9;
					Para_List[1] = ((Phi+15)%360)*radianInv;
					
					Optimize_Water_Pose();
					
					if( E_Dimer_Min < E_Min_MM )	{	// save the water with lowest interaction energy
						E_Min_MM = E_Dimer_Min;
						r_Emin_MM = Para_Best[0];
						SaveDummyAtoms();
						memcpy(Para_Save_for_QM, Para_Best, sizeof(double)*N_PARA);
						memcpy(Para_List, Para_Best, sizeof(double)*N_PARA);
						theta_Save = theta;
					}
					
					char szName[256];
					sprintf(szName, "theta-%d-cluster-%d.xyz", theta, Phi/60 + 1);
					Dimer.WriteXYZ(szName);
				}
				
				if( E_Min_MM > -1.0 )	{	// to generate water along a perp orientation
					for(theta=-30; theta<=30; theta+=60)	{	// only -30 and 30
						Rotate_Dummy_Atoms(i, 1.0*theta);	// generate dummy atoms 
						
						for(Phi=0; Phi<360; Phi+=60)	{
							Para_List[0] = 1.9;
							Para_List[1] = ((Phi+15)%360)*radianInv;
							
							Optimize_Water_Pose();
							
							if( E_Dimer_Min < E_Min_MM )	{	// save the water with lowest interaction energy
								E_Min_MM = E_Dimer_Min;
								r_Emin_MM = Para_Best[0];
								SaveDummyAtoms();
								memcpy(Para_Save_for_QM, Para_Best, sizeof(double)*N_PARA);
								memcpy(Para_List, Para_Best, sizeof(double)*N_PARA);
								theta_Save = theta;
							}
							
							char szName[256];
							sprintf(szName, "theta-%d-cluster-%d.xyz", theta, Phi/60 + 1);
							Dimer.WriteXYZ(szName);
							
						}				
					}
				}
				
				RestoreDummyAtoms();
				fprintf(fFile_Run_Log, "mol_H_donor.cpp> %2d H %2d donor: %2d atom %3s MM, E_min = %7.4lf r_min = %7.4lf  ", 
					Idx_Donor+1, ActiveHAtom+1, i+1, szElemName[i], E_Min_MM, r_Emin_MM);
				fprintf(fFile_Run_Log, "theta = %3d  r = %7.3lf phi = %7.3lf\n", 
					theta_Save, Para_Save_for_QM[0], Para_Save_for_QM[1]);
				
				if(E_Min_MM > -2.0)	{	// skip this acceptor if the interaction is too weak
          fprintf(fFile_Run_Log, "mol_H_donor.cpp> Interaction too weak\n" );
					Idx_Donor++;
					continue;
				}
				
				
				fOut = fopen("mol-waters.pdb", "w");
				for(j=0; j<n_Atom_In_Mol; j++)	{
					fprintf(fOut, "ATOM%7d  %-3s %-3s A%4d    %8.3lf%8.3lf%8.3lf  1.00  0.00\n", 
						j+1, Dimer.AtomName[j], "MOL", 1, Dimer.x[j], Dimer.y[j], Dimer.z[j]);
				}
				fflush(fOut);
				
				sprintf(szWorkDir, "%s/dat-donor-%d-%d/QM", szCurDir, Idx_Donor+1, ActiveHAtom+1);
        printf("Making [%s]\n", szWorkDir );
				mkdir( szWorkDir, 0700 );
				chdir(szWorkDir);
				
				szDirList[nJob] = strdup( szWorkDir);
				DonorH_List[nJob] = ActiveHAtom;
				HeavyAtom_List[nJob] = ActiveAtom;
				Theta_List[nJob] = theta_Save;
				Para_MM_Save[nJob][0] = Para_Save_for_QM[0];
				Para_MM_Save[nJob][1] = Para_Save_for_QM[1];
				DonorIdx_List[nJob] = Idx_Donor;
				
				Output_Optimize_Water_Pose_QM();	// to write input file for Gaussian 
				nJob++;
				
				Idx_Donor++;
			}
		}
		
	}
	
//	OutputCRD();
//	Generate_Gaussian_Input();


//	RunAllJobs();

  if( !complete ){ 
		fclose(fFile_Run_Log);
		exit(0); 
	}
//  RunAllJobsSynchronously( nJob, szExe_G09, (char**) szInputList, (char**) szOutputList, (char**) szDirList );


	chdir(szCurDir);

	E_Dimer_Far = Get_E_Mol_Water_Far_QM();
	Extract_All_E_Dimer();

	fclose(fFile_Run_Log);

	return 0;
}

void Extract_All_E_Dimer(void)
{
	int i, ReadItem, AtomNum, iTmp;
	char szName_Result[256], szLine[256], szStr_Energy[256], *ReadLine;
	FILE *fIn, *fEnergy, *fWater;

	for(i=0; i<nJob; i++)	{
		chdir(szDirList[i]);

		//start	to extract the optimized geometry of water
		strcpy(szName_Result, szOutputList[i]);
		fIn = fopen(szName_Result, "r");
		if(fIn == NULL)	{
			Quit_With_Error_Msg("Fail to open the output of the optimization of dimer.\nQuit\n");
		}
		if( ! Try_To_Find_Tag(fIn, szName_Result, " Optimization completed", szLine) ) {
			fseek( fIn, 0, SEEK_SET );
			To_Find_Tag( fIn,  szName_Result, " Optimization stopped", szLine); 
		}
		To_Find_Tag(fIn, szName_Result, " orientation:", szLine);
		Skip_N_Line(fIn, szName_Result, 4+n_Atom_In_Mol);	// including the coordinates of X1 and X2.
		
		while(1)	{
			if(feof(fIn))	{
				Quit_With_Error_Msg("Reaching the end of the file. Error in extracting the coordinates of water in Acceptor.\n");
			}
			ReadLine = fgets(szLine, 256, fIn);
			ReadItem = sscanf(szLine, "%d %d", &iTmp, &AtomNum);
			if( ReadItem !=2 )	{
				Quit_With_Error_Msg("Error in extracting the coordinates of water in Acceptor.\n");
			}
			if(AtomNum == -1)	{	// -1 means dummy atom !!!!
			}
			else	{
				ReadItem = sscanf(szLine+34, "%lf %lf %lf", &O_Wat_x, &O_Wat_y, &O_Wat_z);
				if(ReadItem != 3)	{
					Quit_With_Error_Msg("Fail to get the coordinate of OH2 from the output of Gaussian for restrained optimization of dimer.\nQuit\n");
				}
				else	{
					break;
				}
			}
		}
		
		while(1)	{
			if(feof(fIn))	{
				Quit_With_Error_Msg("Reaching the end of the file. Error in extracting the coordinates of water in Acceptor.\n");
			}
			ReadLine = fgets(szLine, 256, fIn);
			ReadItem = sscanf(szLine, "%d %d", &iTmp, &AtomNum);
			if( ReadItem !=2 )	{
				Quit_With_Error_Msg("Error in extracting the coordinates of water in Acceptor.\n");
			}
			if(AtomNum == -1)	{	// -1 means dummy atom !!!!
			}
			else	{
				ReadItem = sscanf(szLine+34, "%lf %lf %lf", &H1_Wat_x, &H1_Wat_y, &H1_Wat_z);
				if(ReadItem != 3)	{
					Quit_With_Error_Msg("Fail to get the coordinate of H1 from the output of Gaussian for restrained optimization of dimer.\nQuit\n");
				}
				else	{
					break;
				}
			}
		}
		
		while(1)	{
			if(feof(fIn))	{
				Quit_With_Error_Msg("Reaching the end of the file. Error in extracting the coordinates of water in Acceptor.\n");
			}
			ReadLine = fgets(szLine, 256, fIn);
			ReadItem = sscanf(szLine, "%d %d", &iTmp, &AtomNum);
			if( ReadItem !=2 )	{
				Quit_With_Error_Msg("Error in extracting the coordinates of water in Acceptor.\n");
			}
			if(AtomNum == -1)	{	// -1 means dummy atom !!!!
			}
			else	{
				ReadItem = sscanf(szLine+34, "%lf %lf %lf", &H2_Wat_x, &H2_Wat_y, &H2_Wat_z);
				if(ReadItem != 3)	{
					Quit_With_Error_Msg("Fail to get the coordinate of H2 from the output of Gaussian for restrained optimization of dimer.\nQuit\n");
				}
				else	{
					break;
				}
			}
		}
				
		//end	to extract the optimized geometry of water 
		
		x_Opt_QM[n_Atom_In_Mol+0] = O_Wat_x;	y_Opt_QM[n_Atom_In_Mol+0] = O_Wat_y;	z_Opt_QM[n_Atom_In_Mol+0] = O_Wat_z;
		x_Opt_QM[n_Atom_In_Mol+1] = H1_Wat_x;	y_Opt_QM[n_Atom_In_Mol+1] = H1_Wat_y;	z_Opt_QM[n_Atom_In_Mol+1] = H1_Wat_z;
		x_Opt_QM[n_Atom_In_Mol+2] = H2_Wat_x;	y_Opt_QM[n_Atom_In_Mol+2] = H2_Wat_y;	z_Opt_QM[n_Atom_In_Mol+2] = H2_Wat_z;
		
		//start	to extract the optimized geometry of water
		E_Min = 1.0E100;
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
				if(FindString(szLine, "SCF Done:  E(RHF)") >= 0)	{
					sscanf(szLine+20, "%s", szStr_Energy);
					Replace_D_E(szStr_Energy);
					ReadItem = sscanf(szStr_Energy, "%lf", &E_Min);
				}
			}
			
		}
		//start	to extract the optimized geometry of water
		
		fclose(fIn);
		
		E_Min_QM = (E_Min-E_Dimer_Far)*HARTREE_To_KCAL;
    printf("QM: %d : [%s] E_Min_QM is %f\n", i, szDirList[i],  E_Min_QM );


		memcpy(Dimer.x+n_Atom_In_Mol, x_Opt_QM+n_Atom_In_Mol, sizeof(double)*3);	// update the water configuration with lowest E_int
		memcpy(Dimer.y+n_Atom_In_Mol, y_Opt_QM+n_Atom_In_Mol, sizeof(double)*3);
		memcpy(Dimer.z+n_Atom_In_Mol, z_Opt_QM+n_Atom_In_Mol, sizeof(double)*3);
		
		memcpy(x_Wat_QM, x_Opt_QM+n_Atom_In_Mol, sizeof(double)*3);	// update the water configuration with lowest E_int
		memcpy(y_Wat_QM, y_Opt_QM+n_Atom_In_Mol, sizeof(double)*3);
		memcpy(z_Wat_QM, z_Opt_QM+n_Atom_In_Mol, sizeof(double)*3);
		
		r_Emin_QM = Cal_Dist_Two_Atoms(DonorH_List[i], n_Atom_In_Mol);	// H - OH2
	
		fWater = fopen("../mol-waters.pdb", "a+");
		fseek(fWater, 0, SEEK_END);
		fprintf(fWater, "ATOM%7d  %-3s %-3s A%4d    %8.3lf%8.3lf%8.3lf  1.00  0.00\n", 
			n_Atom_In_Mol+1, "OH2", "WAT", 2, Dimer.x[n_Atom_In_Mol], Dimer.y[n_Atom_In_Mol], Dimer.z[n_Atom_In_Mol]);
		fprintf(fWater, "ATOM%7d  %-3s %-3s A%4d    %8.3lf%8.3lf%8.3lf  1.00  0.00\n", 
			n_Atom_In_Mol+2, "H1", "WAT", 2, Dimer.x[n_Atom_In_Mol+1], Dimer.y[n_Atom_In_Mol+1], Dimer.z[n_Atom_In_Mol+1]);
		fprintf(fWater, "ATOM%7d  %-3s %-3s A%4d    %8.3lf%8.3lf%8.3lf  1.00  0.00\n", 
			n_Atom_In_Mol+3, "H2", "WAT", 2, Dimer.x[n_Atom_In_Mol+2], Dimer.y[n_Atom_In_Mol+2], Dimer.z[n_Atom_In_Mol+2]);
		fprintf(fWater, "END \n");
		fclose(fWater);
		
		fprintf(fFile_Run_Log, "mol_H_donor.cpp> %2d H %2d donor: %2d atom %3s QM, E_min = %7.4lf r_min = %7.4lf\n", 
			DonorIdx_List[i]+1, DonorH_List[i]+1, HeavyAtom_List[i]+1, szElemName[HeavyAtom_List[i]], E_Min_QM, r_Emin_QM);
		fflush(fFile_Run_Log);
		
		chdir(szCurDir);
		fEnergy = fopen(szMolWat_E_Int, "a+");
				
		fseek(fEnergy, SEEK_END, 0);
		fprintf(fEnergy, "#Energy    %3d donor %2d atom %3s QM ", DonorIdx_List[i]+1, DonorH_List[i]+1, szElemName[HeavyAtom_List[i]]);
		fprintf(fEnergy, "%12.6lf %9.4lf ", E_Min_QM, r_Emin_QM);
		fprintf(fEnergy, "\n");
		
		fprintf(fEnergy, "Donor %5d %12.6lf %9.4lf \n", DonorH_List[i]+1, E_Min_QM, r_Emin_QM);
		fprintf(fEnergy, "Saved_MM_Parameters %3d %22.16lf %22.16lf\n", 
			Theta_List[i], Para_MM_Save[i][0], Para_MM_Save[i][1]);
		fclose(fEnergy);
	}

				
//#ifndef _WIN32
//	lockf(fileno(fOut),F_LOCK,0L);
//#endif
				
//#ifndef _WIN32
//	lockf(fileno(fOut),F_ULOCK,0L);
//#endif
}

int Is_N_Planar(int Idx, double& x, double& y, double& z)
{
	int i_O, i_A, i_B, i_C;
	double vx_OA, vy_OA, vz_OA;
	double vx_OB, vy_OB, vz_OB;
	double vx_OC, vy_OC, vz_OC;
	double vx_Norm, vy_Norm, vz_Norm, Projection;

	i_O = Idx;
	i_A = BondList[i_O][0];
	i_B = BondList[i_O][1];
	i_C = BondList[i_O][2];

	Neighbor = i_A;

	vx_OA = x_Mol[i_A] - x_Mol[i_O];
	vy_OA = y_Mol[i_A] - y_Mol[i_O];
	vz_OA = z_Mol[i_A] - z_Mol[i_O];
	Normalizex(vx_OA, vy_OA, vz_OA);

	vx_OB = x_Mol[i_B] - x_Mol[i_O];
	vy_OB = y_Mol[i_B] - y_Mol[i_O];
	vz_OB = z_Mol[i_B] - z_Mol[i_O];
	Normalizex(vx_OB, vy_OB, vz_OB);

	vx_OC = x_Mol[i_C] - x_Mol[i_O];
	vy_OC = y_Mol[i_C] - y_Mol[i_O];
	vz_OC = z_Mol[i_C] - z_Mol[i_O];
	Normalizex(vx_OC, vy_OC, vz_OC);

	vx_Norm = vy_OB*vz_OC - vy_OC*vz_OB;
	vy_Norm = vz_OB*vx_OC - vz_OC*vx_OB;
	vz_Norm = vx_OB*vy_OC - vx_OC*vy_OB;
	Normalizex(vx_Norm, vy_Norm, vz_Norm);

	Projection = vx_OA*vx_Norm + vy_OA*vy_Norm + vz_OA*vz_Norm;

	if(fabs(Projection) < 0.25)	{	// planar. Only the direction along Norm.
		return 1;
	}
	else	{
		return 0;
	}
}


int Is_Three_Atoms_Colinear(int ia, int ib, int ic)
{
	double v1_x, v1_y, v1_z;
	double v2_x, v2_y, v2_z;
	double dot_v1_v2;

	v1_x = x_Mol[ib] - x_Mol[ia];
	v1_y = y_Mol[ib] - y_Mol[ia];
	v1_z = z_Mol[ib] - z_Mol[ia];

	v2_x = x_Mol[ic] - x_Mol[ib];
	v2_y = y_Mol[ic] - y_Mol[ib];
	v2_z = z_Mol[ic] - z_Mol[ib];

	Normalizex(v1_x, v1_y, v1_z);
	Normalizex(v2_x, v2_y, v2_z);

	dot_v1_v2 = fabs(Dot_Product(v1_x, v1_y, v1_z, v2_x, v2_y, v2_z));
	if(dot_v1_v2 > cos(15.0*radianInv))	{	// smaller than 10 degree
		return 1;
	}
	else	{
		return 0;
	}
}

void Generate_Dummy_Atoms(int Idx)
{
	int i, v_Aux_Sel, Is_Colinear;
	double v_Aux[2][3]={{1, 1, 0}, {0, 1, 1}}, dot_1, dot_2;
	double v_x_H_X2, v_y_H_X2, v_z_H_X2;

	x_Dummy_1 = x_Mol[ActiveAtom];	// dummy atom 1
	y_Dummy_1 = y_Mol[ActiveAtom];	// dummy atom 1
	z_Dummy_1 = z_Mol[ActiveAtom];	// dummy atom 1

	v_x_X1_H = x_Mol[ActiveHAtom] - x_Dummy_1;
	v_y_X1_H = y_Mol[ActiveHAtom] - y_Dummy_1;
	v_z_X1_H = z_Mol[ActiveHAtom] - z_Dummy_1;
	Normalizex(v_x_X1_H, v_y_X1_H, v_z_X1_H);
	
	if(Mol.nAtom == 2)	{
		dot_1 = fabs(Dot_Product(v_x_X1_H, v_y_X1_H, v_z_X1_H, v_Aux[0][0], v_Aux[0][1], v_Aux[0][2]));
		dot_2 = fabs(Dot_Product(v_x_X1_H, v_y_X1_H, v_z_X1_H, v_Aux[1][0], v_Aux[1][1], v_Aux[1][2]));
		
		if(dot_1 <= dot_2)	{
			v_Aux_Sel = 0;
		}
		else	{
			v_Aux_Sel = 1;
		}
		
		Cross_Product(v_x_X1_H, v_y_X1_H, v_z_X1_H, v_Aux[v_Aux_Sel][0], v_Aux[v_Aux_Sel][1], v_Aux[v_Aux_Sel][2], v_x_H_X2, v_y_H_X2, v_z_H_X2);
	}
	else	{
		if(Bond_Count[ActiveAtom] == 0)	{
			GetBondList(ActiveAtom);
		}
		for(i=0; i<Bond_Count[ActiveAtom]; i++)	{
			if(BondList[ActiveAtom][i] != ActiveHAtom)	{	// 2nd neighbor atom
				break;
			}
		}
		if(i >= Bond_Count[ActiveAtom])	{
			Quit_With_Error_Msg("Donor> Error in Generate_Dummy_Atoms_Bond_1().\nQuit\n");
		}
		i = BondList[ActiveAtom][i];

		//start	to check whether three atoms (ActiveHAtom, ActiveAtom, i) are colinear or not
		Is_Colinear = Is_Three_Atoms_Colinear(ActiveHAtom, ActiveAtom, i);
		if(Is_Colinear)	{
			dot_1 = fabs(Dot_Product(v_x_X1_H, v_y_X1_H, v_z_X1_H, v_Aux[0][0], v_Aux[0][1], v_Aux[0][2]));
			dot_2 = fabs(Dot_Product(v_x_X1_H, v_y_X1_H, v_z_X1_H, v_Aux[1][0], v_Aux[1][1], v_Aux[1][2]));
			
			if(dot_1 <= dot_2)	{
				v_Aux_Sel = 0;
			}
			else	{
				v_Aux_Sel = 1;
			}
			
			Cross_Product(v_x_X1_H, v_y_X1_H, v_z_X1_H, v_Aux[v_Aux_Sel][0], v_Aux[v_Aux_Sel][1], v_Aux[v_Aux_Sel][2], v_x_H_X2, v_y_H_X2, v_z_H_X2);
		}
		else	{
			Cross_Product(v_x_X1_H, v_y_X1_H, v_z_X1_H, x_Mol[ActiveAtom]-x_Mol[i], y_Mol[ActiveAtom]-y_Mol[i], z_Mol[ActiveAtom]-z_Mol[i], v_x_H_X2, v_y_H_X2, v_z_H_X2);
		}
		//end	to check whether three atoms (Idx, i_A, i) are colinear or not
	}

	Normalizex(v_x_H_X2, v_y_H_X2, v_z_H_X2);
	
	x_Dummy_2 = x_Mol[ActiveHAtom] + R_DUMMY_2 * v_x_H_X2;
	y_Dummy_2 = y_Mol[ActiveHAtom] + R_DUMMY_2 * v_y_H_X2;
	z_Dummy_2 = z_Mol[ActiveHAtom] + R_DUMMY_2 * v_z_H_X2;

	x_Dummy_3 = x_Mol[ActiveHAtom] + R_DUMMY_3 * v_x_X1_H;
	y_Dummy_3 = y_Mol[ActiveHAtom] + R_DUMMY_3 * v_y_X1_H;
	z_Dummy_3 = z_Mol[ActiveHAtom] + R_DUMMY_3 * v_z_X1_H;
	
}

void Rotate_Dummy_Atoms(int Idx, double theta)
{
	double v1_New_x, v1_New_y, v1_New_z;
	double v2_x, v2_y, v2_z, v2_New_x, v2_New_y, v2_New_z;
	double v3_x, v3_y, v3_z;
//	double v3_New_x, v3_New_y, v3_New_z;
	double RotM[3][3], cos_theta, sin_theta;

	theta *= radianInv;
	cos_theta = cos(theta);
	sin_theta = sin(theta);

	v_x_X1_H = x_Mol[ActiveHAtom] - x_Dummy_1_Org;	// X vector
	v_y_X1_H = y_Mol[ActiveHAtom] - y_Dummy_1_Org;
	v_z_X1_H = z_Mol[ActiveHAtom] - z_Dummy_1_Org;
	Normalizex(v_x_X1_H, v_y_X1_H, v_z_X1_H);

	v2_x = x_Mol[ActiveHAtom] - x_Dummy_2_Org;
	v2_y = y_Mol[ActiveHAtom] - y_Dummy_2_Org;
	v2_z = z_Mol[ActiveHAtom] - z_Dummy_2_Org;
	Normalizex(v2_x, v2_y, v2_z);

	Cross_Product(v_x_X1_H, v_y_X1_H, v_z_X1_H, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z);
	Normalizex(v3_x, v3_y, v3_z);

	RotM[0][0] = v3_x*v3_x + (1-v3_x*v3_x)*cos_theta;
	RotM[0][1] = v3_x*v3_y*(1.0-cos_theta) - v3_z*sin_theta;
	RotM[0][2] = v3_x*v3_z*(1.0-cos_theta) + v3_y*sin_theta;
	
	RotM[1][0] = v3_x*v3_y*(1.0-cos_theta) + v3_z*sin_theta;
	RotM[1][1] = v3_y*v3_y + (1-v3_y*v3_y)*cos_theta;
	RotM[1][2] = v3_y*v3_z*(1.0-cos_theta) - v3_x*sin_theta;
	
	RotM[2][0] = v3_x*v3_z*(1.0-cos_theta) - v3_y*sin_theta;
	RotM[2][1] = v3_y*v3_z*(1.0-cos_theta) + v3_x*sin_theta;
	RotM[2][2] = v3_z*v3_z + (1-v3_z*v3_z)*cos_theta;

	v1_New_x = RotM[0][0] * v_x_X1_H + RotM[0][1] * v_y_X1_H + RotM[0][2] * v_z_X1_H;
	v1_New_y = RotM[1][0] * v_x_X1_H + RotM[1][1] * v_y_X1_H + RotM[1][2] * v_z_X1_H;
	v1_New_z = RotM[2][0] * v_x_X1_H + RotM[2][1] * v_y_X1_H + RotM[2][2] * v_z_X1_H;

	v2_New_x = RotM[0][0] * v2_x + RotM[0][1] * v2_y + RotM[0][2] * v2_z;
	v2_New_y = RotM[1][0] * v2_x + RotM[1][1] * v2_y + RotM[1][2] * v2_z;
	v2_New_z = RotM[2][0] * v2_x + RotM[2][1] * v2_y + RotM[2][2] * v2_z;

//	v3_New_x = RotM[0][0] * v3_x + RotM[0][1] * v3_y + RotM[0][2] * v3_z;	// not changed
//	v3_New_y = RotM[1][0] * v3_x + RotM[1][1] * v3_y + RotM[1][2] * v3_z;	// not changed
//	v3_New_z = RotM[2][0] * v3_x + RotM[2][1] * v3_y + RotM[2][2] * v3_z;	// not changed

	x_Dummy_1 = x_Mol[ActiveHAtom] - R_DUMMY_1*v1_New_x;
	y_Dummy_1 = y_Mol[ActiveHAtom] - R_DUMMY_1*v1_New_y;
	z_Dummy_1 = z_Mol[ActiveHAtom] - R_DUMMY_1*v1_New_z;

	x_Dummy_2 = x_Mol[ActiveHAtom] - R_DUMMY_2*v2_New_x;
	y_Dummy_2 = y_Mol[ActiveHAtom] - R_DUMMY_2*v2_New_y;
	z_Dummy_2 = z_Mol[ActiveHAtom] - R_DUMMY_2*v2_New_z;

	v_x_X1_H = v1_New_x;
	v_y_X1_H = v1_New_y;
	v_z_X1_H = v1_New_z;

	x_Dummy_3 = x_Mol[ActiveHAtom] + R_DUMMY_3 * v_x_X1_H;
	y_Dummy_3 = y_Mol[ActiveHAtom] + R_DUMMY_3 * v_y_X1_H;
	z_Dummy_3 = z_Mol[ActiveHAtom] + R_DUMMY_3 * v_z_X1_H;
}

void Make_List_H_Acceptor(void)
{
	FILE *fIn;
	int i;
	char ErrorMsg[256];

	fIn = fopen(szFileElem, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file %s\nQuit\n", szFileElem);
		Quit_With_Error_Msg(ErrorMsg);
	}

	memset(Is_N, 0, sizeof(int)*n_Atom_In_Mol);
	memset(Is_O, 0, sizeof(int)*n_Atom_In_Mol);
	memset(Is_F, 0, sizeof(int)*n_Atom_In_Mol);
	memset(Is_S, 0, sizeof(int)*n_Atom_In_Mol);

	for(i=0; i<n_Atom_In_Mol; i++)	{
		fscanf(fIn, "%s", szElemName[i]);
	}
	fclose(fIn);

	for(i=0; i<n_Atom_In_Mol; i++)	{
		if(strcmp(szElemName[i], "N")==0)	{
			Is_N[i] = 1;
			GetBondList(i);
		}
		else if(strcmp(szElemName[i], "O")==0)	{
			Is_O[i] = 1;
			GetBondList(i);
		}
		else if(strcmp(szElemName[i], "F")==0)	{
			Is_F[i] = 1;
			GetBondList(i);
		}
		else if(strcmp(szElemName[i], "S")==0)	{
			Is_S[i] = 1;
			GetBondList(i);
		}
	}
}

void GetBondList(int Idx)
{
	int i, nBond, iPos, *pBondList, Atom_2;

	nBond = Mol.nBond;
	pBondList = Mol.BondList;

	Bond_Count[Idx] = 0;
	nCount_H[Idx] = 0;

	for(i=0; i<nBond; i++)	{
		iPos = 2*i;
		if(pBondList[iPos] == Idx)	{
			Atom_2 = pBondList[iPos+1];
			BondList[Idx][Bond_Count[Idx]] = Atom_2;
			Bond_Count[Idx]++;

			if( strcmp(szElemName[Atom_2], "H")==0 )	{
				BondWithH[Idx][nCount_H[Idx]] = Atom_2;
				nCount_H[Idx]++;
			}
		}
		if(pBondList[iPos+1] == Idx)	{
			Atom_2 = pBondList[iPos];
			BondList[Idx][Bond_Count[Idx]] = Atom_2;
			Bond_Count[Idx]++;
			if( strcmp(szElemName[Atom_2], "H")==0 )	{
				BondWithH[Idx][nCount_H[Idx]] = Atom_2;
				nCount_H[Idx]++;
			}
		}
	}
	return;
}


//void Optimize_Water_Pose_QM(void)
void Output_Optimize_Water_Pose_QM(void)
{
	FILE *fOut;
	int i;
	char szName_Result[256]="qm.out", szDonor[256], szDonorH[256];

	//start	to generate the Gaussian script for optimization with ModRedundant
	fOut = fopen(szName_Gaussian, "w");
	fprintf(fOut, "%s", szQM_Level_Int_E_Dimer_Opt_Cmd);
#if 0
	sprintf(szDonor, "%s99", szElemName[ActiveAtom]);
	sprintf(szDonorH, "%s99", szElemName[ActiveHAtom]);
	for(i=0; i<n_Atom_In_Mol; i++)	{
		if(i==ActiveAtom)	{
			fprintf(fOut, "%s  %10.5lf%10.5lf%10.5lf\n", szDonor, x_ReadIn[i], y_ReadIn[i], z_ReadIn[i]);
		}
		else if(i==ActiveHAtom)	{
			fprintf(fOut, "%s  %10.5lf%10.5lf%10.5lf\n", szDonorH, x_ReadIn[i], y_ReadIn[i], z_ReadIn[i]);
		}
		else	{
			fprintf(fOut, "%s    %10.5lf%10.5lf%10.5lf\n", szElemName[i], x_ReadIn[i], y_ReadIn[i], z_ReadIn[i]);
		}
	}
	fprintf(fOut, "%s   %10.5lf%10.5lf%10.5lf\n", "X1", x_Dummy_1, y_Dummy_1, z_Dummy_1);
	fprintf(fOut, "%s   %10.5lf%10.5lf%10.5lf\n", "X2", x_Dummy_2, y_Dummy_2, z_Dummy_2);
	fprintf(fOut, "%s   %10.5lf%10.5lf%10.5lf\n", "X3", x_Dummy_3, y_Dummy_3, z_Dummy_3);

	fprintf(fOut, "O1W  %s  roh X2 90.0 %s 180.0 \n", szDonorH, szDonor);
	fprintf(fOut, "H1W  O1W %10.4lf %s  %10.4lf X2 phi \n", b0_H_O_HF, szDonorH, 180.0-0.5*theta0_H_O_H_HF);
	fprintf(fOut, "H2W  O1W %10.4lf H1W %10.4lf X3 0.0 \n", b0_H_O_HF, theta0_H_O_H_HF);

	fprintf(fOut, "\nroh %10.4lf\nphi %lf\n", Para_Save_for_QM[0]+0.15, Para_Save_for_QM[1]*radian);	// shifted
#else
	for(i=0; i<n_Atom_In_Mol; i++)	{
			fprintf(fOut, "%s    %10.5lf%10.5lf%10.5lf\n", szElemName[i], x_ReadIn[i], y_ReadIn[i], z_ReadIn[i]);
	}
	
	int x1 = n_Atom_In_Mol + 1;
	int x2 = n_Atom_In_Mol + 2;
	int x3 = n_Atom_In_Mol + 3;
	int o1w= n_Atom_In_Mol + 4;
	int h1w= n_Atom_In_Mol + 5;
	int h2w= n_Atom_In_Mol + 6;

	fprintf(fOut, "%s   %10.5lf%10.5lf%10.5lf\n", "X", x_Dummy_1, y_Dummy_1, z_Dummy_1);
	fprintf(fOut, "%s   %10.5lf%10.5lf%10.5lf\n", "X", x_Dummy_2, y_Dummy_2, z_Dummy_2);
	fprintf(fOut, "%s   %10.5lf%10.5lf%10.5lf\n", "X", x_Dummy_3, y_Dummy_3, z_Dummy_3);

	fprintf(fOut, "O  %d  roh %d 90.0 %d 180.0 \n", ActiveHAtom+1, x2, ActiveAtom+1);
	fprintf(fOut, "H  %d %10.4lf %d  %10.4lf %d phi \n", o1w, b0_H_O_HF, ActiveHAtom+1, 180.0-0.5*theta0_H_O_H_HF, x2);
	fprintf(fOut, "H  %d %10.4lf %d %10.4lf %d 0.0 \n", o1w, b0_H_O_HF, h1w, theta0_H_O_H_HF, x3);

	fprintf(fOut, "\nroh = %10.4lf\nphi = %lf\n", Para_Save_for_QM[0]+0.15, Para_Save_for_QM[1]*radian);	// shifted

#endif
	
	fclose(fOut);
	//end	to generate the Gaussian script for optimization with ModRedundant

	szInputList[nJob] = strdup( szName_Gaussian );
	szOutputList[nJob]= strdup( szName_Result );
	
//	sprintf(szCmd, "%s < %s > %s", szExe_G09, szName_Gaussian, szName_Result);
//	system(szCmd);

/*	
	//start	to extract the optimized geometry of water
	fIn = fopen(szName_Result, "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open the output of the optimization of dimer.\nQuit\n");
	}
	To_Find_Tag(fIn, szName_Result, " Optimization completed", szLine);
	To_Find_Tag(fIn, szName_Result, " orientation:", szLine);
	Skip_N_Line(fIn, szName_Result, 4+n_Atom_In_Mol);	// including the coordinates of X1 and X2.

	while(1)	{
		if(feof(fIn))	{
			Quit_With_Error_Msg("Reaching the end of the file. Error in extracting the coordinates of water in Acceptor.\n");
		}
		ReadLine = fgets(szLine, 256, fIn);
		ReadItem = sscanf(szLine, "%d %d", &iTmp, &AtomNum);
		if( ReadItem !=2 )	{
			Quit_With_Error_Msg("Error in extracting the coordinates of water in Acceptor.\n");
		}
		if(AtomNum == -1)	{	// -1 means dummy atom !!!!
		}
		else	{
			ReadItem = sscanf(szLine+34, "%lf %lf %lf", &O_Wat_x, &O_Wat_y, &O_Wat_z);
			if(ReadItem != 3)	{
				Quit_With_Error_Msg("Fail to get the coordinate of OH2 from the output of Gaussian for restrained optimization of dimer.\nQuit\n");
			}
			else	{
				break;
			}
		}
	}

	while(1)	{
		if(feof(fIn))	{
			Quit_With_Error_Msg("Reaching the end of the file. Error in extracting the coordinates of water in Acceptor.\n");
		}
		ReadLine = fgets(szLine, 256, fIn);
		ReadItem = sscanf(szLine, "%d %d", &iTmp, &AtomNum);
		if( ReadItem !=2 )	{
			Quit_With_Error_Msg("Error in extracting the coordinates of water in Acceptor.\n");
		}
		if(AtomNum == -1)	{	// -1 means dummy atom !!!!
		}
		else	{
			ReadItem = sscanf(szLine+34, "%lf %lf %lf", &H1_Wat_x, &H1_Wat_y, &H1_Wat_z);
			if(ReadItem != 3)	{
				Quit_With_Error_Msg("Fail to get the coordinate of H1 from the output of Gaussian for restrained optimization of dimer.\nQuit\n");
			}
			else	{
				break;
			}
		}
	}

	while(1)	{
		if(feof(fIn))	{
			Quit_With_Error_Msg("Reaching the end of the file. Error in extracting the coordinates of water in Acceptor.\n");
		}
		ReadLine = fgets(szLine, 256, fIn);
		ReadItem = sscanf(szLine, "%d %d", &iTmp, &AtomNum);
		if( ReadItem !=2 )	{
			Quit_With_Error_Msg("Error in extracting the coordinates of water in Acceptor.\n");
		}
		if(AtomNum == -1)	{	// -1 means dummy atom !!!!
		}
		else	{
			ReadItem = sscanf(szLine+34, "%lf %lf %lf", &H2_Wat_x, &H2_Wat_y, &H2_Wat_z);
			if(ReadItem != 3)	{
				Quit_With_Error_Msg("Fail to get the coordinate of H2 from the output of Gaussian for restrained optimization of dimer.\nQuit\n");
			}
			else	{
				break;
			}
		}
	}


	
	//end	to extract the optimized geometry of water 

	x_Opt_QM[n_Atom_In_Mol+0] = O_Wat_x;	y_Opt_QM[n_Atom_In_Mol+0] = O_Wat_y;	z_Opt_QM[n_Atom_In_Mol+0] = O_Wat_z;
	x_Opt_QM[n_Atom_In_Mol+1] = H1_Wat_x;	y_Opt_QM[n_Atom_In_Mol+1] = H1_Wat_y;	z_Opt_QM[n_Atom_In_Mol+1] = H1_Wat_z;
	x_Opt_QM[n_Atom_In_Mol+2] = H2_Wat_x;	y_Opt_QM[n_Atom_In_Mol+2] = H2_Wat_y;	z_Opt_QM[n_Atom_In_Mol+2] = H2_Wat_z;

	//start	to extract the optimized geometry of water
	E_Min = 1.0E100;
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
			if(FindString(szLine, "SCF Done:  E(RHF)") >= 0)	{
				sscanf(szLine+20, "%s", szStr_Energy);
				Replace_D_E(szStr_Energy);
				ReadItem = sscanf(szStr_Energy, "%lf", &E_Min);
			}
		}

	}
	//start	to extract the optimized geometry of water

	fclose(fIn);

	return E_Min;
*/
}


double Optimize_Water_Pose(void)
{
	nlopt_opt opt;

	Low_Bound[0] = 1.6;
	Up_Bound[0] = 2.8;
	Low_Bound[1] = -10.0;
	Up_Bound[1] = 10.0;

	opt = nlopt_create(NLOPT_LD_LBFGS, N_PARA); /* algorithm and dimensionality */
//	opt = nlopt_create(NLOPT_LD_LBFGS_NOCEDAL, N_PARA); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, Low_Bound);
	nlopt_set_upper_bounds(opt, Up_Bound);
	nlopt_set_min_objective(opt, Callback_Eval_Gradient, NULL);
	nlopt_set_xtol_rel(opt, 1.0E-5);

	E_Min=1.0E20;
	CalObjectiveFunction();
//	printf("Iteration %4d  E = %10.5lf\n", Iteration, E_Dimer_Int);

	FailCount = 0;
	Iteration = 0;

	memcpy(Para_Best, Para_List, sizeof(double)*N_PARA);
	if (nlopt_optimize(opt, Para_List, &E_Dimer_Min) < 0) {
		printf("nlopt failed!\n");
	}
	else {
		printf("Obj_Min = %16.10lf\n", E_Dimer_Min);
	}

	nlopt_destroy(opt);

	if(Para_Best[1] < -PI2)	{
		Para_Best[1] += PI2;
	}
	else if(Para_Best[1] > PI2)	{
		Para_Best[1] -= PI2;
	}

	memcpy(Para_List, Para_Best, sizeof(double)*N_PARA);
	E_Dimer_Min = CalObjectiveFunction();


	return E_Dimer_Min;
}

double Callback_Eval_Gradient(unsigned n, const double *x, double *grad, void *my_func_data)
{
	int i, n_Local;
//	FILE *fOut;

	n_Local = (int)n;
	for(i=0; i<n_Local; i++)	{
		Para_List[i] = x[i];
	}

	if(grad)	{
		CalGradient();	//gradient will be assigned into objGrad automatically
		Iteration++;
		for(i=0; i<n_Local; i++)	{
			grad[i] = Grad_List[i];
		}
//		printf("Iteration %4d  E = %16.8lf\n", Iteration, E_Dimer_Int);

		if(E_Dimer_Int < E_Min)	{
			E_Min = E_Dimer_Int;
			FailCount = 0;

			memcpy(Para_Best, x, sizeof(double)*n);
//			sprintf(szCmd, "cp %s opt.gjf", szName_Gaussian);
//			system(szCmd);
//			fOut = fopen("E_min.txt", "w");
//			fprintf(fOut, "%20.12lf\n%20.12lf\n%20.12lf\n%20.12lf\n", Para_List[0], Para_List[1], Para_List[2], Para_List[3]);
//			fprintf(fOut, "E = %15.8lf\n", E_Min);
//			fclose(fOut);
		}
		else	{
			if(FailCount > 4)      {       // cannot converge due to limited accuracy of numerical gradient
				printf("Too much failed tries.\nI have tried my best.\nQuit.\n");
//				printf("E_Min = %12.6lf\n", E_Min);
//				exit(0);
				memset(grad, 0, sizeof(double)*n);
			}
			FailCount++;
		}
	}
	else	{	// cal object function only
		CalObjectiveFunction();
	}

    return E_Dimer_Int;
}

double CalObjectiveFunction(void)	// only used for the optimization in MM
{
	GenerateWater();
//	Generate_Gaussian_Input();
	E_Dimer_Int = Dimer.Cal_E_Int();

	return E_Dimer_Int;
}

double Get_E_Mol_Water_Far_QM()
{
	FILE *fOut, *fIn;
	char szFileEnergy[]="energy_sum.txt";
	int i, ReadItem;
	double E_Mol_Water=1e100;
	fIn = fopen(szFileEnergy, "r");
	if(fIn != NULL)	{
		ReadItem = fscanf(fIn, "%lf", &E_Mol_Water);
		fclose(fIn);
		if(ReadItem == 1)	{
      printf("Saved dimer energy %lf\n", E_Mol_Water );
			return E_Mol_Water;
		}
	}
	printf("energy_sum.txt doesn't exist: don_acc_separate needs to be run\n" );
	exit(1);
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


double Extract_E_Gaussian_Output(char szFileName[], char szTag[])
{
	FILE *fIn;
	char *ReadLine, szLine[256], szStr_Energy[256], ErrorMsg[256];
	int ReadItem;
	double E_Dimer;

	fIn = fopen(szFileName, "r");
	if(FindString(szTag, "SCF Done:")>=0)	{
		while(1)	{
			if(feof(fIn))	{
				break;
			}
			ReadLine = fgets(szLine, 256, fIn);
			if(ReadLine)	{
				if(FindString(szLine, " SCF Done:  E(RHF)") >= 0)	{
					ReadItem = sscanf(szLine+20, "%s", szStr_Energy);
					Replace_D_E(szStr_Energy);
					ReadItem = sscanf(szStr_Energy, "%lf", &E_Dimer);
					if(ReadItem == 1)	{
						fclose(fIn);
						return E_Dimer;
					}
				}
			}
		}
	}
	else if(FindString(szTag, "Counterpoise: corrected energy =")>=0)	{
		while(1)	{
			if(feof(fIn))	{
				break;
			}
			ReadLine = fgets(szLine, 256, fIn);
			if(ReadLine)	{
				if(FindString(szLine, "Counterpoise: corrected energy =") >= 0)	{
					ReadItem = sscanf(szLine+33, "%s", szStr_Energy);
					Replace_D_E(szStr_Energy);
					ReadItem = sscanf(szStr_Energy, "%lf", &E_Dimer);
					if(ReadItem == 1)	{
						fclose(fIn);
						return E_Dimer;
					}
				}
			}
		}
	}
	else if(FindString(szTag, "MM energy: ")>=0)	{
		while(1)	{
			if(feof(fIn))	{
				break;
			}
			ReadLine = fgets(szLine, 256, fIn);
			if(ReadLine)	{
				if(FindString(szLine, "MM energy: ") >= 0)	{
					ReadItem = sscanf(szLine+11, "%s", szStr_Energy);
					Replace_D_E(szStr_Energy);
					ReadItem = sscanf(szStr_Energy, "%lf", &E_Dimer);
					if(ReadItem == 1)	{
						fclose(fIn);
						return E_Dimer;
					}
				}
			}
		}
	}
	else	{
		sprintf(ErrorMsg, "Extract_E_Gaussian_Output> Unsupported tags for extracting energies from Gaussian output. %s\nQuit\n", szTag);
		fclose(fIn);
		Quit_With_Error_Msg(ErrorMsg);
	}
	return 0.0;
}

#define DELTA_PARA	(1.0E-4)
//#define DELTA_PARA	(5.0E-3)
double CalGradient(void)
{
	int i;
	double Para_Save[N_PARA];
	double f_Left, f_Right;
	
	for(i=0; i<N_PARA; i++)	{
		Para_Save[i] = Para_List[i];	//save current parameters
	}

	for(i=0; i<N_PARA; i++)	{
		Para_List[i] = Para_Save[i] - DELTA_PARA;
		f_Left = CalObjectiveFunction();
		
		Para_List[i] = Para_Save[i] + DELTA_PARA;
		f_Right = CalObjectiveFunction();
		
		Grad_List[i] = (f_Right-f_Left)/(2.0*DELTA_PARA);
//		printf("%10.5lf\n", Grad_List[i]);
		
		Para_List[i] = Para_Save[i];	//restore save (correct) parameter
	}
	
	CalObjectiveFunction();

	return E_Dimer_Int;
}

#undef DELTA_PARA

void ReadMolCRD(void)
{
	FILE *fIn;
	int ReadItem, i;
	char szLine[256], *ReadLine, ErrorMsg[256];


	fIn = fopen(szNameMolCRD, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file %s.\nQuit\n", szNameMolCRD);
		Quit_With_Error_Msg(ErrorMsg);
	}

	while(1)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Error to find the line which contains the number of atoms.\nQuit.\n");
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}
		fgets(szLine, 256, fIn);
		ReadItem = sscanf(szLine, "%d", &n_Atom_In_Mol);
		if(ReadItem == 1)	{
			break;
		}
	}

	for(i=0; i<n_Atom_In_Mol; i++)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Error in reading the coodinates of Mol.\nFile is not complete. \nQuit.\n");
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			sprintf(ErrorMsg, "Error in reading the coodinates of Mol.\nFile is not complete. \nQuit.\n");
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}

		ReadItem = sscanf(szLine+6, "%d", &(ResIDMol[i]));
		ReadItem = sscanf(szLine+11, "%s %s", szResNameMol[i], szAtomNameMol[i]);
		ReadItem = sscanf(szLine+20, "%lf%lf%lf", &(x_Mol[i]), &(y_Mol[i]), &(z_Mol[i]));
		if(ReadItem != 3)	{
			sprintf(ErrorMsg, "Error in reading the coodinates of Mol.\nFile is not complete. \nQuit.\n");
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}
	}

	fclose(fIn);
}

void Gen_xyz(double x_A, double y_A, double z_A, 
			 double x_B, double y_B, double z_B, 
			 double x_C, double y_C, double z_C, 
			 double r, double theta, double phi, 
			 double& x_D, double& y_D, double& z_D)
{
	double vz_x, vz_y, vz_z, BA_x, BA_y, BA_z;
	double vy_x, vy_y, vy_z, vx_x, vx_y, vx_z;
	double cos_phi, sin_phi, cos_theta, sin_theta;
	double r_x, r_y, r_z;

	cos_phi = cos(phi);
	sin_phi = sin(phi);

	vz_x = x_C - x_B;
	vz_y = y_C - y_B;
	vz_z = z_C - z_B;
	Normalizex(vz_x, vz_y, vz_z);	

	BA_x = x_A - x_B;
	BA_y = y_A - y_B;
	BA_z = z_A - z_B;
	Normalizex(BA_x, BA_y, BA_z);	

	vy_x = vz_y*BA_z - vz_z*BA_y;
	vy_y = vz_z*BA_x - vz_x*BA_z;
	vy_z = vz_x*BA_y - vz_y*BA_x;
	Normalizex(vy_x, vy_y, vy_z);

	vx_x = vy_y*vz_z - vy_z*vz_y;
	vx_y = vy_z*vz_x - vy_x*vz_z;
	vx_z = vy_x*vz_y - vy_y*vz_x;

	cos_theta = cos(PI-theta);
	sin_theta = sin(PI-theta);
	
	r_z = r*cos_theta;
	r_x = r*sin_theta*cos_phi;
	r_y = r*sin_theta*sin_phi;
	
	x_D = x_C + vx_x*r_x + vy_x*r_y + vz_x*r_z;
	y_D = y_C + vx_y*r_x + vy_y*r_y + vz_y*r_z;
	z_D = z_C + vx_z*r_x + vy_z*r_y + vz_z*r_z;

	return;
}


void GenerateWater(void)	// with six degree of freedom for water molecule
{
	double r;

	r = Para_List[0];
	O_Wat_x = x_Mol[ActiveHAtom] + r*v_x_X1_H;
	O_Wat_y = y_Mol[ActiveHAtom] + r*v_y_X1_H;
	O_Wat_z = z_Mol[ActiveHAtom] + r*v_z_X1_H;

	Gen_xyz(x_Dummy_2, y_Dummy_2, z_Dummy_2, x_Mol[ActiveHAtom], y_Mol[ActiveHAtom], z_Mol[ActiveHAtom], 
		O_Wat_x, O_Wat_y, O_Wat_z, b0_H_O_EXP, PI-0.5*theta0_H_O_H_EXP*radianInv, Para_List[1], H1_Wat_x, H1_Wat_y, H1_Wat_z);

	Gen_xyz(x_Dummy_3, y_Dummy_3, z_Dummy_3, H1_Wat_x, H1_Wat_y, H1_Wat_z, 
		O_Wat_x, O_Wat_y, O_Wat_z, b0_H_O_EXP, theta0_H_O_H_EXP*radianInv, 0.0, H2_Wat_x, H2_Wat_y, H2_Wat_z);

	Dimer.x[Idx_O_Wat] = O_Wat_x;	Dimer.y[Idx_O_Wat] = O_Wat_y;	Dimer.z[Idx_O_Wat] = O_Wat_z;
	Dimer.x[Idx_H1_Wat] = H1_Wat_x;	Dimer.y[Idx_H1_Wat] = H1_Wat_y;	Dimer.z[Idx_H1_Wat] = H1_Wat_z;
	Dimer.x[Idx_H2_Wat] = H2_Wat_x;	Dimer.y[Idx_H2_Wat] = H2_Wat_y;	Dimer.z[Idx_H2_Wat] = H2_Wat_z;
}

inline double Dot_Product(double x_a, double y_a, double z_a, double x_b, double y_b, double z_b)
{
	return (x_a*x_b + y_a*y_b + z_a*z_b);
}

inline void Cross_Product(double x_a, double y_a, double z_a, double x_b, double y_b, double z_b, double& x_Product, double& y_Product, double& z_Product)
{
	x_Product = y_a*z_b - y_b*z_a;
	y_Product = z_a*x_b - z_b*x_a;
	z_Product = x_a*y_b - x_b*y_a;
}

double Normalizex(double& vect_x, double& vect_y, double& vect_z)
{
	double dist, dist_Inv;

	dist = sqrt(vect_x*vect_x + vect_y*vect_y + vect_z*vect_z);
	if(dist < 1.0E-100)	{
		printf("");
	}
	dist_Inv = 1.0/dist;
	vect_x *= dist_Inv;
	vect_y *= dist_Inv;
	vect_z *= dist_Inv;
	return dist;
}

void OutputCRD(void)
{
	FILE *fOut;
	int i, ResBase, Water_Count=0, nAtom;

	fOut = fopen("solvated-mol.crd", "w");
	for(i=0; i<n_Atom_In_Mol; i++)	{
		fprintf(fOut, "%5d%5d %-5s%-4s%10.5lf%10.5lf%10.5lf %-4s %-5d  0.00000\n", 
			i+1, ResIDMol[i], szResNameMol[i], szAtomNameMol[i], x_Mol[i], y_Mol[i], z_Mol[i], szResNameMol[i], ResIDMol[i]);
	}
	ResBase = ResIDMol[n_Atom_In_Mol-1]+1;
	nAtom = n_Atom_In_Mol;



	for(Water_Count=0; Water_Count<6; Water_Count++)	{
//		GenerateWater(H1_Wat_x, H1_Wat_y, H1_Wat_z, 60.0*Water_Count);
		GenerateWater();

		fprintf(fOut, "%5d%5d SWM4 OH2 %10.5lf%10.5lf%10.5lf WAT  %-5d  0.00000\n", 
			nAtom+1, ResBase, O_Wat_x, O_Wat_y, O_Wat_z, Water_Count+1);
		fprintf(fOut, "%5d%5d SWM4 H1  %10.5lf%10.5lf%10.5lf WAT  %-5d  0.00000\n", 
			nAtom+2, ResBase, H1_Wat_x, H1_Wat_y, H1_Wat_z, Water_Count+1);
		fprintf(fOut, "%5d%5d SWM4 H2  %10.5lf%10.5lf%10.5lf WAT  %-5d  0.00000\n", 
			nAtom+3, ResBase, H2_Wat_x, H2_Wat_y, H2_Wat_z, Water_Count+1);
		nAtom += 3;
		ResBase++;
	}

	fclose(fOut);
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


void Generate_Gaussian_Input(void)
{
	FILE *fOut;
	int i;

	fOut = fopen(szName_Gaussian, "w");
	fprintf(fOut, "%s", szQM_Level_Int_E_Dimer_Cmd);
	for(i=0; i<n_Atom_In_Mol; i++)	{
		fprintf(fOut, "%s  %10.5lf%10.5lf%10.5lf  1\n", szElemName[i], x_Mol[i], y_Mol[i], z_Mol[i]);
	}

	fprintf(fOut, "O  %10.5lf%10.5lf%10.5lf  2\n", O_Wat_x, O_Wat_y, O_Wat_z);
	fprintf(fOut, "H  %10.5lf%10.5lf%10.5lf  2\n", H1_Wat_x, H1_Wat_y, H1_Wat_z);
	fprintf(fOut, "H  %10.5lf%10.5lf%10.5lf  2\n\n", H2_Wat_x, H2_Wat_y, H2_Wat_z);
	
	fclose(fOut);
}

void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in mol_H_donor.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}


void Setup_QM_Level(void)
{
	FILE *fIn;
	char *ReadLine, szLine[256], szSubStr_1[256], szSubStr_2[256], ErrorMsg[256];
	char szKey_QM_MEM[256], szKey_QM_NPROC[256], szQM_Level_Opt[256], szQM_Level_Dimer_Opt[256], szQM_Level_ESP[256], szQM_Level_E_Dimer[256], szQM_Level_E_Monomer[256];
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
      strcpy(szKey_QM_MEM, getenv("MEMORY" ));
      strcat( szKey_QM_MEM, "GB" );
  }
  // MJH stop
  //


	sprintf(szQM_Level_Int_E_Dimer_Cmd, "%%mem=%s\n%%nproc=%s\n%sMol-Water Dimer\n\n%d,1\n",
		szKey_QM_MEM, szKey_QM_NPROC, szQM_Level_E_Dimer, netcharge_mol);
	sprintf(szQM_Level_Int_E_Dimer_Opt_Cmd, "%%mem=%s\n%%nproc=%s\n%sMol-Water Dimer\n\n%d,1 \n",
		szKey_QM_MEM, szKey_QM_NPROC, szQM_Level_Dimer_Opt, netcharge_mol);
	sprintf(szQM_Level_Int_E_Single_Mol_Cmd, "%%mem=%s\n%%nproc=%s\n%s\n\nMol\n\n%d,1\n",
		szKey_QM_MEM, szKey_QM_NPROC, szQM_Level_E_Monomer, netcharge_mol);
	sprintf(szQM_Level_Int_E_Single_Wat_Cmd, "%%mem=%s\n%%nproc=%s\n%s\n\nWater\n\n0,1\n",
		szKey_QM_MEM, szKey_QM_NPROC, szQM_Level_E_Monomer);
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

int Try_To_Find_Tag(FILE *fIn, char szFileName[], char szTag[], char szLine[])
{
	char *ReadLine, ErrorMsg[256];

	while(1)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Fail to find the tag: %s in file %s\nQuit\n", szTag, szFileName);
      return 0;
		}

		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			sprintf(ErrorMsg, "Fail to find the tag: %s in file %s\nQuit\n", szTag, szFileName);
      return 0;
		}
		else	{
			if(FindString(szLine, szTag) >= 0)	{
				return 1;
			}
		}
	}

	return 0;
}



int To_Find_Tag(FILE *fIn, char szFileName[], char szTag[], char szLine[])
{
	char *ReadLine, ErrorMsg[256];

	while(1)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Fail to find the tag: %s in file %s\nQuit\n", szTag, szFileName);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}

		ReadLine = fgets(szLine, 256, fIn);
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

void ReadPreviousCG(void)
{
	FILE *fIn;
	int nAtom, i;

	fIn = fopen(szCGList, "r");
	if(fIn == NULL)	{
		return;
	}
	nAtom = Dimer.nAtom - 3;
	for(i=0; i<nAtom; i++)	{
		fscanf(fIn, "%lf", &(Dimer.CG[i]));
		Mol.CG[i] = Dimer.CG[i];	// assign same charges
	}
	fclose(fIn);

	Mol.Setup_NonBondParameters();
	Dimer.Setup_NonBondParameters();
	
	printf("Read previous charges.\n");
}



double CalRMSD(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum)
{
	double RMSD, midxa, midya, midza, midxb, midyb, midzb;
	double xa_tmp[MAX_ATOM], ya_tmp[MAX_ATOM], za_tmp[MAX_ATOM];
	double xb_tmp[MAX_ATOM], yb_tmp[MAX_ATOM], zb_tmp[MAX_ATOM];
	int i, iMax;

	memcpy(xa_tmp+1, xa, sizeof(double)*AtomNum);
	memcpy(ya_tmp+1, ya, sizeof(double)*AtomNum);
	memcpy(za_tmp+1, za, sizeof(double)*AtomNum);

	memcpy(xb_tmp+1, xb, sizeof(double)*AtomNum+3);	// water is included
	memcpy(yb_tmp+1, yb, sizeof(double)*AtomNum+3);	// water is included
	memcpy(zb_tmp+1, zb, sizeof(double)*AtomNum+3);	// water is included


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
	}
	iMax = AtomNum+3;
	for(i=1; i<=iMax; i++)	{
		xb_tmp[i]-=midxb;
		yb_tmp[i]-=midyb;
		zb_tmp[i]-=midzb;
	}

	quatfit(xa_tmp, ya_tmp, za_tmp, xb_tmp, yb_tmp, zb_tmp, AtomNum);
	RMSD=rmsfit(xa_tmp, ya_tmp, za_tmp, xb_tmp, yb_tmp, zb_tmp, AtomNum);

	for(i=1; i<=iMax; i++)	{	// generate the rotated coodinates
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
	int ia, ia_Max;

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

	ia_Max = n+3;	// include water
//	for(ia=1; ia<=n; ia++)	{
	for(ia=1; ia<=ia_Max; ia++)	{
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

double Cal_Dist_Two_Atoms(int ia, int ib)
{
	double dx, dy, dz, r;

	dx = Dimer.x[ia] - Dimer.x[ib];
	dy = Dimer.y[ia] - Dimer.y[ib];
	dz = Dimer.z[ia] - Dimer.z[ib];
	r = sqrt(dx*dx + dy*dy + dz*dz);
	return r;
}

void SaveDummyAtoms(void)
{
	x_Dummy_1_Save = x_Dummy_1;
	y_Dummy_1_Save = y_Dummy_1;
	z_Dummy_1_Save = z_Dummy_1;
	
	x_Dummy_2_Save = x_Dummy_2;
	y_Dummy_2_Save = y_Dummy_2;
	z_Dummy_2_Save = z_Dummy_2;

	x_Dummy_3_Save = x_Dummy_3;
	y_Dummy_3_Save = y_Dummy_3;
	z_Dummy_3_Save = z_Dummy_3;
	
	v_x_X1_H_Save = v_x_X1_H;
	v_y_X1_H_Save = v_y_X1_H;
	v_z_X1_H_Save = v_z_X1_H;
}

void RestoreDummyAtoms(void)
{
	x_Dummy_1 = x_Dummy_1_Save;
	y_Dummy_1 = y_Dummy_1_Save;
	z_Dummy_1 = z_Dummy_1_Save;
	                          
	x_Dummy_2 = x_Dummy_2_Save;
	y_Dummy_2 = y_Dummy_2_Save;
	z_Dummy_2 = z_Dummy_2_Save;
	
	x_Dummy_3 = x_Dummy_3_Save;
	y_Dummy_3 = y_Dummy_3_Save;
	z_Dummy_3 = z_Dummy_3_Save;
	
	v_x_X1_H = v_x_X1_H_Save;
	v_y_X1_H = v_y_X1_H_Save;
	v_z_X1_H = v_z_X1_H_Save;
}



