/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

// This is to be used in the 150-water-fitting step. It sets up and processes the QM calculations
// for the separate water and molecule QM calcs, which previously acceptor and donor duplicated 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

class CWORKER
{
public:
	int JobID, Status, nCore;
	char szInput[256], szOutput[256];

	void Start_Job(int CoreNum, int ID, char szGjf[], char szOut[]);	// start running a job
	int Get_Job_Status(void);
	void Update_Core_Number_Input_File(void);
};

CWORKER *pWorker;
int nJob=0, nJobDone=0, nWorker=0, nCore_Per_Node=1;
int JobStatus[MAX_JOB], Acceptor_List[MAX_JOB], AcceptorIdx_List[MAX_JOB], Theta_List[MAX_JOB];
char *szDirList[MAX_JOB], *szInputList[MAX_JOB], *szOutputList[MAX_JOB];
double Para_MM_Save[MAX_JOB][2];
const char szGAUSS_SCRDIR[]="GAUSS_SCRDIR";
char szGAUSS_SCRDIR_Base[256];

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

#include <limits.h>

char szExe_G09[PATH_MAX];
char szExe_Cale_MM[PATH_MAX];
//char szMyPath[]="../mypath.txt";
char szMolWat_E_Int[]="E-mol-wat.txt";
char szCGList[]="cg-list.txt";
char szCurDir[PATH_MAX];
char szWorkDir[PATH_MAX];

double x_Wat_MM[3], y_Wat_MM[3], z_Wat_MM[3];	// to store the optimized MM water configurations
double x_Wat_QM[3], y_Wat_QM[3], z_Wat_QM[3];	// to store the optimized MM water configurations

int Idx_Acceptor, ActiveAtom, Neighbor;

#define N_PARA	(2)
double Para_List[N_PARA], Grad_List[N_PARA], Para_Best[N_PARA], Low_Bound[N_PARA], Up_Bound[N_PARA];
char szElemName[MAX_ATOM_MOL][8];

int n_Atom_In_Mol;
int ResIDMol[MAX_ATOM_MOL];
int Idx_O_Wat, Idx_H1_Wat, Idx_H2_Wat;

double x_ReadIn[MAX_ATOM_MOL], y_ReadIn[MAX_ATOM_MOL], z_ReadIn[MAX_ATOM_MOL];

double x_Opt_QM[MAX_ATOM_MOL], y_Opt_QM[MAX_ATOM_MOL], z_Opt_QM[MAX_ATOM_MOL];
double r_Emin_MM=0.0, r_Emin_QM=0.0, E_Min_MM=0.0, E_Min_QM=0.0;
//double E_Min_MM_Shifted=0.0;

char szNameMolCRD[256];
char szAtomNameMol[MAX_ATOM_MOL][8], szResNameMol[MAX_ATOM_MOL][8];
//double x_Mol[MAX_ATOM_MOL], y_Mol[MAX_ATOM_MOL], z_Mol[MAX_ATOM_MOL];
double *x_Mol, *y_Mol, *z_Mol;
double E_Min=1.0E20;

double H1_Wat_x, H1_Wat_y, H1_Wat_z;
double O_Wat_x, O_Wat_y, O_Wat_z;
double H2_Wat_x, H2_Wat_y, H2_Wat_z;

double E_Dimer_Int, E_Dimer, E_Dimer_Min, E_Dimer_Far;

#define R_DUMMY_1	(1.5)
#define R_DUMMY_2	(1.5)
double x_Dummy_1, y_Dummy_1, z_Dummy_1, x_Dummy_1_Save, y_Dummy_1_Save, z_Dummy_1_Save, x_Dummy_1_Org, y_Dummy_1_Org, z_Dummy_1_Org;
double x_Dummy_2, y_Dummy_2, z_Dummy_2, x_Dummy_2_Save, y_Dummy_2_Save, z_Dummy_2_Save, x_Dummy_2_Org, y_Dummy_2_Org, z_Dummy_2_Org;
double x_Dummy_3, y_Dummy_3, z_Dummy_3;
double v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor, v_x_X1_Acceptor_Save, v_y_X1_Acceptor_Save, v_z_X1_Acceptor_Save;
double Para_Save_for_QM[N_PARA];
int theta, theta_Save;

CForceField ForceField;
CMol Mol, Dimer;

int netcharge_mol;
void Get_Netcharge_From_Xpsf(void);

char szName_Gaussian[]="dimer/qm.gjf";
char szFileElem[]="elem-list.txt";
char szFilePsf[]="mol.xpsf";
char szFileDimerPsf[]="mol-wat.xpsf";
char szFileCrd[]="mol-opt.xyz";

int Is_N[MAX_ATOM_MOL], Is_O[MAX_ATOM_MOL], Is_F[MAX_ATOM_MOL], Is_S[MAX_ATOM_MOL];
int Bond_Count[MAX_ATOM_MOL], BondList[MAX_ATOM_MOL][4];


void ReadMolCRD(void);
void GenerateWater(void);
double Normalizex(double& vect_x, double& vect_y, double& vect_z);
void OutputCRD(void);
void Generate_Gaussian_Input(void);
void Generate_Ini_Configuration(void);
double Optimize_Water_Pose(void);
void Output_Optimize_Water_Pose_QM(void);

void Quit_With_Error_Msg(char szMsg[]);
void Make_List_H_Acceptor(void);
void GetBondList(int Idx);
int Is_N_Planar(int Idx, double& x, double& y, double& z);
int Is_Three_Atoms_Colinear(int ia, int ib, int ic);

void SaveDummyAtoms(void);
void RestoreDummyAtoms(void);

void Generate_Dummy_Atoms_Bond_1(int Idx);
void Generate_Dummy_Atoms_Bond_2(int Idx);
void Generate_Dummy_Atoms_Bond_3(int Idx);
void Generate_Dummy_Atoms_Bond_Perp(int Idx);	// Sign = 0 / 1, positive or negative
void Generate_Dummy_Atoms_Bond_Perp_Flipped(int Idx);	// Sign = 0 / 1, positive or negative
void Rotate_Dummy_Atoms(int Idx, double theta);

double Extract_E_Gaussian_Output(char szFileName[], char szTag[]);
void Get_E_Mol_Water_Far_QM_prepare();
double Get_E_Mol_Water_Far_QM_complete();
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
	int i, j, Phi=0, AcceptorSelect = -1;
	int ToGenerate_Water;
	FILE *fOut , *fIn;
	char szCmd[256], ErrorMsg[256], *szEnv;

  timebomb();

	getcwd(szCurDir, 256);

	fFile_Run_Log = fopen(szLogName, "a+");




	fseek(fFile_Run_Log, SEEK_END, 0);
//	Get_EXE_Path("G09_EXE_PATH", szExe_G09);
	Get_EXE_Path("BIN_G09", szExe_G09);
	
	ForceField.ReadForceField(szMMForceField);
	Dimer.ReadPSF(szFileDimerPsf, 1);
	Dimer.AssignForceFieldParameters(&ForceField);
	
	Mol.ReadPSF(szFilePsf, 0);
	Mol.AssignForceFieldParameters(&ForceField);

	printf("Expecting input in [%s]", szFileCrd );

	Mol.ReadXYZ(szFileCrd);

	for(i=0; i<Mol.nAtom; i++)	{	// Readin QM optimized structure
		x_ReadIn[i] = Mol.x[i];
		y_ReadIn[i] = Mol.y[i];
		z_ReadIn[i] = Mol.z[i];
	}
	
	n_Atom_In_Mol = Mol.nAtom;

  fIn = fopen(szFileElem, "r");
  if(fIn == NULL) {
    sprintf(ErrorMsg, "Fail to open file %s\nQuit\n", szFileElem);
    Quit_With_Error_Msg(ErrorMsg);
  }
  for(i=0; i<n_Atom_In_Mol; i++)  {
    fscanf(fIn, "%s", szElemName[i]);
  }
  fclose(fIn);

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


	printf( "Running QM on water and molecule separately.. \n" );
	if( argc==2 && ! strcmp( argv[1], "--prepare" ) ) {
		Get_E_Mol_Water_Far_QM_prepare();
	}
	else if( argc==2 && ! strcmp( argv[1], "--complete" ) ) {
		Get_E_Mol_Water_Far_QM_complete();
	}
	else {
		printf("Neither --prepare nor --complete  passed on command line\n");
		exit(1);
	}
	return 0;
}


double Get_E_Mol_Water_Far_QM_complete()
{
	FILE *fOut, *fIn;
	char szFileEnergy[]="energy_sum.txt";
	int i, ReadItem;
	char szName_Gaussian_Single_Mol[]="single-mol/qm.gjf", szName_Gaussian_water[]="single-wat/qm.gjf";
	char szResult_Single_Mol[]="single-mol/qm.out", szResult_Single_Wat[]="single-wat/qm.out";
	double E_Mol, E_Water, E_Mol_Water;

	fIn = fopen(szFileEnergy, "r");
	if(fIn != NULL)	{
		ReadItem = fscanf(fIn, "%lf", &E_Mol_Water);
		fclose(fIn);
		if(ReadItem == 1)	{
			return E_Mol_Water;
		}
	}

	E_Mol = Extract_E_Gaussian_Output(szResult_Single_Mol, " SCF Done:");

	E_Water = Extract_E_Gaussian_Output(szResult_Single_Wat, " SCF Done:");
	
	fOut = fopen(szFileEnergy, "w");
	fprintf(fOut, "%18.14E\n", E_Mol+E_Water);
	fclose(fOut);

	return (E_Mol+E_Water);
}

void Get_E_Mol_Water_Far_QM_prepare()
{
	FILE *fOut, *fIn;
	char szFileEnergy[]="energy_sum.txt";
	int i, ReadItem;
	char szName_Gaussian_Single_Mol[]="single-mol/qm.gjf", szName_Gaussian_water[]="single-wat/qm.gjf";
	char szResult_Single_Mol[]="single-mol/qm.out", szResult_Single_Wat[]="single-wat/qm.out";
	double E_Mol, E_Water, E_Mol_Water;

	unlink( szFileEnergy );

	make_directory( "single-mol", 0700 );
	make_directory( "single-wat", 0700 );

	fOut = fopen(szName_Gaussian_Single_Mol, "w");
	fprintf(fOut, "%s", szQM_Level_Int_E_Single_Mol_Cmd);
	for(i=0; i<n_Atom_In_Mol; i++)	{
		fprintf(fOut, "%s  %10.5lf%10.5lf%10.5lf\n", szElemName[i], x_ReadIn[i], y_ReadIn[i], z_ReadIn[i]);
	}
	fprintf(fOut, "\n\n");
	fclose(fOut);


	fOut = fopen(szName_Gaussian_water, "w");
	fprintf(fOut, "%s", szQM_Level_Int_E_Single_Wat_Cmd);
	fprintf(fOut, "O  %10.5lf%10.5lf%10.5lf\n", 0.0, 0.0, 0.0);
	fprintf(fOut, "H  %10.5lf%10.5lf%10.5lf\n", b0_H_O_HF, 0.0, 0.0);
	fprintf(fOut, "H  %10.5lf%10.5lf%10.5lf\n\n", b0_H_O_HF*cos(theta0_H_O_H_HF*radianInv), b0_H_O_HF*sin(theta0_H_O_H_HF*radianInv), 0.0);
	fprintf(fOut, "\n\n");
	fclose(fOut);
	
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

  if( !fIn ) {
		printf("Could not open file [%s]\n", szFileName );
		exit(1);
	}

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
		fflush(fFile_Run_Log);
		fclose(fIn);
		Quit_With_Error_Msg(ErrorMsg);
	}
	return 0.0;
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
	H1_Wat_x = x_Mol[ActiveAtom] + r*v_x_X1_Acceptor;
	H1_Wat_y = y_Mol[ActiveAtom] + r*v_y_X1_Acceptor;
	H1_Wat_z = z_Mol[ActiveAtom] + r*v_z_X1_Acceptor;

	x_Dummy_3 = H1_Wat_x + (x_Dummy_2 - x_Mol[ActiveAtom]);
	y_Dummy_3 = H1_Wat_y + (y_Dummy_2 - y_Mol[ActiveAtom]);
	z_Dummy_3 = H1_Wat_z + (z_Dummy_2 - z_Mol[ActiveAtom]);

	Gen_xyz(x_Dummy_2, y_Dummy_2, z_Dummy_2, x_Dummy_3, y_Dummy_3, z_Dummy_3, 
		H1_Wat_x, H1_Wat_y, H1_Wat_z, b0_H_O_EXP, PI*0.5, PI, O_Wat_x, O_Wat_y, O_Wat_z);

	Gen_xyz(x_Dummy_3, y_Dummy_3, z_Dummy_3, H1_Wat_x, H1_Wat_y, H1_Wat_z, 
		O_Wat_x, O_Wat_y, O_Wat_z, b0_H_O_EXP, theta0_H_O_H_EXP*radianInv, Para_List[1], H2_Wat_x, H2_Wat_y, H2_Wat_z);

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
	fprintf(fOut, "Error in mol_H_acceptor.cpp\n");
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
	
	v_x_X1_Acceptor_Save = v_x_X1_Acceptor;
	v_y_X1_Acceptor_Save = v_y_X1_Acceptor;
	v_z_X1_Acceptor_Save = v_z_X1_Acceptor;
}

void RestoreDummyAtoms(void)
{
	x_Dummy_1 = x_Dummy_1_Save;
	y_Dummy_1 = y_Dummy_1_Save;
	z_Dummy_1 = z_Dummy_1_Save;
	                          
	x_Dummy_2 = x_Dummy_2_Save;
	y_Dummy_2 = y_Dummy_2_Save;
	z_Dummy_2 = z_Dummy_2_Save;
	
	v_x_X1_Acceptor = v_x_X1_Acceptor_Save;
	v_y_X1_Acceptor = v_y_X1_Acceptor_Save;
	v_z_X1_Acceptor = v_z_X1_Acceptor_Save;
}




