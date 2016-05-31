/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
//#include <mpi.h> eliot

#ifdef _WIN32
#include <direct.h>	//for win
#else
#include <unistd.h>	//for linux
#endif

#include "ff.h"
#include "nlopt.h"

#define radian	(57.29577951308232088)
#define radianInv	(0.017453292519943295)
#define PI	(3.14159265358979323846)
#define PI2	(6.28318530717958647692)

#define MAX_NUM_ATOM	(256)
#define MAX_ROTAMER		(4096)

//start	to data related with MPI and parallelization
//#define MY_MPI_ROOT	(0)
//#define MAX_N_PROC	(512)
int	ProgID=0, nProc=1;
//end	to data related with MPI and parallelization

//start	to data used for rotamer
double w_rotamer_List[MAX_ROTAMER], w_rotamer_Sum;
double x_Rotamer[MAX_ROTAMER][MAX_NUM_ATOM], y_Rotamer[MAX_ROTAMER][MAX_NUM_ATOM], z_Rotamer[MAX_ROTAMER][MAX_NUM_ATOM];
double x_Rotamer_Opt[MAX_ROTAMER][MAX_NUM_ATOM], y_Rotamer_Opt[MAX_ROTAMER][MAX_NUM_ATOM], z_Rotamer_Opt[MAX_ROTAMER][MAX_NUM_ATOM];
double E_Rotamer_QM[MAX_ROTAMER], E_Min_Rotamer, E_Rotamer_MM[MAX_ROTAMER], rmsd_Rotamer[MAX_ROTAMER];
int nRotamer=-1;
void Read_QM_Rotamer_Data(void);
void Cal_E_MM_Rotamer(void);
void Setup_Weight_Rotamer(void);
//end	to data used for rotamer


#ifndef MINAB 
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif


CMol Mol_ESP;
CForceField ForceField;
FILE *fPara;
FILE *fFitting;

double netcharge_mol=0.0;
double w_Rotamer=0.05;


#define N_MAX_DIH	(32)
#define N_PARA_MAX	(N_MAX_DIH*10+1)

int n_Para;
int Is_Para_Fixed[N_PARA_MAX];
double Para_List[N_PARA_MAX], Para_Best[N_PARA_MAX];
double Chi_SQ, Chi_SQ_1D, Chi_SQ_Rotamer, Chi_SQ_Min, Chi_SQ_Min_Global=1.0E100;

double CalObjectiveFunction(double Para_List[], int Get_E_Rotamer);

char szPhiToScan[]="soft-dih-list.txt";
int n_Phi=0;
int DihList[N_MAX_DIH][4], IdxDihSelect[N_MAX_DIH];
double Phi_To_Set[N_MAX_DIH];
int Read_Soft_DihedralList(void);
void To_Setup_All_Dihedral_Constraint(void);

void Setup_Phi_Constraint_List(void);
//int Is_Soft_Dihedral(int ib, int ic);
int Is_Soft_Dihedral(int ia, int ib, int ic, int id);

//start	to data and functions related with constrained optimization
double Geoometry_Optimization_With_Constraint(CMol *pMol);
double func_Geo_Opt_Constraint(unsigned n, const double *x, double *grad, void *my_func_data);
double Geo_Constraint_Bond(unsigned n, const double *x, double *grad, void *data);
double Geo_Constraint_Angle(unsigned n, const double *x, double *grad, void *data);
double Geo_Constraint_Dihedral(unsigned n, const double *x, double *grad, void *data);
//end	to data and functions related with constrained optimization

//start	data related with torsion fitting
#define MAX_SCAN	(72)
#define HARTREE_To_KCAL	(627.509469)
#define N_TOR_PARA	(11)

//int nScan;
//double Phi_Scaned[MAX_SCAN], E_QM_Scaned[MAX_SCAN], E_MM_Scaned[MAX_SCAN], E_QM_MM_Diff[MAX_SCAN], E_QM_MM_Diff_Fit[MAX_SCAN], dE[MAX_SCAN], w_Point[MAX_SCAN];
//int IsPara_Fixed[N_TOR_PARA], FailCount, Iteration;
double Para_Tor_List[N_MAX_DIH][N_TOR_PARA];
int n_Tor_Para_List[N_MAX_DIH];
double Ini_Tor_Para[N_MAX_DIH][32][N_TOR_PARA];

int nScan_List[N_MAX_DIH];
double w_Scan_List[N_MAX_DIH][MAX_SCAN], w_Scan_Sum;
double Phi_Scan_List[N_MAX_DIH][MAX_SCAN], x_Scan_list[N_MAX_DIH][MAX_SCAN][MAX_NUM_ATOM], y_Scan_list[N_MAX_DIH][MAX_SCAN][MAX_NUM_ATOM], z_Scan_list[N_MAX_DIH][MAX_SCAN][MAX_NUM_ATOM];
double E_Scan_MM[N_MAX_DIH][MAX_SCAN], E_Scan_MM_Mean, E_Scan_QM[N_MAX_DIH][MAX_SCAN], E_Scan_QM_Mean, E_Shift_0;

nlopt_opt opt;


void Read_1D_Scan_QM_Data(void);

double Callback_Eval_Gradient(unsigned n, const double *x, double *grad, void *my_func_data);
void Read_Tor_Para_1D_Fitting(void);
void Assign_Torsion_Parameters(void);
void Collect_Torsion_Parameters(void);
void Distribute_Torsion_Parameters(void);
void Cal_E_MM_Scaned(void);
double Optimize_Torsion_Parameters(void);
double Fitting_Torsion_Parameters(void);
void Enumerate_Ini_Torsion_Parameters(int IdxDih);	// depth
void Save_Parameters(void);
void Read_Parameters(char szName[]);
void Output_1D_QM_MM(void);
void Setup_Weight_1D_Scan(void);
void Output_Rotamer_QM_MM(void);

double Do_Geometry_Optimization(CMol *pMol);
double func_Geo_Optimization(unsigned n, const double *x, double *grad, void *my_func_data);


void Output_1D_QM_MM(void)
{
	FILE *fOut;
	char szName[256];
	int i, j, nAtom;
	double E_Min_1D_QM, E_Min_1D_MM_Mean, E_Min_1D_QM_Mean, E_Min_1D_QM_MM_Shift, Weight_Sum;

	if(ProgID != 0)	{
		return;
	}

	nAtom = Mol_ESP.nAtom;
	for(i=0; i<n_Phi; i++)	{
		sprintf(szName, "1d-qm-mm-%d.dat", i+1);
		fOut = fopen(szName, "w");

		E_Min_1D_QM = 1.0E100;
		for(j=0; j<nScan_List[i]; j++)	{
			memcpy(Mol_ESP.x, x_Scan_list[i][j], sizeof(double)*nAtom);
			memcpy(Mol_ESP.y, y_Scan_list[i][j], sizeof(double)*nAtom);
			memcpy(Mol_ESP.z, z_Scan_list[i][j], sizeof(double)*nAtom);
			E_Scan_MM[i][j] = Mol_ESP.Cal_E(0);

			if(E_Min_1D_QM > E_Scan_QM[i][j])	{
				E_Min_1D_QM = E_Scan_QM[i][j];
			}
		}

		E_Min_1D_MM_Mean = E_Min_1D_QM_Mean = 0.0;
		Weight_Sum = 0.0;
		for(j=0; j<nScan_List[i]; j++)	{
			Weight_Sum += w_Scan_List[i][j];
			E_Min_1D_MM_Mean += (w_Scan_List[i][j] * E_Scan_MM[i][j] );
			E_Min_1D_QM_Mean += (w_Scan_List[i][j] * E_Scan_QM[i][j] );
		}
		E_Min_1D_QM_MM_Shift = (E_Min_1D_QM_Mean - E_Min_1D_MM_Mean)/Weight_Sum;


		for(j=0; j<nScan_List[i]; j++)	{
			fprintf(fOut, "%.1lf %.5lf %.5lf\n", Phi_Scan_List[i][j], E_Scan_QM[i][j]-E_Min_1D_QM, E_Scan_MM[i][j]+E_Min_1D_QM_MM_Shift-E_Min_1D_QM);	// the relative 1D energy
		}

		fclose(fOut);
	}
}

void Output_Rotamer_QM_MM(void)
{
	int i;
	FILE *fOut;

	CalObjectiveFunction(Para_List, 1);

	if(ProgID != 0)	{
		return;
	}

	fOut = fopen("rotamer-E.dat", "w");
	for(i=0; i<nRotamer; i++)	{
		fprintf(fOut, "%4d %.5lf %.5lf\n", i+1, E_Rotamer_QM[i], E_Rotamer_MM[i] + Para_List[10*n_Phi]);
	}
	fclose(fOut);
}

void Read_Parameters(char szName[])
{
	FILE *fIn;

	fIn = fopen(szName, "r");
	for(int i=0; i<n_Para; i++)	{
		fscanf(fIn, "%lf", &(Para_List[i]));
	}
	fclose(fIn);

	Distribute_Torsion_Parameters();
	Assign_Torsion_Parameters();

	Cal_E_MM_Scaned();
	Cal_E_MM_Rotamer();

	Output_1D_QM_MM();
}

void Save_Parameters(void)
{
	FILE *fOut;
	int i;

	if(ProgID == 0)	{
		fOut = fopen("saved-para.dat", "w");
		for(i=0; i<n_Para; i++)	{
			fprintf(fOut, "%.15lf\n", Para_List[i]);
		}
		fprintf(fOut, "\nChi_SQ = %.12lf\n", Chi_SQ);
		fclose(fOut);
	}

	Output_1D_QM_MM();
	Output_Rotamer_QM_MM();
}

double Fitting_Torsion_Parameters(void)
{
	Enumerate_Ini_Torsion_Parameters(0);

	return Chi_SQ_Min;
}

int Run_Opt=1;
void Enumerate_Ini_Torsion_Parameters(int IdxDih)	// depth
{
	int iState;

	if(IdxDih < n_Phi)	{
		for(iState=0; iState<n_Tor_Para_List[IdxDih]; iState++)	{
			memcpy(Para_Tor_List[IdxDih], Ini_Tor_Para[IdxDih][iState], sizeof(double)*10);	// 10 parameters for one dihedral
			Enumerate_Ini_Torsion_Parameters(IdxDih+1);
		}
	}
	else	{	// to do parameter optimization
		fprintf(fFitting, "Starting optimization with %6d initial parameter:\n", Run_Opt);
		fflush(fFitting);
		Run_Opt++;

		Optimize_Torsion_Parameters();
	}
}

int FailCount=0;
double Optimize_Torsion_Parameters(void)
{
	double lb[N_PARA_MAX], ub[N_PARA_MAX];
	double Chi_SQ;
	int i, j, Pos;

	n_Para = 10*n_Phi + 1;

	memset(Is_Para_Fixed, 0, sizeof(int)*n_Para);
//	memset(lb, 0, sizeof(double)*n_Para);
	lb[n_Phi*10] = -HUGE_VAL;
	for(i=0; i<n_Phi; i++)	{
		Pos = 10*i;
		for(j=0; j<10; j+=2)	{	// ten parameters
			lb[Pos+j] = -15.0;
			ub[Pos+j] = 15.0;
		}
		for(j=1; j<10; j+=2)	{	// ten parameters
			lb[Pos+j] = -M_PI;
			ub[Pos+j] = +M_PI;
		}

	}
	ub[n_Phi*10] = HUGE_VAL;


	opt = nlopt_create(NLOPT_LN_COBYLA, n_Para);

	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, ub);
	nlopt_set_min_objective(opt, Callback_Eval_Gradient, NULL);
	nlopt_set_xtol_rel(opt, 1e-6);


	Collect_Torsion_Parameters();
	Assign_Torsion_Parameters();
	Cal_E_MM_Scaned();
	Cal_E_MM_Rotamer();
	Para_List[n_Phi*10] = E_Scan_QM_Mean-E_Scan_MM_Mean;	// energy shift
	E_Shift_0 = fabs(E_Scan_QM_Mean-E_Scan_MM_Mean);
	for(i=0; i<n_Phi; i++)	{
		Pos = 10*i;

		for(j=1; j<10; j+=2)	{	// the phase is fixed
			Is_Para_Fixed[Pos+j] = 1;
		}
	}

	Read_Parameters("saved-para.dat");

	Chi_SQ_Min = 1.0E100;
	memcpy(Para_Best, Para_List, sizeof(double)*n_Para);

	Chi_SQ = CalObjectiveFunction(Para_List, 0);
	if(ProgID == 0)	{
		fprintf(fFitting, "Iteration %4d  Chi^2 = %.8lf  Chi^2(1D) = %.8lf  Chi^2(rotamer) = %.8lf\n", 
			0, Chi_SQ, Chi_SQ_1D, Chi_SQ_Rotamer);
		fflush(fFitting);
	}

	FailCount = 0;
  int ret= nlopt_optimize(opt, Para_List, &Chi_SQ);
	if (ret < 0) {
		printf("nlopt failed! :%d\n", ret);
	}
	else {
		printf("Chi_SQ_min = %lf  Chi_SQ = %lf\n", Chi_SQ_Min, Chi_SQ);
	}


	memcpy(Para_List, Para_Best, sizeof(double)*n_Para);
	CalObjectiveFunction(Para_List, 1);

	nlopt_destroy(opt);

	return Chi_SQ_Min;
}

//end	data related with torsion fitting


double rmsfit(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);
void quatfit(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);
void jacobi(int n, int np, double a[5][5], double d[5], double v[5][5], double b[5], double z[5]);
double CalRMSD(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);





FILE *fFile_Conf;
FILE *fFile_Run_Log;	// will be shared by other source code
char szName_Run_Log[]="fitting-1d-rotamer-log.txt";
char szName_Conf_File[256];
char szName_XPSF[256];
char szName_CRD[256];
char szName_Force_Field[512];
int Flag_Target_Int_Dimer=0, Flag_Target_Int_Water=0;
int IsTestMode = 0;
char szFinalPara[]="final-para.txt";

extern int FindString(char szBuff[], char szTag[]);

//start	to data used to save the read in of configuration file
#define MAX_LINE_CONF		(256)
#define MAX_LINE_CONF_LEN	(512)
char szReadConf_File[MAX_LINE_CONF][MAX_LINE_CONF_LEN];
int n_Line_Conf_File=0;
int ReadConfFile(char szName[]);
int FindItemInConfFile(char szItem[], int LineStart);
//end	to data used to save the read in of configuration file

//start	to data and subroutines related with update the rtf file
#define MAX_LINE_RTF	(1000)
#define MAX_LEN_RTF		(256)
char szName_Old_Rtf[]="mol.rtf";
char szName_New_Rtf[]="mol-esp.rtf";
char szTxt_Rtf[MAX_LINE_RTF][MAX_LEN_RTF];

int To_Export_New_Rtf, nLine_Rtf, Idx_Res_Mol, Idx_Mol_Atom_Start, Idx_Mol_Atom_End;
void Read_Rtf_File(void);
void Update_Rtf_File(void);
//end	to data and subroutines related with update the rtf file

double w_H_Donor_Acceptor=1.0, w_charge=1.0, w_water_int_E=1.0;


void Quit_With_Error_Msg(char szMsg[]);
int FindTags_As_First_String(char szBuff[], char szTag[]);


#define DELTA_PARA	(1.0E-4)
void CalGradient(double Grad_List[])
{
	int i, nAtom;
	double Para_Save[N_PARA_MAX], Delta_Para_List[N_PARA_MAX];
	double f_Left, f_Right;
	
	nAtom = Mol_ESP.nAtom;

	for(i=0; i<n_Para; i++)	{
		Para_Save[i] = Para_List[i];	//save current parameters
		Delta_Para_List[i] = DELTA_PARA;
	}
	Delta_Para_List[n_Para - 1] *= E_Shift_0;

	Distribute_Torsion_Parameters();
	Assign_Torsion_Parameters();
	Cal_E_MM_Rotamer();	// do full optimization from the original coordinates once


	for(i=0; i<n_Para; i++)	{
		if(Is_Para_Fixed[i] == 0)	{
			Para_List[i] = Para_Save[i] - Delta_Para_List[i];
			f_Left = CalObjectiveFunction(Para_List, 0);
			
			Para_List[i] = Para_Save[i] + Delta_Para_List[i];
			f_Right = CalObjectiveFunction(Para_List, 0);
			
			Grad_List[i] = (f_Right-f_Left)/(2.0*Delta_Para_List[i]);
			
			Para_List[i] = Para_Save[i];	//restore save (correct) parameter
		}
		else	{
			Grad_List[i] = 0.0;
		}
	}
	
	CalObjectiveFunction(Para_List, 0);
	
	return;
}
#undef	DELTA_PARA

double CalObjectiveFunction(double Para_List[], int Get_E_Rotamer)
{
	int i, j, nAtom, time_1, time_2;
	double d_E, E_Shift, w_1D;
	double Chi_SQ_Rotamer_Local;
	double E_Rotamer_MM_Local[MAX_ROTAMER];
	double E_Min_1D_MM_Mean, E_Min_1D_QM_Mean, E_Min_1D_QM_MM_Shift, Weight_Sum;

	Chi_SQ = Chi_SQ_1D = Chi_SQ_Rotamer = 0.0;
	Chi_SQ_Rotamer_Local = 0.0;

	w_1D = 1.0 - w_Rotamer;

	time_1 = time(NULL);

	Distribute_Torsion_Parameters();
	Assign_Torsion_Parameters();

	E_Shift = Para_List[10*n_Phi];

	nAtom = Mol_ESP.nAtom;
	for(i=0; i<n_Phi; i++)	{
		for(j=0; j<nScan_List[i]; j++)	{
			memcpy(Mol_ESP.x, x_Scan_list[i][j], sizeof(double)*nAtom);
			memcpy(Mol_ESP.y, y_Scan_list[i][j], sizeof(double)*nAtom);
			memcpy(Mol_ESP.z, z_Scan_list[i][j], sizeof(double)*nAtom);
			E_Scan_MM[i][j] = Mol_ESP.Cal_E(0);
		}

		E_Min_1D_MM_Mean = E_Min_1D_QM_Mean = 0.0;
		Weight_Sum = 0.0;
		for(j=0; j<nScan_List[i]; j++)	{
			Weight_Sum += w_Scan_List[i][j];
			E_Min_1D_MM_Mean += (w_Scan_List[i][j] * E_Scan_MM[i][j] );
			E_Min_1D_QM_Mean += (w_Scan_List[i][j] * E_Scan_QM[i][j] );
		}
		E_Min_1D_QM_MM_Shift = (E_Min_1D_QM_Mean - E_Min_1D_MM_Mean)/Weight_Sum;

		for(j=0; j<nScan_List[i]; j++)	{
			d_E = (E_Scan_QM[i][j] - E_Scan_MM[i][j] - E_Min_1D_QM_MM_Shift);	// fit the relative 1D energy profile
			Chi_SQ_1D += (d_E*d_E*w_Scan_List[i][j]);
		}
	}

	
	Chi_SQ_1D *= (w_1D/w_Scan_Sum);


	for(i=0; i<nRotamer; i++)	{
		E_Rotamer_MM[i] = E_Rotamer_MM_Local[i] = 0.0;
		if(i%nProc != ProgID)	{	//only proceed when i%nProc == ProgID, do job decomposition
			continue;
		}
		memcpy(Mol_ESP.x, x_Rotamer_Opt[i], sizeof(double)*nAtom);
		memcpy(Mol_ESP.y, y_Rotamer_Opt[i], sizeof(double)*nAtom);
		memcpy(Mol_ESP.z, z_Rotamer_Opt[i], sizeof(double)*nAtom);
		Mol_ESP.FullGeometryOptimization_LBFGS_step(1);
		E_Rotamer_MM_Local[i] = Mol_ESP.Cal_E(0);
		d_E = E_Rotamer_QM[i] - E_Rotamer_MM_Local[i] - E_Shift;
		Chi_SQ_Rotamer_Local += (d_E*d_E*w_rotamer_List[i]);
	}
	if(Get_E_Rotamer)	{
//eliot (from Matt Harvey)
	for( int i = 0 ; i < nRotamer; i++ ) { E_Rotamer_MM[i] += E_Rotamer_MM_Local[i]; }
                Chi_SQ_Rotamer = Chi_SQ_Rotamer_Local;

//eliot		MPI_Allreduce(E_Rotamer_MM_Local, E_Rotamer_MM, nRotamer, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}
//eliot	MPI_Allreduce(&Chi_SQ_Rotamer_Local, &Chi_SQ_Rotamer, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//	MPI_Bcast(&Chi_SQ_Rotamer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	printf(" w_rotamer_Sum = %f\n", w_rotamer_Sum );

	Chi_SQ_Rotamer *= (w_Rotamer/w_rotamer_Sum);

	Chi_SQ = Chi_SQ_1D + Chi_SQ_Rotamer;

//	MPI_Barrier(MPI_COMM_WORLD);
//	printf("%lf\n", Chi_SQ);
//	MPI_Finalize();

	time_2 = time(NULL);
/*
	if(ProgID == 0)	{
		FILE *fOut;
		fOut = fopen("time-rec.txt", "a+");
		fseek(fOut, 0, SEEK_END);
		fprintf(fOut, "%d\n", time_2-time_1);
		fclose(fOut);
	}
*/
	printf("Chi_SQ=%f\n", Chi_SQ );
	return Chi_SQ;
}

int Iteration;
double Callback_Eval_Gradient(unsigned n, const double *x, double *grad, void *my_func_data)
{
	memcpy(Para_List, x, sizeof(double)*n);

	if(grad)	{
		CalGradient(grad);	//gradient will be assigned into objGrad automatically
	}
	else {
		Distribute_Torsion_Parameters();
		Assign_Torsion_Parameters();
		Cal_E_MM_Rotamer();

		CalObjectiveFunction((double*)x, 0);
	}

		Iteration++;
		if(ProgID==0)	{
			fprintf(fFitting, "Iteration %4d  Chi^2 = %.8lf  Chi^2(1D) = %.8lf  Chi^2(rotamer) = %.8lf\n", 
				Iteration, Chi_SQ, Chi_SQ_1D, Chi_SQ_Rotamer);
			fflush(fFitting);
		}

		if(Chi_SQ < Chi_SQ_Min_Global)	{
			Save_Parameters();
			memcpy(Para_Best, Para_List, sizeof(double)*n);
			Chi_SQ_Min_Global = Chi_SQ;
      printf("Iteration %d : ", Iteration );
      for(int i=0; i<n; i++ ) { printf("%f ", Para_Best[i]); }
      printf(" : %f\n", Chi_SQ );
		}

		if(Chi_SQ < Chi_SQ_Min)	{
			Chi_SQ_Min = Chi_SQ;
			
			if(Chi_SQ + 2.0E-5 > Chi_SQ_Min)	{
				FailCount++;
			}
			else	{
				FailCount = 0;
			}
		}
		else	{
			FailCount++;
			if(FailCount > 6)      {       // cannot make further progress
        printf("Terminating..\n");
				nlopt_force_stop(opt);	// terminate the optimization
			}
		}

    return Chi_SQ;
}


int main(int argc, char **argv)
{
	char ErrorMsg[256];

  timebomb();

	fFile_Run_Log = fopen(szName_Run_Log, "w");
	if(fFile_Run_Log==NULL)	{
		sprintf(ErrorMsg, "Fail to create the log file.\n\n");
		Quit_With_Error_Msg(ErrorMsg);
	}

	strcpy(szName_Conf_File, argv[1]);
	ReadConfFile(szName_Conf_File);

//eliot	MPI_Init(&argc, &argv);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProgID);
//	MPI_Comm_size(MPI_COMM_WORLD, &nProc);
	

	ForceField.ReadForceField(szName_Force_Field);
	
	Mol_ESP.ReadPSF(szName_XPSF, 0);
	Read_Rtf_File();
	
	strcpy(Mol_ESP.MolName, Mol_ESP.ResName[0]);
	Mol_ESP.AssignForceFieldParameters(&ForceField);

	Read_Soft_DihedralList();

	Mol_ESP.Is_Phi_Psi_Constrained = 0;
	Mol_ESP.E_CMap_On = 0;

	Read_QM_Rotamer_Data();

  if( nRotamer == 0 ) {
    printf("No rotamers found. The should be at least one in all-rotamer.dat. Something has gone wrong in the previous step\n" );
		exit(0);
  }

	Read_Tor_Para_1D_Fitting();
	Assign_Torsion_Parameters();

	Read_1D_Scan_QM_Data();

//	Cal_E_MM_Scaned();
//	Cal_E_MM_Rotamer();

	fFitting = fopen("fitting.dat", "w");
	Fitting_Torsion_Parameters();
	fclose(fFitting);

//	FILE *fOut;
//	fOut = fopen("rotamer-E.dat", "w");
//	for(i=0; i<nRotamer; i++)	{
//		fprintf(fOut, "%d %lf %lf %lf\n", i+1, E_Rotamer_QM[i], E_Rotamer_MM[i], rmsd_Rotamer[i]);
//	}
//	fclose(fOut);

//	Cal_E_MM_QM_Diff(7);
	
//	for(i=0; i<n_Phi; i++)	{
//		Cal_E_MM_QM_Diff(i);
//	}

//	Fit_Torsion_Parameters(7);

//	fPara = fopen("torsion-para.dat", "w");
//	for(i=0; i<n_Phi; i++)	{
//		Fit_Torsion_Parameters(i);
//	}
//	fclose(fPara);


//	Geoometry_Optimization_With_Constraint(&Mol_ESP);

	fflush(fFile_Run_Log);
	

	fclose(fFile_Run_Log);

//	MPI_Finalize();
	
	return 0;
}

void To_Setup_All_Dihedral_Constraint(void)
{
	int i;

	Mol_ESP.nPhi_Fixed = 0;

	for(i=0; i<n_Phi; i++)	{
		Mol_ESP.To_Fix_Dihedral(DihList[i][0], DihList[i][1], DihList[i][2], DihList[i][3], 
			Mol_ESP.Query_Dihedral(DihList[i][0], DihList[i][1], DihList[i][2], DihList[i][3], 0, NULL));
	}
}

void GetSecondItem(char szLine[], char szSecond[])
{
	char szTmp[256], ErrorMsg[256];
	int ReadItem;

	ReadItem = sscanf(szLine, "%s %s", szTmp, szSecond);
	if(ReadItem != 2)	{
		sprintf(ErrorMsg, "Fail to extract the second string from, \n%s\nQuit\n", szLine);
		Quit_With_Error_Msg(ErrorMsg);
	}
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

int ReadConfFile(char szName[])
{
	int LineLen, i;
	FILE *fIn;
	char ErrorMsg[256];
	char szValue[32];

	n_Line_Conf_File = 0;
	
	fIn = fopen(szName, "r");
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		fgets(szReadConf_File[n_Line_Conf_File], MAX_LINE_CONF_LEN, fIn);
		LineLen = strlen(szReadConf_File[n_Line_Conf_File]);
		if( LineLen > 1)	{	// a valid line
			n_Line_Conf_File ++;
			if(n_Line_Conf_File > MAX_LINE_CONF)	{
				sprintf(ErrorMsg, "Error in ReadConfFile(). n_Line_Conf_File > MAX_LINE_CONF\nQuit\n");
				Quit_With_Error_Msg(ErrorMsg);
			}
		}
		if(LineLen >= MAX_LINE_CONF_LEN)	{
			sprintf(ErrorMsg, "Error in ReadConfFile(). LineLen >= MAX_LINE_CONF_LEN\nQuit\n");
			Quit_With_Error_Msg(ErrorMsg);
		}
	}
	fclose(fIn);

	//start	to assign some parameters 
	for(i=0; i<n_Line_Conf_File; i++)	{
		if(FindTags_As_First_String(szReadConf_File[i], "FILE_FORCE_FIELD"))	{
			GetSecondItem(szReadConf_File[i], szName_Force_Field);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "FILE_PSF"))	{
			GetSecondItem(szReadConf_File[i], szName_XPSF);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "FILE_XYZ"))	{
			GetSecondItem(szReadConf_File[i], szName_CRD);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "w_Rotamer"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			w_Rotamer = atof(szValue);
			fprintf(fFile_Run_Log, "w_Rotamer = %.3lf\n", w_Rotamer);
			fflush(fFile_Run_Log);
		}
		
	}
	//end	to assign some parameters 

	return n_Line_Conf_File;
}

int FindItemInConfFile(char szItem[], int LineStart)
{
	int Idx;

	for(Idx=LineStart; Idx<n_Line_Conf_File; Idx++)	{
		if(FindTags_As_First_String(szReadConf_File[Idx], szItem))	{
			return Idx;
		}
	}

	return -1;
}

void Read_Rtf_File(void)
{
	FILE *fIn;
	char *ReadLine;
	int i;

	fIn = fopen(szName_Old_Rtf, "r");
	if(fIn == NULL)	{
		fprintf(fFile_Run_Log, "Warning: Fail to open file %s\nThe new rtf will not be saved.\n", szName_Old_Rtf);
		To_Export_New_Rtf = 0;
		return;
	}
	else	{
		To_Export_New_Rtf = 1;
	}

	nLine_Rtf = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szTxt_Rtf[nLine_Rtf], 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		else	{
			nLine_Rtf++;
		}
	}

	fclose(fIn);

	for(i=0; i<nLine_Rtf; i++)	{
		if(FindString(szTxt_Rtf[i], "RESI MOL")>=0)	{
			sprintf(szTxt_Rtf[i], "RESI MOL %6.3lf\n", netcharge_mol);
			Idx_Res_Mol = i;
		}
		if(strncmp(szTxt_Rtf[i], "ATOM", 4)==0)	{
			break;
		}
	}
	Idx_Mol_Atom_Start = i;
	for(; i<nLine_Rtf; i++)	{
		if(strncmp(szTxt_Rtf[i], "ATOM", 4)!=0)	{
			break;
		}
	}
	Idx_Mol_Atom_End = i - 1;
	return;
}

void Update_Rtf_File(void)
{
	FILE *fOut;
	int i;

	if(To_Export_New_Rtf == 0)	{
		return;
	}

	fOut = fopen(szName_New_Rtf, "w");

	for(i=0; i<nLine_Rtf; i++)	{
		fprintf(fOut, "%s", szTxt_Rtf[i]);
	}
	fclose(fOut);
}

int Read_Soft_DihedralList(void)
{
	FILE *fIn;
	int ReadItem;
	char szLine[256], *ReadLine, ErrorMsg[256];

	n_Phi = 0;

	fIn = fopen(szPhiToScan, "r");

	while(1)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Fail to open file: %s\nQuit\n", szPhiToScan);
			Quit_With_Error_Msg(ErrorMsg);
		}
		ReadLine = fgets(szLine, 128, fIn);
		if(ReadLine == NULL)	{
			break;
		}

		ReadItem = sscanf(szLine, "%d %d %d %d", 
			&(DihList[n_Phi][0]), &(DihList[n_Phi][1]), &(DihList[n_Phi][2]), &(DihList[n_Phi][3]));

		if(ReadItem == 4)	{
			DihList[n_Phi][0]--;
			DihList[n_Phi][1]--;
			DihList[n_Phi][2]--;
			DihList[n_Phi][3]--;
			IdxDihSelect[n_Phi] = Mol_ESP.Query_Dihedral_Index(DihList[n_Phi][0], DihList[n_Phi][1], DihList[n_Phi][2], DihList[n_Phi][3]);

			if(IdxDihSelect[n_Phi] < 0)	{
				Quit_With_Error_Msg("Fail to identify the index of one soft dihedral.\n");
			}
			n_Phi++;
		}
		else	{
			break;
		}
	}

	fclose(fIn);

	Setup_Phi_Constraint_List();

	return n_Phi;
}


void Setup_Phi_Constraint_List(void)
{
	int i, iPos, ia, ib, ic, id, nDihedral, *DihedralList;

	nDihedral = Mol_ESP.nDihedral;
	DihedralList = Mol_ESP.DihedralList;

	// to constrain soft dihedrals
	for(i=0; i<nDihedral; i++)	{
		iPos = 4*i;

		ia = DihedralList[iPos+0];
		ib = DihedralList[iPos+1];
		ic = DihedralList[iPos+2];
		id = DihedralList[iPos+3];

//		if(Is_Soft_Dihedral(ib, ic))	{	// in the list of soft dihedrals
		if(Is_Soft_Dihedral(ia, ib, ic, id))	{	// in the list of soft dihedrals
			Mol_ESP.Is_Phi_Constrained[i] = 0;
		}
		else	{	// only constrain those dihedrals within rings / methyls
			Mol_ESP.Is_Phi_Constrained[i] = 1;
		}
	}
	for(i=0; i<n_Phi; i++)	{	// to constrain soft dihedrals
		Mol_ESP.Is_Phi_Constrained[IdxDihSelect[i]] = 1;
	}
}

int Is_Soft_Dihedral(int ia, int ib, int ic, int id)
{
	int i;

	for(i=0; i<n_Phi; i++)	{
//		if( ( (ib == DihList[i][1]) && (ic == DihList[i][2]) ) || ( (ib == DihList[i][2]) && (ic == DihList[i][1]) ) )	{
		if( ( (ia == DihList[i][0]) && (ib == DihList[i][1]) && (ic == DihList[i][2]) && (id == DihList[i][3]) ) || 
			( (ia == DihList[i][3]) && (ib == DihList[i][2]) && (ic == DihList[i][1]) && (id == DihList[i][0]) ) )	{
			return 1;
		}
	}
	return 0;
}


void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in 1D-rotamer-fitting.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
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

		n_Tor_Para_List[Idx_Phi] = 0;
		while(1)	{
			for(j=0; j<(N_TOR_PARA-1); j++)	{
//				ReadItem = fscanf(fIn, "%lf", &(Para_Tor_List[i][j]));
				ReadItem = fscanf(fIn, "%lf", &(Ini_Tor_Para[Idx_Phi][n_Tor_Para_List[Idx_Phi]][j]));
				if(ReadItem != 1)	{
					break;
				}
			}
			if(j==(N_TOR_PARA-1))	{	// a valid record
				n_Tor_Para_List[Idx_Phi]++;
			}
			else	{
				break;
			}
		}


		fclose(fIn);


		n_Tor_Para_List[Idx_Phi] = min(n_Tor_Para_List[Idx_Phi], 1);	// to control the number of optimizations !!!
	}
}

void Collect_Torsion_Parameters(void)
{
	int i, Pos;

	for(i=0; i<n_Phi; i++)	{
		Pos = i * 10;
		memcpy(Para_List+Pos, Para_Tor_List[i], sizeof(double)*10);
	}
}

void Distribute_Torsion_Parameters(void)
{
	int i, Pos;

	for(i=0; i<n_Phi; i++)	{
		Pos = i * 10;
		memcpy(Para_Tor_List[i], &(Para_List[Pos]), sizeof(double)*10);
	}
}

void Assign_Torsion_Parameters(void)
{
	int i_Phi, Idx;
	double *pPara_k_Dih, *pPara_phi;

	for(i_Phi=0; i_Phi<n_Phi; i_Phi++)	{
		Idx = IdxDihSelect[i_Phi];
		pPara_k_Dih = Mol_ESP.Para_k_Dih[Idx];
		pPara_phi = Mol_ESP.Para_phi[Idx];

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

void Read_QM_Rotamer_Data(void)
{
	FILE *fIn;
	char szLine[256], *ReadLine, szTmp[256], ErrorMsg[256];
	int ReadItem, i, nAtom;

	nAtom = Mol_ESP.nAtom;

	fIn = fopen("all-rotamer.dat", "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file: all-rotamer.dat\nQuit\n");
		Quit_With_Error_Msg(ErrorMsg);
	}
	nRotamer = 0;
	E_Min_Rotamer = 1.0E100;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}		
		if(strncmp(szLine, "E_Rotamer", 9)==0)	{
			ReadItem = sscanf(szLine, "%s %lf", szTmp, &(E_Rotamer_QM[nRotamer]));
			if(ReadItem == 2)	{
				fgets(szLine, 256, fIn);
				if(strncmp(szLine, "Coordinate", 10)==0)	{
					for(i=0; i<nAtom; i++)	{
						fgets(szLine, 256, fIn);
						sscanf(szLine, "%lf %lf %lf", &(x_Rotamer[nRotamer][i]), &(y_Rotamer[nRotamer][i]), &(z_Rotamer[nRotamer][i]));
					}
					if(E_Min_Rotamer > E_Rotamer_QM[nRotamer])	{
						E_Min_Rotamer = E_Rotamer_QM[nRotamer];
					}
					nRotamer++;
				}
			}
		}
	}
	fclose(fIn);

	for(i=0; i<nRotamer; i++)	{
		E_Rotamer_QM[i] = (E_Rotamer_QM[i] - E_Min_Rotamer)*HARTREE_To_KCAL;
	}

}

void Read_1D_Scan_QM_Data(void)
{
	int i, j, ReadItem, nAtom, Count;
	FILE *fIn;
	char szName[256], *ReadLine, szLine[256], szTmp[256], ErrorMsg[256];
	double fTmp;

	nAtom = Mol_ESP.nAtom;

	Count = 0;
	E_Scan_QM_Mean = 0.0;
	for(i=0; i<n_Phi; i++)	{
		sprintf(szName, "mm-tor-1D-idx-%d.dat", i+1);
		fIn = fopen(szName, "r");
		if(fIn == NULL)	{
			sprintf(ErrorMsg, "Fail to open file: %s\nQuit\n", szName);
			Quit_With_Error_Msg(ErrorMsg);
		}

		nScan_List[i] = 0;
		while(1)	{
			if(feof(fIn))	{
				break;
			}
			ReadLine = fgets(szLine, 256, fIn);
			if(ReadLine == NULL)	{
				break;
			}
			if(strncmp(szLine, "E_Scan", 6) == 0)	{
				ReadItem = sscanf(szLine, "%s %lf %lf %s %lf", szTmp, &(E_Scan_QM[i][nScan_List[i]]), &fTmp, szTmp, &(Phi_Scan_List[i][nScan_List[i]]));
				if(ReadItem == 5)	{
					ReadLine = fgets(szLine, 256, fIn);
					if(strncmp(szLine, "Coordinate", 10)==0)	{
						for(j=0; j<nAtom; j++)	{
							fgets(szLine, 256, fIn);
							sscanf(szLine, "%lf %lf %lf", &(x_Scan_list[i][nScan_List[i]][j]), &(y_Scan_list[i][nScan_List[i]][j]), &(z_Scan_list[i][nScan_List[i]][j])); 
						}
					}
				}
				E_Scan_QM[i][nScan_List[i]] = (E_Scan_QM[i][nScan_List[i]] - E_Min_Rotamer) * HARTREE_To_KCAL;
//				E_Scan_QM_Mean += E_Scan_QM[i][nScan_List[i]];
				Count++;
				nScan_List[i]++;
			}
		}

		fclose(fIn);
	}
//	E_Scan_QM_Mean /= Count;

	Setup_Weight_1D_Scan();
}

double CalRMSD(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum)
{
	double RMSD, midxa, midya, midza, midxb, midyb, midzb;
	double xa_tmp[MAX_NUM_ATOM], ya_tmp[MAX_NUM_ATOM], za_tmp[MAX_NUM_ATOM];
	double xb_tmp[MAX_NUM_ATOM], yb_tmp[MAX_NUM_ATOM], zb_tmp[MAX_NUM_ATOM];
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

//	for(i=1; i<=iMax; i++)	{	// generate the rotated coodinates
//		xb[i-1]=midxb+xb_tmp[i];
//		yb[i-1]=midyb+yb_tmp[i];
//		zb[i-1]=midzb+zb_tmp[i];
//	}


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

void Cal_E_MM_Scaned(void)
{
	int i, j, nAtom, Count=0;
	double E_Min_Local_QM;

	E_Scan_MM_Mean = 0.0;
	E_Scan_QM_Mean = 0.0;
	nAtom = Mol_ESP.nAtom;

	for(i=0; i<n_Phi; i++)	{
		E_Min_Local_QM = 1.0E100;
		for(j=0; j<nScan_List[i]; j++)	{
			if(E_Scan_QM[i][j] < E_Min_Local_QM)	{	// to find the lowest QM energy
				E_Min_Local_QM = E_Scan_QM[i][j];
			}
		}

		for(j=0; j<nScan_List[i]; j++)	{
			memcpy(Mol_ESP.x, x_Scan_list[i][j], sizeof(double)*nAtom);
			memcpy(Mol_ESP.y, y_Scan_list[i][j], sizeof(double)*nAtom);
			memcpy(Mol_ESP.z, z_Scan_list[i][j], sizeof(double)*nAtom);
			E_Scan_MM[i][j] = Mol_ESP.Cal_E(0);

			if( E_Scan_QM[i][j] < (E_Min_Local_QM+6.0) )	{	// skip configuration with much higher energies
				E_Scan_QM_Mean += E_Scan_QM[i][j];
				E_Scan_MM_Mean += E_Scan_MM[i][j];
				Count ++;
			}
		}
	}
	E_Scan_QM_Mean /= Count;
	E_Scan_MM_Mean /= Count;
}


double x_Best[MAX_NUM_ATOM*3];
nlopt_opt opt_Geo;
int Iter_Geo_Opt;
double Do_Geometry_Optimization(CMol *pMol)
{
	int i, nAtom, nDim, iPos, Return_Opt;
	double x[MAX_NUM_ATOM*3], x_Ini[MAX_NUM_ATOM], y_Ini[MAX_NUM_ATOM], z_Ini[MAX_NUM_ATOM], E_min=1.0E100;

//	pMol->Restrain_All_Torsions(1);

	nAtom = pMol->nAtom;
	nDim = 3*nAtom;

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
		pMol->Restrain_All_Torsions(0);
		return (pMol->Cal_E(0));
	}
	else	{	// try FullGeometryOptimization_LBFGS()
		memcpy(pMol->x, x_Ini, sizeof(double)*nAtom);	// restore the structure
		memcpy(pMol->y, y_Ini, sizeof(double)*nAtom);
		memcpy(pMol->z, z_Ini, sizeof(double)*nAtom);
		
		pMol->FullGeometryOptimization_LBFGS_step(0);
		pMol->Restrain_All_Torsions(0);
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


void Cal_E_MM_Rotamer(void)
{
	int i, nAtom;
	double E_Rotamer_MM_Local[MAX_ROTAMER];

//	fOut = fopen("coord-opt-mm.dat", "w");
	nAtom = Mol_ESP.nAtom;
	for(i=0; i<nRotamer; i++)	{
		E_Rotamer_MM[i] = E_Rotamer_MM_Local[i] = 0.0;
		if(i%nProc != ProgID)	{	//only proceed when i%nProc == ProgID, do job decomposition
			continue;
		}
		memcpy(Mol_ESP.x, x_Rotamer[i], sizeof(double)*nAtom);
		memcpy(Mol_ESP.y, y_Rotamer[i], sizeof(double)*nAtom);
		memcpy(Mol_ESP.z, z_Rotamer[i], sizeof(double)*nAtom);

		E_Rotamer_MM_Local[i] = Do_Geometry_Optimization(&Mol_ESP);

//		Mol_ESP.FullGeometryOptimization_LBFGS(0);
//		E_Rotamer_MM_Local[i] = Mol_ESP.Cal_E(0);
		
		memcpy(x_Rotamer_Opt[i], Mol_ESP.x, sizeof(double)*nAtom);
		memcpy(y_Rotamer_Opt[i], Mol_ESP.y, sizeof(double)*nAtom);
		memcpy(z_Rotamer_Opt[i], Mol_ESP.z, sizeof(double)*nAtom);
//		rmsd_Rotamer[i] = CalRMSD(Mol_ESP.x, Mol_ESP.y, Mol_ESP.z, x_Rotamer[i], y_Rotamer[i], z_Rotamer[i], nAtom);

//		for(j=0; j<nAtom; j++)	{
//			fprintf(fOut, "%.6lf %.6lf %.6lf\n", Mol_ESP.x[j], Mol_ESP.y[j], Mol_ESP.z[j])
//		}
	}
//	fclose(fOut);

//	MPI_Allreduce(E_Rotamer_MM_Local, E_Rotamer_MM, nRotamer, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
// eliot , taken from Matt Harvey
        for( int i = 0; i < nRotamer; i++ ) {
                E_Rotamer_MM[i] = E_Rotamer_MM_Local[i];
        }

	Setup_Weight_Rotamer();
}

void Setup_Weight_Rotamer(void)
{
	int i, nAtom;
	double E_Min_QM, E_Min_MM, kbT=2.5;

	nAtom = Mol_ESP.nAtom;
	E_Min_QM = E_Min_MM = 1.0E100;
	for(i=0; i<nRotamer; i++)	{
		if(E_Min_QM > E_Rotamer_QM[i])	{
			E_Min_QM = E_Rotamer_QM[i];
		}
		if(E_Min_MM > E_Rotamer_MM[i])	{
			E_Min_MM = E_Rotamer_MM[i];
		}
	}

	w_rotamer_Sum = 0.0;
	printf("nRotamer = %d\n", nRotamer );

	for(i=0; i<nRotamer; i++)	{
//		w_rotamer_List[i] = exp( -(E_Rotamer_QM[i]-E_Min_QM)/kbT) + exp( -(E_Rotamer_MM[i]-E_Min_MM)/kbT);
//		w_rotamer_List[i] = exp( -(E_Rotamer_QM[i]-E_Min_QM)/kbT) + min(exp( -min((E_Rotamer_MM[i]-E_Min_QM),0.0)/2.0), 12.0);
//		w_rotamer_List[i] = exp( -(E_Rotamer_QM[i]-E_Min_QM)/kbT) + min(exp( -(E_Rotamer_MM[i]-E_Min_QM)/2.0), 12.0);
		w_rotamer_List[i] = exp( -(E_Rotamer_QM[i]-E_Min_QM)/kbT) + min(exp( -(E_Rotamer_MM[i]-E_Min_MM - 3.0)/1.0), 8.0);
		w_rotamer_Sum += w_rotamer_List[i];
		printf("%d %f\n", i, w_rotamer_List[i] );
	}

}

void Setup_Weight_1D_Scan(void)
{
	int i, j, nAtom;
	double E_Min;

	w_Scan_Sum = 0.0;
	nAtom = Mol_ESP.nAtom;
	for(i=0; i<n_Phi; i++)	{
		E_Min = 1.0E100;
		for(j=0; j<nScan_List[i]; j++)	{
			if(E_Scan_QM[i][j] < E_Min)	{
				E_Min = E_Scan_QM[i][j];
			}
		}
		for(j=0; j<nScan_List[i]; j++)	{
			if( (E_Scan_QM[i][j] - E_Min) < 20.0)	{	// cutoff 20kcal/mol
				w_Scan_List[i][j] = 1.0;
				w_Scan_Sum += 1.0;
			}
			else	{
				w_Scan_List[i][j] = 0.0;
			}
		}
	}
}

