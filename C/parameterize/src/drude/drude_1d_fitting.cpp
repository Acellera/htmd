/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include <mpi.h>

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

#ifndef MINAB
#define MINAB 
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif

int iseed = -127;	//the initial seed for random number generator

CMol Mol_ESP;
CForceField ForceField;


double netcharge_mol=0.0;


#define N_MAX_DIH	(16)

char szPhiToScan[]="soft-dih-list.txt";
int n_Phi=0;
int DihList[N_MAX_DIH][4], IdxDihSelect[N_MAX_DIH];
double Phi_To_Set[N_MAX_DIH];
int Read_Soft_DihedralList(void);
void To_Setup_All_Dihedral_Constraint(void);
void Quit_With_Error_Msg(char szMsg[]);

//start	to data and functions related with constrained optimization
double Geoometry_Optimization_With_Constraint(CMol *pMol);
double func_Geo_Opt_Constraint(unsigned n, const double *x, double *grad, void *my_func_data);
double Geo_Constraint_Bond(unsigned n, const double *x, double *grad, void *data);
double Geo_Constraint_Angle(unsigned n, const double *x, double *grad, void *data);
double Geo_Constraint_Dihedral(unsigned n, const double *x, double *grad, void *data);
//end	to data and functions related with constrained optimization

void Setup_Phi_Constraint_List(void);
//int Is_Soft_Dihedral(int ib, int ic);
int Is_Soft_Dihedral(int ia, int ib, int ic, int id);

//start	data related with torsion fitting
#define MAX_SCAN	(96)
#define HARTREE_To_KCAL	(627.509469)
#define N_TOR_PARA	(11)

int nScan;
double Phi_Scaned[MAX_SCAN], E_QM_Scaned[MAX_SCAN], E_MM_Scaned[MAX_SCAN], E_QM_MM_Diff[MAX_SCAN], E_QM_MM_Diff_Fit[MAX_SCAN], dE[MAX_SCAN], w_Point[MAX_SCAN];
double Para_Tor[N_TOR_PARA], Para_Best[N_TOR_PARA], Chi_SQ_Min;
int IsPara_Fixed[N_TOR_PARA], FailCount, Iteration;

int ActiveRun;	// the active run among 32 runs
double All_Para_Save[32][N_TOR_PARA], Chi_SQ_Min_Local[32], Chi_SQ_Min_Global;	// to save all results and optimal parameters for 32 runs


double Callback_Eval_Gradient(unsigned n, const double *x, double *grad, void *my_func_data);
void Cal_E_MM_QM_Diff(int Idx);
void Fit_Torsion_Parameters(int Idx);
double Optimize_Torsion_Parameters(void);
double CalGradient_Tor(double x[], double grad[]);
double CalObjectiveFunction_Tor(double x[]);

//end	data related with torsion fitting

int Gen_iseed(void);
double rand(int &iseed);

double Callback_Eval_Gradient(unsigned n, const double *x, double *grad, void *my_func_data)
{
	int n_Local;
	double Chi_SQ;

	n_Local = (int)n;

	if(grad)	{
		Chi_SQ = CalGradient_Tor((double*)x, grad);	//gradient will be assigned into objGrad automatically
		Iteration++;

		if(Chi_SQ < Chi_SQ_Min)	{
			Chi_SQ_Min = Chi_SQ;

			memcpy(Para_Best, x, sizeof(double)*n);

			FailCount = 0;
		}
		else	{
			FailCount++;
		}
		if(Chi_SQ < Chi_SQ_Min_Local[ActiveRun])	{
			Chi_SQ_Min_Local[ActiveRun] = Chi_SQ;
			memcpy(All_Para_Save[ActiveRun], x, sizeof(double)*n);
		}
	}
	else	{	// cal object function only
		Chi_SQ = CalObjectiveFunction_Tor((double*)x);
	}

    return Chi_SQ;
}


// including 1/2/3/4/6, E0 . Total 11 parameters
double CalObjectiveFunction_Tor(double x[])	// x - Para_Tor
{
	int i;
	double f, Phi;

	f = 0.0;
	for(i=0; i<nScan; i++)	{
		Phi = Phi_Scaned[i] * radianInv;
		E_QM_MM_Diff_Fit[i] = x[0] * (1.0 + cos(1.0*Phi - x[1]) ) 
			+ x[2] * (1.0 + cos(2.0*Phi - x[3]) ) 
			+ x[4] * (1.0 + cos(3.0*Phi - x[5]) ) 
			+ x[6] * (1.0 + cos(4.0*Phi - x[7]) ) 
			+ x[8] * (1.0 + cos(6.0*Phi - x[9]) ) 
			+ x[10];
		dE[i] = E_QM_MM_Diff_Fit[i] - E_QM_MM_Diff[i];
		f += (dE[i]*dE[i] * w_Point[i]);
	}

	return f;
}

double CalGradient_Tor(double x[], double grad[])	// x - Para_Tor
{
	int i;
	double f, Phi;

	f = CalObjectiveFunction_Tor(x);
	memset(grad, 0, sizeof(double)*N_TOR_PARA);

	for(i=0; i<nScan; i++)	{
		Phi = Phi_Scaned[i] * radianInv;

		if(!IsPara_Fixed[0])	{
			grad[0] += ( dE[i] * w_Point[i] * ( 1.0 + cos(1.0*Phi - x[1]) ) );
		}
		if(!IsPara_Fixed[2])	{
			grad[2] += ( dE[i] * w_Point[i] * ( 1.0 + cos(2.0*Phi - x[3]) ) );
		}
		if(!IsPara_Fixed[4])	{
			grad[4] += ( dE[i] * w_Point[i] * ( 1.0 + cos(3.0*Phi - x[5]) ) );
		}
		if(!IsPara_Fixed[6])	{
			grad[6] += ( dE[i] * w_Point[i] * ( 1.0 + cos(4.0*Phi - x[7]) ) );
		}
		if(!IsPara_Fixed[8])	{
			grad[8] += ( dE[i] * w_Point[i] * ( 1.0 + cos(6.0*Phi - x[9]) ) );
		}
		if(!IsPara_Fixed[10])	{
			grad[10] += (dE[i] * w_Point[i]);
		}

	}
	for(i=0; i<N_TOR_PARA; i++)	{
		grad[i] *= 2.0;
	}

	return f;
}

double Optimize_Torsion_Parameters(void)
{
	double lb[N_TOR_PARA] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -100.0 }; /* lower bounds */
	double ub[N_TOR_PARA] = { 12.0, 8.0, 12.0, 8.0, 12.0, 8.0, 12.0, 8.0, 12.0, 8.0,  100.0 }; /* lower bounds */
	double Phase[32][5]={{0, 0, 0, 0, 0}, {0, 0, 0, 0, PI}, {0, 0, 0, PI, 0}, {0, 0, 0, PI, PI}, {0, 0, PI, 0, 0}, {0, 0, PI, 0, PI}, {0, 0, PI, PI, 0}, {0, 0, PI, PI, PI}, {0, PI, 0, 0, 0}, {0, PI, 0, 0, PI}, 
						{0, PI, 0, PI, 0}, {0, PI, 0, PI, PI}, {0, PI, PI, 0, 0}, {0, PI, PI, 0, PI}, {0, PI, PI, PI, 0}, {0, PI, PI, PI, PI}, {PI, 0, 0, 0, 0}, {PI, 0, 0, 0, PI}, {PI, 0, 0, PI, 0}, {PI, 0, 0, PI, PI},
						{PI, 0, PI, 0, 0}, {PI, 0, PI, 0, PI}, {PI, 0, PI, PI, 0}, {PI, 0, PI, PI, PI}, {PI, PI, 0, 0, 0}, {PI, PI, 0, 0, PI}, {PI, PI, 0, PI, 0}, {PI, PI, 0, PI, PI}, {PI, PI, PI, 0, 0}, {PI, PI, PI, 0, PI}, 
						{PI, PI, PI, PI, 0}, {PI, PI, PI, PI, PI}};
	nlopt_opt opt;
	double Chi_SQ;
	int i, j;

	Chi_SQ_Min = 1.0E100;
	opt = nlopt_create(NLOPT_LD_LBFGS, N_TOR_PARA);

	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, ub);
	nlopt_set_min_objective(opt, Callback_Eval_Gradient, NULL);
	nlopt_set_xtol_rel(opt, 1e-8);

	iseed = time(NULL);
	
	for(i=0; i<32; i++)	{
		ActiveRun = i;
		Chi_SQ_Min_Local[i] = 1.0E200;
		FailCount = Iteration = 0;

//		for(j=0; j<5; j++)	{
//			Para_Tor[2*j] = 0.1+1.0*rand(iseed);
//		}
		
		Para_Tor[0] = Para_Tor[2] = Para_Tor[4] = Para_Tor[6] = Para_Tor[8] = 1.0;


		for(j=0; j<5; j++)	{
			Para_Tor[2*j+1] = Phase[i][j];
		}
//		Para_Tor[1] = Para_Tor[3] = Para_Tor[5] = Para_Tor[7] = Para_Tor[9] = 0.0;
		Para_Tor[10] = 0.0;
		
		IsPara_Fixed[1] = IsPara_Fixed[3] = IsPara_Fixed[5] = IsPara_Fixed[7] = IsPara_Fixed[9] = 1;

		//for test
//		IsPara_Fixed[6] = 1;	Para_Tor[6] = 0.0;
		
		if (nlopt_optimize(opt, Para_Tor, &Chi_SQ) < 0) {
			printf("nlopt failed!\n");
		}
		else {
			printf("Chi_SQ_min = %lf  Chi_SQ = %lf\n", Chi_SQ_Min, Chi_SQ);
		}
	}

	memcpy(Para_Tor, Para_Best, sizeof(double)*N_TOR_PARA);
	CalObjectiveFunction_Tor(Para_Tor);


	nlopt_destroy(opt);

	return Chi_SQ_Min;
}

void Fit_Torsion_Parameters(int Idx)
{
	FILE *fOut;
	char szName[256];
	int i, j, iTmp, i_Run;
	double E_QM_MM_Diff_Min=1.0E100;
/*
	sprintf(szName, "qm-mm-diff-%d.dat", Idx+1);
	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file: %s\nQuit\n", szName);
		Quit_With_Error_Msg(ErrorMsg);
	}

	nScan = 0;
	while(1)	{
		ReadItem = fscanf(fIn, "%lf %lf %lf", &(Phi_Scaned[nScan]), &(E_QM_Scaned[nScan]), &(E_QM_MM_Diff[nScan]));
		if(ReadItem == 3)	{
			if(E_QM_MM_Diff_Min > E_QM_MM_Diff[nScan])	{
				E_QM_MM_Diff_Min = E_QM_MM_Diff[nScan];
			}
			nScan++;
		}
		else	{
			break;
		}
	}
	fclose(fIn);
*/	
	for(i=0; i<nScan; i++)	{
		if( (E_QM_MM_Diff_Min > E_QM_MM_Diff[i]) && (w_Point[i] > 0.0) )	{
			E_QM_MM_Diff_Min = E_QM_MM_Diff[i];
		}
	}

	for(i=0; i<nScan; i++)	{
		E_QM_MM_Diff[i] -= E_QM_MM_Diff_Min;

		if(E_QM_Scaned[i] > 20.0)	{	// !!!!
			w_Point[i] = 0.0;;
		}
	}

	Optimize_Torsion_Parameters();

	sprintf(szName, "fitting-1d-%d.dat", Idx+1);
	fOut = fopen(szName, "w");
	for(i=0; i<nScan; i++)	{
		fprintf(fOut, "%.1lf %lf %lf %lf %lf\n", Phi_Scaned[i], E_QM_Scaned[i], E_QM_Scaned[i] - E_QM_MM_Diff[i] + E_QM_MM_Diff_Fit[i], E_QM_MM_Diff[i], E_QM_MM_Diff_Fit[i]);
//		fprintf(fOut, "%.1lf %lf %lf\n", Phi_Scaned[i], E_QM_MM_Diff[i], E_QM_MM_Diff_Fit[i]);
	}
	fclose(fOut);


	Chi_SQ_Min_Global = 1.0E200;
	for(i_Run=0; i_Run<32; i_Run++)	{
		if(Chi_SQ_Min_Global > Chi_SQ_Min_Local[i_Run])	{
			Chi_SQ_Min_Global = Chi_SQ_Min_Local[i_Run];
		}
	}

	int Sort_List[32];

	for(i=0; i<32; i++)	{
		Sort_List[i] = i;
	}
	for(i=0; i<32; i++)	{	// sort Chi_SQ_Min_Local[]
		for(j=i+1; j<32; j++)	{
			if(Chi_SQ_Min_Local[Sort_List[i]] > Chi_SQ_Min_Local[Sort_List[j]])	{
				iTmp = Sort_List[i];
				Sort_List[i] = Sort_List[j];
				Sort_List[j] = iTmp;
			}
		}
	}

	sprintf(szName, "torsion-para-%d.dat", Idx+1);
	fOut = fopen(szName, "w");
	for(i_Run=0; i_Run<32; i_Run++)	{
		if(fabs(Chi_SQ_Min_Global - Chi_SQ_Min_Local[Sort_List[i_Run]]) < 0.01)	{	// save all parameters giving comparable chi^2
			for(i=0; i<(N_TOR_PARA-1); i++)	{
				fprintf(fOut, "%.12lf ", All_Para_Save[Sort_List[i_Run]][i]);
			}
			fprintf(fOut, "\n");
		}
	}
	fclose(fOut);
}



void Cal_E_MM_QM_Diff(int Idx)
{
	FILE *fIn, *fOut;
	char szName[256], szLine[256], *ReadLine, szTmp[256], ErrorMsg[256];
	double E_QM_Min, E_MM_Min, E_MM_Read, E_MM_Cal;
	double Para_k_Dih_Save[6], Para_phi_Dih_Save[6];
	double x_Save[MAX_NUM_ATOM], y_Save[MAX_NUM_ATOM], z_Save[MAX_NUM_ATOM];
	int i, j, ReadItem, nAtom, nRealAtom, *pRealAtom_List, Idx_Atom;

	memcpy(Para_k_Dih_Save, Mol_ESP.Para_k_Dih[IdxDihSelect[Idx]]+1, sizeof(double)*6);
	memcpy(Para_phi_Dih_Save, Mol_ESP.Para_phi[IdxDihSelect[Idx]]+1, sizeof(double)*6);
	memset(Mol_ESP.Para_k_Dih[IdxDihSelect[Idx]]+1, 0, sizeof(double)*6);	// turn of this dihedral
	memset(Mol_ESP.Para_phi[IdxDihSelect[Idx]]+1, 0, sizeof(double)*6);

	sprintf(szName, "tor-1D-idx-%d.dat", Idx+1);
	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		sprintf(ErrorMsg, "Fail to open file: %s\nQuit\n", szName);
		printf("%s", ErrorMsg);
		Quit_With_Error_Msg(ErrorMsg);
	}

	sprintf(szName, "mm-tor-1D-idx-%d.dat", Idx+1);	// the one with optimized geometry in MM force field
	fOut = fopen(szName, "w");

	nRealAtom = Mol_ESP.nRealAtom;
	pRealAtom_List = Mol_ESP.RealAtom_List;
	nAtom = Mol_ESP.nAtom;
	nScan = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		if(strncmp(szLine, "E_Scan", 6)==0)	{
			ReadItem = sscanf(szLine, "%s %lf %lf %s %lf", szTmp, &(E_QM_Scaned[nScan]), &E_MM_Read, szTmp, &(Phi_Scaned[nScan]));
			if(ReadItem != 5)	{
				sprintf(ErrorMsg, "Error in reading file: %s\nQuit\n", szName);
				Quit_With_Error_Msg(ErrorMsg);
			}
			ReadLine = fgets(szLine, 256, fIn);	// skip one line, "Coordinate"

			for(i=0; i<nRealAtom; i++)	{	// read one snapshot
				ReadLine = fgets(szLine, 256, fIn);
				Idx_Atom = pRealAtom_List[i];
				ReadItem = sscanf(szLine, "%lf %lf %lf", &(Mol_ESP.x[Idx_Atom]), &(Mol_ESP.y[Idx_Atom]), &(Mol_ESP.z[Idx_Atom]));
				if(Mol_ESP.mass[Idx_Atom] > 1.5)	{	// a heavy atom. To set up the coordinate for the drude particle
					Mol_ESP.x[Idx_Atom+1] = Mol_ESP.x[Idx_Atom];
					Mol_ESP.y[Idx_Atom+1] = Mol_ESP.y[Idx_Atom];
					Mol_ESP.z[Idx_Atom+1] = Mol_ESP.z[Idx_Atom];
				}
			}
			Mol_ESP.Position_LonePair();	// to set up the coordinates of lone pairs
			E_MM_Cal = Mol_ESP.Cal_E(0);	// the nergy without optimization and dihedral term is turned off

//			Mol_ESP.FullGeometryOptimization_LBFGS();

			To_Setup_All_Dihedral_Constraint();
			E_MM_Scaned[nScan] = Geoometry_Optimization_With_Constraint(&Mol_ESP);

			if(E_MM_Scaned[nScan] > (1.0E3+E_MM_Cal))	{	// not a valid configuration. Something wrong in the constrained optimization
				continue;
			}

			memcpy(x_Save, Mol_ESP.x, sizeof(double)*nAtom);
			memcpy(y_Save, Mol_ESP.y, sizeof(double)*nAtom);
			memcpy(z_Save, Mol_ESP.z, sizeof(double)*nAtom);
			
			fprintf(fOut, "E_Scan %.13E %.13E Phi %.1lf\n", E_QM_Scaned[nScan], E_MM_Scaned[nScan], Phi_Scaned[nScan]);
			fprintf(fOut, "Coordinate\n");
			for(j=0; j<nAtom; j++)	{
				fprintf(fOut, "%14.6lf%14.6lf%14.6lf\n", x_Save[j], y_Save[j], z_Save[j]);
			}
			//		fprintf(fOut, "%.1lf %.6lf %.6lf %.6lf\n", Phi_Scaned[i], E_QM_Scaned[i], E_MM_Scaned[i], E_QM_Scaned[i] - E_MM_Scaned[i]);
			

			nScan++;
		}

	}
	fclose(fOut);

	fclose(fIn);

	E_QM_Min = E_MM_Min = 1.0E100;
	for(i=0; i<nScan; i++)	{
		if(E_QM_Min > E_QM_Scaned[i])	{
			E_QM_Min = E_QM_Scaned[i];
		}
		if(E_MM_Min > E_MM_Scaned[i])	{
			E_MM_Min = E_MM_Scaned[i];
		}
	}

	for(i=0; i<nScan; i++)	{
		E_QM_Scaned[i] = (E_QM_Scaned[i]-E_QM_Min)*HARTREE_To_KCAL;
		E_MM_Scaned[i] -= E_MM_Min;
		E_QM_MM_Diff[i] = E_QM_Scaned[i]-E_MM_Scaned[i];

		if(E_MM_Scaned[i] > 60.0)	{	// an unusual large MM energy
			w_Point[i] = 0.0;
		}
		else	{
			w_Point[i] = 1.0;
		}
	}


/*
	sprintf(szName, "qm-mm-diff-%d.dat", Idx+1);
	fOut = fopen(szName, "w");
	for(i=0; i<nScan; i++)	{
		E_QM_Scaned[i] = (E_QM_Scaned[i]-E_QM_Min)*HARTREE_To_KCAL;
		E_MM_Scaned[i] -= E_MM_Min;
		fprintf(fOut, "%.1lf %.6lf %.6lf\n", Phi_Scaned[i], E_QM_Scaned[i], E_QM_Scaned[i]-E_MM_Scaned[i]);
//		fprintf(fOut, "%.1lf %.12e %.6lf %.6lf\n", Phi_Scaned[i], E_QM_Scaned[i]);
//		fprintf(fOut, "%.1lf %.6lf %.6lf %.6lf\n", Phi_Scaned[i], E_QM_Scaned[i], E_MM_Scaned[i], E_QM_Scaned[i] - E_MM_Scaned[i]);
	}
	fclose(fOut);
*/


	memcpy(Mol_ESP.Para_k_Dih[IdxDihSelect[Idx]]+1, Para_k_Dih_Save, sizeof(double)*6);
	memcpy(Mol_ESP.Para_phi[IdxDihSelect[Idx]]+1, Para_phi_Dih_Save, sizeof(double)*6);
}

int Iter_Opt_Constrained;
double func_Geo_Opt_Constraint(unsigned n, const double *x, double *grad, void *my_func_data)
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
		Iter_Opt_Constrained ++;
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



double Geoometry_Optimization_With_Constraint(CMol *pMol)
{
	nlopt_opt opt, opt_AugLag;
	int Return_Opt, nAtom, nDim, i, iPos;
	double E_min, x[MAX_NUM_ATOM*3];

	Mol_ESP.Restrain_All_Torsions(1);	// apply soft restraints on all soft dihedrals

	nAtom = pMol->nAtom;
	nDim = 3 * nAtom;

	printf("\nPhi before MM optimization:\n");
	for(i=0; i<n_Phi; i++)	{
		printf("%3d %3d %3d %3d %8.4lf\n", DihList[i][0], DihList[i][1], DihList[i][2], DihList[i][3], 
			Mol_ESP.Query_Dihedral(DihList[i][0], DihList[i][1], DihList[i][2], DihList[i][3], 0, NULL));
	}
// ELIOT
//	AUG_LAG_ICM_TOL = 1.0E-17;
	opt_AugLag = nlopt_create(NLOPT_AUGLAG_EQ, nDim);
	nlopt_set_min_objective(opt_AugLag, func_Geo_Opt_Constraint, pMol);

	for(i=0; i<pMol->nBond_Fixed; i++)	{
		pMol->Geo_Fix_r0[i].pMol = pMol;
		nlopt_add_equality_constraint(opt_AugLag, Geo_Constraint_Bond, &(pMol->Geo_Fix_r0[i]), 1.0E-6);
	}

	for(i=0; i<pMol->nAngle_Fixed; i++)	{
		pMol->Geo_Fix_theta0[i].pMol = pMol;
		nlopt_add_equality_constraint(opt_AugLag, Geo_Constraint_Angle, &(pMol->Geo_Fix_theta0[i]), 1.0E-6);
	}

	for(i=0; i<pMol->nPhi_Fixed; i++)	{
		pMol->Geo_Fix_phi0[i].pMol = pMol;
		nlopt_add_equality_constraint(opt_AugLag, Geo_Constraint_Dihedral, &(pMol->Geo_Fix_phi0[i]), 1.0E-6);
	}

	nlopt_set_xtol_rel(opt_AugLag, 2E-8);
	
	opt = nlopt_create(NLOPT_LD_LBFGS_GEO, nDim);
	nlopt_set_min_objective(opt, func_Geo_Opt_Constraint, pMol);
	nlopt_set_xtol_rel(opt, 1E-8);

	nlopt_set_local_optimizer(opt_AugLag, opt);
	nlopt_destroy(opt);

	for(i=0; i<nAtom; i++)	{
		iPos = 3*i;
		x[iPos  ] = pMol->x[i];
		x[iPos+1] = pMol->y[i];
		x[iPos+2] = pMol->z[i];
	}
	
	Iter_Opt_Constrained = 0;
	Return_Opt = nlopt_optimize(opt_AugLag, x, &E_min);
	if (Return_Opt < 0) {
		printf("nlopt failed!\n");
	}
	else {
		printf("After constrained optimization, E = %12.6lf\n", E_min);
	}

	nlopt_destroy(opt_AugLag);


	printf("Phi after MM optimization:\n");
	for(i=0; i<n_Phi; i++)	{
		printf("%3d %3d %3d %3d %8.4lf\n", DihList[i][0], DihList[i][1], DihList[i][2], DihList[i][3], 
			Mol_ESP.Query_Dihedral(DihList[i][0], DihList[i][1], DihList[i][2], DihList[i][3], 0, NULL));
	}

	Mol_ESP.Restrain_All_Torsions(0);	// lift soft restraints on all soft dihedrals

	return (pMol->Cal_E(0));
}

double Geo_Constraint_Bond(unsigned n, const double *x, double *grad, void *data)
{
	int ia, ib, iPos;
	double r, r0, g[2][3];
	CMol* pMol;
	GEO_FIX_BOND *Geo_Fix_r;

	Geo_Fix_r = (GEO_FIX_BOND*)data;
	pMol = (CMol*)(Geo_Fix_r->pMol);
	ia = Geo_Fix_r->ia;
	ib = Geo_Fix_r->ib;
	r0 = Geo_Fix_r->r0;

	if(grad)	{
		r = pMol->Query_Distance(ia, ib, 1, g);
		memset(grad, 0, sizeof(double)*n);

		iPos = 3 * ia;
		grad[iPos  ] = g[0][0];
		grad[iPos+1] = g[0][1];
		grad[iPos+2] = g[0][2];

		iPos = 3 * ib;
		grad[iPos  ] = g[1][0];
		grad[iPos+1] = g[1][1];
		grad[iPos+2] = g[1][2];
	}
	else	{
		r = pMol->Query_Distance(ia, ib, 0, g);
	}

	return (r - r0);
}

double Geo_Constraint_Angle(unsigned n, const double *x, double *grad, void *data)
{
	int ia, ib, ic, iPos;
	double theta, theta0, g[3][3];
	CMol* pMol;
	GEO_FIX_ANGLE *Geo_Fix_theta;

	Geo_Fix_theta = (GEO_FIX_ANGLE*)data;
	pMol = (CMol*)(Geo_Fix_theta->pMol);
	ia = Geo_Fix_theta->ia;
	ib = Geo_Fix_theta->ib;
	ic = Geo_Fix_theta->ic;
	theta0 = Geo_Fix_theta->theta0;

	if(grad)	{
		theta = pMol->Query_Angle(ia, ib, ic, 1, g);
		memset(grad, 0, sizeof(double)*n);

		iPos = 3 * ia;
		grad[iPos  ] = g[0][0];
		grad[iPos+1] = g[0][1];
		grad[iPos+2] = g[0][2];

		iPos = 3 * ib;
		grad[iPos  ] = g[1][0];
		grad[iPos+1] = g[1][1];
		grad[iPos+2] = g[1][2];

		iPos = 3 * ic;
		grad[iPos  ] = g[2][0];
		grad[iPos+1] = g[2][1];
		grad[iPos+2] = g[2][2];
	}
	else	{
		theta = pMol->Query_Angle(ia, ib, ic, 0, g);
	}

	return (theta - theta0);
}

double Geo_Constraint_Dihedral(unsigned n, const double *x, double *grad, void *data)
{
	int ia, ib, ic, id, iPos;
	double phi, phi0, d_phi, g[4][3];
	CMol* pMol;
	GEO_FIX_DIHEDRAL *Geo_Fix_phi;

	Geo_Fix_phi = (GEO_FIX_DIHEDRAL*)data;
	pMol = (CMol*)(Geo_Fix_phi->pMol);
	ia = Geo_Fix_phi->ia;
	ib = Geo_Fix_phi->ib;
	ic = Geo_Fix_phi->ic;
	id = Geo_Fix_phi->id;
	phi0 = Geo_Fix_phi->phi0;

	if(grad)	{
		phi = pMol->Query_Dihedral(ia, ib, ic, id, 1, g);
		memset(grad, 0, sizeof(double)*n);

		iPos = 3 * ia;
		grad[iPos  ] = g[0][0];
		grad[iPos+1] = g[0][1];
		grad[iPos+2] = g[0][2];

		iPos = 3 * ib;
		grad[iPos  ] = g[1][0];
		grad[iPos+1] = g[1][1];
		grad[iPos+2] = g[1][2];

		iPos = 3 * ic;
		grad[iPos  ] = g[2][0];
		grad[iPos+1] = g[2][1];
		grad[iPos+2] = g[2][2];

		iPos = 3 * id;
		grad[iPos  ] = g[3][0];
		grad[iPos+1] = g[3][1];
		grad[iPos+2] = g[3][2];
	}
	else	{
		phi = pMol->Query_Dihedral(ia, ib, ic, id, 0, g);
	}

	d_phi = phi - phi0;
	if(d_phi < -PI)	{
		d_phi += PI2;
	}
	else if(d_phi > PI)	{
		d_phi -= PI2;
	}

	return d_phi;
}


FILE *fFile_Conf;
FILE *fFile_Run_Log;	// will be shared by other source code
char szName_Run_Log[]="1d-fitting-log.txt";
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

int main(int argc, char **argv)
{
	int i;

  timebomb();

	fFile_Run_Log = fopen(szName_Run_Log, "w");
	if(fFile_Run_Log==NULL)	{
		printf("Fail to create the log file.\n\n");
		exit(1);
	}

	ForceField.ReadForceField("drude-mol.prm");
	
	Mol_ESP.drude_ReadPSF("drude-mol.xpsf", 0, 1);	// delete the test charge as the last atom in xpsf
//	Read_Rtf_File();
	
	strcpy(Mol_ESP.MolName, Mol_ESP.ResName[0]);
	Mol_ESP.AssignForceFieldParameters(&ForceField);

	Read_Soft_DihedralList();

	Mol_ESP.Is_Phi_Psi_Constrained = 0;
	Mol_ESP.E_CMap_On = 0;

	for(i=0; i<n_Phi; i++)	{
		Cal_E_MM_QM_Diff(i);
		Fit_Torsion_Parameters(i);
	}

	fflush(fFile_Run_Log);
	

	fclose(fFile_Run_Log);
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

#define IADD   453806245
#define IMUL   314159269
#define MASK   2147483647
#define SCALE  0.4656612873e-9

double rand(int &iseed)
{
	iseed = (iseed * IMUL + IADD) & MASK;
	return (iseed * SCALE);
}
#undef SCALE
#undef MASK
#undef IMUL
#undef IADD


int Gen_iseed(void)
{
	int i, j, iseed, nLen, Sum, LenUse=7;
	char szDir[1024], szString[256], sz_iseed[256];


//	_getcwd(szDir, 1024);
	getcwd(szDir, 1024);

	nLen = strlen(szDir);

	Sum = 0;
	for(i=0; i<nLen; i++)	{
		Sum += (int)szDir[i];
	}

	sprintf(szString, "%18.13e", log(1.0*Sum));

	i = strlen(szString)-1;
	while(1)	{
		if( (szString[i]=='E') || (szString[i]=='e') )	{
			break;
		}
		i--;
	}
	i-=2;

	for(j=0; j<LenUse; j++)	{
		sz_iseed[j] = szString[i-j];
	}
	sz_iseed[j] = 0;
	iseed = atoi(sz_iseed);
//	printf("%s\n", sz_iseed);
	fprintf(fFile_Run_Log, "%s\n", sz_iseed);
	fflush(fFile_Run_Log);

	return iseed;
}

void GetSecondItem(char szLine[], char szSecond[])
{
	char szTmp[256];
	int ReadItem;

	ReadItem = sscanf(szLine, "%s %s", szTmp, szSecond);
	if(ReadItem != 2)	{
		fprintf(fFile_Run_Log, "Fail to extract the second string from, \n%s\nQuit\n", szLine);
		fflush(fFile_Run_Log);
		exit(1);
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
				fprintf(fFile_Run_Log, "Error in ReadConfFile(). n_Line_Conf_File > MAX_LINE_CONF\nQuit\n");
				fflush(fFile_Run_Log);
				exit(1);
			}
		}
		if(LineLen >= MAX_LINE_CONF_LEN)	{
			fprintf(fFile_Run_Log, "Error in ReadConfFile(). LineLen >= MAX_LINE_CONF_LEN\nQuit\n");
			fflush(fFile_Run_Log);
			exit(1);
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

		ReadItem = sscanf(szLine, "%d %d %d %d", 
			&(DihList[n_Phi][0]), &(DihList[n_Phi][1]), &(DihList[n_Phi][2]), &(DihList[n_Phi][3]));

		if(ReadItem == 4)	{
			DihList[n_Phi][0] = Mol_ESP.RealAtom_List[DihList[n_Phi][0]-1];
			DihList[n_Phi][1] = Mol_ESP.RealAtom_List[DihList[n_Phi][1]-1];
			DihList[n_Phi][2] = Mol_ESP.RealAtom_List[DihList[n_Phi][2]-1];
			DihList[n_Phi][3] = Mol_ESP.RealAtom_List[DihList[n_Phi][3]-1];
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
	fprintf(fOut, "Error in 1D-fitting.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}

