/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include <mpi.h> eliot
#include "ff.h"
#include "nlopt.h"

char szForceFiled[]="mol.prm";
char szXpsfFile[]="mol.xpsf";
char szCrdFile[]="mol-opt.xyz";
char szPhiToScan[]="soft-dih-list.txt";

#define N_MAX_DIH	(32)
#define N_STATE		(4)
#define MAX_ATOM_LOCAL	(256)
#define MAX_CONF_CHECK	(480)

CMol Mol;
CForceField ForceField;
int nDihedral, n_Conf;
double Ratio_Pick=1.0;

FILE *fFile_Run_Log, *fResult;
void Quit_With_Error_Msg(char szMsg[]);

int iseed = 7643;
int n_Phi=0, Counter=0;
int DihList[N_MAX_DIH][4], IdxDihSelect[N_MAX_DIH];
double Phi_Set[N_STATE]={-180.0, -90.0, 0.0, 90.0}, Phi_To_Set[N_MAX_DIH];
double x_Save[MAX_ATOM_LOCAL], y_Save[MAX_ATOM_LOCAL], z_Save[MAX_ATOM_LOCAL];


int	ProgID=0, nProc=1;


int Read_Soft_DihedralList(void);
int Get_Dihedral_Index(int ia, int ib, int ic, int id);
void Enumerate_Dihedrals(int IdxDih);
void BackupCoordinates(void);
void RestoreCoordinates(void);
double rand(int &iseed);

int Iter_Geo_Opt;
double func_Geo_Optimization(unsigned n, const double *x, double *grad, void *my_func_data);
double Do_Geometry_Optimization(CMol *pMol);
void Random_Pick_Configuration_Optimization(void);

void Enumerate_Dihedrals(int IdxDih)	// depth
{
	int iState, i, j;
	char szName[256], szIniName[256];
	FILE *fIn;
	double r;

	if(IdxDih < n_Phi)	{
		for(iState=0; iState<N_STATE; iState++)	{
			Phi_To_Set[IdxDih] = Phi_Set[iState];
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
//			printf("%.0lf ", Phi_To_Set[i]);
		}

		Do_Geometry_Optimization(&Mol);
//		Mol.FullGeometryOptimization_LBFGS();
		Mol.WriteXYZ(szName);

		fprintf(fResult, "%6d %.7E ", Counter, Mol.E_Total);
		for(j=0; j<n_Phi; j++)	{
			fprintf(fResult, "%.0lf ", Mol.QueryDihedral(IdxDihSelect[j]));
		}
		fprintf(fResult, "\n");


		return;
	}
}

void Random_Pick_Configuration_Optimization(void)
{
	int j, IdxPhi, IdxPhi_2, IdxConf, iState, Max_Conf, Counter=0;
	char szName[256];

	iseed += (1001*ProgID + 101*ProgID*ProgID);
	Max_Conf = (int)(1.0*MAX_CONF_CHECK/nProc);

	for(IdxPhi=0; IdxPhi<n_Phi; IdxPhi++)	{	// to select a dihedral
		for(iState=0; iState<N_STATE; iState++)	{	// total four states for selected dihedral
			for(IdxConf=0; IdxConf<Max_Conf; IdxConf++)	{
				RestoreCoordinates();

				for(IdxPhi_2=0; IdxPhi_2<n_Phi; IdxPhi_2++)	{
					Mol.QueryDihedral(IdxDihSelect[IdxPhi_2]);
					if(IdxPhi_2 != IdxPhi)	{
						Mol.Edit_Dihedral(IdxDihSelect[IdxPhi_2], Phi_Set[(int)(rand(iseed)*N_STATE)]);
					}
					else	{
						Mol.Edit_Dihedral(IdxDihSelect[IdxPhi], Phi_Set[iState]);	// not random
					}
				}
				
				
				Do_Geometry_Optimization(&Mol);
				Counter++;
				
				fprintf(fResult, "%6d %.7E ", Counter, Mol.E_Total);
				for(j=0; j<n_Phi; j++)	{
					fprintf(fResult, "%.0lf ", Mol.QueryDihedral(IdxDihSelect[j]));
				}
				fprintf(fResult, "\n");
				
				sprintf(szName, "./opt-pdb/opt-id=%d-%d.xyz", ProgID, Counter);
				Mol.WriteXYZ(szName);
			}
		}
	}

}

int main(int argc, char **argv)
{
  timebomb();

	char szName[256];

	fFile_Run_Log = fopen("mm-pes.log", "w");
/*	eliot
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProgID);
	MPI_Comm_size(MPI_COMM_WORLD, &nProc);
*/
	if(ProgID == 0)	{
		system("make_directory opt-pdb");
	}

	sprintf(szName, "mm-pes-id-%d.dat", ProgID);
	fResult = fopen(szName, "w");

	ForceField.ReadForceField(szForceFiled);
	Mol.ReadPSF(szXpsfFile, 0);
	Mol.AssignForceFieldParameters(&ForceField);
	Mol.ReadXYZ(szCrdFile);
	BackupCoordinates();

	Read_Soft_DihedralList();

//	Enumerate_Dihedrals(0);

	Random_Pick_Configuration_Optimization();

	fclose(fFile_Run_Log);
	fclose(fResult);
/*
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
*/
	return 0;
}


double x_Best[MAX_ATOM_LOCAL*3];
nlopt_opt opt_Geo;
double Do_Geometry_Optimization(CMol *pMol)
{
	int i, nAtom, nDim, iPos, Return_Opt, IterStep;
	double x[MAX_ATOM_LOCAL*3], x_Ini[MAX_ATOM_LOCAL], y_Ini[MAX_ATOM_LOCAL], z_Ini[MAX_ATOM_LOCAL], E_min=1.0E100;

	nAtom = pMol->nAtom;
	nDim = 3*nAtom;

	if(n_Phi <= 5)	{	// this step is not so expensive
		IterStep = pMol->FullGeometryOptimization_LBFGS_niter();
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
		
		IterStep = pMol->FullGeometryOptimization_LBFGS_niter();
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


void Quit_With_Error_Msg(char szMsg[])
{
	fprintf(fFile_Run_Log, "%s", szMsg);
	fflush(fFile_Run_Log);

	FILE *fOut;
	fOut = fopen("error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in mm_pes.cpp\n");
	fclose(fOut);


	exit(1);
}

int Read_Soft_DihedralList(void)
{
	FILE *fIn;
	int ReadItem;

	n_Phi = 0;

	fIn = fopen(szPhiToScan, "r");

	while(1)	{
		ReadItem = fscanf(fIn, "%d %d %d %d", &(DihList[n_Phi][0]), &(DihList[n_Phi][1]), &(DihList[n_Phi][2]), &(DihList[n_Phi][3]));
		if(ReadItem == 4)	{
			DihList[n_Phi][0]--;
			DihList[n_Phi][1]--;
			DihList[n_Phi][2]--;
			DihList[n_Phi][3]--;
			IdxDihSelect[n_Phi] = Mol.Query_Dihedral_Index(DihList[n_Phi][0], DihList[n_Phi][1], DihList[n_Phi][2], DihList[n_Phi][3]);
			if(IdxDihSelect[n_Phi] < 0)	{
				Quit_With_Error_Msg("Fail to identify the index of one soft dihedral.\n");
			}
			
			Mol.BuildSegmentList_Dihedrals(IdxDihSelect[n_Phi]);

			n_Phi++;
		}
		else	{
			break;
		}
	}

	fclose(fIn);

	n_Conf = pow(N_STATE*1.0, n_Phi*1.0);
	if(n_Conf >= MAX_CONF_CHECK)	{	// too many possible combinations, just randomly pick a part of them
		Ratio_Pick = 1.0 * MAX_CONF_CHECK / n_Conf;
	}
	else	{
		Ratio_Pick = 1.0;
	}

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


#define IADD   453806245
#define IMUL   314159269
#define MASK   2147483647
#define SCALE  0.4656612873e-9

double rand(int &iseed)
{
	iseed = (iseed * IMUL + IADD) & MASK;
	return (iseed * SCALE);
}
