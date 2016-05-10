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

#ifndef PI
#define PI	(3.14159265358979323846)
#define PI2	(6.28318530717958647692)
#define radian	(57.29577951308232088)
#define radianInv	(0.017453292519943295)
#endif

#define MAX_ATOM_TYPE	(256)	//This may be a little small for ala10
#define MAX_NUM_EQUIV	(32)	// the maximum of number of atoms sharing the same type
#define MAX_N_ATOM_NEUTRAL	(256)	//the total number of atoms used in neutral charge constrain
#define MAX_NUM_PARA	(256)
#define MAX_PARA_SAVE	(256)
#define MAX_MOL			(32)

#define N_CONF_DIMER	(12)

//start	to define the lower and upper boundaries of parameters
#define ATM_TYPE_HEAVY	(1)
#define ATM_TYPE_DRUDE	(2)
#define ATM_TYPE_H		(3)
#define ATM_TYPE_LP		(4)

#define LOW_BOUND_CG_HEAVY_ATM	(0.1)
#define UP_BOUND_CG_HEAVY_ATM	(4.0)
#define LOW_BOUND_CG_DRUDE	(-3.6)
#define UP_BOUND_CG_DRUDE	(-0.2)
//#define LOW_BOUND_CG_H	(0.04)
//#define LOW_BOUND_CG_H	(0.01)
//#define LOW_BOUND_CG_H	(0.033) // Changed Dec. 5th, 2011
//#define LOW_BOUND_CG_H	(0.043) // Changed Dec. 29th, 2011
//#define LOW_BOUND_CG_H	(0.03) // Changed Jan. 1st, 2012
#define LOW_BOUND_CG_H	(0.00)
#define UP_BOUND_CG_H	(1.00)
#define LOW_BOUND_CG_LP	(-0.70)

//#define UP_BOUND_CG_LP	(-0.10)
#define UP_BOUND_CG_LP	(-0.05)	// Changed Dec. 5th, 2011

//#define LOW_BOUND_THOLE	(0.25)
//#define UP_BOUND_THOLE	(2.00)

#define LOW_BOUND_THOLE	(0.05)
#define UP_BOUND_THOLE	(1.40)

#define LOW_BOUND_ANISO	(0.30)
#define UP_BOUND_ANISO	(1.50)

#define LOW_BOUND_LJ_E	(-2.00)
#define UP_BOUND_LJ_E	(-5.0E-3)
#define LOW_BOUND_LJ_r	(0.00)
#define UP_BOUND_LJ_r	(5.5)

#define LOW_BOUND_NBFIX_E	(-2.00)
#define UP_BOUND_NBFIX_E	(-5.0E-3)
#define LOW_BOUND_NBFIX_r	(1.60)
#define UP_BOUND_NBFIX_r	(5.5)

//#define LOW_BOUND_LP_r		(0.12)
#define LOW_BOUND_LP_r		(0.1)
#define UP_BOUND_LP_r		(1.20)

#define LOW_BOUND_LP_r_NEG		(-1.20)
//#define UP_BOUND_LP_r_NEG		(-0.12)
#define UP_BOUND_LP_r_NEG		(-0.01)

//#define LOW_BOUND_LP_THETA		(1.0)
#define LOW_BOUND_LP_THETA		(0.1)
#define UP_BOUND_LP_THETA		(PI)

#define D_CG_HEAVY_ATOM	(0.05)
#define D_CG_H_ATOM	(0.02)
#define D_CG_LP	(0.05)
#define D_THOLE	(0.10)
#define D_ANISO	(0.06)


//end	to define the lower and upper boundaries of parameters

#ifndef MINAB
#define MINAB
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif

//start	data related with reading in fitted torsion parameters
#define N_MAX_DIH	(32)
int n_Phi=0;
int DihList[N_MAX_DIH][4], IdxDihSelect[N_MAX_DIH];
char szPhiToScan[]="soft-dih-list.txt";
int Read_Soft_DihedralList(void);
void Read_Fitted_Torsion_Parameters(void);
//end	data related with reading in fitted torsion parameters



int iseed = 627;	//the initial seed for random number generator

static double Chi_SQ_Min=1.0E20;

//start	data used for picking stable parameters
int Is_Saved_Para[MAX_PARA_SAVE], iActiveSave=0;	// the flag for saved parameters, the pointer to the next available space
double Para_Set_Save[MAX_PARA_SAVE][MAX_NUM_PARA];
double Para_Set_Saved_Chi2[MAX_PARA_SAVE];
//end	data used for picking stable parameters

double netcharge_mol=0.0;
int SameAtom_Count;

//start	to data related with MPI and parallelization
#define MY_MPI_ROOT	(0)
#define MAX_N_PROC	(64)
int	ProgID=0, nProc=1;
//end	to data related with MPI and parallelization

FILE *fFile_Conf;
FILE *fFile_Run_Log;	// will be shared by other source code
char szName_Run_Log[]="fitting-log.txt";
char szName_Conf_File[256];
char szName_XPSF[256];
//char szName_CRD[256];
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
#define MAX_ANISOTROPY	(256)
char szName_Old_Rtf[]="drude-mol.rtf";
char szName_New_Rtf[]="drude-esp-mol.rtf";
char szTxt_Rtf[MAX_LINE_RTF][MAX_LEN_RTF];

int To_Export_New_Rtf, nLine_Rtf, Idx_Res_Mol, Idx_Mol_Atom_Start, Idx_Mol_Atom_End, nAnisotropy, Aniso_List[MAX_ANISOTROPY];
void Read_Rtf_File(void);
void Update_Rtf_File(void);
//end	to data and subroutines related with update the rtf file

int To_Use_Miller_Alpha=0;
double w_H_Donor_Acceptor=1.0, w_charge=2.0, w_alpha=0.5, w_thole=0.12, w_anisotropy=0.1, w_water_E_min=0.5, w_water_R_min=8.0, Scale_Alpha=1.0, w_water_int_E=1.0;
#define SCALE_DOWN_ALPHA	(0.90)
double Scale_TargetAlpha=1.0, AutoScale_Alpha=1.0;
//double w_H_Donor_Acceptor=1.0, w_charge=1.0, w_alpha=1.0, w_thole=1.0, w_water_int_E=1.0, Scale_Alpha=1.0;

double rmsfit(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);
void quatfit(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n);
void jacobi(int n, int np, double a[5][5], double d[5], double v[5][5], double b[5], double z[5]);
double CalRMSD_Drude(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int n, CMol *pMol);

double Callback_Eval_Gradient(unsigned n, const double *x, double *grad, void *my_func_data);
int Gen_iseed(void);
double rand(int &iseed);

void Quit_With_Error_Msg(char szMsg[]);
int FindTags_As_First_String(char szBuff[], char szTag[]);

void Get_Netcharge_From_Xpsf(void);

int SplitString(char szBuff[], char ItemList[][32]);
void Para_Aniso_K_A(double K11, double K22, double K33, double& A11, double& A22);
void Para_Aniso_A_K(double& K11, double& K22, double& K33, double A11, double A22);

int Not_Fit_Elec_Para=0;
int To_Fit_Aniso=0;	// the flag for fitting anisotropy parameters or not
int To_Fit_LP_Para=0;	// the flag for fitting lone pair parameters or not. distance and angle
int To_Fit_LJ=0;	// the flag for fitting LJ parameters or not
int Target_E_Int_Dimer=0;	// the flag for fitting QM energies for dimer interactions or not
int Target_E_Int_Water=0;	// the flag for fitting QM energies for compound-water interactions or not

double QM_Dipole_x=0.0, QM_Dipole_y=0.0, QM_Dipole_z=0.0, QM_Dipole=0.0;
double MM_Dipole_x, MM_Dipole_y, MM_Dipole_z, MM_Dipole;


//start	to data for ESP fitting
//#define NEIGHBOR_CUT	(3.6)
#define NEIGHBOR_CUT	(3.0)
#define N_MAX_GRID_ESP	(5000)	// the maximum number of grid points with ESP
#define N_MAX_CONF		(500)	// the maximum number of configurations with various positions of test charge
int n_Conf=0, n_Grid=0;
double ESP_Grid_x[N_MAX_GRID_ESP], ESP_Grid_y[N_MAX_GRID_ESP], ESP_Grid_z[N_MAX_GRID_ESP];
double ESP_Coord_x[N_MAX_CONF][MAX_NUM_PARA], ESP_Coord_y[N_MAX_CONF][MAX_NUM_PARA], ESP_Coord_z[N_MAX_CONF][MAX_NUM_PARA];	// coordinates of all perturbed configurations
double ESP_Coord_Opt_x[N_MAX_CONF][MAX_NUM_PARA], ESP_Coord_Opt_y[N_MAX_CONF][MAX_NUM_PARA], ESP_Coord_Opt_z[N_MAX_CONF][MAX_NUM_PARA];		// pre-optimized coordinates of all perturbed configurations
double ESP_Coord_Neutral_Opt_x[MAX_NUM_PARA], ESP_Coord_Neutral_Opt_y[MAX_NUM_PARA], ESP_Coord_Neutral_Opt_z[MAX_NUM_PARA];	// pre-optimized coordinates used for unperturbed ESP
double ESP_Grid_QM[N_MAX_CONF][N_MAX_GRID_ESP], ESP_Grid_QM_0[N_MAX_GRID_ESP];
double ESP_Grid_MM[N_MAX_CONF][N_MAX_GRID_ESP], ESP_Grid_MM_0[N_MAX_GRID_ESP];

double dist_All[N_MAX_GRID_ESP][MAX_ATOM_TYPE];
int N_HEAVY_ATM;
char szFile_pot0[1024]="./QM_ESP/nma.pot0";	// the unperturbed ESP and coordinates for grid points
char szFile_qpos[1024]="./QM_ESP/nma.qpos.all";		// the coodinates of all configurations
char szFile_pot[1024]="./QM_ESP/nma.pot";		// the perturbed ESP
//end	to data for ESP fitting

double w_ESP_Grid[N_MAX_GRID_ESP], w_ESP_Grid_Sum;
int H_Donor_Acceptor[MAX_ATOM_TYPE];

nlopt_opt opt;


//start	to data used for fitting LJ parameters by force matching
#define MAX_N_PERT_MOL_AR	(800)
int nConf_Mol_Ar=0, i_Ref_Conf;
int Used_Ar_Pert[MAX_N_PERT_MOL_AR];
double f_Pert_x_QM[MAX_N_PERT_MOL_AR][MAX_NUM_PARA];
double f_Pert_y_QM[MAX_N_PERT_MOL_AR][MAX_NUM_PARA];
double f_Pert_z_QM[MAX_N_PERT_MOL_AR][MAX_NUM_PARA];
double Ar_Pert_x[MAX_N_PERT_MOL_AR][MAX_NUM_PARA];
double Ar_Pert_y[MAX_N_PERT_MOL_AR][MAX_NUM_PARA];
double Ar_Pert_z[MAX_N_PERT_MOL_AR][MAX_NUM_PARA];
double E_Ar_Pert_QM[MAX_N_PERT_MOL_AR], E_Ar_Pert_MM[MAX_N_PERT_MOL_AR];
char szFile_MolArForce[256]="", szFile_MolArEnergy[256]="", szFile_MolAr_PSF[256]="";
//end	to data used for fitting LJ parameters by force matching

//start	to data used for fitting LJ parameters by force matching
#define MAX_N_PERT_MOL_WAT	(500)
int nConf_Mol_Wat=0;
int Used_Mol_Wat_Int[MAX_N_PERT_MOL_WAT], nUsed_Conf_Mol_Wat_Int;
double Coord_Dimer_Wat_x[MAX_N_PERT_MOL_WAT][MAX_NUM_PARA];
double Coord_Dimer_Wat_y[MAX_N_PERT_MOL_WAT][MAX_NUM_PARA];
double Coord_Dimer_Wat_z[MAX_N_PERT_MOL_WAT][MAX_NUM_PARA];
double E_Int_Wat_QM[MAX_N_PERT_MOL_WAT], E_Int_Wat_MM[MAX_N_PERT_MOL_WAT];
char szFile_MolWaterEnergy[256]="", szFile_MolWater_PSF[256]="";
//end	to data used for fitting LJ parameters by force matching

double FAR_DRUDE_R=0.10;	// the cutoff used to detect strange configurations to decide whether down-scaled polarizabilities are needed or not


//start	to data used for fitting E_Min and R_Min
#define N_PARA_WAT	(2)
#define R_DUMMY_1	(1.5)
#define R_DUMMY_2	(1.5)
#define R_DUMMY_3	(5.0)

#define MAX_DONOR	(100)
#define MAX_ACCEPTOR	(100)

#define b0_H_O_EXP	(0.9572)	// tip3 and exp
#define theta0_H_O_H_EXP	(104.52)	// tip3, exp?
#define Lower_Bound_R_H_Bond	(1.50)
#define Upper_Bound_R_H_Bond	(2.8)

double SCALE_QM_E_MIN=1.0;
double SHIFT_QM_R_MIN=0.0;
double SHIFT_QM_R_MIN_CHARGED=0.0;


int N_Donor, N_Acceptor;
int Donor_List_Active[MAX_DONOR];	// the donor list with QM data
int Acceptor_List_Acceptor[MAX_ACCEPTOR];
int theta_Acceptor[MAX_ACCEPTOR], theta_Donor[MAX_DONOR];
int theta_Acceptor_Save[MAX_ACCEPTOR], theta_Donor_Save[MAX_DONOR];
double E_Min_QM_Donor[MAX_DONOR], E_Min_QM_Acceptor[MAX_ACCEPTOR];
double r_Min_QM_Donor[MAX_DONOR], r_Min_QM_Acceptor[MAX_ACCEPTOR];
double E_Min_MM_Donor[MAX_DONOR], E_Min_MM_Acceptor[MAX_ACCEPTOR];
double r_Min_MM_Donor[MAX_DONOR], r_Min_MM_Acceptor[MAX_ACCEPTOR];
double Para_Save_Donor[MAX_DONOR][N_PARA_WAT], Para_Save_Donor_Org[MAX_DONOR][N_PARA_WAT];
double x_Wat_Donor[MAX_DONOR][3], y_Wat_Donor[MAX_DONOR][3], z_Wat_Donor[MAX_DONOR][3];
double Para_Save_Acceptor[MAX_DONOR][N_PARA_WAT], Para_Save_Acceptor_Org[MAX_DONOR][N_PARA_WAT];
double x_Wat_Acceptor[MAX_DONOR][3], y_Wat_Acceptor[MAX_DONOR][3], z_Wat_Acceptor[MAX_DONOR][3];
int ToFit_Donor[MAX_DONOR], ToFit_Acceptor[MAX_DONOR], nFit_Donor, nFit_Acceptor;
double x_Mol_Full_Opt[MAX_NUM_PARA], y_Mol_Full_Opt[MAX_NUM_PARA], z_Mol_Full_Opt[MAX_NUM_PARA];
double Normalize_drude(double& vect_x, double& vect_y, double& vect_z);
double Dot_Product(double x_a, double y_a, double z_a, double x_b, double y_b, double z_b);
void Cross_Product(double x_a, double y_a, double z_a, double x_b, double y_b, double z_b, double& x_Product, double& y_Product, double& z_Product);
//end	to data used for fitting E_Min and R_Min


#define N_ATM_DIMER		(38)
double x_Dimer[N_CONF_DIMER][N_ATM_DIMER], y_Dimer[N_CONF_DIMER][N_ATM_DIMER], z_Dimer[N_CONF_DIMER][N_ATM_DIMER];
double E_Int_QM_Dimer[N_CONF_DIMER], E_Int_MM_Dimer[N_CONF_DIMER];


typedef struct	{
	CMol* pMol;	//the pointer to the molecule
	int AtomIdx;	//the index of atom in this molecule which has specific type
}ATOMPARA;

typedef struct	{
	CMol* pMol;
	char MolName[8];	//Each molecule should have a unique name
	char ChemName[16];	//Chemical name/type, used to query LJ parameter
	int AtomIdx;	//the index of atom in this molecule
	int Atom_Count;	//the total number of atoms which have this atom type
	int CG_Represent;	//the pointer to the representative atom type with same charge
	int Thole_Represent;	//the pointer to the representative atom type with same thole
	int IsRedundant;
	int NotFitCharge;	// if true, the charge will be calculated from user defined sum
	int NotFitThole;	// if true, the charge will be calculated from user defined sum
	int NotFitLJ;
	int TypeID;
	
	int IsDrudeHost;	//the flag for host of drude / heavy atoms
	double CG, alpha, thole;	//only heavy atoms have alpha and thole
//	double LJ_E_Min, LJ_r_Min;	// LJ parameter
	double a_Miller;	// the alpha values from Miller
	double mass;
	
	int n_Atom_Neutral;	//the lnegth of atoms used to calculate the charge of this atom
	double weight[MAX_N_ATOM_NEUTRAL];
	int Cal_CG_List[MAX_N_ATOM_NEUTRAL];	//the list of atoms used to calculate the charge of this atom
}ATOMTYPE;


int Mol_Active_Atom;

class COptimizeWaterDimer
{
public:
	CMol *pDimer;
	double Low_Bound[N_PARA_WAT], Up_Bound[N_PARA_WAT], Para_Best[N_PARA_WAT], Para_List[N_PARA_WAT], Grad_List[N_PARA_WAT];
	double E_Dimer_Int, E_Min;

	double *x, *y, *z;
	double x_Dummy_1, y_Dummy_1, z_Dummy_1;
	double x_Dummy_2, y_Dummy_2, z_Dummy_2;
	double x_Dummy_3, y_Dummy_3, z_Dummy_3;
	double v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor;
	double v_x_X1_H, v_y_X1_H, v_z_X1_H;
	int ActiveAtom, ActiveHAtom, IsAcceptor, nAtom;

	int Bond_Count[MAX_ATOM_TYPE], BondList[MAX_ATOM_TYPE][4], Bonded_Atom[MAX_ATOM_TYPE], Is_H[MAX_ATOM_TYPE];

	void Init(CMol* pMol);
	void GetBondList(int Idx);
	void Generate_Dummy_Atoms(int Idx, int bAcceptor, int theta);	// generate three dummy atoms
	void Generate_Dummy_Atoms_Acceptor_1(void);
	void Generate_Dummy_Atoms_Acceptor_2(void);
	void Generate_Dummy_Atoms_Acceptor_3(void);
	void Generate_Dummy_Atoms_Donor(void);
	void Rotate_Dummy_Atoms(int theta);						// to generate the coordinates of dummy atoms with rotation of theta degree
	void GenerateWater(void);							// to generate water configuration based on dummy atoms
	void Gen_xyz(double x_A, double y_A, double z_A, double x_B, double y_B, double z_B, double x_C, double y_C, double z_C, double r, double theta, double phi, double& x_D, double& y_D, double& z_D);

	double CalGradient_Wat_Opt(void);
	double CalObjectiveFunction_Mol_Wat_Dimer_E_Int(void);	// to calculate the interaction energy betweem compound and water
	double Optimize_Water_Pose(double Para_In[], double& E_Int_Min, double& r_Min);	// to optimize water positions with two degrees of freedom
	int IsGradOptimized(double ftol);
	int Is_Three_Atoms_Colinear(int ia, int ib, int ic);

};



class CFitting
{
public:
	CForceField *myForceField;
	int n_Mol, n_Type, n_Para, n_Para_CG, n_Para_Thole, n_Para_Aniso, n_Para_LP, n_Para_LJ, n_LJ_Fit, n_Para_NBFix, n_NBFix_Fit, n_Aniso_Fit, n_LP_Fit, n_CG_Cal;
	int Last_CG_To_Cal;
	CMol* pMol_List[MAX_MOL];
	CMol *pMol_ESP, *pNMA_Single, *pNMA_Dimer, *pMol_Ar, *pMol_Wat;
	int Para_CG[MAX_NUM_PARA], Para_Thole[MAX_NUM_PARA], CG_Eq_List[MAX_NUM_PARA], Para_Aniso[MAX_NUM_PARA], Para_LP[MAX_NUM_PARA], Para_Aniso_n_Equiv[MAX_NUM_PARA], Aniso_Equiv[MAX_NUM_PARA][16];
	int Para_LJ[MAX_NUM_PARA], Para_NBFix[MAX_NUM_PARA];
	double CG_Var_Range[MAX_NUM_PARA];	// The range used to determine the lower and upper boundaries
	double Para_List[MAX_NUM_PARA], *pCG, *pThole, *pAniso, *pLP, *pLJ, *pNBFix;
	double LowBound[MAX_NUM_PARA], UpBound[MAX_NUM_PARA];
	double NLOpt_LowBound[MAX_NUM_PARA], NLOpt_UpBound[MAX_NUM_PARA];
	int Is_Para_Fixed[MAX_NUM_PARA];
	double NLOpt_Para[MAX_NUM_PARA];
	int n_Para_Free, Idx_NLOpt_Org[MAX_NUM_PARA], Idx_Org_NLOpt[MAX_NUM_PARA];
	double Grad_List[MAX_NUM_PARA];
	ATOMTYPE AtomType[MAX_ATOM_TYPE];
	int n_Heavy_Atm, HeavyAtmList[MAX_ATOM_TYPE];
	double Alpha_0[MAX_ATOM_TYPE];	// target alpha value. Miller's is used. 

	int LP_Rep[MAX_ATOM_TYPE];

	int n_CG_Atom, CG_Atom_Group_Size[MAX_ATOM_TYPE], CG_Atom_Group_List[MAX_ATOM_TYPE][4];
	double cg_0[MAX_ATOM_TYPE];	// target charges

	int n_User_CG, User_CG_Count[MAX_ATOM_TYPE], User_CG_Group_List[MAX_ATOM_TYPE][4];	// data used for user-defined CG constraints
	double User_CG_q0[MAX_ATOM_TYPE], User_CG_k[MAX_ATOM_TYPE];

	double Chi_SQ, Chi_SQ_ESP, Chi_SQ_RSTR_CG, Chi_SQ_RSTR_Alpha, Chi_SQ_RSTR_Thole, Chi_SQ_RSTR_Anisotropy, Chi_SQ_E_int, Chi_SQ_Wat_Emin, Chi_SQ_Wat_Rmin, Chi_SQ_f_Pert, Chi_SQ_User;
	
	
	CFitting(void);
	void RegisterAtomType(CMol* pMol, int n);	//n - is the number of atom will be precesses
	void DetermineTypeID(void);
	void ProcessAtomTypeConstraint(void);
	int IndexAtomType(char szMolName[], int Idx);	// >=0 found this type; -1 - this is a new type
	void DistributeData(void);
	int FindTheRepresentative(int Idx);	//if Idx is an redundant atom type, find its representative; otherwise, just return Idx
	void PrintInfo(void);
	double CalObjectiveFunction(void);
	void CalGradient(void);
	
	void GenIniRandomParameters(void);
	void Optimization_NLOpt(void);
	
	double Cal_Chi_SQ_ESP(void);
	double Cal_Chi_SQ_ChargeRestraint(void);
	double Cal_Chi_SQ_TholeRestraint(void);
	double Cal_Chi_SQ_AlphaRestraint(void);
	double Cal_Chi_SQ_AnisotropyRestraint(void);
	double Cal_Chi_SQ_E_Int_Dimer(void);
	double Cal_Chi_SQ_User(void);
	double Cal_Chi_SQ_Force_Pert_Ar(void);
	double Cal_Chi_SQ_Int_Energy_Mol_Water(int OutputFlag);
	double Cal_Chi_SQ_Int_Energy_Mol_Ar(void);
	void UpdateMolParameters(void);
	void UpdateMolParameters_NMA_Single(void);
	void UpdateMolParameters_NMA_Dimer(void);
	void UpdateMolParameters_Mol_Ar(void);
	void UpdateMolParameters_Mol_Water(void);
	void ReadESPData(void);
	void ReadDimerData(void);
	void ReadArPertForceData(void);
	void ReadArPertEnergyData(void);
	void ReadMolWaterIntEnergyData(void);
	void ReadMolWater_Emin_Rmin(void);
	void PreOptimizeConfigurations(void);
	int FindAtomsinMol(char szAtmName[], int IdxList[]);
	void SetupMillerPolarizability(void);
	void SetupTargetChargesFromIniCharges(void);
	void AssignCoordinate(CMol *pMol, double x[], double y[], double z[], int n);
	void CollectCoordinate(CMol *pMol, double x[], double y[], double z[], int n);
	void MergeDrudeForces(CMol *pMol, double fx[], double fy[], double fz[]);
	void BuildRealAtomList(CMol *pMol, int& LenList, int *List);

	void Generate_Water_Confs_Donor(void);
	void Generate_Water_Confs_Acceptor(void);
	void To_Detect_Invalid_Acceptor_Donor(void);

	int To_Find_The_First_LP(int AtomIdx);
	
	void SaveParameters(void);
	void SaveParameters_CHARMM(void);
	void ReadParameters(char szName[]);
	int IsThereFarDrudes(double& R_Drude_Max);
	int To_Find_Stable_Parameters(void);

};

CFitting fitting;
COptimizeWaterDimer OptWater;


CFitting::CFitting(void)
{
	
	n_Mol = 0;
	pMol_ESP = pMol_Wat = pMol_Ar = NULL;
	n_Para = n_Para_CG = n_Para_Thole = n_Para_Aniso = n_Para_LP = n_Para_LJ = n_LJ_Fit = n_NBFix_Fit = n_Para_NBFix = n_Aniso_Fit = n_LP_Fit = n_CG_Cal = 0;
	
	memset(AtomType, 'x', sizeof(ATOMTYPE)*MAX_ATOM_TYPE);	//set with strange symbols used for debug convenience
	myForceField = NULL;
}


void CFitting::RegisterAtomType(CMol* pMol, int n)	//when there exists a test charge, n = nAtom = 1; otherwise n = nAtom
{
	int i;
	char szMolName[8];
	
	n_Type = 0;
	strcpy(szMolName, pMol->MolName);
	for(i=0; i<n; i++)	{	//all atoms will be added into the list
		AtomType[n_Type].pMol = pMol;
		strcpy(AtomType[n_Type].MolName, szMolName);
		strcpy(AtomType[n_Type].ChemName, pMol->ChemName[i]);
		AtomType[n_Type].AtomIdx = i;
		AtomType[n_Type].Atom_Count = 1;	// all are assumed unique
		AtomType[n_Type].IsRedundant=0;
		AtomType[n_Type].IsDrudeHost = 0;
		AtomType[n_Type].NotFitCharge=0;
		AtomType[n_Type].NotFitThole=1;		//default: do not fit thole
		AtomType[n_Type].NotFitLJ=0;	// default: fit LJ parameter
		AtomType[n_Type].n_Atom_Neutral=0;
//		AtomType[n_Type].LJ_E_Min = pMol->Para_LJ_Epsilon[i];
//		AtomType[n_Type].LJ_r_Min = pMol->Para_LJ_Sigma[i];
		AtomType[n_Type].CG = pMol->CG[i];
		AtomType[n_Type].alpha = pMol->alpha[i];
		AtomType[n_Type].thole = pMol->thole[i];
		AtomType[n_Type].mass = pMol->mass[i];
		if(pMol->IsDrude[i])	{
			AtomType[n_Type-1].IsDrudeHost = 1;	//the previous one is a host of this drude
			AtomType[n_Type-1].NotFitThole = 0;	//if Drude host
		}
		n_Type++;
	}
	if(strcmp(pMol->AtomName[n-1], "CAL")==0)	{	//the test charge
//		printf("Warning!!!\nThe test charge is added as a variable atom type.\n");
		fprintf(fFile_Run_Log, "Warning!!!\nThe test charge is added as a variable atom type.\n");
		fflush(fFile_Run_Log);
	}

	DetermineTypeID();	// determine atom types for all atoms in mol, pMol_ESP

	pMol_List[n_Mol] = pMol;
	n_Mol++;
}

void CFitting::DetermineTypeID(void)
{
	int i;

	for(i=0; i<n_Type; i++)	{
//		if( strncmp(pMol_ESP->AtomName[i], "LP", 2) == 0 )	{	// This is a lone pair
		if( pMol_ESP->mass[i] < 1.0E-10 )	{	// This is a lone pair
			AtomType[i].TypeID = ATM_TYPE_LP;
		}
		else if(pMol_ESP->AtomName[i][0] == 'H')	{	// H atom
			AtomType[i].TypeID = ATM_TYPE_H;
		}
		else if(pMol_ESP->IsDrude[i])	{	// drude
			AtomType[i].TypeID = ATM_TYPE_DRUDE;
		}
		else if( ((i+1)<n_Type) && (pMol_ESP->IsDrude[i+1]) )	{	// heavy atom
			AtomType[i].TypeID = ATM_TYPE_HEAVY;
		}
		else	{
//			printf("Unknow atom type: %s\nQuit\n", pMol_ESP->AtomName[i]);
			fprintf(fFile_Run_Log, "Unknow atom type: %s\nQuit\n", pMol_ESP->AtomName[i]);
			fflush(fFile_Run_Log);
			exit(1);
		}
	}
}

#define MAX_LEN_LINE	(1024)
#define MAX_ITEM		(512)
void CFitting::ProcessAtomTypeConstraint(void)
{
	char szLine[MAX_LEN_LINE], ItemList[MAX_ITEM][32], szDrudeName[32], szAtomName[32], szTmp[256];
	int n_Real_Atom, RealAtom_List[MAX_NUM_PARA], Atom_CG_Constrained;
	int Idx_in_Conf, IdxNBFix_Rec, *IsNBFix_To_Be_Fitted, IdxLJ_Rec, *IsLJ_To_Be_Fitted;
	int n_Item, Idx_Rep, Idx_Redundant, Idx_CG_Not_Fit, i, j, AtomIdx, nAtm_From_Name, IdxList[512], Idx_Item, Para_Offset;
	int nAtm_in_Neutral_Frag, IdxList_Frag[MAX_N_ATOM_NEUTRAL], Idx_Sel, Count_Idx_Sel, Counter;
	CMol *pMol;


	SetupTargetChargesFromIniCharges();
	BuildRealAtomList(pMol_ESP, n_Real_Atom, RealAtom_List);	// include the Ar probe atom
	
	//start	to read the constraint rules
	Idx_in_Conf = 0;
	IsNBFix_To_Be_Fitted = new int[myForceField->n_Rec_NBFix];
	memset(IsNBFix_To_Be_Fitted, 0, sizeof(int)*myForceField->n_Rec_NBFix);
	IsLJ_To_Be_Fitted = new int[myForceField->n_Rec_LJ];
	memset(IsLJ_To_Be_Fitted, 0, sizeof(int)*myForceField->n_Rec_LJ);

	while(1)	{
		Idx_in_Conf = FindItemInConfFile("equivalent", Idx_in_Conf);
		if(Idx_in_Conf < 0)	{	// have checked all records in configuration file
			break;
		}
		strcpy(szLine, szReadConf_File[Idx_in_Conf]);
		
		n_Item = SplitString(szLine, ItemList);
		
		for(Idx_Item=1; Idx_Item<n_Item; Idx_Item++)	{
			nAtm_From_Name = 1;	// since index of atom is the input now!
			IdxList[0] = RealAtom_List[atoi(ItemList[Idx_Item]) - 1];
			//				nAtm_From_Name = FindAtomsinMol(ItemList[Idx_Item], IdxList);
			if(Idx_Item==1)	{
				Idx_Rep = IndexAtomType(pMol_ESP->MolName, IdxList[0]);
				
				for(i=1; i<nAtm_From_Name; i++)	{
					Idx_Redundant = IndexAtomType(pMol_ESP->MolName, IdxList[i]);	//index starting from 0
					
					AtomType[Idx_Redundant].IsRedundant = 1;
					AtomType[Idx_Redundant].CG_Represent = Idx_Rep;
					if(AtomType[Idx_Redundant].IsDrudeHost)	{
						AtomType[Idx_Redundant].Thole_Represent = Idx_Rep;
					}
					AtomType[Idx_Rep].Atom_Count++;
				}
			}
			else	{
				for(i=0; i<nAtm_From_Name; i++)	{
					Idx_Redundant = IndexAtomType(pMol_ESP->MolName, IdxList[i]);	//index starting from 0
					
					AtomType[Idx_Redundant].IsRedundant = 1;
					AtomType[Idx_Redundant].CG_Represent = Idx_Rep;
					if(AtomType[Idx_Redundant].IsDrudeHost)	{
						AtomType[Idx_Redundant].Thole_Represent = Idx_Rep;

						AtomType[Idx_Redundant+1].IsRedundant = 1;	// the equivalent drude
						AtomType[Idx_Redundant+1].CG_Represent = Idx_Rep+1;
						AtomType[Idx_Rep+1].Atom_Count++;
					}
					AtomType[Idx_Rep].Atom_Count++;
				}
			}
		}
		
		Idx_in_Conf++;
	}

	if(Not_Fit_Elec_Para == 0)	{	// to fit electrostatical parameters
		Idx_in_Conf = 0;
		while(1)	{
			Idx_in_Conf = FindItemInConfFile("neutral", Idx_in_Conf);
			if(Idx_in_Conf < 0)	{	// have checked all records in configuration file
				break;
			}
			strcpy(szLine, szReadConf_File[Idx_in_Conf]);
			
			n_Item = SplitString(szLine, ItemList);
			
			nAtm_in_Neutral_Frag = 0;
			
			for(Idx_Item=1; Idx_Item<n_Item; Idx_Item++)	{
				nAtm_From_Name = 1;	// since index of atom is the input now!
				IdxList[0] = RealAtom_List[atoi(ItemList[Idx_Item]) - 1];
//				nAtm_From_Name = FindAtomsinMol(ItemList[Idx_Item], IdxList);
				for(i=0; i<nAtm_From_Name; i++)	{
					IdxList_Frag[nAtm_in_Neutral_Frag] = FindTheRepresentative(IdxList[i]);
					nAtm_in_Neutral_Frag++;

					//start	to add the associated drude
					if(pMol_ESP->mass[IdxList_Frag[nAtm_in_Neutral_Frag-1]] > 1.50)	{	// it is a heavy atom
						IdxList_Frag[nAtm_in_Neutral_Frag] = IdxList_Frag[nAtm_in_Neutral_Frag-1] + 1;	// the index of the drude is next to the heavy atom
						nAtm_in_Neutral_Frag++;
					}
					//end	to add the associated drude
				}
			}
			for(i=0; i<nAtm_in_Neutral_Frag; i++)	{
				
				if((AtomType[IdxList_Frag[i]].NotFitCharge == 0) && (pMol_ESP->mass[AtomType[IdxList_Frag[i]].AtomIdx] > 1.0))	{	// a non-drude, non-LP particle (real atom)
//				if(AtomType[IdxList_Frag[i]].IsDrudeHost == 1)	{	// the host of a heavy atom
					break;
				}
			}
			if(i>=nAtm_in_Neutral_Frag)	{	// fail to find a heavy atom within the fragment. Sowmthing wrong?
				fprintf(fFile_Run_Log, "Fail to find a heavy atom within the fragment in string.\n%s\nQuit.\n", szLine);
				fflush(fFile_Run_Log);
				exit(1);
			}
			Count_Idx_Sel = 0;
			Idx_Sel = IdxList_Frag[i];
			
			for(i=0; i<nAtm_in_Neutral_Frag; i++)	{
				if(IdxList_Frag[i] == Idx_Sel)	{
					Count_Idx_Sel++;
				}
			}
			Idx_CG_Not_Fit = Idx_Sel;
			AtomType[Idx_CG_Not_Fit].NotFitCharge = 1;	//the charge will be calculated by the user defined Eq.
			
			Counter = 0;
			AtomType[Idx_CG_Not_Fit].n_Atom_Neutral = 0;
			for(i=0; i<nAtm_in_Neutral_Frag; i++)	{
				if(IdxList_Frag[i] != Idx_CG_Not_Fit)	{
					AtomType[Idx_CG_Not_Fit].weight[Counter] = -1.0/Count_Idx_Sel;
					AtomType[Idx_CG_Not_Fit].Cal_CG_List[Counter] = IdxList_Frag[i];
					Counter++;
				}
			}
			AtomType[Idx_CG_Not_Fit].n_Atom_Neutral = Counter;
			
			Idx_in_Conf++;
		}
		
		
		Idx_in_Conf = 0;
		while(1)	{
			Idx_in_Conf = FindItemInConfFile("Not_Fit_Thole", Idx_in_Conf);
			if(Idx_in_Conf < 0)	{	// have checked all records in configuration file
				break;
			}
			strcpy(szLine, szReadConf_File[Idx_in_Conf]);
			n_Item = SplitString(szLine, ItemList);
			
			for(Idx_Item=1; Idx_Item<n_Item; Idx_Item++)	{
				nAtm_From_Name = 1;	// since index of atom is the input now!
				IdxList[0] = RealAtom_List[atoi(ItemList[Idx_Item]) - 1];
//				nAtm_From_Name = FindAtomsinMol(ItemList[Idx_Item], IdxList);
				
				for(i=0; i<nAtm_From_Name; i++)	{
					Idx_Redundant = IndexAtomType(pMol_ESP->MolName, IdxList[i]);	//index starting from 0
					
					if(AtomType[Idx_Redundant].IsDrudeHost)	{
						AtomType[Idx_Redundant].NotFitThole = 1;	// just do not fit this thole parameter
						AtomType[Idx_Redundant].Thole_Represent = Idx_Redundant;	// pointing to itself
						AtomType[Idx_Redundant].thole = AtomType[Idx_Redundant].pMol->thole[Idx_Redundant];	// the default thole parameter
					}
				}
			}
			Idx_in_Conf++;
		}
		
		
		
		/*
		else if(strncmp(szLine, "same_thole", 10) == 0)	{	// to process entries defining list same thole for two non-equivalent atoms. obsolete. 
		n_Item = SplitString(szLine, ItemList);
		
		  Idx_Rep = IndexAtomType(ItemList[1], atoi(ItemList[2])-1);	//index starting from 0
		  if(AtomType[Idx_Rep].IsDrudeHost != 1)	{	//not a drude host
		  //				printf("Atom %s in molecule %s is not a host of drude! \n", 
		  fprintf(fFile_Run_Log, "Atom %s in molecule %s is not a host of drude! \n", 
		  ItemList[2], ItemList[1]);
		  fflush(fFile_Run_Log);
		  }
		  for(i=3; i<n_Item; i+=2)	{
		  Idx_Redundant = IndexAtomType(ItemList[i], atoi(ItemList[i+1])-1);	//index starting from 0
		  if(AtomType[Idx_Redundant].IsDrudeHost != 1)	{	//not a drude host
		  //					printf("Atom %s in molecule %s is not a host of drude! \n", 
		  fprintf(fFile_Run_Log, "Atom %s in molecule %s is not a host of drude! \n", 
		  ItemList[i+1], ItemList[i]);
		  fflush(fFile_Run_Log);
		  }
		  
			AtomType[Idx_Redundant].NotFitThole = 1;
			AtomType[Idx_Redundant].Thole_Represent = Idx_Rep;
			}
			}
		*/
		
		//start	to apply the constraint, the whole molecule is neutral
		nAtm_in_Neutral_Frag = 0;
		for(i=0; i<n_Type; i++)	{
			IdxList_Frag[i] = FindTheRepresentative(i);
		}
		nAtm_in_Neutral_Frag = n_Type;
		
		for(i=0; i<nAtm_in_Neutral_Frag; i++)	{	// Need improvement!!!
//			if(AtomType[IdxList_Frag[i]].IsDrudeHost == 1)	{	// the host of a heavy atom
			if((AtomType[IdxList_Frag[i]].NotFitCharge == 0) && (pMol_ESP->IsDrude[AtomType[IdxList_Frag[i]].AtomIdx]==0))	{	// not a drude
				break;
			}
		}
		if(i>=nAtm_in_Neutral_Frag)	{	// fail to find a valid atom within the molecule. Sowmthing wrong?
			Quit_With_Error_Msg("Fail to find a heavy atom within the molecule.\nQuit.\n");
		}
		Count_Idx_Sel = 0;
		Idx_Sel = IdxList_Frag[i];
		
		for(i=0; i<nAtm_in_Neutral_Frag; i++)	{
			if(IdxList_Frag[i] == Idx_Sel)	{
				Count_Idx_Sel++;
			}
		}
		Idx_CG_Not_Fit = Idx_Sel;
		AtomType[Idx_CG_Not_Fit].NotFitCharge = 1;	//the charge will be calculated by the user defined Eq.
		Last_CG_To_Cal = Idx_CG_Not_Fit;
		SameAtom_Count = Count_Idx_Sel;

		
		Counter = 0;
		AtomType[Idx_CG_Not_Fit].n_Atom_Neutral = 0;
		for(i=0; i<nAtm_in_Neutral_Frag; i++)	{
			if(IdxList_Frag[i] != Idx_CG_Not_Fit)	{
				AtomType[Idx_CG_Not_Fit].weight[Counter] = -1.0/Count_Idx_Sel;
				AtomType[Idx_CG_Not_Fit].Cal_CG_List[Counter] = IdxList_Frag[i];
				Counter++;
			}
		}
		AtomType[Idx_CG_Not_Fit].n_Atom_Neutral = Counter;
		//end	to apply the constraint, the whole molecule is neutral
		//start	to read the constraint rules
		
		//start	to parse the user-defined CG constraints
		n_User_CG = 0;
		memset(User_CG_Count, 0, sizeof(int)*MAX_ATOM_TYPE);

		Idx_in_Conf = 0;
		while(1)	{
			Idx_in_Conf = FindItemInConfFile("constrain_CG", Idx_in_Conf);
			if(Idx_in_Conf < 0)	{	// have checked all records in configuration file
				break;
			}
			strcpy(szLine, szReadConf_File[Idx_in_Conf]);
			n_Item = sscanf(szLine, "%s %d %lf", szTmp, &Atom_CG_Constrained, &(User_CG_k[n_User_CG]));
			if(n_Item == 3)	{
				nAtm_From_Name = 1;	// since index of atom is the input now!
				IdxList[0] = RealAtom_List[Atom_CG_Constrained - 1];
//				nAtm_From_Name = FindAtomsinMol(szAtomName, IdxList);
				if(nAtm_From_Name >= 1)	{
					AtomIdx = IdxList[0];	// only consider one atom for constraint

					if( (pMol_ESP->mass[AtomIdx] < 1.0E-10) || ( (pMol_ESP->mass[AtomIdx] > 1.0) && (pMol_ESP->mass[AtomIdx] < 1.01) ) )	{	// lone pair or H atom
						User_CG_Count[n_User_CG] = 1;
						User_CG_Group_List[n_User_CG][0] = AtomIdx;
						User_CG_q0[n_User_CG] = pMol_ESP->CG[AtomIdx];
						fprintf(fFile_Run_Log, "Add user-defined CG constraint for atom %s with q0= %6.3lf and k = %6.3lf\n", pMol_ESP->AtomName[IdxList[0]], User_CG_q0[n_User_CG], User_CG_k[n_User_CG]);
						n_User_CG++;
					}
					else if(pMol_ESP->mass[AtomIdx] > 2.0)	{	// a heavy atom. The sum of charges of host and drude is constrained. 
						User_CG_Count[n_User_CG] = 2;
						User_CG_Group_List[n_User_CG][0] = AtomIdx;		// the heavy atom
						User_CG_Group_List[n_User_CG][1] = AtomIdx+1;	// the drude
						User_CG_q0[n_User_CG] = pMol_ESP->CG[AtomIdx] + pMol_ESP->CG[AtomIdx+1];
						fprintf(fFile_Run_Log, "Add user-defined CG constraint for atom %s with q0= %6.3lf and k = %6.3lf\n", pMol_ESP->AtomName[IdxList[0]], User_CG_q0[n_User_CG], User_CG_k[n_User_CG]);
						n_User_CG++;
					}
				}
			}

			Idx_in_Conf++;
		}
		fflush(fFile_Run_Log);
		//start	to parse the user-defined CG constraints

		
		PrintInfo();
		
		
		//start	to compile a parameter list and setup the value
		n_Para_CG = 0;
		n_Para_Thole = 0;
		n_CG_Cal = 0;
		for(i=0; i<n_Type; i++)	{
			if(AtomType[i].IsRedundant == 0)	{	//parameters will be fitted
				if(AtomType[i].NotFitThole == 0)	{	//drude host, thole will be fitted
					AtomType[i].thole = AtomType[i].pMol->thole[AtomType[i].AtomIdx];
					Para_Thole[n_Para_Thole] = i;
					n_Para_Thole++;
				}
				
				if(AtomType[i].NotFitCharge == 0)	{	//charge will be fitted
					AtomType[i].CG = AtomType[i].pMol->CG[AtomType[i].AtomIdx];
					Para_CG[n_Para_CG] = i;
					n_Para_CG++;
				}
				
				if(AtomType[i].n_Atom_Neutral > 0)	{
					CG_Eq_List[n_CG_Cal] = i;
					n_CG_Cal++;
				}
			}
		}
		
		//start	to set up para for Aniso
		n_Aniso_Fit = 0;
		if(To_Fit_Aniso)	{
			memset(Para_Aniso, 0, sizeof(int)*MAX_NUM_PARA);
			memset(Para_Aniso_n_Equiv, 0, sizeof(int)*MAX_NUM_PARA);
			memset(Aniso_Equiv, 0, sizeof(int)*MAX_NUM_PARA*16);
			
			int Atm_O, Atm_O_Master;	// the index of O and the representative O atom
			for(i=0; i<pMol_ESP->nAniso; i++)	{
				Atm_O = pMol_ESP->AnisoList[i][0];	// O atom
				Atm_O_Master = FindTheRepresentative(Atm_O);
				
				if(Atm_O_Master == Atm_O)	{	// a new entry
					Para_Aniso[n_Aniso_Fit] = i;
					n_Aniso_Fit++;
				}
				else	{	// to find the entry
					for(j=0; j<n_Aniso_Fit; j++)	{
						if(Atm_O_Master == pMol_ESP->AnisoList[Para_Aniso[j]][0])	{
							break;
						}
					}
					if(j>=n_Aniso_Fit)	{
						//					printf("Fail to find the equivalent anisotropy term for %d entry.\nQuit\n", i+1);
						fprintf(fFile_Run_Log, "Fail to find the equivalent anisotropy term for %d entry.\nQuit\n", i+1);
						fflush(fFile_Run_Log);
						exit(1);
					}
					else	{
						Aniso_Equiv[j][Para_Aniso_n_Equiv[j]] = i;
						Para_Aniso_n_Equiv[j] ++;
					}
				}
			}
		}
		n_Para_Aniso = n_Aniso_Fit * 2;
		//end	to set up para for Aniso
	
		//start	to set up para for lone pairs
		n_LP_Fit = 0;

		if(To_Fit_LP_Para)	{
			memset(Para_LP, 0, sizeof(int)*MAX_NUM_PARA);

			for(i=0; i<pMol_ESP->nLP; i++)	{
				LP_Rep[i] = To_Find_The_First_LP(FindTheRepresentative(pMol_ESP->LP_PosHost[i][0]));
				if( (pMol_ESP->LP_Dist[i] > UP_BOUND_LP_r_NEG) && (pMol_ESP->LP_Dist[i] < 0.0) )	{	// a dummy LP
					LP_Rep[i] = -1;		// a dummy LP, not to fit
				}
				else if(LP_Rep[i] == i)	{	// it is a representative LP. To add two parameters to fit
					Para_LP[n_LP_Fit] = i;
					n_LP_Fit++;
				}
			}
		}
		n_Para_LP = n_LP_Fit * 2;
		//end	to set up para for lone pairs

	}
	else	{
		n_Para_CG = n_Para_Thole = n_Para_Aniso = n_Para_LP = 0;
	}


	//start	to set up para for LJ
	if(To_Fit_LJ)	{
		n_LJ_Fit = 0;

		for(i=0; i<n_Type; i++)	{
			IdxLJ_Rec = pMol_ESP->LJ_Para_Rec[AtomType[i].AtomIdx];
			if( (IdxLJ_Rec >= 0) && (IsLJ_To_Be_Fitted[IdxLJ_Rec]==0) )	{
				IsLJ_To_Be_Fitted[IdxLJ_Rec] = 1;
				Para_LJ[n_LJ_Fit] = IdxLJ_Rec;	// pointing to the LJ record in ForceFiled
				n_LJ_Fit++;
//				printf("LJ_para %3d\n", i+1);
			}
		}
		n_Para_LJ = n_LJ_Fit*2;
	}
	//end	to set up para for LJ

	//start	to determine the number of NBFix parameters to fit
	n_NBFix_Fit = 0;
	Idx_in_Conf = 0;
	while(1)	{
		Idx_in_Conf = FindItemInConfFile("Fit_Para_NBFIX", Idx_in_Conf);
		if(Idx_in_Conf < 0)	{	// have checked all records in configuration file
			break;
		}
		strcpy(szLine, szReadConf_File[Idx_in_Conf]);
		n_Item = SplitString(szLine, ItemList);

		if(n_Item >= 3)	{	// maybe a valid line 
			IdxNBFix_Rec = myForceField->QueryNBFix(ItemList[1], ItemList[2]);

			if(IdxNBFix_Rec>=0)	{	// a valid record
				if(IsNBFix_To_Be_Fitted[IdxNBFix_Rec]==0)	{	// only add to the list of fitting one time
					IsNBFix_To_Be_Fitted[IdxNBFix_Rec] = 1;
					Para_NBFix[n_NBFix_Fit] = IdxNBFix_Rec;
					n_NBFix_Fit ++;
				}
			}
		}
		Idx_in_Conf++;
	}
	n_Para_NBFix = n_NBFix_Fit * 2; 
	//end	to determine the number of NBFix parameters to fit


	n_Para = n_Para_CG + n_Para_Thole + n_Para_Aniso + n_Para_LP + n_Para_LJ + n_Para_NBFix;
	
	pCG = Para_List;
	for(i=0; i<n_Para_CG; i++)	{
		pCG[i] = AtomType[Para_CG[i]].CG;
	}

	pThole = pCG+n_Para_CG;
	for(i=0; i<n_Para_Thole; i++)	{
		pThole[i] = AtomType[Para_Thole[i]].thole;
	}

	pAniso = pThole+n_Para_Thole;
	for(i=0; i<n_Aniso_Fit; i++)	{
		double K11, K22, K33;
		pMol = pMol_ESP;
		K11 = pMol->Para_Aniso[Para_Aniso[i]][0];
		K22 = pMol->Para_Aniso[Para_Aniso[i]][1];
		K33 = pMol->Para_Aniso[Para_Aniso[i]][2];
		Para_Aniso_K_A(K11, K22, K33, pAniso[i*2], pAniso[i*2+1]);
	}

	pLP = pAniso+n_Para_Aniso;
	for(i=0; i<n_LP_Fit; i++)	{
		pLP[i*2  ] = pMol_ESP->LP_Dist[Para_LP[i]];
		pLP[i*2+1] = pMol_ESP->LP_Theta[Para_LP[i]] * radianInv;	// save the angle in radian
	}

	pLJ = pLP+n_Para_LP;
	for(i=0; i<n_LJ_Fit; i++)	{	// parameters for 14 is not fitted here !!!!
		IdxLJ_Rec = Para_LJ[i];
		pLJ[i*2  ] = myForceField->LJ_Rec[IdxLJ_Rec].para[1];	// E_Min
		pLJ[i*2+1] = myForceField->LJ_Rec[IdxLJ_Rec].para[2];	// r_Min
	}

	pNBFix = pLJ+n_Para_LJ;
	for(i=0; i<n_NBFix_Fit; i++)	{
		IdxNBFix_Rec = Para_NBFix[i];
		pNBFix[i*2  ] = myForceField->NBFix_Rec[IdxNBFix_Rec].para[0];
		pNBFix[i*2+1] = myForceField->NBFix_Rec[IdxNBFix_Rec].para[1];
	}
	//end	to compile a parameter list and setup the value


	//start	to print the info about those parameters to fit
//	printf("Info about parameters to fit: \n");
	fprintf(fFile_Run_Log, "Info about parameters to fit: \n");
	fflush(fFile_Run_Log);
	for(i=0; i<n_Para_CG; i++)	{
		pMol = AtomType[Para_CG[i]].pMol;
		AtomIdx = AtomType[Para_CG[i]].AtomIdx;

		switch(AtomType[Para_CG[i]].TypeID)	{
		case ATM_TYPE_HEAVY:
			LowBound[i] = LOW_BOUND_CG_HEAVY_ATM;
			UpBound[i] = UP_BOUND_CG_HEAVY_ATM;
			CG_Var_Range[i] = D_CG_HEAVY_ATOM;
			break;
		case ATM_TYPE_DRUDE:
			LowBound[i] = LOW_BOUND_CG_DRUDE;
			UpBound[i] = UP_BOUND_CG_DRUDE;
			CG_Var_Range[i] = D_CG_HEAVY_ATOM;
			break;
		case ATM_TYPE_H:
			LowBound[i] = LOW_BOUND_CG_H;
			UpBound[i] = UP_BOUND_CG_H;
			CG_Var_Range[i] = D_CG_H_ATOM;
			break;
		case ATM_TYPE_LP:
			LowBound[i] = LOW_BOUND_CG_LP;
			UpBound[i] = UP_BOUND_CG_LP;
			CG_Var_Range[i] = D_CG_LP;
			break;
		default:
//			printf("Wrong type ID for atom: %s\nQuit\n", pMol->AtomName[AtomIdx]);
			fprintf(fFile_Run_Log, "Wrong type ID for atom: %s\nQuit\n", pMol->AtomName[AtomIdx]);
			fflush(fFile_Run_Log);
			exit(1);
			break;
		}
		
//		printf("Parameter %3d,  CG    of atom %4s in mol %5s. Idx = %d Bound (%6.3lf,%6.3lf)\n", 
		fprintf(fFile_Run_Log, "Parameter %2d,   CG       of atom %4s in mol %5s. Idx = %d Bound (%6.3lf,%6.3lf)\n", 
			i+1, pMol->AtomName[AtomIdx],pMol->MolName, AtomIdx+1, LowBound[i], UpBound[i]);
		fflush(fFile_Run_Log);
	}
	
	Para_Offset = n_Para_CG;
	for(i=0; i<n_Para_Thole; i++)	{
		pMol = AtomType[Para_Thole[i]].pMol;
		AtomIdx = AtomType[Para_Thole[i]].AtomIdx;

		LowBound[i+Para_Offset] = LOW_BOUND_THOLE;
		UpBound[i+Para_Offset] = UP_BOUND_THOLE;
		CG_Var_Range[i+Para_Offset] = D_THOLE;
		
		fprintf(fFile_Run_Log, "Parameter %2d,  Thole     of atom %4s in mol %5s. Idx = %d Bound (%6.3lf,%6.3lf)\n", 
			i+1+Para_Offset, pMol->AtomName[AtomIdx],pMol->MolName, AtomIdx+1, LowBound[i+n_Para_CG], UpBound[i+n_Para_CG]);
		fflush(fFile_Run_Log);
	}

	Para_Offset += n_Para_Thole;
	for(i=0; i<n_Aniso_Fit; i++)	{
		LowBound[i*2+Para_Offset] = LOW_BOUND_ANISO;
		UpBound[i*2+Para_Offset] = UP_BOUND_ANISO;
		LowBound[i*2+1+Para_Offset] = LOW_BOUND_ANISO;
		UpBound[i*2+1+Para_Offset] = UP_BOUND_ANISO;
		CG_Var_Range[i*2+Para_Offset] = D_ANISO;
		CG_Var_Range[i*2+1+Para_Offset] = D_ANISO;
		AtomIdx = pMol_ESP->AnisoList[Para_Aniso[i]][0];

		fprintf(fFile_Run_Log, "Parameter %2d,  aniso_A11 of atom %4s in mol %5s. Bound (%6.3lf,%6.3lf)\n", 
			i*2+1+Para_Offset, pMol->AtomName[AtomIdx], pMol->MolName, LowBound[i*2+Para_Offset], UpBound[i*2+Para_Offset]);
		fprintf(fFile_Run_Log, "Parameter %2d,  aniso_A22 of atom %4s in mol %5s. Bound (%6.3lf,%6.3lf)\n", 
			i*2+2+Para_Offset, pMol->AtomName[AtomIdx], pMol->MolName, LowBound[i*2+1+Para_Offset], UpBound[i*2+1+Para_Offset]);
		fflush(fFile_Run_Log);
	}

	Para_Offset += n_Para_Aniso;
	for(i=0; i<n_LP_Fit; i++)	{
		if(Para_List[i*2+Para_Offset] < 0.0)	{	// negative, "bisector"
			LowBound[i*2+Para_Offset] = LOW_BOUND_LP_r_NEG;
			UpBound[i*2+Para_Offset] = UP_BOUND_LP_r_NEG;
		}
		else	{
			LowBound[i*2+Para_Offset] = LOW_BOUND_LP_r;
			UpBound[i*2+Para_Offset] = UP_BOUND_LP_r;
		}
		LowBound[i*2+1+Para_Offset] = LOW_BOUND_LP_THETA;
		UpBound[i*2+1+Para_Offset] = UP_BOUND_LP_THETA;

		fprintf(fFile_Run_Log, "Parameter %2d, distance   of lone pair %6s.                Bound (%6.3lf,%6.3lf)\n", 
			i*2+1+Para_Offset, pMol_ESP->AtomName[Para_LP[i]], LowBound[i*2+Para_Offset], UpBound[i*2+Para_Offset]);
		fprintf(fFile_Run_Log, "Parameter %2d, angle      of lone pair %6s.                Bound (%6.3lf,%6.3lf)\n", 
			i*2+2+Para_Offset, pMol_ESP->AtomName[Para_LP[i]], LowBound[i*2+1+Para_Offset], UpBound[i*2+1+Para_Offset]);
		fflush(fFile_Run_Log);
	}

	Para_Offset += n_Para_LP;
	for(i=0; i<n_LJ_Fit; i++)	{
		LowBound[i*2+Para_Offset] = LOW_BOUND_LJ_E;
		UpBound[i*2+Para_Offset] = UP_BOUND_LJ_E;
		LowBound[i*2+1+Para_Offset] = LOW_BOUND_LJ_r;
		UpBound[i*2+1+Para_Offset] = UP_BOUND_LJ_r;

		fprintf(fFile_Run_Log, "Parameter %2d, LJ_E_min   of Chem name %6s.       Bound (%6.3lf,%6.3lf)\n", 
			i*2+1+Para_Offset, myForceField->LJ_Rec[Para_LJ[i]].Chem, LowBound[i*2+Para_Offset], UpBound[i*2+Para_Offset]);
		fprintf(fFile_Run_Log, "Parameter %2d, LJ_r_min   of Chem name %6s.       Bound (%6.3lf,%6.3lf)\n", 
			i*2+2+Para_Offset, myForceField->LJ_Rec[Para_LJ[i]].Chem, LowBound[i*2+1+Para_Offset], UpBound[i*2+1+Para_Offset]);
		fflush(fFile_Run_Log);
	}

	Para_Offset += n_Para_LJ;
	for(i=0; i<n_NBFix_Fit; i++)	{
		LowBound[i*2+Para_Offset] = LOW_BOUND_NBFIX_E;
		UpBound[i*2+Para_Offset] = UP_BOUND_NBFIX_E;
		LowBound[i*2+1+Para_Offset] = LOW_BOUND_NBFIX_r;
		UpBound[i*2+1+Para_Offset] = UP_BOUND_NBFIX_r;

		IdxNBFix_Rec = Para_NBFix[i];
		fprintf(fFile_Run_Log, "Parameter %2d, NBFix_E_min of atom %6s and  %6s.   Bound (%6.3lf,%6.3lf)\n", 
			i*2+1+Para_Offset, myForceField->NBFix_Rec[IdxNBFix_Rec].Chem_1,myForceField->NBFix_Rec[IdxNBFix_Rec].Chem_2, LowBound[i*2+Para_Offset], UpBound[i*2+Para_Offset]);
		fprintf(fFile_Run_Log, "Parameter %2d, NBFix_r_min of atom %6s and  %6s.   Bound (%6.3lf,%6.3lf)\n", 
			i*2+2+Para_Offset, myForceField->NBFix_Rec[IdxNBFix_Rec].Chem_1,myForceField->NBFix_Rec[IdxNBFix_Rec].Chem_2, LowBound[i*2+1+Para_Offset], UpBound[i*2+1+Para_Offset]);
		fflush(fFile_Run_Log);
	}
	//end	to print the info about those parameters to fit


	SetupMillerPolarizability();	// setup Miller's polarizability used for constraints

	delete []IsNBFix_To_Be_Fitted;
	delete []IsLJ_To_Be_Fitted;
}

#define ALPHA_TYPE_C	(6)
#define ALPHA_TYPE_N	(7)
#define ALPHA_TYPE_O	(8)
#define ALPHA_TYPE_F	(9)
#define ALPHA_TYPE_P	(15)
#define ALPHA_TYPE_S	(16)
#define ALPHA_TYPE_CL	(17)
#define ALPHA_TYPE_AR	(18)
#define ALPHA_TYPE_BR	(35)
#define ALPHA_TYPE_I	(53)

//#define ALPHA_H	(0.387)
//#define ALPHA_H	(0.271)	// about 70% of 0.387
#define ALPHA_H	(0.224)	// about 0.387*0.6

void CFitting::SetupMillerPolarizability(void)
{
	int i, Is_LP_DRUDE[MAX_ATOM_TYPE], Is_H[MAX_ATOM_TYPE], H_Count[MAX_ATOM_TYPE], C_Count[MAX_ATOM_TYPE], S_Count[MAX_ATOM_TYPE], TypeID[MAX_ATOM_TYPE];
	int Bond_List[MAX_ATOM_TYPE][4], Bond_Count[MAX_ATOM_TYPE];
	int SP2_O_Count[MAX_ATOM_TYPE], SP2_Flag[MAX_ATOM_TYPE], Bonded_Sp2_Atom[MAX_ATOM_TYPE], Bonded_Sp2_Atom_Restricted[MAX_ATOM_TYPE];
	int Idx, Atom_1, Atom_2, nBond;
	double mass;
	
	memset(Is_LP_DRUDE, 0, sizeof(int)*MAX_ATOM_TYPE);
	memset(Is_H, 0, sizeof(int)*MAX_ATOM_TYPE);
	memset(H_Count, 0, sizeof(int)*MAX_ATOM_TYPE);
	memset(C_Count, 0, sizeof(int)*MAX_ATOM_TYPE);
	memset(S_Count, 0, sizeof(int)*MAX_ATOM_TYPE);
	memset(Bond_Count, 0, sizeof(int)*MAX_ATOM_TYPE);
	memset(Bond_List, 0, sizeof(int)*MAX_ATOM_TYPE*4);
	memset(TypeID, 0, sizeof(int)*MAX_ATOM_TYPE);
	memset(SP2_Flag, 0, sizeof(int)*MAX_ATOM_TYPE);
	memset(SP2_O_Count, 0, sizeof(int)*MAX_ATOM_TYPE);
	memset(Bonded_Sp2_Atom, 0, sizeof(int)*MAX_ATOM_TYPE);
	memset(Bonded_Sp2_Atom_Restricted, 0, sizeof(int)*MAX_ATOM_TYPE);

	memset(H_Donor_Acceptor, 0, sizeof(int)*MAX_ATOM_TYPE);

	for(i=0; i<n_Type; i++)	{
		mass = pMol_ESP->mass[i];

		if( (mass > 1.00) && (mass < 1.10) )	{	// atom mass
			Is_H[i] = 1;
		}
		else	{
			Is_H[i] = 0;
		}

		if( mass < 0.8 )	{	// drude or lone pair
			Is_LP_DRUDE[i] = 1;
		}

		if(AtomType[i].IsDrudeHost)	{
			mass += pMol_ESP->mass[i+1];	// including the drude

			if( (mass > 11.5) && (mass < 12.5) )	{	// C
				TypeID[i] = ALPHA_TYPE_C;
			}
			else if( (mass > 13.5) && (mass < 14.5) )	{	// N
				TypeID[i] = ALPHA_TYPE_N;
				H_Donor_Acceptor[i] = 1;
			}
			else if( (mass > 15.5) && (mass < 16.5) )	{	// O
				TypeID[i] = ALPHA_TYPE_O;
				H_Donor_Acceptor[i] = 1;
			}
			else if( (mass > 18.5) && (mass < 19.5) )	{	// F
				TypeID[i] = ALPHA_TYPE_F;
				H_Donor_Acceptor[i] = 1;
			}
			else if( (mass > 30.5) && (mass < 31.5) )	{	// P
				TypeID[i] = ALPHA_TYPE_P;
			}
			else if( (mass > 31.5) && (mass < 32.5) )	{	// S
				TypeID[i] = ALPHA_TYPE_S;
			}
			else if( (mass > 35.0) && (mass < 36.0) )	{	// CL
				TypeID[i] = ALPHA_TYPE_CL;
			}
			else if( (mass > 39.9) && (mass < 40.0) )	{	// AR
				TypeID[i] = ALPHA_TYPE_AR;
			}
			else if( (mass > 79.5) && (mass < 80.5) )	{	// BR
				TypeID[i] = ALPHA_TYPE_BR;
			}
			else if( (mass > 126.5) && (mass < 127.5) )	{	// I
				TypeID[i] = ALPHA_TYPE_I;
			}
			else	{
//				printf("There is %s atom in this molecule.\nThere is no alpha value for it. \nQuit\n", pMol_ESP->AtomName[i]);
				fprintf(fFile_Run_Log, "Warning> There is %s atom in this molecule.\nThere is no Miller's alpha value for it. \n", pMol_ESP->AtomName[i]);
				fflush(fFile_Run_Log);
//				exit(1);
			}
		}
	}

	//start	to count H atoms connected for each heavy atom
	nBond = pMol_ESP->nBond;
	for(i=0; i<nBond; i++)	{
		Atom_1 = pMol_ESP->BondList[2*i];
		Atom_2 = pMol_ESP->BondList[2*i+1];

		if( (Is_LP_DRUDE[Atom_1]) || (Is_LP_DRUDE[Atom_2]) )	{	// not a bond between two real atoms
			continue;
		}

		Bond_List[Atom_1][Bond_Count[Atom_1]] = Atom_2;
		Bond_Count[Atom_1] ++;
		Bond_List[Atom_2][Bond_Count[Atom_2]] = Atom_1;
		Bond_Count[Atom_2] ++;

		if(Is_H[Atom_1])	{
			H_Count[Atom_2]++;
		}
		if(Is_H[Atom_2])	{
			H_Count[Atom_1]++;
		}
		if(TypeID[Atom_1] == ALPHA_TYPE_C)	{
			C_Count[Atom_2]++;
		}
		if(TypeID[Atom_2] == ALPHA_TYPE_C)	{
			C_Count[Atom_1]++;
		}
		if(TypeID[Atom_1] == ALPHA_TYPE_S)	{
			S_Count[Atom_2]++;
		}
		if(TypeID[Atom_2] == ALPHA_TYPE_S)	{
			S_Count[Atom_1]++;
		}
	}
	//end	to count H atoms connected for each heavy atom

	for(i=0; i<n_Type; i++)	{
		if( Is_H[i] )	{	// H
			if(H_Donor_Acceptor[Bond_List[i][0]])	{	// bonded with a H acceptor, N/O/F
				H_Donor_Acceptor[i] = 1;
			}
		}
	}


	//start	to indentify those restrict sp2 atoms
	for(i=0; i<n_Type; i++)	{
		if(AtomType[i].IsDrudeHost)	{
			if(TypeID[i] == ALPHA_TYPE_C)	{	// Carbon
				if(Bond_Count[i] == 3)	{
					SP2_Flag[i] = 1;
				}
			}
			else if(TypeID[i] == ALPHA_TYPE_N)	{
				if(Bond_Count[i] == 2)	{
					SP2_Flag[i] = 1;
				}
			}
			else if(TypeID[i] == ALPHA_TYPE_O)	{
				if(Bond_Count[i] == 1)	{
					SP2_Flag[i] = 1;
					SP2_O_Count[Bond_List[i][0]] ++;
				}
			}
			else if(TypeID[i] == ALPHA_TYPE_S)	{
				if(Bond_Count[i] == 1)	{
					SP2_Flag[i] = 1;
				}
			}
		}
	}
	//end	to indentify those restrict sp2 atoms

	for(Idx=0; Idx<n_Type; Idx++)	{
		Bonded_Sp2_Atom_Restricted[Idx] = 0;
		for(i=0; i<Bond_Count[Idx]; i++)	{
			if( SP2_Flag[Bond_List[Idx][i]] == 1)	{	// real SP2 atom
				Bonded_Sp2_Atom_Restricted[Idx]++;
			}
		}
	}

	for(i=0; i<n_Type; i++)	{
		if(TypeID[i] == ALPHA_TYPE_N)	{
			if(Bond_Count[i] == 3)	{
				if(Bonded_Sp2_Atom_Restricted[i] >= 1)	{	// directly connected with at least one sp2 atom
					SP2_Flag[i] = 2;	// promoted sp2 atom
				}
			}
		}
		else if(TypeID[i] == ALPHA_TYPE_O)	{
			if(Bond_Count[i] == 2)	{
				if(Bonded_Sp2_Atom_Restricted[i] == 2)	{	// directly connected with two sp2 atoms
					SP2_Flag[i] = 2;	// promoted sp2 atom
				}
			}
		}
		else if(TypeID[i] == ALPHA_TYPE_S)	{
			if(Bond_Count[i] == 2)	{
				if(Bonded_Sp2_Atom_Restricted[i] == 2)	{	// directly connected with two sp2 atoms
					SP2_Flag[i] = 2;	// promoted sp2 atom
				}
			}
		}
	}

	for(Idx=0; Idx<n_Type; Idx++)	{
		Bonded_Sp2_Atom[Idx] = 0;
		for(i=0; i<Bond_Count[Idx]; i++)	{
			if( SP2_Flag[Bond_List[Idx][i]] >= 1)	{	// real SP2 atom + promoted SP2 atom
				Bonded_Sp2_Atom[Idx]++;
			}
		}
	}

//	printf("Miller's alpha values:\n");
	if(To_Use_Miller_Alpha)	{
		fprintf(fFile_Run_Log, "Miller's alpha values (the following values are scaled by %.2lf in fitting):\n", Scale_TargetAlpha);
		fflush(fFile_Run_Log);
	}
	n_Heavy_Atm = 0;

	for(i=0; i<n_Type; i++)	{
		if(TypeID[i])	{	// a heavy atom
			switch(TypeID[i])	{
			case ALPHA_TYPE_F:
//				AtomType[i].a_Miller = -0.296;
				AtomType[i].a_Miller = -0.266;	// 0.296 * 0.90
				break;
			case ALPHA_TYPE_CL:
//				AtomType[i].a_Miller = -2.315;
				AtomType[i].a_Miller = -1.968;	// 2.315 * 0.85
				break;
			case ALPHA_TYPE_BR:
//				AtomType[i].a_Miller = -3.013;
				AtomType[i].a_Miller = -2.410;	// 3.013 * 0.80
				break;
			case ALPHA_TYPE_I:
//				AtomType[i].a_Miller = -5.415;
				AtomType[i].a_Miller = -2.7;
				break;
			case ALPHA_TYPE_P:
//				AtomType[i].a_Miller = -2.063;
				AtomType[i].a_Miller = -1.650;	// 2.063 * 0.8
				break;
			case ALPHA_TYPE_C:
				if(Bond_Count[i] == 4)	{	// sp3
//					AtomType[i].a_Miller = -1.061 - ALPHA_H*H_Count[i];
					AtomType[i].a_Miller = -0.848 - ALPHA_H*H_Count[i];	// 1.061 * 0.80
				}
				else if(Bond_Count[i] == 3)	{	// sp2
					if(Bonded_Sp2_Atom[i]==3)	{
//						AtomType[i].a_Miller = -1.896;
						AtomType[i].a_Miller = -1.138;	// 1.896 * 0.60
					}
					else	{
//						AtomType[i].a_Miller = -1.352 - ALPHA_H*H_Count[i];
						AtomType[i].a_Miller = -0.946 - ALPHA_H*H_Count[i];	// 1.352 * 0.70
					}
				}
				else if(Bond_Count[i] == 2)	{	// sp
//					AtomType[i].a_Miller = -1.283 - ALPHA_H*H_Count[i];
					AtomType[i].a_Miller = -0.962 - ALPHA_H*H_Count[i];	// 1.283 * 0.75
				}
				else	{
					fprintf(fFile_Run_Log, "Warning> No Miller's alpha value for this C atom. Index = %d %s\n", i+1, pMol_ESP->AtomName[i]);
					fflush(fFile_Run_Log);
					AtomType[i].a_Miller = pMol_ESP->alpha[i];	// using the alpha assigned in xpsf
//					exit(1);
				}

				break;
			case ALPHA_TYPE_N:
				if(Bond_Count[i] == 4)	{
					AtomType[i].a_Miller = -1.0;	// Added by Lei
				}
				else if(Bond_Count[i] == 3)	{
					if(SP2_Flag[i] > 0)	{
//						AtomType[i].a_Miller = -1.090 - ALPHA_H*H_Count[i];
						AtomType[i].a_Miller = -0.763 - ALPHA_H*H_Count[i];	// 1.090 * 0.7
					}
					else	{
//						AtomType[i].a_Miller = -0.964 - ALPHA_H*H_Count[i];
						AtomType[i].a_Miller = -0.674 - ALPHA_H*H_Count[i];	// 0.964 * 0.7
					}
				}
				else if(Bond_Count[i] == 2)	{	// sp2
//					AtomType[i].a_Miller = -1.030 - ALPHA_H*H_Count[i];
					AtomType[i].a_Miller = -0.721 - ALPHA_H*H_Count[i];	// 1.030 * 0.70
				}
				else if(Bond_Count[i] == 1)	{	// sp
//					AtomType[i].a_Miller = -0.956 - ALPHA_H*H_Count[i];
					AtomType[i].a_Miller = -0.669 - ALPHA_H*H_Count[i];	// 0.956 * 0.70
				}
				else	{
					fprintf(fFile_Run_Log, "Warning> No Miller's alpha value for this N atom. Index = %d %s\n", i+1, pMol_ESP->AtomName[i]);
					fflush(fFile_Run_Log);
					AtomType[i].a_Miller = pMol_ESP->alpha[i];	// using the alpha assigned in xpsf
//					exit(1);
				}
				break;
			case ALPHA_TYPE_O:
				if(Bond_Count[i] == 2)	{
					if(SP2_Flag[i] > 0)	{
						AtomType[i].a_Miller = -0.274;
					}
					else	{
//						AtomType[i].a_Miller = -0.637 - ALPHA_H*H_Count[i];
						AtomType[i].a_Miller = -0.478 - ALPHA_H*H_Count[i];	// 0.637 * 0.75
					}
				}
				else if(Bond_Count[i] == 1)	{	// sp2
					if(SP2_O_Count[Bond_List[i][0]] == 2)	{	// anionic oxygen
//						AtomType[i].a_Miller = -0.858;
						AtomType[i].a_Miller = -0.644;	// 0.858 * 0.75
					}
					else	{	// regular C=O 
//						AtomType[i].a_Miller = -0.569;
						AtomType[i].a_Miller = -0.427;	// 0.569 * 0.75
					}
				}
				else	{
					fprintf(fFile_Run_Log, "Warning> No Miller's alpha value for this O atom. Index = %d %s\n", i+1, pMol_ESP->AtomName[i]);
					fflush(fFile_Run_Log);
					AtomType[i].a_Miller = pMol_ESP->alpha[i];	// using the alpha assigned in xpsf
//					exit(1);
				}
				break;
			case ALPHA_TYPE_S:
				if(Bond_Count[i] == 2)	{
					if( H_Count[i] == 2 )	{	// H2S, m231
						AtomType[i].a_Miller = -1.700;
					}
					else if( (H_Count[i] == 1) && (C_Count[i] == 1) )	{ // H-S-C
						AtomType[i].a_Miller = -1.700;
					}
					else if( C_Count[i] == 2 )	{	// C-S-C
						AtomType[i].a_Miller = -1.700;
					}
					else if( (C_Count[i] == 1) && (S_Count[i] == 1) )	{	// C-S-S
						AtomType[i].a_Miller = -1.700;
					}
					else	{
						AtomType[i].a_Miller = -1.700;
					}
				}
				else if(Bond_Count[i] == 1)	{	// sp2
					AtomType[i].a_Miller = -1.700;
				}
				else	{
					fprintf(fFile_Run_Log, "Warning> No Miller's alpha value for this S atom. Index = %d %s\n", i+1, pMol_ESP->AtomName[i]);
					fflush(fFile_Run_Log);
					AtomType[i].a_Miller = pMol_ESP->alpha[i];	// using the alpha assigned in xpsf
//					exit(1);
				}
				break;
			case ALPHA_TYPE_AR:
				AtomType[i].a_Miller = -1.000;	// Only for test !!!
				break;
			default:
				fprintf(fFile_Run_Log, "Warning> No Miller's alpha value for this atom. Index = %d %s\n", i+1, pMol_ESP->AtomName[i]);
				fflush(fFile_Run_Log);
				AtomType[i].a_Miller = pMol_ESP->alpha[i];	// using the alpha assigned in xpsf
//				exit(1);
				break;
			}
			
			HeavyAtmList[n_Heavy_Atm] = i;
			if(To_Use_Miller_Alpha)	{
				fprintf(fFile_Run_Log, "%3d atom %s, %7.3lf\n", i+1, pMol_ESP->AtomName[i], AtomType[i].a_Miller);
				fflush(fFile_Run_Log);
				Alpha_0[n_Heavy_Atm] = AtomType[i].a_Miller * Scale_TargetAlpha;	// only scale Miller's alpha used in our automatic procedure
			}
			else	{
				Alpha_0[n_Heavy_Atm] = pMol_ESP->alpha[i];	// NOT scaled
			}
			n_Heavy_Atm++;
		}
	}	


	
}
#undef ALPHA_H

#undef ALPHA_TYPE_C	
#undef ALPHA_TYPE_N	
#undef ALPHA_TYPE_O	
#undef ALPHA_TYPE_F	
#undef ALPHA_TYPE_P	
#undef ALPHA_TYPE_S	
#undef ALPHA_TYPE_CL
#undef ALPHA_TYPE_AR
#undef ALPHA_TYPE_BR
#undef ALPHA_TYPE_I	

void CFitting::DistributeData(void)
{
	int i, j, Idx_CG_Cal, nElem, Idx_Dependent, IdxHost;
	double CG;
	
//	if(nProc > 1)	{	//run paralle, syncronize all parameters
	//	MPI_Bcast(Para_List, n_Para, MPI_DOUBLE, MY_MPI_ROOT,MPI_COMM_WORLD);
//	}
	
	//start	to set up parameters for those key eliments
	for(i=0; i<n_Para_CG; i++)	{
		AtomType[Para_CG[i]].CG = pCG[i];
	}
	for(i=0; i<n_Para_Thole; i++)	{
		AtomType[Para_Thole[i]].thole = pThole[i];
	}
	for(i=0; i<n_LJ_Fit; i++)	{
		myForceField->LJ_Rec[Para_LJ[i]].para[1] = pLJ[2*i  ];	// E_Min
		myForceField->LJ_Rec[Para_LJ[i]].para[2] = pLJ[2*i+1];	// r_Min
	}
	for(i=0; i<n_NBFix_Fit; i++)	{
		myForceField->NBFix_Rec[Para_NBFix[i]].para[0] = pNBFix[2*i  ];
		myForceField->NBFix_Rec[Para_NBFix[i]].para[1] = pNBFix[2*i+1];
	}
	//end	to set up parameters for those key eliments
	
	//start	to calculte the CGs for those atoms with Eqs defined
	for(i=0; i<n_CG_Cal; i++)	{
		if(i == Last_CG_To_Cal)	{	// this one should be calculated at last
			continue;
		}
		Idx_CG_Cal = CG_Eq_List[i];
		nElem = AtomType[Idx_CG_Cal].n_Atom_Neutral;
		CG = 0.0;
		
		for(j=0; j<nElem; j++)	{
			Idx_Dependent = AtomType[Idx_CG_Cal].Cal_CG_List[j];
			CG += (AtomType[Idx_CG_Cal].weight[j] * AtomType[Idx_Dependent].CG);
		}
		AtomType[Idx_CG_Cal].CG = CG;
	}

	//start	to calculate the CG that makes the overall molecule neutral
	if(Not_Fit_Elec_Para == 0)	{
		Idx_CG_Cal = Last_CG_To_Cal;
//		i = Last_CG_To_Cal;
//		Idx_CG_Cal = CG_Eq_List[i];
		nElem = AtomType[Idx_CG_Cal].n_Atom_Neutral;
//		CG = 0.0;
		CG = netcharge_mol/SameAtom_Count;
		
		for(j=0; j<nElem; j++)	{
			Idx_Dependent = AtomType[Idx_CG_Cal].Cal_CG_List[j];
			CG += (AtomType[Idx_CG_Cal].weight[j] * AtomType[Idx_Dependent].CG);
		}
		AtomType[Idx_CG_Cal].CG = CG;
	}
	//end	to calculate the CG that makes the overall molecule neutral
	
	//end	to calculte the CGs for those atoms with Eqs defined
	
	//start	to assign CG and thole for redundant atom types
	for(i=0; i<n_Type; i++)	{
		if(AtomType[i].IsRedundant)	{
			IdxHost = AtomType[i].CG_Represent;	//the representative atom
			
			AtomType[i].CG = AtomType[IdxHost].CG;

			if(AtomType[i].IsDrudeHost)	{
				AtomType[i].thole = AtomType[IdxHost].thole;
			}
		}
		else	{
			if( (AtomType[i].NotFitThole == 1) && (AtomType[i].IsDrudeHost) )	{
				AtomType[i].thole = AtomType[AtomType[i].Thole_Represent].thole;
			}
		}
	}
	//end	to assign CG and thole for redundant atom types
	
	UpdateMolParameters();
}

void CFitting::UpdateMolParameters(void)
{
	int i, nAtom;
	CMol* pMol;
	double CG_Drude, CG_Total;

	pMol = pMol_ESP;
	//start	to update CG[] and thole[]
	for(i=0; i<n_Type; i++)	{
		pMol->CG[i] = AtomType[i].CG;
		if(AtomType[i].IsDrudeHost)	{
			pMol->thole[i] = AtomType[i].thole;
		}
	}
	//end	to update CG[] and thole[]
	
	//start	to calculate alpha from CG[] for drudes' hosts
	nAtom = pMol->nAtom;
	CG_Total = 0.0;
	
	for(i=0; i<nAtom; i++)	{
		CG_Total += pMol->CG[i];
		if(pMol->IsDrude[i])	{	//find a drude. Calculate alpha for its host
			CG_Drude = pMol->CG[i];
			pMol->alpha[i-1] = CG_Drude*CG_Drude*MD_COULOMB/(2.0*K_DRUDE);	// alpha
			if(CG_Drude < 0.0)	{
				pMol->alpha[i-1] = -(pMol->alpha[i-1]);
			}
		}
	}
	//end	to calculate alpha from CG[] for drudes' hosts

	pMol->Setup_NonBondParameters();
	pMol->Setup_TholePairParameters();

	//start	to update anisotropy parameters
	int Idx_Aniso, Idx_Aniso_Same, Idx_Para, Idx_Equiv;
	for(i=0; i<n_Aniso_Fit; i++)	{
		Idx_Aniso = Para_Aniso[i];
		Idx_Para = 2*i;

		Para_Aniso_A_K(pMol->Para_Aniso[Idx_Aniso][0],pMol->Para_Aniso[Idx_Aniso][1],pMol->Para_Aniso[Idx_Aniso][2], pAniso[Idx_Para], pAniso[Idx_Para+1]);	//aussme there is only one anisotropy term

		for(Idx_Equiv=0; Idx_Equiv<Para_Aniso_n_Equiv[i]; Idx_Equiv++)	{
			Idx_Aniso_Same = Aniso_Equiv[i][Idx_Equiv];

			pMol->Para_Aniso[Idx_Aniso_Same][0] = pMol->Para_Aniso[Idx_Aniso][0];	// copy the parameters for all equivalent anisotropy terms
			pMol->Para_Aniso[Idx_Aniso_Same][1] = pMol->Para_Aniso[Idx_Aniso][1];
			pMol->Para_Aniso[Idx_Aniso_Same][2] = pMol->Para_Aniso[Idx_Aniso][2];
		}
	}
	//end	to update anisotropy parameters

	//start	to update LP parameters
	int LPIdx, nLP;
	for(i=0; i<n_LP_Fit; i++)	{
		LPIdx = Para_LP[i];
		pMol->LP_Dist[LPIdx] = pLP[i*2  ];
		pMol->LP_Theta[LPIdx] = pLP[i*2+1];	// in radian. Degree is used in xpsf

		pMol->LP_Sin_Theta[LPIdx] = sin(pMol->LP_Theta[LPIdx]);
		pMol->LP_Cos_Theta[LPIdx] = cos(pMol->LP_Theta[LPIdx]);
	}
	if(n_LP_Fit > 0)	{
		nLP = pMol->nLP;
		for(i=0; i<nLP; i++)	{
			if(LP_Rep[i] >= 0)	{	// not a dummy LP
				if(LP_Rep[i] != i)	{	// a redundant LP, to copy the parameters from the representative LP
					pMol->LP_Dist[i] = pMol->LP_Dist[LP_Rep[i]];
					pMol->LP_Theta[i] = pMol->LP_Theta[LP_Rep[i]];
					pMol->LP_Sin_Theta[i] = pMol->LP_Sin_Theta[LP_Rep[i]];
					pMol->LP_Cos_Theta[i] = pMol->LP_Cos_Theta[LP_Rep[i]];
				}
			}
		}
	}
	//end	to update LP parameters

	if(Target_E_Int_Dimer)	{
		UpdateMolParameters_NMA_Single();
		UpdateMolParameters_NMA_Dimer();
	}

	if(Target_E_Int_Water)	{
		UpdateMolParameters_Mol_Water();
	}

	if(To_Fit_LJ)	{
		UpdateMolParameters_Mol_Ar();
	}
}

void CFitting::UpdateMolParameters_Mol_Ar(void)
{
	int i, nAtom;

	if(pMol_Ar == NULL)	{
		Quit_With_Error_Msg("From CFitting::UpdateMolParameters_Mol_Ar> pMol_Ar is not valid!\nQuit\n");
	}

	nAtom = pMol_Ar->nAtom - 2;	// substract the last two particles, Ar and DAr
	for(i=0; i<nAtom; i++)	{
		pMol_Ar->CG[i] = pMol_ESP->CG[i];
		pMol_Ar->Para_LJ_Epsilon[i] = pMol_ESP->Para_LJ_Epsilon[i];
		pMol_Ar->Para_LJ_Sigma[i] = pMol_ESP->Para_LJ_Sigma[i];
		if(pMol_Ar->IsDrude[i])	{	//copy parameters for the host atoms, alpha and thole
			pMol_Ar->alpha[i-1] = pMol_ESP->alpha[i-1];
			pMol_Ar->thole[i-1] = pMol_ESP->thole[i-1];
		}
	}

/*
	//start	to only used for Ar-Ar fitting !!!
	for(i=0; i<2; i++)	{
		pMol_Ar->CG[i+2] = pMol_ESP->CG[i];
		pMol_Ar->Para_LJ_Epsilon[i+2] = pMol_ESP->Para_LJ_Epsilon[i];
		pMol_Ar->Para_LJ_Sigma[i+2] = pMol_ESP->Para_LJ_Sigma[i];
		if(pMol_Ar->IsDrude[i+2])	{	//copy parameters for the host atoms, alpha and thole
			pMol_Ar->alpha[i-1+2] = pMol_ESP->alpha[i-1];
			pMol_Ar->thole[i-1+2] = pMol_ESP->thole[i-1];
		}
	}
	//end	to only used for Ar-Ar fitting !!!
*/

	pMol_Ar->Setup_NonBondParameters();
	pMol_Ar->Setup_TholePairParameters();

}

void CFitting::UpdateMolParameters_Mol_Water(void)
{
	int i, nAtom;
	CMol *pMol;

	pMol = pMol_Wat;
	if(pMol == NULL)	{
		Quit_With_Error_Msg("From CFitting::UpdateMolParameters_Mol_Water> pMol_Wat is not valid!\nQuit\n");
	}

	nAtom = pMol_ESP->nAtom - 1;	// skip the test charge
	for(i=0; i<nAtom; i++)	{
		pMol->CG[i] = pMol_ESP->CG[i];
		pMol->Para_LJ_Epsilon[i] = pMol_ESP->Para_LJ_Epsilon[i];
		pMol->Para_LJ_Sigma[i] = pMol_ESP->Para_LJ_Sigma[i];
		if(pMol->IsDrude[i])	{	//copy parameters for the host atoms, alpha and thole
			pMol->alpha[i-1] = pMol_ESP->alpha[i-1];
			pMol->thole[i-1] = pMol_ESP->thole[i-1];
		}
	}

	pMol->Setup_NonBondParameters();
	pMol->Setup_TholePairParameters();


	//start	to update anisotropy parameters
	int n_Aniso;
	n_Aniso = pMol_ESP->nAniso;
	for(i=0; i<n_Aniso; i++)	{		// the two PSF have to be same!!!
		pMol->Para_Aniso[i][0] = pMol_ESP->Para_Aniso[i][0];
		pMol->Para_Aniso[i][1] = pMol_ESP->Para_Aniso[i][1];
		pMol->Para_Aniso[i][2] = pMol_ESP->Para_Aniso[i][2];
	}
	//end	to update anisotropy parameters

	//start	to update LP parameters
	int nLP;
	nLP = pMol_ESP->nLP;
	for(i=0; i<nLP; i++)	{
		if(LP_Rep[i] >= 0)	{	// not a dummy LP
				pMol->LP_Dist[i] = pMol_ESP->LP_Dist[i];
				pMol->LP_Theta[i] = pMol_ESP->LP_Theta[i];
				pMol->LP_Sin_Theta[i] = pMol_ESP->LP_Sin_Theta[i];
				pMol->LP_Cos_Theta[i] = pMol_ESP->LP_Cos_Theta[i];
		}
	}
	//end	to update LP parameters

}


void CFitting::UpdateMolParameters_NMA_Single(void)
{
	int i, nAtom;
//	int n_NB_Pair, n_Thole_Pair, iPair, Atom_1, Atom_2;
//	double *pCG, *pThole, *pAlpha;

	if(pMol_ESP == NULL)	{
		Quit_With_Error_Msg("From CFitting::UpdateMolParameters_NMA_Single> pMol_ESP is not valid!\nQuit\n");
	}
	if(pNMA_Single == NULL)	{
		Quit_With_Error_Msg("From CFitting::UpdateMolParameters_NMA_Single> pNMA_Single is not valid!\nQuit\n");
	}
	
	nAtom = pNMA_Single->nAtom;
	
	for(i=0; i<nAtom; i++)	{
		pNMA_Single->CG[i] = pMol_ESP->CG[i];
		pNMA_Single->Para_LJ_Epsilon[i] = pMol_ESP->Para_LJ_Epsilon[i];
		pNMA_Single->Para_LJ_Sigma[i] = pMol_ESP->Para_LJ_Sigma[i];
		if(pNMA_Single->IsDrude[i])	{	//copy parameters for the host atoms, alpha and thole
			pNMA_Single->alpha[i-1] = pMol_ESP->alpha[i-1];
			pNMA_Single->thole[i-1] = pMol_ESP->thole[i-1];
		}
	}
	pNMA_Single->Setup_NonBondParameters();
	pNMA_Single->Setup_TholePairParameters();


	//start	to update anisotropy parameters
	int n_Aniso;
	n_Aniso = pMol_ESP->nAniso;
	for(i=0; i<n_Aniso; i++)	{		// the two PSF have to be same!!!
		pNMA_Single->Para_Aniso[i][0] = pMol_ESP->Para_Aniso[i][0];
		pNMA_Single->Para_Aniso[i][1] = pMol_ESP->Para_Aniso[i][1];
		pNMA_Single->Para_Aniso[i][2] = pMol_ESP->Para_Aniso[i][2];
	}
	//end	to update anisotropy parameters

}

#define MAX_ATM_DIMER	(1024)
void CFitting::UpdateMolParameters_NMA_Dimer(void)
{
	int i, j, nAtom_NMA, nAtom_Dimer, nAniso;
//	int n_NB_Pair, n_Thole_Pair, iPair, Atom_1, Atom_2;
//	double *pCG, *pThole, *pAlpha;
	static int FirstRun=1, IdxList[MAX_ATM_DIMER], AnisoList[MAX_ATM_DIMER];
	
	nAtom_Dimer = pNMA_Dimer->nAtom;
	nAtom_NMA = pNMA_Single->nAtom;

	if(nAtom_Dimer > MAX_ATM_DIMER)	{
//		printf("nAtom_Dimer > MAX_ATM_DIMER\nQuit\n");
		fprintf(fFile_Run_Log, "nAtom_Dimer > MAX_ATM_DIMER\nQuit\n");
		fflush(fFile_Run_Log);
		exit(1);
	}

	//start	to assign the atom type list
	if(FirstRun == 1)	{
		FirstRun = 0;

		for(i=0; i<nAtom_Dimer; i++)	{
			for(j=0; j<nAtom_NMA; j++)	{
				if(strcmp(pNMA_Dimer->AtomName[i], pNMA_Single->AtomName[j]) == 0)	{	// the same atom name. We assume same atom names mean same atom type, parameters
					break;
				}
			}
			if(j>=nAtom_NMA)	{	// fail to find the atom type in alad
//				printf("Fail to find the %3d atom type %s in monomer\nQuit\n", i+1, pNMA_Dimer->AtomName[i]);
				fprintf(fFile_Run_Log, "Fail to find the %3d atom type %s in monomer\nQuit\n", i+1, pNMA_Dimer->AtomName[i]);
				fflush(fFile_Run_Log);
				exit(1);
			}
			IdxList[i] = j;
		}

		for(i=0; i<pNMA_Dimer->nAniso; i++)	{
			for(j=0; j<pNMA_Single->nAniso; j++)	{
				if(strcmp(pNMA_Dimer->AtomName[pNMA_Dimer->AnisoList[i][0]], pNMA_Single->AtomName[pNMA_Single->AnisoList[j][0]]) == 0)	{	// the same atom name. We assume same atom names mean same atom type, parameters
					break;
				}
			}
			if(j>=pNMA_Single->nAniso)	{	// fail to find the atom type in alad
//				printf("Fail to find the %3d atom type %s in monomer\nQuit\n", pNMA_Dimer->AnisoList[i][0]+1, pNMA_Dimer->AtomName[pNMA_Dimer->AnisoList[i][0]]);
				fprintf(fFile_Run_Log, "Fail to find the %3d atom type %s in monomer\nQuit\n", pNMA_Dimer->AnisoList[i][0]+1, pNMA_Dimer->AtomName[pNMA_Dimer->AnisoList[i][0]]);
				fflush(fFile_Run_Log);
				exit(1);
			}
			AnisoList[i] = j;
		}
	}
	//end	to assign the atom type list

	for(i=0; i<nAtom_Dimer; i++)	{
		pNMA_Dimer->CG[i] = pNMA_Single->CG[IdxList[i]];
		pNMA_Dimer->Para_LJ_Epsilon[i] = pNMA_Single->Para_LJ_Epsilon[IdxList[i]];
		pNMA_Dimer->Para_LJ_Sigma[i] = pNMA_Single->Para_LJ_Sigma[IdxList[i]];
		if(pNMA_Dimer->IsDrude[i])	{	//copy parameters for the host atoms, alpha and thole
			pNMA_Dimer->alpha[i-1] = pNMA_Single->alpha[IdxList[i-1]];
			pNMA_Dimer->thole[i-1] = pNMA_Single->thole[IdxList[i-1]];
		}
//		printf("%4s %10.7lf\n", pMol_Ala5->AtomName[i], pMol_Ala5->CG[i]);
	}

	pNMA_Dimer->Setup_NonBondParameters();
	pNMA_Dimer->Setup_TholePairParameters();

	//start	to update anisotropy parameters
	int IdxAniso;
	nAniso = pNMA_Dimer->nAniso;
	for(i=0; i<nAniso; i++)	{
		IdxAniso = AnisoList[i];
		pNMA_Dimer->Para_Aniso[i][0] = pMol_ESP->Para_Aniso[IdxAniso][0];
		pNMA_Dimer->Para_Aniso[i][1] = pMol_ESP->Para_Aniso[IdxAniso][1];
		pNMA_Dimer->Para_Aniso[i][2] = pMol_ESP->Para_Aniso[IdxAniso][2];
	}
	//end	to update anisotropy parameters
}
#undef MAX_ATM_DIMER




void CFitting::ReadESPData(void)
{
	int i, j, ReadItem, IdxGrid, ReadMore=1, n_Real_Atom, AtomIdx, RealAtom_List[MAX_NUM_PARA];
	double x, y, z, ESP;
	char szName[256], szLine[256], *ReadLine;
	FILE *fIn;
	
	BuildRealAtomList(pMol_ESP, n_Real_Atom, RealAtom_List);	// include the Ar probe atom

	fprintf(fFile_Run_Log, "Entering ReadESPData()\n");
	fflush(fFile_Run_Log);
	if(pMol_ESP == NULL)	{
		fprintf(fFile_Run_Log, "pMol_ESP  == NULL\nQuit\n");
		fflush(fFile_Run_Log);
		exit(1);
	}
	N_HEAVY_ATM = pMol_ESP->nAtom - 1 - pMol_ESP->nDrude - pMol_ESP->nLP;
	//start	to read the coordinates of the grid checked in ESP
	fIn = fopen(szFile_pot0, "r");
	if(fIn == NULL)	{
		fprintf(fFile_Run_Log, "Fail to open file: %s\nQuit.\n", szFile_pot0);
		fflush(fFile_Run_Log);
		exit(1);
	}
	n_Grid = 0;
	while(1)	{
		ReadItem = fscanf(fIn, "%lf %lf %lf %lf", &x, &y, &z, &(ESP_Grid_QM_0[n_Grid]));
		if(ReadItem == 4)	{	// a valid line
			ESP_Grid_x[n_Grid] = x;
			ESP_Grid_y[n_Grid] = y;
			ESP_Grid_z[n_Grid] = z;
			
			n_Grid++;
			if(n_Grid >= N_MAX_GRID_ESP)	{
				fprintf(fFile_Run_Log, "n_Grid >= N_MAX_GRID_ESP\nQuit\n");
				fflush(fFile_Run_Log);
				exit(1);
			}
		}
		else	{
			break;
		}
	}
	fclose(fIn);
	fprintf(fFile_Run_Log, "There are %5d grid points in %s.\n", n_Grid, szFile_pot0);
	fflush(fFile_Run_Log);
	//end	to read the coordinates of the grid checked in ESP
	
	
	//start	to read the coordinates of test charge 
	n_Conf = 0;
	fIn = fopen(szFile_qpos, "r");	// read the coordinates of all configurations including test charge
	if(fIn == NULL)	{
		fprintf(fFile_Run_Log, "Fail to open file: %s\nQuit.\n", szFile_qpos);
		fflush(fFile_Run_Log);
		exit(1);
	}
	ReadMore = 1;
	while(ReadMore)	{
		if(feof(fIn))	{
			break;
		}
		for(i=0; i<n_Real_Atom; i++)	{
			ReadLine = fgets(szLine, 256, fIn);
			if(ReadLine != NULL)	{
				AtomIdx = RealAtom_List[i];
				ReadItem = sscanf(szLine, "%lf %lf %lf", &(ESP_Coord_x[n_Conf][AtomIdx]), &(ESP_Coord_y[n_Conf][AtomIdx]), &(ESP_Coord_z[n_Conf][AtomIdx]));
			}
			else	{
				ReadItem = -1;
			}
			if(ReadItem != 3)	{
				if( i==0 )	{
					ReadMore = 0;
					break;
				}
				else	{
					fclose(fIn);
					Quit_With_Error_Msg("Error in reading data file for coodinates of molecule and test charges.\nThe data for the coordinate are not complete. \n");
				}
			}
			if(pMol_ESP->IsDrude[AtomIdx+1])	{	// a drude host
				ESP_Coord_x[n_Conf][AtomIdx+1] = ESP_Coord_x[n_Conf][AtomIdx];
				ESP_Coord_y[n_Conf][AtomIdx+1] = ESP_Coord_y[n_Conf][AtomIdx];
				ESP_Coord_z[n_Conf][AtomIdx+1] = ESP_Coord_z[n_Conf][AtomIdx];
			}
		}

		if(ReadMore==0)	{
			break;
		}

		n_Conf++;
	}
	fclose(fIn);
	fprintf(fFile_Run_Log, "There are %5d configurations in %s.\n", n_Conf, szFile_qpos);
	fflush(fFile_Run_Log);
	//end	to read the coordinates of test charge 
	
	
	//start	to read QM ESP data for all configurations 
	fIn = fopen(szFile_pot, "r");
	if(fIn == NULL)	{
		fprintf(fFile_Run_Log, "Fail to open file: %s\nQuit.\n", szFile_pot);
		fflush(fFile_Run_Log);
		exit(1);
	}
	
	for(j=0; j<n_Conf; j++)	{
		for(IdxGrid=0; IdxGrid<n_Grid; IdxGrid++)	{
			ReadItem = fscanf(fIn, "%lf", &ESP);
			if(ReadItem == 1)	{	// a valid data point
				ESP_Grid_QM[j][IdxGrid] = ESP;
			}
			else	{
				fprintf(fFile_Run_Log, "Fatal error in reading file: %s\nQuit\n", szFile_pot);
				fflush(fFile_Run_Log);
				exit(1);
			}
		}
	}
	fclose(fIn);
	//end	to read QM ESP data for all configurations 
	
	
	fprintf(fFile_Run_Log, "Leaving ReadESPData()\n");
	fflush(fFile_Run_Log);

	return;
}

#define RATIO_ACTIVE_INACTIVE	(1.2)

double CFitting::Cal_Chi_SQ_ESP(void)
{
	static int FirstRun=1, ToCalculate[MAX_ATOM_TYPE];
	int i, iConf, iGrid, nAtom, iActive, nActive_NonBond, iPair, Atom_1, Atom_2;
	int Count_Grid_H_Nonor_Acceptor;
	double *x, *y, *z;	//coordinate of molecule
	double *cg;	//the charge array of atoms
	double dist_min_Active, dist_min_Inactive;
	double ESP, d_ESP, dx, dy, dz, dist, Test_Charge, Chi_SQ_Unpert=0.0, Chi_SQ_Pert=0.0;
	double Chi_SQ_Pert_Local, Chi_SQ_Pert_Buff[MAX_N_PROC];
	
	Chi_SQ_ESP = 0.0;
	
	x = pMol_ESP->x;
	y = pMol_ESP->y;
	z = pMol_ESP->z;
	cg = pMol_ESP->CG;
	
	nAtom = pMol_ESP->nAtom - 1;	//skip the last one, test charge


	if(FirstRun)	{
		FirstRun = 0;
		memset(dist_All, 0, sizeof(double)*MAX_ATOM_TYPE*N_MAX_GRID_ESP);

		for(i=0; i<nAtom; i++)	{
			if(pMol_ESP->IsDrude[i])	{	// the positions of drudes will be changed
				ToCalculate[i] = 1;
			}
			else if( pMol_ESP->IsLonePair[i] && To_Fit_LP_Para)	{	// the positions of lone pair will be changed too. 
				ToCalculate[i] = 1;
			}
			else	{
				ToCalculate[i] = 0;
			}
		}
		ToCalculate[nAtom] = 1;	// test charge. The coordinate is changed. 

		iConf=0;	// assuming all configurations have same coodinates except the position of test charge
		AssignCoordinate(pMol_ESP, ESP_Coord_x[iConf], ESP_Coord_y[iConf], ESP_Coord_z[iConf], pMol_ESP->nAtom);
		pMol_ESP->Position_LonePair();
		
		for(iGrid=0; iGrid<n_Grid; iGrid++)	{
			ESP = 0.0;
			
			for(i=0; i<nAtom; i++)	{
				if(ToCalculate[i]==0)	{
					dx = x[i] - ESP_Grid_x[iGrid];
					dy = y[i] - ESP_Grid_y[iGrid];
					dz = z[i] - ESP_Grid_z[iGrid];
					dist = sqrt(dx*dx + dy*dy +dz*dz);
					dist_All[iGrid][i] = dist;
				}
			}
		}

		//start	to set up the weight for all grid points
		w_ESP_Grid_Sum = 0.0;
		Count_Grid_H_Nonor_Acceptor = 0;
		for(iGrid=0; iGrid<n_Grid; iGrid++)	{
			dist_min_Active = dist_min_Inactive = 1.0E10;
			
			for(i=0; i<nAtom; i++)	{
				if(H_Donor_Acceptor[i])	{	// it is an atom serving as H-bond donor/acceptor
					if(dist_min_Active > dist_All[iGrid][i])	{
						dist_min_Active = dist_All[iGrid][i];
					}
				}
				else	{
					if( (dist_min_Inactive > dist_All[iGrid][i]) && (pMol_ESP->mass[i] > 0.8) )	{	// a real atom
						dist_min_Inactive = dist_All[iGrid][i];
					}
				}
			}
			
			if( (dist_min_Active <= (dist_min_Inactive*RATIO_ACTIVE_INACTIVE)) && (dist_min_Active < NEIGHBOR_CUT) )	{
				w_ESP_Grid[iGrid] = w_H_Donor_Acceptor;
				Count_Grid_H_Nonor_Acceptor++;
			}
			else	{
				w_ESP_Grid[iGrid] = 1.0;
			}
			w_ESP_Grid_Sum += w_ESP_Grid[iGrid];
		}
		
		fprintf(fFile_Run_Log, "%4d / %4d grids points are identified as the neighbors of H donor/acceptor.\n", 
			Count_Grid_H_Nonor_Acceptor, n_Grid);
		fflush(fFile_Run_Log);
		//end	to set up the weight for all grid points
	}
	

	// !!! Be careful. It is assumed that the configuration of molecules are exactly same. !!!


	//start	to calculate MM ESP for unperturbed configuration
	Test_Charge = pMol_ESP->CG[nAtom];
	pMol_ESP->CG[nAtom] = 0.0;

	nActive_NonBond = pMol_ESP->nActive_NonBond;
	for(iActive=0; iActive<nActive_NonBond; iActive++)	{	// update parameters
		iPair = pMol_ESP->Active_List_NonBond[iActive];
		Atom_1 = pMol_ESP->NB_List_i[iPair];
		Atom_2 = pMol_ESP->NB_List_j[iPair];
		pMol_ESP->Para_Elec_Pair[iPair] = MD_COULOMB * cg[Atom_1] * cg[Atom_2];
	}


//	pMol_ESP->SavePdb("before.pdb");
	AssignCoordinate(pMol_ESP, ESP_Coord_Neutral_Opt_x, ESP_Coord_Neutral_Opt_y, ESP_Coord_Neutral_Opt_z, pMol_ESP->nAtom);
	pMol_ESP->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
	pMol_ESP->Cal_Dipole_drude(MM_Dipole_x, MM_Dipole_y, MM_Dipole_z, MM_Dipole);
//	pMol_ESP->SavePdb("after.pdb");
	
	Chi_SQ_Unpert = 0.0;
	for(iGrid=0; iGrid<n_Grid; iGrid++)	{
		ESP = 0.0;
		
		for(i=0; i<nAtom; i++)	{	// unperturbed configuration, no test charge
			if(ToCalculate[i]==0)	{
				dist = dist_All[iGrid][i];
			}
			else	{
				dx = x[i] - ESP_Grid_x[iGrid];
				dy = y[i] - ESP_Grid_y[iGrid];
				dz = z[i] - ESP_Grid_z[iGrid];
				dist = sqrt(dx*dx + dy*dy +dz*dz);
			}
			
			ESP += (cg[i]/dist);
		}
		ESP_Grid_MM_0[iGrid] = ESP;
		
		d_ESP = ESP - ESP_Grid_QM_0[iGrid];

		Chi_SQ_Unpert += (d_ESP*d_ESP*w_ESP_Grid[iGrid]);
	}
	pMol_ESP->CG[nAtom] = Test_Charge;
	for(iActive=0; iActive<nActive_NonBond; iActive++)	{	// update parameters
		iPair = pMol_ESP->Active_List_NonBond[iActive];
		Atom_1 = pMol_ESP->NB_List_i[iPair];
		Atom_2 = pMol_ESP->NB_List_j[iPair];
		pMol_ESP->Para_Elec_Pair[iPair] = MD_COULOMB * cg[Atom_1] * cg[Atom_2];
	}
	//end	to calculate MM ESP for unperturbed configuration


	Chi_SQ_Pert = 0.0;
	Chi_SQ_Pert_Local = 0.0;

	for(iConf=0; iConf<n_Conf; iConf++)	{
		if(iConf%nProc != ProgID)	{	//only proceed when iConf%nProc == ProgID, do job decomposition
			continue;
		}
		AssignCoordinate(pMol_ESP, ESP_Coord_Opt_x[iConf], ESP_Coord_Opt_y[iConf], ESP_Coord_Opt_z[iConf], pMol_ESP->nAtom);

		pMol_ESP->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
//		pMol_ESP->SavePdb("tmp.pdb");
		
		//start	to accumulate ESP difference between QM and MM
		for(iGrid=0; iGrid<n_Grid; iGrid++)	{
			ESP = 0.0;
			
			for(i=0; i<nAtom; i++)	{
				if(ToCalculate[i]==0)	{
					dist = dist_All[iGrid][i];
				}
				else	{
					dx = x[i] - ESP_Grid_x[iGrid];
					dy = y[i] - ESP_Grid_y[iGrid];
					dz = z[i] - ESP_Grid_z[iGrid];
					dist = sqrt(dx*dx + dy*dy +dz*dz);

				}

				ESP += (cg[i]/dist);
			}
			
			d_ESP = ( (ESP - ESP_Grid_QM[iConf][iGrid]));

			Chi_SQ_Pert_Local += (d_ESP*d_ESP*w_ESP_Grid[iGrid]);
		}
		//end	to accumulate ESP difference between QM and MM
	}

	Chi_SQ_Pert = Chi_SQ_Pert_Local;
/*
	MPI_Gather(&Chi_SQ_Pert_Local,1,MPI_DOUBLE,Chi_SQ_Pert_Buff,1,MPI_DOUBLE,MY_MPI_ROOT,MPI_COMM_WORLD);
	if(ProgID == MY_MPI_ROOT)	{
		Chi_SQ_Pert = 0.0;
		for(i=0; i<nProc; i++)	{
			Chi_SQ_Pert += Chi_SQ_Pert_Buff[i];	//calculate the total Chi_SQ_ESP
		}
	}
	MPI_Bcast(&Chi_SQ_Pert, 1, MPI_DOUBLE, MY_MPI_ROOT,MPI_COMM_WORLD);
*/	

	Chi_SQ_ESP = Chi_SQ_Unpert*n_Conf + Chi_SQ_Pert;

	
	return (Chi_SQ_ESP*10000.0/(w_ESP_Grid_Sum*n_Conf*2.0));
}
#undef	RATIO_ACTIVE_INACTIVE

int CFitting::To_Find_Stable_Parameters(void)
{
	int i, Idx, ExistFarDrude;
	double R_Drude_Max;

	Idx = (iActiveSave-1+MAX_PARA_SAVE)%MAX_PARA_SAVE;	// the last saved one

	for(i=0; i<MAX_PARA_SAVE; i++)	{	// to check all saved parameters
		if(Is_Saved_Para[Idx] == 0)	{	// not a valid parameter set
			break;
		}
		if( ((Para_Set_Saved_Chi2[Idx]-Chi_SQ_Min) / Chi_SQ_Min) > 0.03)	{	// the chi^2 should be close
			break;
		}

		memcpy(Para_List, Para_Set_Save[Idx], sizeof(double)*n_Para);
		DistributeData();
		ExistFarDrude = IsThereFarDrudes(R_Drude_Max);	// to check the current parameter set

		if(ProgID == MY_MPI_ROOT)	{
			fprintf(fFile_Run_Log, "To_Find_Stable_Parameters> Idx = %3d   ExistFarDrude = %3d  chi^2 = %.6lf\n", Idx, ExistFarDrude, Para_Set_Saved_Chi2[Idx]);
		}
		
		if(ExistFarDrude == 0)	{	// a stable parameter set. To output this selected parameter set
			if(i == 0)	{	// checking the current best parameters
				if(ProgID == MY_MPI_ROOT)	{
					fprintf(fFile_Run_Log, "\nThe current parameters are stable, chi^2_min= %.6lf chi^2=%.6lf\nR_Drude_Max = %.4lf\n\n", Chi_SQ_Min, Para_Set_Saved_Chi2[Idx], R_Drude_Max);
				}
				return 1;
			}

			PreOptimizeConfigurations();

			if(Target_E_Int_Water)	{
				Generate_Water_Confs_Acceptor();
				Generate_Water_Confs_Donor();
			}


			CalObjectiveFunction();
			SaveParameters();
			SaveParameters_CHARMM();

			if(ProgID == MY_MPI_ROOT)	{
				fprintf(fFile_Run_Log, "\nReplace the current parameters with a more stable one, chi^2_min= %.6lf chi^2=%.6lf\nR_Drude_Max = %.4lf\n\n", Chi_SQ_Min, Chi_SQ, R_Drude_Max);
			}

			return 1;
		}
		
		Idx = (Idx-1+MAX_PARA_SAVE)%MAX_PARA_SAVE;	// move to the previous best
	}

	if(ProgID == MY_MPI_ROOT)	{
		fprintf(fFile_Run_Log, "\nFail to find a set of stable parameters. The alpha will be scaled down.\n\n");
	}

	return 0;	// to rescale alpha
}


#define X_PERT	(0.09)
int CFitting::IsThereFarDrudes(double& R_SQ_MAX)
{
	int i, nAtom, nActive_NonBond, iPair, Atom_1, Atom_2, iActive, bFarDrude=0, IdxMax;
	double Test_Charge, *cg, dx, dy, dz, R_SQ, R_SQ_CUT;

	R_SQ_MAX = 0.0;
	iseed += (ProgID*ProgID);

	R_SQ_CUT = FAR_DRUDE_R * FAR_DRUDE_R;
	nAtom = pMol_ESP->nAtom;
	cg = pMol_ESP->CG;

	//start	to calculate MM ESP for unperturbed configuration
	Test_Charge = pMol_ESP->CG[nAtom-1];
	pMol_ESP->CG[nAtom-1] = 0.0;

	nActive_NonBond = pMol_ESP->nActive_NonBond;
	for(iActive=0; iActive<nActive_NonBond; iActive++)	{	// update parameters
		iPair = pMol_ESP->Active_List_NonBond[iActive];
		Atom_1 = pMol_ESP->NB_List_i[iPair];
		Atom_2 = pMol_ESP->NB_List_j[iPair];
		pMol_ESP->Para_Elec_Pair[iPair] = MD_COULOMB * cg[Atom_1] * cg[Atom_2];
	}


	pMol_ESP->nAtom--;	// remove the last atom, test charge
	IdxMax = pMol_ESP->nAtom-1;
	for(int Check=0; Check<30; Check++)	{
		AssignCoordinate(pMol_ESP, ESP_Coord_Neutral_Opt_x, ESP_Coord_Neutral_Opt_y, ESP_Coord_Neutral_Opt_z, pMol_ESP->nAtom);

		//start	to perturb the positions of drudes to check the stabilities
		for(Atom_1=0; Atom_1<IdxMax; Atom_1++)	{
			if(pMol_ESP->mass[Atom_1] > 1.2)	{	// a heavy atom
				Atom_2 = Atom_1 + 1;
				pMol_ESP->x[Atom_2] = pMol_ESP->x[Atom_1] + X_PERT*(rand(iseed)-0.5);
				pMol_ESP->y[Atom_2] = pMol_ESP->y[Atom_1] + X_PERT*(rand(iseed)-0.5);
				pMol_ESP->z[Atom_2] = pMol_ESP->z[Atom_1] + X_PERT*(rand(iseed)-0.5);
			}
		}
		//end	to perturb the positions of drudes to check the stabilities

		pMol_ESP->FullGeometryOptimization_LBFGS_drude(0);	//let all particles relax. Only optimize the geometry of the molecule

		//start	to check the bond lengths between heavy atoms and their drudes
		for(Atom_1=0; Atom_1<IdxMax; Atom_1++)	{
			if(pMol_ESP->mass[Atom_1] > 1.2)	{	// a heavy atom
				Atom_2 = Atom_1 + 1;
				dx = pMol_ESP->x[Atom_1] - pMol_ESP->x[Atom_2];
				dy = pMol_ESP->y[Atom_1] - pMol_ESP->y[Atom_2];
				dz = pMol_ESP->z[Atom_1] - pMol_ESP->z[Atom_2];
				R_SQ = dx*dx + dy*dy + dz*dz;
				if(R_SQ > R_SQ_MAX)	{
					R_SQ_MAX = R_SQ;
				}
				if(R_SQ >= R_SQ_CUT)	{
					bFarDrude = 1;
					break;
				}
			}
		}
		//end	to check the bond lengths between heavy atoms and their drudes

		if(bFarDrude)	{
			break;
		}
	}

	fprintf(fFile_Run_Log, "\nR_MAX_Drude = %.4lf  ProcID = %d\n\n", sqrt(R_SQ_MAX), ProgID);
//	if(bFarDrude)	{
//		fprintf(fFile_Run_Log, "The electrostatic parameters will be refitted by scaling down target polarizability!\n\n");
//	}

	pMol_ESP->nAtom++;	// restore the last atom, test charge

	
	pMol_ESP->CG[nAtom-1] = Test_Charge;
	for(iActive=0; iActive<nActive_NonBond; iActive++)	{	// update parameters
		iPair = pMol_ESP->Active_List_NonBond[iActive];
		Atom_1 = pMol_ESP->NB_List_i[iPair];
		Atom_2 = pMol_ESP->NB_List_j[iPair];
		pMol_ESP->Para_Elec_Pair[iPair] = MD_COULOMB * cg[Atom_1] * cg[Atom_2];
	}
	//end	to calculate MM ESP for unperturbed configuration

/*
	int bFarDrude_List[MAX_N_PROC];
	double R_SQ_MAX_List[MAX_N_PROC];
	MPI_Gather(&bFarDrude,1,MPI_INT,bFarDrude_List,1,MPI_INT,MY_MPI_ROOT,MPI_COMM_WORLD);
	MPI_Gather(&R_SQ_MAX,1,MPI_DOUBLE,R_SQ_MAX_List,1,MPI_DOUBLE,MY_MPI_ROOT,MPI_COMM_WORLD);
	if(ProgID == MY_MPI_ROOT)	{
		bFarDrude = 0;
		for(i=0; i<nProc; i++)	{
			bFarDrude += bFarDrude_List[i];	
		}
		R_SQ_MAX = 0.0;
		for(i=0; i<nProc; i++)	{
			if(R_SQ_MAX < R_SQ_MAX_List[i])	{
				R_SQ_MAX = R_SQ_MAX_List[i];
			}
		}
*/
		R_SQ_MAX = sqrt(R_SQ_MAX);


		fprintf(fFile_Run_Log, "\nR_MAX_Drude_All = %.4lf  ProcID = %d\n\n", R_SQ_MAX, ProgID);
//	}
/*
	MPI_Bcast(&bFarDrude, 1, MPI_INT, MY_MPI_ROOT,MPI_COMM_WORLD);
	MPI_Bcast(&R_SQ_MAX, 1, MPI_DOUBLE, MY_MPI_ROOT,MPI_COMM_WORLD);
*/
	return bFarDrude;
}
#undef X_PERT

void CFitting::PreOptimizeConfigurations(void)
{
	double *x, *y, *z, *cg, Test_Charge;
	int nAtom, nActive_NonBond, iActive, iConf, iPair, Atom_1, Atom_2;
	static int First;

	x = pMol_ESP->x;
	y = pMol_ESP->y;
	z = pMol_ESP->z;
	cg = pMol_ESP->CG;
	nAtom = pMol_ESP->nAtom - 1;	//skip the last one, test charge

	// !!! Be careful. It is assumed that the configuration of molecules are exactly same. !!!

	//start	to calculate MM ESP for unperturbed configuration
	AssignCoordinate(pMol_ESP, ESP_Coord_x[0], ESP_Coord_y[0], ESP_Coord_z[0], pMol_ESP->nAtom);
	pMol_ESP->Position_LonePair();
//	pMol_ESP->SavePdb("1.pdb");
	Test_Charge = pMol_ESP->CG[nAtom];
	pMol_ESP->CG[nAtom] = 0.0;

	nActive_NonBond = pMol_ESP->nActive_NonBond;
	for(iActive=0; iActive<nActive_NonBond; iActive++)	{	// update parameters
		iPair = pMol_ESP->Active_List_NonBond[iActive];
		Atom_1 = pMol_ESP->NB_List_i[iPair];
		Atom_2 = pMol_ESP->NB_List_j[iPair];
		pMol_ESP->Para_Elec_Pair[iPair] = MD_COULOMB * cg[Atom_1] * cg[Atom_2];
	}

	if(First)	{
		First = 0;
		pMol_ESP->SavePdb("opt-qm.pdb");
	}

	pMol_ESP->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
	CollectCoordinate(pMol_ESP, ESP_Coord_Neutral_Opt_x, ESP_Coord_Neutral_Opt_y, ESP_Coord_Neutral_Opt_z, pMol_ESP->nAtom);

//	pMol_ESP->SavePdb("2.pdb");

	if(Target_E_Int_Water)	{	// to do fully geometry optimization
		pMol_ESP->Geo_Opt_Drude_Only = 0;

		pMol_ESP->Restrain_All_Torsions(1);
		pMol_ESP->nAtom--;
		pMol_ESP->FullGeometryOptimization_LBFGS_drude(0);
		pMol_ESP->nAtom++;
		pMol_ESP->Restrain_All_Torsions(0);
		CollectCoordinate(pMol_ESP, x_Mol_Full_Opt, y_Mol_Full_Opt, z_Mol_Full_Opt, pMol_ESP->nAtom);	// to be used for mol-water interactions
		CollectCoordinate(pMol_ESP, pMol_Wat->x, pMol_Wat->y, pMol_Wat->z, pMol_ESP->nAtom);	// to be used for mol-water interactions

//		pMol_ESP->SavePdb("3.pdb");


		pMol_ESP->Geo_Opt_Drude_Only = 1;	// to restore the flag
	}

	pMol_ESP->CG[nAtom] = Test_Charge;
	for(iActive=0; iActive<nActive_NonBond; iActive++)	{	// update parameters
		iPair = pMol_ESP->Active_List_NonBond[iActive];
		Atom_1 = pMol_ESP->NB_List_i[iPair];
		Atom_2 = pMol_ESP->NB_List_j[iPair];
		pMol_ESP->Para_Elec_Pair[iPair] = MD_COULOMB * cg[Atom_1] * cg[Atom_2];
	}
	//end	to calculate MM ESP for unperturbed configuration


	for(iConf=0; iConf<n_Conf; iConf++)	{
		if(iConf%nProc != ProgID)	{	//only proceed when iConf%nProc == ProgID, do job decomposition
			continue;
		}
		AssignCoordinate(pMol_ESP, ESP_Coord_x[iConf], ESP_Coord_y[iConf], ESP_Coord_z[iConf], pMol_ESP->nAtom);
		pMol_ESP->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
		CollectCoordinate(pMol_ESP, ESP_Coord_Opt_x[iConf], ESP_Coord_Opt_y[iConf], ESP_Coord_Opt_z[iConf], pMol_ESP->nAtom);
	}

}

void CFitting::CollectCoordinate(CMol *pMol, double x[], double y[], double z[], int n)
{
	int i;
	double *px, *py, *pz;

	px = pMol->x;
	py = pMol->y;
	pz = pMol->z;

	for(i=0; i<n; i++)	{
		x[i] = px[i];
		y[i] = py[i];
		z[i] = pz[i];
	}
}


void CFitting::AssignCoordinate(CMol *pMol, double x[], double y[], double z[], int n)
{
	int i;
	double *px, *py, *pz;

	px = pMol->x;
	py = pMol->y;
	pz = pMol->z;

	for(i=0; i<n; i++)	{
		px[i] = x[i];
		py[i] = y[i];
		pz[i] = z[i];
	}
}

void CFitting::MergeDrudeForces(CMol *pMol, double fx[], double fy[], double fz[])
{
	int i, nAtom, *IsDrude, *IsDrudeHost;
	double *gx, *gy, *gz;

	nAtom = pMol->nAtom;
	gx = pMol->grad_x;
	gy = pMol->grad_y;
	gz = pMol->grad_z;
	IsDrude = pMol->IsDrude;
	IsDrudeHost = pMol->IsDrudeHost;

	for(i=0; i<nAtom; i++)	{
		if(IsDrudeHost[i])	{	// merge the forces on host and drude
			fx[i] = -(gx[i] + gx[i+1]);
			fy[i] = -(gy[i] + gy[i+1]);
			fz[i] = -(gz[i] + gz[i+1]);
		}
		else if(IsDrude[i])	{	// set force as zero
			fx[i] = 0.0;
			fy[i] = 0.0;
			fz[i] = 0.0;
		}
		else	{
			fx[i] = -gx[i];
			fy[i] = -gy[i];
			fz[i] = -gz[i];
		}
	}
}

#define N_ATM_NMA	(19)
double CFitting::Cal_Chi_SQ_E_Int_Dimer(void)
{
	double E_A, E_B, E_A_B, d_E;
	int iConf, i;

	Chi_SQ_E_int = 0.0;

	for(iConf=0; iConf<N_CONF_DIMER; iConf++)	{

		//start	to set up the first NMA coodinates
		for(i=0; i<N_ATM_NMA; i++)	{
			pNMA_Single->x[i] = x_Dimer[iConf][i];
			pNMA_Single->y[i] = y_Dimer[iConf][i];
			pNMA_Single->z[i] = z_Dimer[iConf][i];
		}
		//end	to set up the first NMA coodinates
		pNMA_Single->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
		pNMA_Single->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
		pNMA_Single->Geo_Opt_Drude_Only = 0;	//only drude will be relaxed
		E_A = pNMA_Single->Cal_E(0);

		//start	to set up the second NMA coodinates
		for(i=0; i<N_ATM_NMA; i++)	{
			pNMA_Single->x[i] = x_Dimer[iConf][i+N_ATM_NMA];
			pNMA_Single->y[i] = y_Dimer[iConf][i+N_ATM_NMA];
			pNMA_Single->z[i] = z_Dimer[iConf][i+N_ATM_NMA];
		}
		//end	to set up the second NMA coodinates
		pNMA_Single->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
		pNMA_Single->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
		pNMA_Single->Geo_Opt_Drude_Only = 0;	//only drude will be relaxed
		E_B = pNMA_Single->Cal_E(0);

		//start	to set up the NMA dimer coodinates
		for(i=0; i<N_ATM_DIMER; i++)	{
			pNMA_Dimer->x[i] = x_Dimer[iConf][i];
			pNMA_Dimer->y[i] = y_Dimer[iConf][i];
			pNMA_Dimer->z[i] = z_Dimer[iConf][i];
		}
		//end	to set up the NMA dimer coodinates
		pNMA_Dimer->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
		pNMA_Dimer->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
		pNMA_Dimer->Geo_Opt_Drude_Only = 0;	//only drude will be relaxed
		E_A_B = pNMA_Dimer->Cal_E(0);

		E_Int_MM_Dimer[iConf] = E_A_B - E_A - E_B;
		d_E = E_Int_MM_Dimer[iConf] - E_Int_QM_Dimer[iConf];
		Chi_SQ_E_int += (d_E*d_E);
	}

	return Chi_SQ_E_int;
}
#undef N_ATM_NMA

/*
double CFitting::Cal_Chi_SQ_E_Int_Dimer(void)
{
	double E_A_B_Far, E_A_B, d_E;
	int iConf, i;

	Chi_SQ_E_int = 0.0;

	for(iConf=0; iConf<N_CONF_DIMER; iConf++)	{
		//start	to set up the NMA dimer coodinates
		for(i=0; i<N_ATM_DIMER; i++)	{
			pNMA_Dimer->x[i] = x_Dimer[iConf][i];
			pNMA_Dimer->y[i] = y_Dimer[iConf][i];
			pNMA_Dimer->z[i] = z_Dimer[iConf][i];
		}
		//end	to set up the NMA dimer coodinates
		pNMA_Dimer->BackupCoordinates();

		pNMA_Dimer->TranslateLastSegment();	// move the last segment (the other molecule) far
		pNMA_Dimer->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
		pNMA_Dimer->FullGeometryOptimization_LBFGS_drude();	//let the drudes relax
		pNMA_Dimer->Geo_Opt_Drude_Only = 0;	//only drude will be relaxed
		E_A_B_Far = pNMA_Dimer->Cal_E(0);

		pNMA_Dimer->TranslateLastSegmentBack();	// move the last segment (the other molecule) far
		pNMA_Dimer->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
		pNMA_Dimer->FullGeometryOptimization_LBFGS_drude();	//let the drudes relax
		pNMA_Dimer->Geo_Opt_Drude_Only = 0;	//only drude will be relaxed
		E_A_B = pNMA_Dimer->Cal_E(0);

		E_Int_MM_Dimer[iConf] = E_A_B - E_A_B_Far;
		d_E = E_Int_MM_Dimer[iConf] - E_Int_QM_Dimer[iConf];
		Chi_SQ_E_int += (d_E*d_E);
	}

	return Chi_SQ_E_int;
}
*/


void CFitting::BuildRealAtomList(CMol *pMol, int& LenList, int *List)
{
	int i, nAtom, *IsDrude, *IsLP;

	nAtom = pMol->nAtom;
	IsDrude = pMol->IsDrude;
	IsLP = pMol->IsLonePair;

	LenList = 0;

	for(i=0; i<nAtom; i++)	{
		if(IsDrude[i])	{	// drude
		}
		else if(IsLP[i])	{	// lone pair
		}
		else	{
			List[LenList] = i;
			LenList ++;
		}
	}
}

double CFitting::Cal_Chi_SQ_Int_Energy_Mol_Ar(void)
{
	double E_A_B_Far, E_A_B, d_E;
	int iConf;
	int LET_DRUDE_RELAX = 1;
	
	Chi_SQ_E_int = 0.0;


	for(iConf=0; iConf<nConf_Mol_Ar; iConf++)	{
//	for(iConf=400; iConf<nConf_Mol_Ar; iConf++)	{
		if(Used_Ar_Pert[iConf] == 0)	{
			continue;
		}

		AssignCoordinate(pMol_Ar, Ar_Pert_x[iConf], Ar_Pert_y[iConf], Ar_Pert_z[iConf], pMol_Ar->nAtom);
		pMol_Ar->BackupCoordinates();
		
		pMol_Ar->TranslateLastSegment();	// move the last segment (the other molecule) far
		pMol_Ar->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
		pMol_Ar->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
		pMol_Ar->Geo_Opt_Drude_Only = 0;	//only drude will be relaxed
		E_A_B_Far = pMol_Ar->Cal_E(0);
		
//		pMol_Ar->WriteCRDFile("A_B_far.crd");
		
		
		
		pMol_Ar->TranslateLastSegmentBack();	// move the last segment (the other molecule) far
		pMol_Ar->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
		pMol_Ar->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
		pMol_Ar->Geo_Opt_Drude_Only = 0;	//only drude will be relaxed
		E_A_B = pMol_Ar->Cal_E(0);
		
		
//		pMol_Ar->WriteCRDFile("A_B.crd");
		
		
		
		E_Ar_Pert_MM[iConf] = E_A_B - E_A_B_Far;
		d_E = E_Ar_Pert_MM[iConf] - E_Ar_Pert_QM[iConf];
		Chi_SQ_E_int += (d_E*d_E);
	}
	
	return Chi_SQ_E_int;
}


double CFitting::Cal_Chi_SQ_Int_Energy_Mol_Water(int OutputFlag)
{
	int iDonor, iAcceptor, nAtom;
	int LET_DRUDE_RELAX = 1;
	double dE, dR, E_Mol_Water_Far;

	Chi_SQ_Wat_Emin = Chi_SQ_Wat_Rmin = 0.0;

	nAtom = pMol_ESP->nAtom;

	AssignCoordinate(pMol_Wat, x_Mol_Full_Opt, y_Mol_Full_Opt, z_Mol_Full_Opt, nAtom-1);	// set up optimized geometry. Ignore the effect of parameter change in numerical derivative on the geometry of the compound

	//start	to determine the energy when the compound and water are far
	if(nFit_Acceptor !=0)	{
		iAcceptor=0;
		while(1)	{
			if(ToFit_Acceptor[iAcceptor] != 0)	{
				break;
			}
		}
		Mol_Active_Atom = Acceptor_List_Acceptor[iAcceptor];
		OptWater.IsAcceptor = 1;
		OptWater.Generate_Dummy_Atoms(Mol_Active_Atom, 1, theta_Acceptor[iAcceptor]);
		OptWater.Para_List[0] = 2.5;	OptWater.Para_List[1] = 0.0;	// the orientation is NOT important at all!
		OptWater.GenerateWater();
//		pMol_Wat->SavePdb("tmp.pdb");
	}
	else	{
		iDonor=0;
		while(1)	{
			if(ToFit_Donor[iDonor] != 0)	{
				break;
			}
		}
		Mol_Active_Atom = Donor_List_Active[iDonor];
		OptWater.IsAcceptor = 0;
		OptWater.Generate_Dummy_Atoms(Mol_Active_Atom, 0, theta_Donor[iDonor]);
		OptWater.Para_List[0] = 2.5;	OptWater.Para_List[1] = 0.0;	// the orientation is NOT important at all!
		OptWater.GenerateWater();
//		pMol_Wat->SavePdb("tmp.pdb");
	}

	pMol_Wat->TranslateLastSegment();	// move the last segment (the other molecule) far
//	pMol_Wat->SavePdb("tmp.pdb");
	pMol_Wat->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
	pMol_Wat->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
	pMol_Wat->Geo_Opt_Drude_Only = 0;
	E_Mol_Water_Far = pMol_Wat->Cal_E(0);
	//end	to determine the energy when the compound and water are far


	for(iAcceptor=0; iAcceptor<N_Acceptor; iAcceptor++)	{
		if(ToFit_Acceptor[iAcceptor] == 0)	{
			continue;
		}
		Mol_Active_Atom = Acceptor_List_Acceptor[iAcceptor];
		OptWater.Generate_Dummy_Atoms(Mol_Active_Atom, 1, theta_Acceptor[iAcceptor]);
		OptWater.Optimize_Water_Pose(Para_Save_Acceptor[iAcceptor], E_Min_MM_Acceptor[iAcceptor], r_Min_MM_Acceptor[iAcceptor]);
		E_Min_MM_Acceptor[iAcceptor] -= E_Mol_Water_Far;
		
//		OutputFlag = 1;
		if(OutputFlag)	{
			char szName[256];
			sprintf(szName, "mm-acceptor-%d.pdb", iAcceptor+1);
			pMol_Wat->SavePdb(szName);
		}
//		OutputFlag = 0;

		dE = E_Min_MM_Acceptor[iAcceptor] - E_Min_QM_Acceptor[iAcceptor];
		dR = r_Min_MM_Acceptor[iAcceptor] - r_Min_QM_Acceptor[iAcceptor];

		Chi_SQ_Wat_Emin += (dE*dE);
		Chi_SQ_Wat_Rmin += (dR*dR);
	}

	for(iDonor=0; iDonor<N_Donor; iDonor++)	{
		if(ToFit_Donor[iDonor] == 0)	{
			continue;
		}
		Mol_Active_Atom = Donor_List_Active[iDonor];

		OptWater.Generate_Dummy_Atoms(Mol_Active_Atom, 0, theta_Donor[iDonor]);
		OptWater.Optimize_Water_Pose(Para_Save_Donor[iDonor], E_Min_MM_Donor[iDonor], r_Min_MM_Donor[iDonor]);
		E_Min_MM_Donor[iDonor] -= E_Mol_Water_Far;
		
//		OutputFlag = 1;
		if(OutputFlag)	{
			char szName[256];
			sprintf(szName, "mm-donor-%d.pdb", iDonor+1);
			pMol_Wat->SavePdb(szName);
		}
//		OutputFlag = 0;

		dE = E_Min_MM_Donor[iDonor] - E_Min_QM_Donor[iDonor];
		dR = r_Min_MM_Donor[iDonor] - r_Min_QM_Donor[iDonor];

		Chi_SQ_Wat_Emin += (dE*dE);
		Chi_SQ_Wat_Rmin += (dR*dR);
	}

	return ( (Chi_SQ_Wat_Emin*w_water_E_min + Chi_SQ_Wat_Rmin*w_water_R_min)/(nFit_Acceptor+nFit_Donor));
}


/*
double CFitting::Cal_Chi_SQ_Int_Energy_Mol_Water(void)
{
	double E_A_B_Far, E_A_B, d_E;
	int iConf;
	int LET_DRUDE_RELAX = 1;

	Chi_SQ_E_int = 0.0;

	for(iConf=0; iConf<nConf_Mol_Wat; iConf++)	{
		if(Used_Mol_Wat_Int[iConf] == 0)	{
			continue;
		}

		AssignCoordinate(pMol_Wat, Coord_Dimer_Wat_x[iConf], Coord_Dimer_Wat_y[iConf], Coord_Dimer_Wat_z[iConf], pMol_Wat->nAtom);
		pMol_Wat->BackupCoordinates();

		pMol_Wat->TranslateLastSegment();	// move the last segment (the other molecule) far
		pMol_Wat->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
		pMol_Wat->FullGeometryOptimization_LBFGS_drude();	//let the drudes relax
		pMol_Wat->Geo_Opt_Drude_Only = 0;	//only drude will be relaxed
		E_A_B_Far = pMol_Wat->Cal_E(0);

//		pMol_Wat->WriteCRDFile("A_B_far.crd");

		pMol_Wat->TranslateLastSegmentBack();	// move the last segment (the other molecule) far
		pMol_Wat->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
		pMol_Wat->FullGeometryOptimization_LBFGS_drude();	//let the drudes relax
		pMol_Wat->Geo_Opt_Drude_Only = 0;	//only drude will be relaxed
		E_A_B = pMol_Wat->Cal_E(0);


//		pMol_Wat->WriteCRDFile("A_B.crd");

		E_Int_Wat_MM[iConf] = E_A_B - E_A_B_Far;
		d_E = E_Int_Wat_MM[iConf] - E_Int_Wat_QM[iConf];
		Chi_SQ_E_int += (d_E*d_E);
	}

	return (Chi_SQ_E_int*w_water_int_E/nUsed_Conf_Mol_Wat_Int);
}
*/

double CFitting::Cal_Chi_SQ_Force_Pert_Ar(void)
{
	double f_Ref_x[MAX_NUM_PARA], f_Ref_y[MAX_NUM_PARA], f_Ref_z[MAX_NUM_PARA];
	double f_Merged_x[MAX_NUM_PARA], f_Merged_y[MAX_NUM_PARA], f_Merged_z[MAX_NUM_PARA];
	double d_fx, d_fy, d_fz;
	int nAtom, Real_Atom_List[MAX_NUM_PARA], n_Real_Atom;
	int i, iConf, IdxAtom;
	int LET_DRUDE_RELAX = 0;
	static int First_Run_Force_Pert_Ar=1;

	if(First_Run_Force_Pert_Ar)	{
		First_Run_Force_Pert_Ar = 0;
		nAtom = pMol_Ar->nAtom;

		for(iConf=0; iConf<=i_Ref_Conf; iConf++)	{
			if(Used_Ar_Pert[iConf] == 0)	{	// not used for fitting
				continue;
			}
			
			AssignCoordinate(pMol_Ar, Ar_Pert_x[iConf], Ar_Pert_y[iConf], Ar_Pert_z[iConf], pMol_Ar->nAtom);
			pMol_Ar->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
			pMol_Ar->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
			
			for(i=0; i<nAtom; i++)	{	// save the optimized coordinates
				Ar_Pert_x[iConf][i] = pMol_Ar->x[i];
				Ar_Pert_y[iConf][i] = pMol_Ar->y[i];
				Ar_Pert_z[iConf][i] = pMol_Ar->z[i];
			}
		}
	}

	Chi_SQ_f_Pert = 0.0;

	//start	to build a list of real atoms
	BuildRealAtomList(pMol_Ar, n_Real_Atom, Real_Atom_List);
	n_Real_Atom --;	// to skip the last Ar, probe atom

	//start	to calculate the MM force for the reference configuration
	AssignCoordinate(pMol_Ar, Ar_Pert_x[i_Ref_Conf], Ar_Pert_y[i_Ref_Conf], Ar_Pert_z[i_Ref_Conf], pMol_Ar->nAtom);

	if(LET_DRUDE_RELAX)	{
		pMol_Ar->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
		pMol_Ar->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
	}
	pMol_Ar->Geo_Opt_Drude_Only = 0;	//only drude will be relaxed

	pMol_Ar->Cal_E(0);
	MergeDrudeForces(pMol_Ar, f_Ref_x, f_Ref_y, f_Ref_z);
	//end	to calculate the MM force for the reference configuration

	for(iConf=0; iConf<i_Ref_Conf; iConf++)	{
		if(Used_Ar_Pert[iConf] == 0)	{	// not used for fitting
			continue;
		}

		AssignCoordinate(pMol_Ar, Ar_Pert_x[iConf], Ar_Pert_y[iConf], Ar_Pert_z[iConf], pMol_Ar->nAtom);
		if(LET_DRUDE_RELAX)	{
			pMol_Ar->Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
			pMol_Ar->FullGeometryOptimization_LBFGS_drude(1);	//let the drudes relax
			pMol_Ar->Geo_Opt_Drude_Only = 0;	//only drude will be relaxed
		}

		pMol_Ar->Cal_E(0);
		MergeDrudeForces(pMol_Ar, f_Merged_x, f_Merged_y, f_Merged_z);

		for(i=0; i<n_Real_Atom; i++)	{
			IdxAtom = Real_Atom_List[i];
			d_fx = f_Merged_x[IdxAtom] - f_Ref_x[IdxAtom] - f_Pert_x_QM[iConf][IdxAtom];
			d_fy = f_Merged_y[IdxAtom] - f_Ref_y[IdxAtom] - f_Pert_y_QM[iConf][IdxAtom];
			d_fz = f_Merged_z[IdxAtom] - f_Ref_z[IdxAtom] - f_Pert_z_QM[iConf][IdxAtom];
			Chi_SQ_f_Pert += (d_fx*d_fx + d_fy*d_fy + d_fz*d_fz);
		}
	}


	return Chi_SQ_f_Pert;
}

double CFitting::Cal_Chi_SQ_User(void)	// the user-defined CG constraints
{
	int i, jMax, j;
	double q, dq;

	Chi_SQ_User = 0.0;
	
	for(i=0; i<n_User_CG; i++)	{
		jMax = User_CG_Count[i];
		q = 0.0;
		for(j=0; j<jMax; j++)	{
			q += pMol_ESP->CG[User_CG_Group_List[i][j]];
		}
		dq = q - User_CG_q0[i];
		Chi_SQ_User += (User_CG_k[i] * dq * dq);
	}

	return Chi_SQ_User;
}


#define FLAT_CG	(0.03)

double CFitting::Cal_Chi_SQ_ChargeRestraint(void)	// using ini charge as the target 
{
	int i, j, nAtom, *Is_Drude, IdxMax;
	double dq, abs_dq, *q, cg;

	q = pMol_ESP->CG;
	Is_Drude = pMol_ESP->IsDrude;
	nAtom = pMol_ESP->nAtom - 1;	// skip the test charge

	Chi_SQ_RSTR_CG=0.0;
	IdxMax = pMol_ESP->nAtom - 1;
	
	for(i=0; i<nAtom; i++)	{
		if(CG_Atom_Group_Size[i] > 0)	{
			if(CG_Atom_Group_Size[i] <= 2)	{	// H or regular heavy atoms, no lone pairs
				cg = 0.0;
				for(j=0; j<CG_Atom_Group_Size[i]; j++)	{
					cg += q[CG_Atom_Group_List[i][j]];
				}
				dq = cg - cg_0[i];


				abs_dq = fabs(dq);
				if(abs_dq < FLAT_CG)	{
					abs_dq = 0.0;
				}
				else	{
					abs_dq = abs_dq - FLAT_CG;
				}
				Chi_SQ_RSTR_CG += (abs_dq*abs_dq);
			}
			else	{	// heavy atom with lone pair
				cg = 0.0;
				for(j=2; j<CG_Atom_Group_Size[i]; j++)	{	// the charge on all lone pairs
					cg += q[CG_Atom_Group_List[i][j]];
				}
				dq = cg - cg_0[i];	// the sum of charges on all lone pairs is expected to be close to the reference charge
				abs_dq = fabs(dq);
				if(abs_dq < FLAT_CG)	{
					abs_dq = 0.0;
				}
				else	{
					abs_dq = abs_dq - FLAT_CG;
				}
				Chi_SQ_RSTR_CG += (abs_dq*abs_dq);

				cg = 0.0;
				for(j=0; j<2; j++)	{	// the charge on heavy atom and drude
					cg += q[CG_Atom_Group_List[i][j]];
				}
				Chi_SQ_RSTR_CG += (cg*cg);

			}

		}

	}

	return (Chi_SQ_RSTR_CG*w_charge/n_CG_Atom);
}



void CFitting::SetupTargetChargesFromIniCharges(void)
{
	int i, j, k, nAtom, nBond, IdxMax, *Is_Drude, Atom_1, Atom_2, iTmp, Local_List[8];
	int Bond_List[MAX_ATOM_TYPE][4], Bond_Count[MAX_ATOM_TYPE];
	double *q,*mass;

	q = pMol_ESP->CG;
	mass = pMol_ESP->mass;
	nAtom = pMol_ESP->nAtom-1;	// skip the test charge
	IdxMax = nAtom - 1;
	Is_Drude = pMol_ESP->IsDrude;

	memset(cg_0, 0, sizeof(double)*MAX_ATOM_TYPE);
	memset(Bond_Count, 0, sizeof(int)*MAX_ATOM_TYPE);
	memset(Bond_List, 0, sizeof(int)*MAX_ATOM_TYPE*4);

	n_CG_Atom = 0;
	for(i=0; i<nAtom; i++)	{
		if(mass[i] > 0.8)	{	// not a drude or LP
			CG_Atom_Group_Size[i] = 1;
			CG_Atom_Group_List[i][0] = i;
			n_CG_Atom++;
		}
		else	{
			CG_Atom_Group_Size[i] = 0;
		}
	}

	nBond = pMol_ESP->nBond;
	for(i=0; i<nBond; i++)	{
		Atom_1 = pMol_ESP->BondList[2*i];
		Atom_2 = pMol_ESP->BondList[2*i+1];

		if( (mass[Atom_1] > 0.8) && (mass[Atom_2] < 0.8) )	{
			CG_Atom_Group_List[Atom_1][CG_Atom_Group_Size[Atom_1]] = Atom_2;
			CG_Atom_Group_Size[Atom_1] ++;
		}
		else if((mass[Atom_2] > 0.8) && (mass[Atom_1] < 0.8))	{
			CG_Atom_Group_List[Atom_2][CG_Atom_Group_Size[Atom_2]] = Atom_1;
			CG_Atom_Group_Size[Atom_2] ++;
		}
	}

	//start	to sort atoms in the list according  to mass
	for(i=0; i<nAtom; i++)	{
		if(CG_Atom_Group_Size[i] > 2)	{	// there are host atom, drude and lone pairs
			for(j=1; j<CG_Atom_Group_Size[i]; j++)	{
				if( CG_Atom_Group_List[i][j] == (i+1) )	{	// the drude
					if( j != 1 )	{	// to exchange the position
						CG_Atom_Group_List[i][j] = CG_Atom_Group_List[i][1];
						CG_Atom_Group_List[i][1] = i+1;	// the drude
					}

					break;
				}
			}
		}
	}
	//end	to sort atoms in the list according  to mass
	
	for(i=0; i<nAtom; i++)	{
		cg_0[i] = 0.0;
		if(CG_Atom_Group_Size[i] > 0)	{
			for(j=0; j<CG_Atom_Group_Size[i]; j++)	{
				cg_0[i] += q[CG_Atom_Group_List[i][j]];
			}
		}
	}
}


#define FLAT_THOLE	(0.00)
#define TARGET_THOLE	(1.3)
//#define TARGET_THOLE	(0.20)

double CFitting::Cal_Chi_SQ_TholeRestraint(void)
{
	int i;
	double *thole, d_Thole, abs_d_Thole;

	thole = pMol_ESP->thole;

	Chi_SQ_RSTR_Thole = 0.0;

	for(i=0; i<n_Heavy_Atm; i++)	{
		d_Thole = pMol_ESP->thole[HeavyAtmList[i]] - TARGET_THOLE;
		abs_d_Thole = fabs(d_Thole);
		if(abs_d_Thole > FLAT_THOLE)	{
			abs_d_Thole -= FLAT_THOLE;
			Chi_SQ_RSTR_Thole += (abs_d_Thole*abs_d_Thole);
		}
	}

	return (Chi_SQ_RSTR_Thole*w_thole/n_Heavy_Atm);
}

/*
double Cal_Weight_Anisotropy(double K)
{
	return (exp(-0.1*(K+200.0)));	// > -200
}
*/

#define K_CENTER	(50.0)
#define FLAT_K	(200.0)
double Cal_Weight_Anisotropy(double K)
{
	double d_k, abs_d_k, weight;

	d_k = K - K_CENTER;
	abs_d_k = fabs(d_k);

	if(abs_d_k > FLAT_K)	{
		abs_d_k -= FLAT_K;
	}
	else	{
		abs_d_k = 0.0;
	}

	weight = 0.005 * abs_d_k * abs_d_k;

	return weight;	// [-150, 250]
}


double CFitting::Cal_Chi_SQ_AnisotropyRestraint(void)
{
	int i, Idx, nAnisotropy=0;
	double K11, K22, K33;

	Chi_SQ_RSTR_Anisotropy = 0.0;

	for(i=0; i<pMol_ESP->nActive_Aniso; i++)	{
		Idx = pMol_ESP->Active_List_Aniso[i];	//the index of anisotropy

		K11 = pMol_ESP->Para_Aniso[Idx][0];
		K22 = pMol_ESP->Para_Aniso[Idx][1];
		K33 = pMol_ESP->Para_Aniso[Idx][2];

		Chi_SQ_RSTR_Anisotropy += Cal_Weight_Anisotropy(K11);
		Chi_SQ_RSTR_Anisotropy += Cal_Weight_Anisotropy(K22);
		Chi_SQ_RSTR_Anisotropy += Cal_Weight_Anisotropy(K33);

		nAnisotropy++;
	}

	if(nAnisotropy > 0)	{
		return (Chi_SQ_RSTR_Anisotropy*w_anisotropy/nAnisotropy);
	}
	else	{
		return 0.0;
	}
}

#define FLAT_ALPHA	(0.00)

double CFitting::Cal_Chi_SQ_AlphaRestraint(void)
{
	int i;
	double d_Alpha, abs_d_Alpha;

	Chi_SQ_RSTR_Alpha = 0.0;

	for(i=0; i<n_Heavy_Atm; i++)	{
		d_Alpha = pMol_ESP->alpha[HeavyAtmList[i]] - Alpha_0[i]*AutoScale_Alpha;
		abs_d_Alpha = fabs(d_Alpha);
		if(abs_d_Alpha > FLAT_ALPHA)	{
			abs_d_Alpha -= FLAT_ALPHA;
			Chi_SQ_RSTR_Alpha += (abs_d_Alpha*abs_d_Alpha);
		}
	}

	return (Chi_SQ_RSTR_Alpha*w_alpha/n_Heavy_Atm);
}



double CFitting::CalObjectiveFunction(void)
{
	Chi_SQ = Chi_SQ_ESP = Chi_SQ_RSTR_CG = Chi_SQ_RSTR_Alpha = Chi_SQ_RSTR_Thole = Chi_SQ_RSTR_Anisotropy = Chi_SQ_E_int = Chi_SQ_Wat_Emin = Chi_SQ_Wat_Rmin = Chi_SQ_f_Pert = Chi_SQ_User = 0.0;

	DistributeData();

	if(Not_Fit_Elec_Para == 0)	{
		Chi_SQ += Cal_Chi_SQ_ESP();
		Chi_SQ += Cal_Chi_SQ_ChargeRestraint();
		Chi_SQ += Cal_Chi_SQ_AlphaRestraint();
		Chi_SQ += Cal_Chi_SQ_TholeRestraint();
		Chi_SQ += Cal_Chi_SQ_AnisotropyRestraint();
		Chi_SQ += Cal_Chi_SQ_User();
	}

//	Chi_SQ += Cal_Chi_SQ_Int_Energy_Mol_Ar();

	if(Target_E_Int_Water)	{
		Chi_SQ += Cal_Chi_SQ_Int_Energy_Mol_Water(0);
	}

//	if(Target_E_Int_Dimer)	{
//		Chi_SQ += Cal_Chi_SQ_E_Int_Dimer();
//	}

//	if(To_Fit_LJ)	{
//		Chi_SQ += Cal_Chi_SQ_Force_Pert_Ar();
//		Chi_SQ += Cal_Chi_SQ_Int_Energy_Mol_Mol();
//	}
	
//	MPI_Finalize();
	return Chi_SQ;
}

void CFitting::Generate_Water_Confs_Donor(void)
{
	int IdxDonor, nAtom, theta;
	double r_H_OH2, E_Dimer_Min, r_Donor_H_Cut=2.7, E_Min_Local, Para_Opt[N_PARA_WAT];
	CMol *Dimer;

	Dimer = pMol_Wat;
	nAtom = pMol_ESP->nAtom;
	for(IdxDonor=0; IdxDonor<N_Donor; IdxDonor++)	{
		Mol_Active_Atom = Donor_List_Active[IdxDonor];

		E_Min_Local = 1.0E100;

		if(theta_Donor_Save[IdxDonor] != 0)	{	// the orientation is tilted
			theta=theta_Donor_Save[IdxDonor];
//			for(theta=-30; theta<=30; theta+=30)	{	// to get the optimal poses among three orientations
				memcpy(Para_Save_Donor[IdxDonor], Para_Save_Donor_Org[IdxDonor], sizeof(double)*N_PARA_WAT);
				OptWater.Generate_Dummy_Atoms(Mol_Active_Atom, 0, theta);
				OptWater.Optimize_Water_Pose(Para_Save_Donor[IdxDonor], E_Dimer_Min, r_H_OH2);
				
				if(E_Dimer_Min < E_Min_Local)	{
					E_Min_Local = E_Dimer_Min;
					theta_Donor[IdxDonor] = theta;
					memcpy(Para_Opt, Para_Save_Donor[IdxDonor], sizeof(double)*N_PARA_WAT);
					r_Min_MM_Donor[IdxDonor] = r_H_OH2;
					
					memcpy(x_Wat_Donor[IdxDonor], Dimer->x+nAtom, sizeof(double)*3);
					memcpy(y_Wat_Donor[IdxDonor], Dimer->y+nAtom, sizeof(double)*3);
					memcpy(z_Wat_Donor[IdxDonor], Dimer->z+nAtom, sizeof(double)*3);
				}
//			}
		}
		else	{
			theta = 0;
			memcpy(Para_Save_Donor[IdxDonor], Para_Save_Donor_Org[IdxDonor], sizeof(double)*N_PARA_WAT);
			OptWater.Generate_Dummy_Atoms(Mol_Active_Atom, 0, theta);
			OptWater.Optimize_Water_Pose(Para_Save_Donor[IdxDonor], E_Dimer_Min, r_H_OH2);
			
			if(E_Dimer_Min < E_Min_Local)	{
				E_Min_Local = E_Dimer_Min;
				theta_Donor[IdxDonor] = theta;
				memcpy(Para_Opt, Para_Save_Donor[IdxDonor], sizeof(double)*N_PARA_WAT);
				r_Min_MM_Donor[IdxDonor] = r_H_OH2;
				
				memcpy(x_Wat_Donor[IdxDonor], Dimer->x+nAtom, sizeof(double)*3);
				memcpy(y_Wat_Donor[IdxDonor], Dimer->y+nAtom, sizeof(double)*3);
				memcpy(z_Wat_Donor[IdxDonor], Dimer->z+nAtom, sizeof(double)*3);
			}
/*
			if( (fabs(Upper_Bound_R_H_Bond - r_H_OH2) < 0.02) || (fabs(Lower_Bound_R_H_Bond - r_H_OH2) < 0.02) )	{	// touching the boundaries
				for(theta=-30; theta<=30; theta+=60)	{	// to get the optimal poses among three orientations
					memcpy(Para_Save_Donor[IdxDonor], Para_Save_Donor_Org[IdxDonor], sizeof(double)*N_PARA_WAT);
					OptWater.Generate_Dummy_Atoms(Mol_Active_Atom, 0, theta);
					OptWater.Optimize_Water_Pose(Para_Save_Donor[IdxDonor], E_Dimer_Min, r_H_OH2);
					
					if(E_Dimer_Min < E_Min_Local)	{
						E_Min_Local = E_Dimer_Min;
						theta_Donor[IdxDonor] = theta;
						memcpy(Para_Opt, Para_Save_Donor[IdxDonor], sizeof(double)*N_PARA_WAT);
						r_Min_MM_Donor[IdxDonor] = r_H_OH2;
						
						memcpy(x_Wat_Donor[IdxDonor], Dimer->x+nAtom, sizeof(double)*3);
						memcpy(y_Wat_Donor[IdxDonor], Dimer->y+nAtom, sizeof(double)*3);
						memcpy(z_Wat_Donor[IdxDonor], Dimer->z+nAtom, sizeof(double)*3);
					}
				}
			}
*/
		}

		memcpy(Para_Save_Donor[IdxDonor], Para_Opt, sizeof(double)*N_PARA_WAT);;
	}
}

void CFitting::Generate_Water_Confs_Acceptor(void)
{
	int IdxAcceptor, nAtom, theta;
	double E_Dimer_Min, r_Acceptor_H, E_Min_Local, Para_Opt[N_PARA_WAT];
	CMol *Dimer;

	Dimer = pMol_Wat;
	nAtom = pMol_ESP->nAtom;
	for(IdxAcceptor=0; IdxAcceptor<N_Acceptor; IdxAcceptor++)	{
		Mol_Active_Atom = Acceptor_List_Acceptor[IdxAcceptor];

		E_Min_Local = 1.0E100;

		if(theta_Acceptor_Save[IdxAcceptor] != 0)	{	// the orientation is tilted
			theta=theta_Acceptor_Save[IdxAcceptor];
//			for(theta=-45; theta<=45; theta+=45)	{	// to get the optimal poses among three orientations
				memcpy(Para_Save_Acceptor[IdxAcceptor], Para_Save_Acceptor_Org[IdxAcceptor], sizeof(double)*N_PARA_WAT);
				OptWater.Generate_Dummy_Atoms(Mol_Active_Atom, 1, theta);
				OptWater.Optimize_Water_Pose(Para_Save_Acceptor[IdxAcceptor], E_Dimer_Min, r_Acceptor_H);
				
				if(E_Dimer_Min < E_Min_Local)	{
					E_Min_Local = E_Dimer_Min;
					theta_Acceptor[IdxAcceptor] = theta;
					memcpy(Para_Opt, Para_Save_Acceptor[IdxAcceptor], sizeof(double)*N_PARA_WAT);
					r_Min_MM_Acceptor[IdxAcceptor] = r_Acceptor_H;
					
					memcpy(x_Wat_Acceptor[IdxAcceptor], Dimer->x+nAtom, sizeof(double)*3);
					memcpy(y_Wat_Acceptor[IdxAcceptor], Dimer->y+nAtom, sizeof(double)*3);
					memcpy(z_Wat_Acceptor[IdxAcceptor], Dimer->z+nAtom, sizeof(double)*3);
				}
//			}
		}
		else	{	// colinear in QM
			theta = 0;
			memcpy(Para_Save_Acceptor[IdxAcceptor], Para_Save_Acceptor_Org[IdxAcceptor], sizeof(double)*N_PARA_WAT);
			OptWater.Generate_Dummy_Atoms(Mol_Active_Atom, 1, theta);
			OptWater.Optimize_Water_Pose(Para_Save_Acceptor[IdxAcceptor], E_Dimer_Min, r_Acceptor_H);
			
			if(E_Dimer_Min < E_Min_Local)	{
				E_Min_Local = E_Dimer_Min;
				theta_Acceptor[IdxAcceptor] = theta;
				memcpy(Para_Opt, Para_Save_Acceptor[IdxAcceptor], sizeof(double)*N_PARA_WAT);
				r_Min_MM_Acceptor[IdxAcceptor] = r_Acceptor_H;
				
				memcpy(x_Wat_Acceptor[IdxAcceptor], Dimer->x+nAtom, sizeof(double)*3);
				memcpy(y_Wat_Acceptor[IdxAcceptor], Dimer->y+nAtom, sizeof(double)*3);
				memcpy(z_Wat_Acceptor[IdxAcceptor], Dimer->z+nAtom, sizeof(double)*3);
			}
/*
//			if( (E_Dimer_Min > -3.5) && ( (fabs(Upper_Bound_R_H_Bond - r_Acceptor_H) < 0.01) || (fabs(Lower_Bound_R_H_Bond - r_Acceptor_H) < 0.01) ) )	{	// touching the boundaries
			if( (fabs(Upper_Bound_R_H_Bond - r_Acceptor_H) < 0.02) || (fabs(Lower_Bound_R_H_Bond - r_Acceptor_H) < 0.02) )	{	// touching the boundaries
				for(theta=-45; theta<=45; theta+=90)	{	// to get the optimal poses among three orientations
					memcpy(Para_Save_Acceptor[IdxAcceptor], Para_Save_Acceptor_Org[IdxAcceptor], sizeof(double)*N_PARA_WAT);
					OptWater.Generate_Dummy_Atoms(Mol_Active_Atom, 1, theta);
					OptWater.Optimize_Water_Pose(Para_Save_Acceptor[IdxAcceptor], E_Dimer_Min, r_Acceptor_H);
					
					if(E_Dimer_Min < E_Min_Local)	{
						E_Min_Local = E_Dimer_Min;
						theta_Acceptor[IdxAcceptor] = theta;
						memcpy(Para_Opt, Para_Save_Acceptor[IdxAcceptor], sizeof(double)*N_PARA_WAT);
						r_Min_MM_Acceptor[IdxAcceptor] = r_Acceptor_H;
						
						memcpy(x_Wat_Acceptor[IdxAcceptor], Dimer->x+nAtom, sizeof(double)*3);
						memcpy(y_Wat_Acceptor[IdxAcceptor], Dimer->y+nAtom, sizeof(double)*3);
						memcpy(z_Wat_Acceptor[IdxAcceptor], Dimer->z+nAtom, sizeof(double)*3);
					}
				}
			}
*/
		}

		memcpy(Para_Save_Acceptor[IdxAcceptor], Para_Opt, sizeof(double)*N_PARA_WAT);;
	}
}

void CFitting::To_Detect_Invalid_Acceptor_Donor(void)
{
	int i;

	nFit_Donor = 0;
	for(i=0; i<N_Donor; i++)	{
		if( fabs(r_Min_MM_Donor[i] - r_Min_QM_Donor[i]) > 0.7 )	{	// too large distance
			ToFit_Donor[i] = 0;
			fprintf(fFile_Run_Log, "Warning> The interaction of donor atom %d and water is not going to be fitted.\n", Donor_List_Active[i]+1);
		}
		else if( r_Min_QM_Donor[i] > (Upper_Bound_R_H_Bond + 0.2) )	{
			ToFit_Donor[i] = 0;
			fprintf(fFile_Run_Log, "Warning> The interaction of donor atom %d and water is not going to be fitted.\n", Donor_List_Active[i]+1);
		}
		else	{
			ToFit_Donor[i] = 1;
			nFit_Donor++;
		}
	}
	nFit_Acceptor = 0;
	for(i=0; i<N_Acceptor; i++)	{
		if( fabs(r_Min_MM_Acceptor[i] - r_Min_QM_Acceptor[i]) > 0.7 )	{	// too large distance
			ToFit_Acceptor[i] = 0;
			fprintf(fFile_Run_Log, "Warning> The interaction of donor atom %d and water is not going to be fitted.\n", Acceptor_List_Acceptor[i]+1);
		}
		else if( r_Min_QM_Acceptor[i] > (Upper_Bound_R_H_Bond + 0.2) )	{
			ToFit_Donor[i] = 0;
			fprintf(fFile_Run_Log, "Warning> The interaction of donor atom %d and water is not going to be fitted.\n", Donor_List_Active[i]+1);
		}
		else	{
			ToFit_Acceptor[i] = 1;
			nFit_Acceptor++;
		}
	}
	if( (nFit_Donor + nFit_Acceptor) == 0)	{
		fprintf(fFile_Run_Log, "Warning> nFit_Donor + nFit_Acceptor = 0\nThe interactions with water is not fitted.\n");
		Target_E_Int_Water = 0;
	}
}




int Iteration;
FILE *fLog;
static int FailCount=0;

double Callback_Eval_Gradient(unsigned n, const double *x, double *grad, void *my_func_data)
{
	int i, n_Local;

	n_Local = (int)n;
	for(i=0; i<n_Local; i++)	{
		fitting.NLOpt_Para[i] = x[i];
		fitting.Para_List[fitting.Idx_NLOpt_Org[i]] = x[i];
	}

	if(grad)	{
		fitting.CalGradient();	//gradient will be assigned into objGrad automatically
		Iteration++;
		for(i=0; i<n_Local; i++)	{
			grad[i] = fitting.Grad_List[fitting.Idx_NLOpt_Org[i]];
		}
		fprintf(fLog, "Iteration %4d  Chi^2 = %12.6lf Chi^2(ESP) = %12.6lf      Chi^2(RSTR_CG) = %12.6lf     Chi^2(RSTR_Alpha) = %12.6lf      Chi^2(RSTR_Thole) = %12.6lf      Chi^2(RSTR_Aniso) = %12.6lf      Chi^2(Emin) = %12.6lf      Chi^2(Rmin) = %12.6lf      Chi^2(CG_User) = %12.6lf\n", 
			Iteration, fitting.Chi_SQ, fitting.Chi_SQ_ESP*10000.0/(w_ESP_Grid_Sum*n_Conf*2.0), fitting.Chi_SQ_RSTR_CG*w_charge/fitting.n_CG_Atom, 
			fitting.Chi_SQ_RSTR_Alpha*w_alpha/fitting.n_Heavy_Atm, fitting.Chi_SQ_RSTR_Thole*w_thole/fitting.n_Heavy_Atm, fitting.Chi_SQ_RSTR_Anisotropy*w_anisotropy/fitting.pMol_ESP->nActive_Aniso, 
			(fitting.Chi_SQ_Wat_Emin*w_water_E_min)/(nFit_Acceptor+nFit_Donor), (fitting.Chi_SQ_Wat_Rmin*w_water_R_min)/(nFit_Acceptor+nFit_Donor), fitting.Chi_SQ_User);
		fflush(fLog);

		if(fitting.Chi_SQ < Chi_SQ_Min)	{
			if(fitting.Chi_SQ + 0.0001 > Chi_SQ_Min)	{	// too small progess
				FailCount++;
			}
			else	{
				FailCount = 0;
			}

			Chi_SQ_Min = fitting.Chi_SQ;
			fitting.SaveParameters();
			fitting.SaveParameters_CHARMM();

			memcpy(Para_Set_Save[iActiveSave], fitting.Para_List, sizeof(double)*fitting.n_Para);	// backup the new parameter sets
			Is_Saved_Para[iActiveSave] = 1;
			Para_Set_Saved_Chi2[iActiveSave] = Chi_SQ_Min;
			iActiveSave = (iActiveSave+1)%MAX_PARA_SAVE;
		}
		else	{
			if(FailCount >= 10)      {       // cannot converge due to limited accuracy of numerical gradient
				fprintf(fFile_Run_Log, "Too much failed tries.\nPlease check fitting.dat.\nQuit.\n");
				fflush(fFile_Run_Log);

				nlopt_force_stop(opt);	// terminate the optimization
			}
			FailCount++;
		}
	}
	else	{	// cal object function only
		fitting.CalObjectiveFunction();
	}

    return fitting.Chi_SQ;
}


//#define DELTA_PARA	(2.0E-4)	// Important parameters. Should be not too large or too small
//#define DELTA_PARA	(1.0E-4)	// Important parameters. Should be not too large or too small
//#define DELTA_PARA	(5.0E-5)
//#define DELTA_PARA	(2.0E-5)

#define DELTA_PARA	(1.0E-4)

void CFitting::CalGradient(void)
{
	int i;
	double Para_Save[MAX_NUM_PARA];
	double f_Left, f_Right;
	
	for(i=0; i<n_Para; i++)	{
		Para_Save[i] = Para_List[i];	//save current parameters
	}

	PreOptimizeConfigurations();

	if(Target_E_Int_Water)	{
		Generate_Water_Confs_Acceptor();
		Generate_Water_Confs_Donor();
	}


	for(i=0; i<n_Para; i++)	{
		if(Is_Para_Fixed[i] == 0)	{
			Para_List[i] = Para_Save[i] - DELTA_PARA;
			f_Left = CalObjectiveFunction();
			
			Para_List[i] = Para_Save[i] + DELTA_PARA;
			f_Right = CalObjectiveFunction();
			
			Grad_List[i] = (f_Right-f_Left)/(2.0*DELTA_PARA);
			
			Para_List[i] = Para_Save[i];	//restore save (correct) parameter
		}
		else	{
			Grad_List[i] = 0.0;
		}
	}

	
	CalObjectiveFunction();
	
	return;
}
#undef	DELTA_PARA


void CFitting::GenIniRandomParameters(void)
{
	int i;
	double dx;
	FILE *fOut;

	iseed = Gen_iseed() + time(NULL);	//using various iseed

	fOut = fopen("rand-ini-para.dat", "w");
	for(i=0; i<n_Para; i++)	{
		if(Is_Para_Fixed[i]==0)	{
			dx = UpBound[i] - LowBound[i];

//			if(dx < 0.6)	{
//				Para_List[i] = LowBound[i] + rand(iseed)*dx;
//			}
//			else	{
				Para_List[i] = Para_List[i] + 0.8*(rand(iseed) - 0.5);
//			}
			Para_List[i] = min(max(Para_List[i], LowBound[i]), UpBound[i]);
		}


		fprintf(fOut, "%24.16lf\n", Para_List[i]);

		if(Is_Para_Fixed[i]==0)	{
			NLOpt_Para[Idx_Org_NLOpt[i]] = Para_List[i];
		}
	}
	fclose(fOut);
}

#define BOUND_FIX	(0.00001)
void CFitting::Optimization_NLOpt(void)
{
	int i;
//	char szLine[256];
	double Obj_Min;

	n_Para_Free = 0;
	for(i=0; i<n_Para; i++)	{
		if( LowBound[i]+BOUND_FIX >= UpBound[i] )	{	// fixed parameter
			Is_Para_Fixed[i] = 1;
//			printf("Warning> LowBound >= UpBound for parameter %d\n", i+1);
			Para_List[i] = UpBound[i];
			fprintf(fFile_Run_Log, "Warning> LowBound >= UpBound for parameter %d\n", i+1);
			fflush(fFile_Run_Log);
		}
		else	{	// free parameter
			Is_Para_Fixed[i] = 0;
			Idx_Org_NLOpt[i] = n_Para_Free;
			Idx_NLOpt_Org[n_Para_Free] = i;
			NLOpt_LowBound[n_Para_Free] = LowBound[i];
			NLOpt_UpBound[n_Para_Free] = UpBound[i];
			NLOpt_Para[n_Para_Free] = Para_List[i];
			n_Para_Free++;
		}
	}

//	printf("Boundaries list: \n");
//	for(i=0; i<n_Para; i++)	{
//		printf("%6.2lf %6.2lf  \n", LowBound[i], UpBound[i]);
//	}

//	if(ProgID == MY_MPI_ROOT)	{
//		GenIniRandomParameters();
//	}



	ReadParameters("para-opt-start.dat");	// the parameters will be overridden if there exist para-opt-start.dat



	for(i=0; i<n_Para; i++) {	//make sure the parameter is in the bounds
		if(Is_Para_Fixed[i] == 0)	{	// this parameter is not fixed
			NLOpt_Para[Idx_Org_NLOpt[i]] = Para_List[i];
			if(Para_List[i] < LowBound[i])  {
				Para_List[i] = LowBound[i] + 0.0001;
				NLOpt_Para[Idx_Org_NLOpt[i]] = Para_List[i];
				fprintf(fFile_Run_Log, "Warning. The initial value of parameter %3d is changed to the lower bound.\n", i+1);
			}
			else if(Para_List[i] > UpBound[i])   {
				Para_List[i] = UpBound[i] - 0.0001;
				NLOpt_Para[Idx_Org_NLOpt[i]] = Para_List[i];
				fprintf(fFile_Run_Log, "Warning. The initial value of parameter %3d is changed to the upper bound.\n", i+1);
			}       
		}
	}
	fflush(fFile_Run_Log);

	if(IsTestMode)	{
		fprintf(fFile_Run_Log, "There is no fatal error in test model.\nIt is time to run the code in productive mode.\n");
		fprintf(fFile_Run_Log,"Please wait for the fitting results sent by email. Thanks. \nHave a great day!\n");
		exit(0);
	}

	memset(Is_Saved_Para, 0, sizeof(int)*MAX_PARA_SAVE);
	iActiveSave=0;

	opt = nlopt_create(NLOPT_LD_LBFGS, n_Para_Free); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, NLOpt_LowBound);
	nlopt_set_upper_bounds(opt, NLOpt_UpBound);
	nlopt_set_min_objective(opt, Callback_Eval_Gradient, NULL);
	nlopt_set_xtol_rel(opt, 1.0E-10);
	
	Iteration = 0;
	fLog = fopen("fitting.dat", "w");

	PreOptimizeConfigurations();
	
	Generate_Water_Confs_Acceptor();
	Generate_Water_Confs_Donor();
	To_Detect_Invalid_Acceptor_Donor();
	
	CalObjectiveFunction();

	fprintf(fLog, "Iteration %4d  Chi^2 = %12.6lf Chi^2(ESP) = %12.6lf      Chi^2(RSTR_CG) = %12.6lf     Chi^2(RSTR_Alpha) = %12.6lf      Chi^2(RSTR_Thole) = %12.6lf      Chi^2(RSTR_Aniso) = %12.6lf      Chi^2(Emin) = %12.6lf      Chi^2(Rmin) = %12.6lf      Chi^2(CG_User) = %12.6lf\n", 
			Iteration, Chi_SQ, Chi_SQ_ESP*10000.0/(w_ESP_Grid_Sum*n_Conf*2.0), Chi_SQ_RSTR_CG*w_charge/n_CG_Atom, 
			Chi_SQ_RSTR_Alpha*w_alpha/n_Heavy_Atm, Chi_SQ_RSTR_Thole*w_thole/n_Heavy_Atm, Chi_SQ_RSTR_Anisotropy*w_anisotropy/pMol_ESP->nActive_Aniso, 
			(Chi_SQ_Wat_Emin*w_water_E_min)/(nFit_Acceptor+nFit_Donor), (Chi_SQ_Wat_Rmin*w_water_R_min)/(nFit_Acceptor+nFit_Donor), Chi_SQ_User);
	fflush(fLog);

	if (nlopt_optimize(opt, NLOpt_Para, &Obj_Min) < 0) {
//		printf("nlopt failed!\n");
		fprintf(fFile_Run_Log, "nlopt failed!\n");
		fflush(fFile_Run_Log);
	}
	else {
//		printf("Obj_Min = %16.10lf\n", Obj_Min);
		fprintf(fFile_Run_Log, "Obj_Min = %16.10lf\n", Obj_Min);
		fflush(fFile_Run_Log);
	}


	if(To_Find_Stable_Parameters() == 0)	{	// Check whether a down-scaled polarizability is needed
		memset(Is_Saved_Para, 0, sizeof(int)*MAX_PARA_SAVE);
		iActiveSave=0;

		AutoScale_Alpha *= SCALE_DOWN_ALPHA;
		fprintf(fFile_Run_Log, "Important!!!\nLong bond lengths between heavy atoms and their drudes are deteced.\nTarget alphas will be scaled down by %.3lf.\n", AutoScale_Alpha);
		fflush(fFile_Run_Log);

		Iteration = 0;
		FailCount = 0;
		Chi_SQ_Min = 1.0E10;

		CalObjectiveFunction();
		
		if (nlopt_optimize(opt, NLOpt_Para, &Obj_Min) < 0) {
			fprintf(fFile_Run_Log, "nlopt failed!\n");
			fflush(fFile_Run_Log);
		}
		else {
			fprintf(fFile_Run_Log, "Obj_Min = %16.10lf\n", Obj_Min);
			fflush(fFile_Run_Log);
		}
	}
	

	double R_Drude_Max;
	IsThereFarDrudes(R_Drude_Max);	// check again


	fclose(fLog);
	nlopt_destroy(opt);

}

void CFitting::SaveParameters(void)
{
	FILE *fOut;
	int i;
/*	
	if(ProgID != MY_MPI_ROOT)	{
		return;
	}
*/	
	fOut = fopen("para-opt.dat", "a+");
	fseek(fOut, 0, SEEK_END);
	for(i=0; i<n_Para; i++)	{
		fprintf(fOut, "%24.16lf\n", Para_List[i]);
	}
	fprintf(fOut, "\n\n! Chi^2 = %12.6lf Chi^2(ESP) = %12.6lf      Chi^2(RSTR_CG) = %12.6lf     Chi^2(RSTR_Alpha) = %12.6lf      Chi^2(RSTR_Thole) = %12.6lf      Chi^2(RSTR_Aniso) = %12.6lf      Chi^2(Emin) = %12.6lf      Chi^2(Rmin) = %12.6lf      Chi^2(CG_User) = %12.6lf\n", 
		Chi_SQ, Chi_SQ_ESP*10000.0/(w_ESP_Grid_Sum*n_Conf*2.0), 
		Chi_SQ_RSTR_CG*w_charge/n_CG_Atom, 
		Chi_SQ_RSTR_Alpha*w_alpha/n_Heavy_Atm, 
		Chi_SQ_RSTR_Thole*w_thole/n_Heavy_Atm, 
		Chi_SQ_RSTR_Anisotropy*w_anisotropy/pMol_ESP->nActive_Aniso, 
		(Chi_SQ_Wat_Emin*w_water_E_min)/(nFit_Acceptor+nFit_Donor), (Chi_SQ_Wat_Rmin*w_water_R_min)/(nFit_Acceptor+nFit_Donor), 
		Chi_SQ_User);
	fprintf(fOut, "!                  Std_Error(ESP) = %12.6lf  Std_Error(RSTR_CG) = %12.6lf Std_Error(RSTR_Alpha) = %12.6lf  Std_Error(RSTR_Thole) = %12.6lf                                        Std_Error(Emin) = %12.6lf  Std_Error(Rmin) = %12.6lf\n\n\n\n", 
		sqrt(Chi_SQ_ESP/(w_ESP_Grid_Sum*n_Conf*2.0)), sqrt(Chi_SQ_RSTR_CG/n_CG_Atom), sqrt(Chi_SQ_RSTR_Alpha/n_Heavy_Atm), sqrt(Chi_SQ_RSTR_Thole/n_Heavy_Atm), 
		sqrt(Chi_SQ_Wat_Emin/(nFit_Acceptor+nFit_Donor)), sqrt(Chi_SQ_Wat_Rmin/(nFit_Acceptor+nFit_Donor)));
	
	fclose(fOut);


	fOut = fopen("para-check.dat", "w");	// alway update. Will be used as para-opt-start.dat. 
	fseek(fOut, 0, SEEK_END);
	for(i=0; i<n_Para; i++)	{
		fprintf(fOut, "%24.16lf\n", Para_List[i]);
	}
	fprintf(fOut, "\n\n! Chi^2 = %12.6lf Chi^2(ESP) = %12.6lf      Chi^2(RSTR_CG) = %12.6lf     Chi^2(RSTR_Alpha) = %12.6lf      Chi^2(RSTR_Thole) = %12.6lf      Chi^2(RSTR_Aniso) = %12.6lf      Chi^2(Emin) = %12.6lf      Chi^2(Rmin) = %12.6lf      Chi^2(CG_User) = %12.6lf\n", 
		Chi_SQ, Chi_SQ_ESP*10000.0/(w_ESP_Grid_Sum*n_Conf*2.0), 
		Chi_SQ_RSTR_CG*w_charge/n_CG_Atom, 
		Chi_SQ_RSTR_Alpha*w_alpha/n_Heavy_Atm, 
		Chi_SQ_RSTR_Thole*w_thole/n_Heavy_Atm, 
		Chi_SQ_RSTR_Anisotropy*w_anisotropy/pMol_ESP->nActive_Aniso, 
		(Chi_SQ_Wat_Emin*w_water_E_min)/(nFit_Acceptor+nFit_Donor), (Chi_SQ_Wat_Rmin*w_water_R_min)/(nFit_Acceptor+nFit_Donor), 
		Chi_SQ_User);
	fprintf(fOut, "!                  Std_Error(ESP) = %12.6lf  Std_Error(RSTR_CG) = %12.6lf Std_Error(RSTR_Alpha) = %12.6lf  Std_Error(RSTR_Thole) = %12.6lf                                      Std_Error(Emin) = %12.6lf  Std_Error(Rmin) = %12.6lf\n\n\n\n", 
		sqrt(Chi_SQ_ESP/(w_ESP_Grid_Sum*n_Conf*2.0)), sqrt(Chi_SQ_RSTR_CG/n_CG_Atom), sqrt(Chi_SQ_RSTR_Alpha/n_Heavy_Atm), sqrt(Chi_SQ_RSTR_Thole/n_Heavy_Atm), 
		sqrt(Chi_SQ_Wat_Emin/(nFit_Acceptor+nFit_Donor)), sqrt(Chi_SQ_Wat_Rmin/(nFit_Acceptor+nFit_Donor)));
	
	fclose(fOut);
}

void CFitting::SaveParameters_CHARMM(void)
{
	int i, nAtom, IdxMax;
	FILE *fOut, *fFinal;
	CMol* pMol;
	char szBuff[256];
	double CG_Sum, q[MAX_ATOM_TYPE], alpha[MAX_ATOM_TYPE], thole[MAX_ATOM_TYPE], A11, A22;
	double q_non_polar[MAX_ATOM_TYPE];
	
/*	if(ProgID != MY_MPI_ROOT)	{
		return;
	}
*/	
	if(n_Mol < 1)	{
//		printf("n_Mol < 1\nQuit.\n");
		fprintf(fFile_Run_Log, "n_Mol < 1\nQuit.\n");
		fflush(fFile_Run_Log);
		exit(1);
	}
	pMol = pMol_List[0];
	nAtom = pMol->nAtom;
	
	fOut = fopen("para-opt-charmm.txt", "a+");
	if(fOut == NULL)	{
		fprintf(fFile_Run_Log, "Fail to open file: para-opt-charmm.txt\nQuit\n");
		fflush(fFile_Run_Log);
		exit(1);
	}
	fseek(fOut, 0, SEEK_END);

	fFinal = fopen(szFinalPara, "w");
	if(fFinal == NULL)	{
		fprintf(fFile_Run_Log, "Fail to open file: %s\nQuit\n", szFinalPara);
		fflush(fFile_Run_Log);
		exit(1);
	}

	fprintf(fFinal, "Begin your input parameters\n\n");
	int Idx_in_Conf = 0;
	while(1)	{
		Idx_in_Conf = FindItemInConfFile("equivalent", Idx_in_Conf);
		if(Idx_in_Conf < 0)	{	// have checked all records in configuration file
			break;
		}
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
		Idx_in_Conf++;
	}
	fprintf(fFinal, "\n");
	Idx_in_Conf = 0;
	while(1)	{
		Idx_in_Conf = FindItemInConfFile("constrain_CG", Idx_in_Conf);
		if(Idx_in_Conf < 0)	{	// have checked all records in configuration file
			break;
		}
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
		Idx_in_Conf++;
	}
	Idx_in_Conf = 0;
	while(1)	{
		Idx_in_Conf = FindItemInConfFile("fix_CG", Idx_in_Conf);
		if(Idx_in_Conf < 0)	{	// have checked all records in configuration file
			break;
		}
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
		Idx_in_Conf++;
	}
	Idx_in_Conf = 0;
	while(1)	{
		Idx_in_Conf = FindItemInConfFile("fix", Idx_in_Conf);
		if(Idx_in_Conf < 0)	{	// have checked all records in configuration file
			break;
		}
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
		Idx_in_Conf++;
	}
	fprintf(fFinal, "\n");
	Idx_in_Conf = 0;
	while(1)	{
		Idx_in_Conf = FindItemInConfFile("neutral", Idx_in_Conf);
		if(Idx_in_Conf < 0)	{	// have checked all records in configuration file
			break;
		}
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
		Idx_in_Conf++;
	}
	fprintf(fFinal, "\n");

	fprintf(fFinal, "\nnetcharge  %2.0lf\n", netcharge_mol);


	Idx_in_Conf = FindItemInConfFile("w_H_Donor_Acceptor", 0);
	if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
	}
	Idx_in_Conf = FindItemInConfFile("w_charge", 0);
	if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
	}
	Idx_in_Conf = FindItemInConfFile("w_alpha", 0);
	if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
	}
	Idx_in_Conf = FindItemInConfFile("w_thole", 0);
	if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
	}
	Idx_in_Conf = FindItemInConfFile("w_anisotropy", 0);
	if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
	}
	Idx_in_Conf = FindItemInConfFile("Scale_Alpha", 0);
	if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
	}
	Idx_in_Conf = FindItemInConfFile("Scale_TargetAlpha", 0);
	if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
	}
	Idx_in_Conf = FindItemInConfFile("To_Use_Miller_Alpha", 0);
	if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
		fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
	}

	Idx_in_Conf = FindItemInConfFile("Target_E_Int_Water", 0);
	if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
		fprintf(fFinal, "\nUsing molecule-water interactions in fitting\n");
		Idx_in_Conf = FindItemInConfFile("w_water_E_min", 0);
		if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
			fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
		}
		Idx_in_Conf = FindItemInConfFile("w_water_R_min", 0);
		if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
			fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
		}
		Idx_in_Conf = FindItemInConfFile("SCALE_QM_E_MIN_Drude", 0);
		if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
			fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
		}
		Idx_in_Conf = FindItemInConfFile("SHIFT_QM_R_MIN_Drude", 0);
		if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
			fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
		}
		Idx_in_Conf = FindItemInConfFile("SHIFT_QM_R_MIN_CHARGED_Drude", 0);
		if(Idx_in_Conf >= 0)	{	// have checked all records in configuration file
			fprintf(fFinal, "%s", szReadConf_File[Idx_in_Conf]);
		}
	}


	fprintf(fFinal, "\nEnd   your input parameters\n\n");
	fprintf(fFinal, "\nThe results of charge fitting is attached below.\n\n");




	CG_Sum = 0.0;
	fprintf(fOut, "GROUP\n");
	fprintf(fFinal, "GROUP\n");
	IdxMax = nAtom - 1;	// skip the test charge, the last atom
	for(i=0; i<IdxMax; i++)	{
		if(pMol->IsDrude[i])	{	//skip drude
		}
		else	{
			if( ( (i+1) < nAtom) && (pMol->IsDrude[i+1]) )	{
				q[i] = pMol->CG[i]+pMol->CG[i+1];
				sprintf(szBuff, "%8.3lf", q[i]);
				q[i] = atof(szBuff);

				sprintf(szBuff, "%8.3lf", pMol->alpha[i]);
				alpha[i] = atof(szBuff);

				sprintf(szBuff, "%8.3lf", pMol->thole[i]);
				thole[i] = atof(szBuff);

				CG_Sum += q[i];
			}
			else	{
				q[i] = pMol->CG[i];
				sprintf(szBuff, "%8.3lf", q[i]);
				q[i] = atof(szBuff);
				CG_Sum += q[i];
			}
		}
	}
	q[0] -= (CG_Sum-netcharge_mol);	// to make the whole molecule neutral


	FILE *fNewCG;
	fNewCG = fopen("cg-list.txt", "w");	// used for non-polarizable FF
//	fNewDrudeCG = fopen("drude-cg-list.txt", "w");	// used for polarizable FF
	int iRtf_Output_Count, HostAtom, LP_Count;

	iRtf_Output_Count = 0;
	LP_Count = 0;

	for(i=0; i<IdxMax; i++)	{
		q_non_polar[i] = q[i];
	}
	for(i=0; i<IdxMax; i++)	{
		if(pMol->mass[i] < 0.01)	{	// Lone pair. The charge should be merged into the charge of the host atom
			HostAtom = pMol->LP_Host[LP_Count];
			q_non_polar[HostAtom] += q_non_polar[i];
			q_non_polar[i] = 0.0;	// set the charge of LP as zero
			LP_Count++;
		}
	}


	for(i=0; i<IdxMax; i++)	{
		if(pMol->IsDrude[i])	{	//skip drude
		}
		else	{
			fprintf(fOut, "ATOM %-7s %-8s ", pMol->AtomName[i], pMol->ChemName[i]);
			fprintf(fFinal, "ATOM %-7s %-8s ", pMol->AtomName[i], pMol->ChemName[i]);
			if( ( (i+1) < nAtom) && (pMol->IsDrude[i+1]) )	{
				fprintf(fOut, "%6.3lf  ALPHA %6.3lf  THOLE %5.3lf\n", q[i], alpha[i], thole[i]);
				fprintf(fFinal, "%6.3lf  ALPHA %6.3lf  THOLE %5.3lf\n", q[i], alpha[i], thole[i]);
			}
			else	{
				fprintf(fOut, "%6.3lf\n", q[i]);
				fprintf(fFinal, "%6.3lf\n", q[i]);
			}
	
			//start	to update the rtf file
			if(pMol->mass[i] > 1.5)	{	// a heavy atom
				sprintf(szTxt_Rtf[Idx_Mol_Atom_Start+iRtf_Output_Count], "ATOM %-7s%-8s%7.3lf  ALPHA%7.3lf  THOLE%6.3lf\n", 
					pMol->AtomName[i], pMol->ChemName[i], q[i], pMol->alpha[i], pMol->thole[i]);
				iRtf_Output_Count++;
			}
			else	{	// lonepair or H atom
				sprintf(szTxt_Rtf[Idx_Mol_Atom_Start+iRtf_Output_Count], "ATOM %-7s%-8s%7.3lf\n", pMol->AtomName[i], pMol->ChemName[i], q[i]);
				iRtf_Output_Count++;
			}
			//end	to update the rtf file
	
			if(pMol->mass[i] > 0.9)	{	// all real atoms only
				fprintf(fNewCG, "%18.12lf\n", q_non_polar[i]);	// the charge for each real atom
			}
		}

//		fprintf(fNewDrudeCG, "%18.12lf\n", pMol->CG[i]);	// the charge for each particle
	}

	fclose(fNewCG);
//	fclose(fNewDrudeCG);

	fprintf(fFinal, "\n\nWith alpha scaled by %6.3lf, \n", Scale_Alpha);
	for(i=0; i<IdxMax; i++)	{
		if(pMol->IsDrude[i])	{	//skip drude
		}
		else	{
			fprintf(fFinal, "ATOM %-7s %-7s ", pMol->AtomName[i], pMol->ChemName[i]);
			if( ( (i+1) < nAtom) && (pMol->IsDrude[i+1]) )	{
				fprintf(fFinal, "%6.3lf  ALPHA %6.3lf  THOLE %5.3lf\n", q[i], alpha[i]*Scale_Alpha, thole[i]);
			}
			else	{
				fprintf(fFinal, "%6.3lf\n", q[i]);
			}
		}
	}


	if(n_Aniso_Fit > 0)	{
		fprintf(fOut, "\n! Anisotrpy parameters\n");
		fprintf(fOut, "! In toppar file, \n");
		for(i=0; i<pMol->nAniso; i++)	{
			Para_Aniso_K_A(pMol->Para_Aniso[i][0], pMol->Para_Aniso[i][1], pMol->Para_Aniso[i][2], A11, A22);
			fprintf(fOut, "ANISOTROPY %7s %7s %7s %7s  A11 %8.5lf A22 %8.5lf \n", 
				pMol->AtomName[pMol->AnisoList[i][0]], pMol->AtomName[pMol->AnisoList[i][1]], pMol->AtomName[pMol->AnisoList[i][2]], pMol->AtomName[pMol->AnisoList[i][3]], 
				A11, A22);

			// to update anisotropy parameters in rtf file
			sprintf(szTxt_Rtf[Aniso_List[i]], "ANISOTROPY %7s %7s %7s %7s  A11 %8.5lf A22 %8.5lf \n", 
				pMol->AtomName[pMol->AnisoList[i][0]], pMol->AtomName[pMol->AnisoList[i][1]], pMol->AtomName[pMol->AnisoList[i][2]], pMol->AtomName[pMol->AnisoList[i][3]], A11, A22);
		}
	}
	if(n_LP_Fit > 0)	{
		fprintf(fOut, "\n! Lone pair parameters\n");
		for(i=0; i<pMol->nLP; i++)	{
			if(pMol->LP_Dist[i] < 0)	{	// bisector
				fprintf(fOut, "LONEPAIR bisector %-7s%-7s%-7s%-7s distance %5.3lf angle %6.2lf dihe %6.2lf\n", 
					pMol->AtomName[pMol->LPList[i]], pMol->AtomName[pMol->LP_PosHost[i][0]], 
					pMol->AtomName[pMol->LP_PosHost[i][1]], pMol->AtomName[pMol->LP_PosHost[i][2]], 
					-pMol->LP_Dist[i], pMol->LP_Theta[i]*radian,  pMol->LP_Phi[i]);
			}
			else	{	// relative
				fprintf(fOut, "LONEPAIR relative %-7s%-7s%-7s%-7s distance %5.3lf angle %6.2lf dihe %6.2lf\n", 
					pMol->AtomName[pMol->LPList[i]], pMol->AtomName[pMol->LP_PosHost[i][0]], 
					pMol->AtomName[pMol->LP_PosHost[i][1]], pMol->AtomName[pMol->LP_PosHost[i][2]], 
					pMol->LP_Dist[i], pMol->LP_Theta[i]*radian,  pMol->LP_Phi[i]);
			}
		}
	}
	if(n_LJ_Fit > 0)	{
		fprintf(fOut, "\nLJ     parameters, \n");
		for(i=0; i<n_LJ_Fit; i++)	{
			fprintf(fOut, "%7s  0.0  %7.4lf %7.4lf\n", myForceField->LJ_Rec[Para_LJ[i]].Chem, myForceField->LJ_Rec[Para_LJ[i]].para[1], myForceField->LJ_Rec[Para_LJ[i]].para[2]);	// E_Min, r_Min
		}
	}

	if(n_NBFix_Fit > 0)	{
		fprintf(fOut, "\n\nNBFix parameters, \n");
		for(i=0; i<n_NBFix_Fit; i++)	{
			fprintf(fOut, "%-7s %-7s  %9.5lf %9.5lf\n", myForceField->NBFix_Rec[Para_NBFix[i]].Chem_1, myForceField->NBFix_Rec[Para_NBFix[i]].Chem_2, myForceField->NBFix_Rec[Para_NBFix[i]].para[0], myForceField->NBFix_Rec[Para_NBFix[i]].para[1]);
		}
	}

	fprintf(fOut, "\n\n!  Chi^2 = %7.4lf      Chi^2(ESP) = %7.4lf         Chi^2(RSTR_CG) = %7.4lf        Chi^2(RSTR_Alpha) = %7.4lf         Chi^2(RSTR_Thole) = %7.4lf         Chi^2(RSTR_Anisotropy) = %7.4lf         Chi^2(Emin) = %7.4lf         Chi^2(Rmin) = %7.4lf       Chi^2(CG_User) = %7.4lf\n", 
		Chi_SQ, Chi_SQ_ESP*10000.0/(w_ESP_Grid_Sum*n_Conf*2.0), 
		Chi_SQ_RSTR_CG*w_charge/n_CG_Atom, 
		Chi_SQ_RSTR_Alpha*w_alpha/n_Heavy_Atm, 
		Chi_SQ_RSTR_Thole*w_thole/n_Heavy_Atm, 
		Chi_SQ_RSTR_Anisotropy*w_anisotropy/pMol_ESP->nActive_Aniso, 
		Chi_SQ_Wat_Emin*w_water_E_min/(nFit_Acceptor+nFit_Donor), Chi_SQ_Wat_Rmin*w_water_R_min/(nFit_Acceptor+nFit_Donor), Chi_SQ_User);
	fprintf(fOut, "!                   Std_Error(ESP) = %7.4lf     Std_Error(RSTR_CG) = %7.4lf    Std_Error(RSTR_Alpha) = %7.4lf     Std_Error(RSTR_Thole) = %7.4lf                                         Std_Error(Emin) = %7.4lf     Std_Error(Rmin) = %7.4lf\n\n\n\n", 
		sqrt(Chi_SQ_ESP/(w_ESP_Grid_Sum*n_Conf*2.0)), sqrt(Chi_SQ_RSTR_CG/n_CG_Atom), sqrt(Chi_SQ_RSTR_Alpha/n_Heavy_Atm), sqrt(Chi_SQ_RSTR_Thole/n_Heavy_Atm), sqrt(Chi_SQ_Wat_Emin/(nFit_Acceptor+nFit_Donor)), sqrt(Chi_SQ_Wat_Rmin/(nFit_Acceptor+nFit_Donor)));


	if(n_Aniso_Fit > 0)	{
		fprintf(fFinal, "\n! Anisotrpy parameters\n");
		fprintf(fFinal, "! In toppar file, \n");
		for(i=0; i<pMol->nAniso; i++)	{
			Para_Aniso_K_A(pMol->Para_Aniso[i][0], pMol->Para_Aniso[i][1], pMol->Para_Aniso[i][2], A11, A22);
			fprintf(fFinal, "ANISOTROPY %7s %7s %7s %7s  A11 %8.5lf A22 %8.5lf \n", 
				pMol->AtomName[pMol->AnisoList[i][0]], pMol->AtomName[pMol->AnisoList[i][1]], pMol->AtomName[pMol->AnisoList[i][2]], pMol->AtomName[pMol->AnisoList[i][3]], 
				A11, A22);
		}
	}
	if(n_LP_Fit > 0)	{
		fprintf(fFinal, "\n! Lone pair parameters\n");
		for(i=0; i<pMol->nLP; i++)	{
			if(pMol->LP_Dist[i] < 0)	{	// bisector
				fprintf(fFinal, "LONEPAIR bisector %-7s%-7s%-7s%-7s distance %5.3lf angle %6.2lf dihe %6.2lf\n", 
					pMol->AtomName[pMol->LPList[i]], pMol->AtomName[pMol->LP_PosHost[i][0]], 
					pMol->AtomName[pMol->LP_PosHost[i][1]], pMol->AtomName[pMol->LP_PosHost[i][2]], 
					pMol->LP_Dist[i], pMol->LP_Theta[i]*radian,  pMol->LP_Phi[i]);
			}
			else	{	// relative
				fprintf(fFinal, "LONEPAIR relative %-7s%-7s%-7s%-7s distance %5.3lf angle %6.2lf dihe %6.2lf\n", 
					pMol->AtomName[pMol->LPList[i]], pMol->AtomName[pMol->LP_PosHost[i][0]], 
					pMol->AtomName[pMol->LP_PosHost[i][1]], pMol->AtomName[pMol->LP_PosHost[i][2]], 
					pMol->LP_Dist[i], pMol->LP_Theta[i]*radian,  pMol->LP_Phi[i]);
			}
		}
	}

	if(n_LJ_Fit > 0)	{
		fprintf(fFinal, "\nLJ     parameters, \n");
		for(i=0; i<n_LJ_Fit; i++)	{
			fprintf(fFinal, "%7s  0.0  %7.4lf %7.4lf\n", myForceField->LJ_Rec[Para_LJ[i]].Chem, myForceField->LJ_Rec[Para_LJ[i]].para[1], myForceField->LJ_Rec[Para_LJ[i]].para[2]);
		}
	}
	if(n_NBFix_Fit > 0)	{
		fprintf(fFinal, "\n\nNBFix parameters, \n");
		for(i=0; i<n_NBFix_Fit; i++)	{
			fprintf(fFinal, "%-7s %-7s  %9.5lf %9.5lf\n", myForceField->NBFix_Rec[Para_NBFix[i]].Chem_1, myForceField->NBFix_Rec[Para_NBFix[i]].Chem_2, myForceField->NBFix_Rec[Para_NBFix[i]].para[0], myForceField->NBFix_Rec[Para_NBFix[i]].para[1]);
		}
	}

	fprintf(fFinal, "\n\n!  Chi^2 = %7.4lf      Chi^2(ESP) = %7.4lf         Chi^2(RSTR_CG) = %7.4lf        Chi^2(RSTR_Alpha) = %7.4lf         Chi^2(RSTR_Thole) = %7.4lf         Chi^2(RSTR_Aniso) = %7.4lf         Chi^2(Emin) = %7.4lf         Chi^2(Rmin) = %7.4lf       Chi^2(CG_User) = %7.4lf\n", 
		Chi_SQ, Chi_SQ_ESP*10000.0/(w_ESP_Grid_Sum*n_Conf*2.0), 
		Chi_SQ_RSTR_CG*w_charge/n_CG_Atom, 
		Chi_SQ_RSTR_Alpha*w_alpha/n_Heavy_Atm, 
		Chi_SQ_RSTR_Thole*w_thole/n_Heavy_Atm, 
		Chi_SQ_RSTR_Anisotropy*w_anisotropy/pMol_ESP->nActive_Aniso, 
		Chi_SQ_Wat_Emin*w_water_E_min/(nFit_Acceptor+nFit_Donor), Chi_SQ_Wat_Rmin*w_water_R_min/(nFit_Acceptor+nFit_Donor), Chi_SQ_User);
	fprintf(fFinal, "!                   Std_Error(ESP) = %7.4lf     Std_Error(RSTR_CG) = %7.4lf    Std_Error(RSTR_Alpha) = %7.4lf     Std_Error(RSTR_Thole) = %7.4lf                                         Std_Error(Emin) = %7.4lf     Std_Error(Rmin) = %7.4lf\n\n\n\n", 
		sqrt(Chi_SQ_ESP/(w_ESP_Grid_Sum*n_Conf*2.0)), sqrt(Chi_SQ_RSTR_CG/n_CG_Atom), sqrt(Chi_SQ_RSTR_Alpha/n_Heavy_Atm), sqrt(Chi_SQ_RSTR_Thole/n_Heavy_Atm), sqrt(Chi_SQ_Wat_Emin/(nFit_Acceptor+nFit_Donor)), sqrt(Chi_SQ_Wat_Rmin/(nFit_Acceptor+nFit_Donor)));
	
	fclose(fOut);

	fprintf(fFinal, "\n\nAutoScale_Alpha = %.4lf\nFAR_DRUDE_R = %.4lf\n", AutoScale_Alpha, FAR_DRUDE_R);

	fclose(fFinal);

	Update_Rtf_File();



	if(QM_Dipole > 1.0E-10)	{
		fFinal = fopen(szFinalPara, "a+");
		fseek(fFinal, 0, SEEK_END);
		fprintf(fFinal, "\n\nCompare QM and MM dipole: \n");
		fprintf(fFinal, "QM dipole: (%7.4lf, %7.4lf, %7.4lf) %7.4lf\n", QM_Dipole_x, QM_Dipole_y, QM_Dipole_z, QM_Dipole);
		fprintf(fFinal, "MM dipole: (%7.4lf, %7.4lf, %7.4lf) %7.4lf\n\n", MM_Dipole_x, MM_Dipole_y, MM_Dipole_z, MM_Dipole);

		fclose(fFinal);
	}



	if(Target_E_Int_Water)	{
		PreOptimizeConfigurations();
		Cal_Chi_SQ_Int_Energy_Mol_Water(1);

		fFinal = fopen(szFinalPara, "a+");
		fseek(fFinal, 0, SEEK_END);
		fprintf(fFinal, "\n\nList of E_Min, Rmin in MM and QM: \n");
		
		for(int iDonor=0; iDonor<N_Donor; iDonor++)	{
			if(ToFit_Donor[iDonor] == 0)	{
				fprintf(fFinal, "%2d Donor, the H atom %3d is NOT used for fitting due to the big difference between QM and MM.\n", 
					iDonor+1, Donor_List_Active[iDonor]+1);
			}
			else	{
				fprintf(fFinal, "%2d Donor, the H atom %3d, MM E_min = %7.3lf, QM E_min = %7.3lf, MM R_min = %5.3lf, QM R_min = %5.3lf\n", 
					iDonor+1, Donor_List_Active[iDonor]+1, E_Min_MM_Donor[iDonor], E_Min_QM_Donor[iDonor], r_Min_MM_Donor[iDonor], r_Min_QM_Donor[iDonor]);
			}
		}

		fprintf(fFinal, "\n");
		for(int iAcceptor=0; iAcceptor<N_Acceptor; iAcceptor++)	{
			if(ToFit_Acceptor[iAcceptor] == 0)	{
				fprintf(fFinal, "%2d Acceptor, atom %3d  is NOT used for fitting due to the big difference between QM and MM.\n", 
					iAcceptor+1, Acceptor_List_Acceptor[iAcceptor]+1);
			}
			else	{
				fprintf(fFinal, "%2d Acceptor, atom %3d, MM E_min = %7.3lf, QM E_min = %7.3lf, MM R_min = %5.3lf, QM R_min = %5.3lf\n", 
					iAcceptor+1, Acceptor_List_Acceptor[iAcceptor]+1, E_Min_MM_Acceptor[iAcceptor], E_Min_QM_Acceptor[iAcceptor], r_Min_MM_Acceptor[iAcceptor], r_Min_QM_Acceptor[iAcceptor]);
			}
		}

		double rmsd;

		rmsd = CalRMSD_Drude(x_Mol_Full_Opt, y_Mol_Full_Opt, z_Mol_Full_Opt, 
			ESP_Coord_x[0], ESP_Coord_y[0], ESP_Coord_z[0], pMol->nAtom-1, pMol);
		fprintf(fFinal,"Note: RMSD between optimized structures in MM and QM is %.3lf Angstrom.\n", rmsd);


//		fprintf(fFinal,"Note: QM E_min is scaled by a constant %4.2lf for non-charged molecules. QM R_min is shifted by %5.2lf Angstrom for non-charged molecules and %5.2lf Angstrom for charged molecules. \n", 
//			SCALE_QM_E_MIN, SHIFT_QM_R_MIN, SHIFT_QM_R_MIN_CHARGED);
		fclose(fFinal);
	}
}

#define R_CUTOF_AR_PROBE	(5.0)
#define HIGH_FORCE_CUTOFF	(0.30)
#define HIGH_ENERGY_CUTOFF	(0.50)
#define HIGH_FORCE_CUTOFF_SQ	(0.09)
void CFitting::ReadArPertForceData(void)
{
	FILE *fIn;
	char szLine[256];
	int i, iConf, nAtom, n_Real_Atom, Idx_Real_Atom_Max, ReadItem, Idx_Ar, RealAtom_List[MAX_NUM_PARA], AtomIdx;
	double dx, dy, dz, r, dist_Min=1.0E10, f_SQ, xSum, ySum, zSum;

	memset(Used_Ar_Pert, 0, sizeof(int)*MAX_N_PERT_MOL_AR);	// defaut: not used for fitting
	memset(f_Pert_x_QM, 0, sizeof(double)*MAX_N_PERT_MOL_AR*MAX_NUM_PARA);
	memset(f_Pert_y_QM, 0, sizeof(double)*MAX_N_PERT_MOL_AR*MAX_NUM_PARA);
	memset(f_Pert_z_QM, 0, sizeof(double)*MAX_N_PERT_MOL_AR*MAX_NUM_PARA);

	nAtom = pMol_Ar->nAtom;

	fIn = fopen(szFile_MolArForce, "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Error in opening the data file for Ar perturbed forces.\nQuit\n");
	}
	BuildRealAtomList(pMol_Ar, n_Real_Atom, RealAtom_List);	// include the Ar probe atom

	nConf_Mol_Ar = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		fgets(szLine, 256, fIn);
		if(strncmp(szLine, "Coordinate", 10) == 0)	{
			for(i=0; i<n_Real_Atom; i++)	{
				fgets(szLine, 256, fIn);
				AtomIdx = RealAtom_List[i];
				ReadItem = sscanf(szLine, "%lf %lf %lf", &(Ar_Pert_x[nConf_Mol_Ar][AtomIdx]), &(Ar_Pert_y[nConf_Mol_Ar][AtomIdx]), &(Ar_Pert_z[nConf_Mol_Ar][AtomIdx]));
				if(ReadItem != 3)	{
					fclose(fIn);
					Quit_With_Error_Msg("Error in reading data file for Ar perturbed forces.\nThe data for the coordinate are not complete. \n");
				}
				if(pMol_Ar->IsDrude[AtomIdx+1])	{	// a drude host
					Ar_Pert_x[nConf_Mol_Ar][AtomIdx+1] = Ar_Pert_x[nConf_Mol_Ar][AtomIdx];
					Ar_Pert_y[nConf_Mol_Ar][AtomIdx+1] = Ar_Pert_y[nConf_Mol_Ar][AtomIdx];
					Ar_Pert_z[nConf_Mol_Ar][AtomIdx+1] = Ar_Pert_z[nConf_Mol_Ar][AtomIdx];
				}
			}
			fgets(szLine, 256, fIn);
			if(strncmp(szLine, "Force", 5) != 0)	{	// wrong format? 
				fclose(fIn);
				Quit_With_Error_Msg("Error in reading data file for Ar perturbed forces.\nIt is supposed to be the entry for forces. \n");
			}
			for(i=0; i<n_Real_Atom; i++)	{
				fgets(szLine, 256, fIn);
				AtomIdx = RealAtom_List[i];
				ReadItem = sscanf(szLine, "%lf %lf %lf", &(f_Pert_x_QM[nConf_Mol_Ar][AtomIdx]), &(f_Pert_y_QM[nConf_Mol_Ar][AtomIdx]), &(f_Pert_z_QM[nConf_Mol_Ar][AtomIdx]));
				if(ReadItem != 3)	{
					fclose(fIn);
					Quit_With_Error_Msg("Error in reading data file for Ar perturbed forces.\nThe data for the force are not complete. \n");
				}
			}
			nConf_Mol_Ar ++;

			if(nConf_Mol_Ar >= MAX_N_PERT_MOL_AR)	{
				fclose(fIn);
				Quit_With_Error_Msg("nConf_Mol_Ar >= MAX_N_PERT_MOL_AR\nPlease increase MAX_N_PERT_MOL_AR.\nQuit\n");
			}
		}
	}

	fclose(fIn);
	i_Ref_Conf = nConf_Mol_Ar - 1;	// the last frame is used as the reference. Ar is supposed to be very far from the molecule.
	
	//start	to check the distance between Ar from all atoms in the molecule
	Idx_Ar = nAtom - 2;	// the Ar atom
	for(i=0; i<Idx_Ar; i++)	{
		if(pMol_Ar->IsLonePair[i])	{	// skip lonepair since their positions are not correct
			continue;
		}
		dx = Ar_Pert_x[i_Ref_Conf][i] - Ar_Pert_x[i_Ref_Conf][Idx_Ar];
		dy = Ar_Pert_y[i_Ref_Conf][i] - Ar_Pert_y[i_Ref_Conf][Idx_Ar];
		dz = Ar_Pert_z[i_Ref_Conf][i] - Ar_Pert_z[i_Ref_Conf][Idx_Ar];
		r = sqrt(dx*dx + dy*dy + dz*dz);
		if(r < dist_Min)	{
			dist_Min = r;
		}
	}
	if(dist_Min < R_CUTOF_AR_PROBE)	{
		Quit_With_Error_Msg("In the reference configuration of your force data, Ar atom is not far enough from the molecule.\nQuit\n");
	}
	//end	to check the distance between Ar from all atoms in the molecule

	for(iConf=0; iConf<nConf_Mol_Ar; iConf++)	{
		for(i=0; i<nAtom; i++)	{
			f_Pert_x_QM[iConf][i] -= f_Pert_x_QM[i_Ref_Conf][i];
			f_Pert_y_QM[iConf][i] -= f_Pert_y_QM[i_Ref_Conf][i];
			f_Pert_z_QM[iConf][i] -= f_Pert_z_QM[i_Ref_Conf][i];
			f_SQ = f_Pert_x_QM[iConf][i]*f_Pert_x_QM[iConf][i] + f_Pert_y_QM[iConf][i]*f_Pert_y_QM[iConf][i] + f_Pert_z_QM[iConf][i]*f_Pert_z_QM[iConf][i];
			if(f_SQ > HIGH_FORCE_CUTOFF_SQ)	{
				break;
			}
		}
		if(i<nAtom)	{
			Used_Ar_Pert[iConf] = 0;	// not used due to large forces
		}
		else	{
			Used_Ar_Pert[iConf] = 1;
		}
	}

	//start	to translate the molecule to make the CM at origin

	Idx_Real_Atom_Max = n_Real_Atom - 1;
	for(iConf=0; iConf<nConf_Mol_Ar; iConf++)	{
		xSum = ySum = zSum = 0.0;
		for(i=0; i<Idx_Real_Atom_Max; i++)	{
			AtomIdx = RealAtom_List[i];
			xSum += Ar_Pert_x[iConf][AtomIdx];
			ySum += Ar_Pert_y[iConf][AtomIdx];
			zSum += Ar_Pert_z[iConf][AtomIdx];
		}
		xSum/=Idx_Real_Atom_Max;
		ySum/=Idx_Real_Atom_Max;
		zSum/=Idx_Real_Atom_Max;
		for(i=0; i<nAtom; i++)	{
			Ar_Pert_x[iConf][i] -= xSum;
			Ar_Pert_y[iConf][i] -= ySum;
			Ar_Pert_z[iConf][i] -= zSum;
		}
	}

	//start	to for test only, the force projected along atom 3 -> 5 (H1->Ar)
//	int I_A=3, I_B=5;
	int I_A=3, I_B=5;
	double vec_x, vec_y, vec_z, norm, dist, f_proj;

	for(iConf=0; iConf<nConf_Mol_Ar; iConf++)	{
		vec_x = Ar_Pert_x[iConf][I_B] - Ar_Pert_x[iConf][I_A];
		vec_y = Ar_Pert_y[iConf][I_B] - Ar_Pert_y[iConf][I_A];
		vec_z = Ar_Pert_z[iConf][I_B] - Ar_Pert_z[iConf][I_A];
		dist = sqrt(vec_x*vec_x + vec_y*vec_y + vec_z*vec_z);
		norm = 1.0/dist;
		vec_x *= norm;
		vec_y *= norm;
		vec_z *= norm;

		f_proj = f_Pert_x_QM[iConf][I_A]*vec_x + f_Pert_y_QM[iConf][I_A]*vec_y + f_Pert_z_QM[iConf][I_A]*vec_z;

		printf("%6.4lf %6.4lf %6.4lf %10.7lf %10.7lf %10.7lf %10.7lf %10.7lf\n", vec_x, vec_y, vec_z, dist, f_proj, f_Pert_x_QM[iConf][I_A], f_Pert_y_QM[iConf][I_A], f_Pert_z_QM[iConf][I_A]);
	}

	//end	to for test only, the force projected along atom 3 -> 5 (H1->Ar)
}

void CFitting::ReadArPertEnergyData(void)
{
	FILE *fIn;
	char szLine[256];
	int i, iConf, nAtom, n_Real_Atom, ReadItem, Idx_Ar, RealAtom_List[MAX_NUM_PARA], AtomIdx;
	double dx, dy, dz, r, dist_Min=1.0E10;

	memset(Used_Ar_Pert, 0, sizeof(int)*MAX_N_PERT_MOL_AR);	// defaut: not used for fitting
	nAtom = pMol_Ar->nAtom;

	fIn = fopen(szFile_MolArEnergy, "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Error in opening the data file for Ar perturbed forces.\nQuit\n");
	}
	BuildRealAtomList(pMol_Ar, n_Real_Atom, RealAtom_List);	// include the Ar probe atom

	nConf_Mol_Ar = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		fgets(szLine, 256, fIn);
		if(strncmp(szLine, "Coordinate", 10) == 0)	{
			for(i=0; i<n_Real_Atom; i++)	{
				fgets(szLine, 256, fIn);
				AtomIdx = RealAtom_List[i];
				ReadItem = sscanf(szLine, "%lf %lf %lf", &(Ar_Pert_x[nConf_Mol_Ar][AtomIdx]), &(Ar_Pert_y[nConf_Mol_Ar][AtomIdx]), &(Ar_Pert_z[nConf_Mol_Ar][AtomIdx]));
				if(ReadItem != 3)	{
					fclose(fIn);
					Quit_With_Error_Msg("Error in reading data file for Ar perturbed forces.\nThe data for the coordinate are not complete. \n");
				}
				if(pMol_Ar->IsDrude[AtomIdx+1])	{	// a drude host
					Ar_Pert_x[nConf_Mol_Ar][AtomIdx+1] = Ar_Pert_x[nConf_Mol_Ar][AtomIdx];
					Ar_Pert_y[nConf_Mol_Ar][AtomIdx+1] = Ar_Pert_y[nConf_Mol_Ar][AtomIdx];
					Ar_Pert_z[nConf_Mol_Ar][AtomIdx+1] = Ar_Pert_z[nConf_Mol_Ar][AtomIdx];
				}
			}
			fgets(szLine, 256, fIn);
			if(strncmp(szLine, "Energy", 5) != 0)	{	// wrong format? 
				fclose(fIn);
				Quit_With_Error_Msg("Error in reading data file for Ar perturbed forces.\nIt is supposed to be the entry for forces. \n");
			}
			fgets(szLine, 256, fIn);
			ReadItem = sscanf(szLine, "%lf", &(E_Ar_Pert_QM[nConf_Mol_Ar]));
//			printf("%lf\n", E_Ar_Pert_QM[nConf_Mol_Ar]);
			if(ReadItem != 1)	{
				fclose(fIn);
				Quit_With_Error_Msg("Error in reading data file for Ar perturbed forces.\nThe data for the force are not complete. \n");
			}
			nConf_Mol_Ar ++;

			if(nConf_Mol_Ar >= MAX_N_PERT_MOL_AR)	{
				fclose(fIn);
				Quit_With_Error_Msg("nConf_Mol_Ar >= MAX_N_PERT_MOL_AR\nPlease increase MAX_N_PERT_MOL_AR.\nQuit\n");
			}
		}
	}

	fclose(fIn);
	i_Ref_Conf = nConf_Mol_Ar - 1;	// the last frame is used as the reference. Ar is supposed to be very far from the molecule.
	
	//start	to check the distance between Ar from all atoms in the molecule
	Idx_Ar = nAtom - 2;	// the Ar atom
	for(i=0; i<Idx_Ar; i++)	{
		if(pMol_Ar->IsLonePair[i])	{	// skip lonepair since their positions are not correct
			continue;
		}
		dx = Ar_Pert_x[i_Ref_Conf][i] - Ar_Pert_x[i_Ref_Conf][Idx_Ar];
		dy = Ar_Pert_y[i_Ref_Conf][i] - Ar_Pert_y[i_Ref_Conf][Idx_Ar];
		dz = Ar_Pert_z[i_Ref_Conf][i] - Ar_Pert_z[i_Ref_Conf][Idx_Ar];
		r = sqrt(dx*dx + dy*dy + dz*dz);
		if(r < dist_Min)	{
			dist_Min = r;
		}
	}
	if(dist_Min < R_CUTOF_AR_PROBE)	{
		Quit_With_Error_Msg("In the reference configuration of your force data, Ar atom is not far enough from the molecule.\nQuit\n");
	}
	//end	to check the distance between Ar from all atoms in the molecule

	for(iConf=0; iConf<nConf_Mol_Ar; iConf++)	{
		E_Ar_Pert_QM[iConf] -= E_Ar_Pert_QM[i_Ref_Conf];
		if(E_Ar_Pert_QM[iConf] > HIGH_ENERGY_CUTOFF)	{
			Used_Ar_Pert[iConf] = 0;	// not used due to large energy
		}
		else	{
			Used_Ar_Pert[iConf] = 1;
		}
	}
}

#define HIGH_INT_ENERGY_WATER_CUTOFF	(0.0)
void CFitting::ReadMolWaterIntEnergyData(void)
{
	FILE *fIn;
	char szLine[256];
	CMol *pMol;
	int i, iConf, nAtom, n_Real_Atom, ReadItem, RealAtom_List[MAX_NUM_PARA], AtomIdx;
	double dist_Min=1.0E10;

	pMol = pMol_Wat;
	memset(Used_Mol_Wat_Int, 0, sizeof(int)*MAX_N_PERT_MOL_WAT);	// defaut: not used for fitting
	nAtom = pMol->nAtom;

	fIn = fopen(szFile_MolWaterEnergy, "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Error in opening the data file for molecule-water coordinates and energies.\nQuit\n");
	}
	BuildRealAtomList(pMol, n_Real_Atom, RealAtom_List);	// include the water molecule

	nConf_Mol_Wat = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		fgets(szLine, 256, fIn);
		if(FindTags_As_First_String(szLine, "Coordinate"))	{
			for(i=0; i<n_Real_Atom; i++)	{
				fgets(szLine, 256, fIn);
				AtomIdx = RealAtom_List[i];
				ReadItem = sscanf(szLine, "%lf %lf %lf", &(Coord_Dimer_Wat_x[nConf_Mol_Wat][AtomIdx]), &(Coord_Dimer_Wat_y[nConf_Mol_Wat][AtomIdx]), &(Coord_Dimer_Wat_z[nConf_Mol_Wat][AtomIdx]));
				if(ReadItem != 3)	{
					fclose(fIn);
					Quit_With_Error_Msg("Error in reading data file for the dimer of molecule & water.\nThe data for the coordinate are not complete. \n");
				}
				if(pMol->mass[AtomIdx] > 2.0)	{	// It is a heavy atom, then it is a drude host
					Coord_Dimer_Wat_x[nConf_Mol_Wat][AtomIdx+1] = Coord_Dimer_Wat_x[nConf_Mol_Wat][AtomIdx];
					Coord_Dimer_Wat_y[nConf_Mol_Wat][AtomIdx+1] = Coord_Dimer_Wat_y[nConf_Mol_Wat][AtomIdx];
					Coord_Dimer_Wat_z[nConf_Mol_Wat][AtomIdx+1] = Coord_Dimer_Wat_z[nConf_Mol_Wat][AtomIdx];
				}
			}
			fgets(szLine, 256, fIn);
			if(FindTags_As_First_String(szLine, "Energy")==0)	{	// wrong format? 
				fclose(fIn);
				Quit_With_Error_Msg("Error in reading data file for molecule-water interaction energies.\nIt is supposed to be the entry for energy. \n");
			}
			fgets(szLine, 256, fIn);
			ReadItem = sscanf(szLine, "%lf", &(E_Int_Wat_QM[nConf_Mol_Wat]));
			if(ReadItem != 1)	{
				fclose(fIn);
				Quit_With_Error_Msg("Error in reading data file for molecule-water interaction energies.\nThe data for energy are not complete. \n");
			}
			nConf_Mol_Wat ++;

			if(nConf_Mol_Wat >= MAX_N_PERT_MOL_WAT)	{
				fclose(fIn);
				Quit_With_Error_Msg("nConf_Mol_Wat >= MAX_N_PERT_MOL_WAT\nPlease increase MAX_N_PERT_MOL_WAT.\nQuit\n");
			}
		}
	}

	fclose(fIn);

	for(iConf=0; iConf<nConf_Mol_Wat; iConf++)	{
		if(E_Int_Wat_QM[iConf] > HIGH_INT_ENERGY_WATER_CUTOFF)	{
			Used_Mol_Wat_Int[iConf] = 0;	// not used due to large energy
		}
		else	{
			Used_Mol_Wat_Int[iConf] = 1;
			nUsed_Conf_Mol_Wat_Int++;
		}
	}

}
#undef HIGH_INT_ENERGY_WATER_CUTOFF
#undef	R_CUTOF_AR_PROBE
#undef	HIGH_FORCE_CUTOFF
#undef	HIGH_ENERGY_CUTOFF
#undef	HIGH_FORCE_CUTOFF_SQ


//#define SCALE_QM_E_MIN	(1.16)
//#define SHIFT_QM_R_MIN	(-0.20)
//#define SHIFT_QM_R_MIN_CHARGED	(-0.10)
void CFitting::ReadMolWater_Emin_Rmin(void)
{
	FILE *fIn;
	char szLine[256], szTmp[256], *ReadLine, ErrorMsg[256];
	CMol *pMol;
	int i, nAtom, n_Real_Atom, RealAtom_List[MAX_NUM_PARA], Acceptor_Set[MAX_NUM_PARA], Donor_Set[MAX_NUM_PARA], ReadItem;
	double dist_Min=1.0E10;

	memset(Acceptor_Set, 0, sizeof(int)*MAX_NUM_PARA);
	memset(Donor_Set, 0, sizeof(int)*MAX_NUM_PARA);

	pMol =  pMol_ESP;
	nAtom = pMol->nAtom;

	fIn = fopen(szFile_MolWaterEnergy, "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Error in opening the data file for molecule-water coordinates and energies.\nQuit\n");
	}
	BuildRealAtomList(pMol, n_Real_Atom, RealAtom_List);	// include the water molecule
	n_Real_Atom--;	// skip the test charge

	N_Donor = N_Acceptor = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine != NULL)	{
			if(strncmp(szLine, "Acceptor", 8)==0)	{	// the entry for acceptor
				sscanf(szLine+8, "%d %lf %lf", &(Acceptor_List_Acceptor[N_Acceptor]), &(E_Min_QM_Acceptor[N_Acceptor]), &(r_Min_QM_Acceptor[N_Acceptor]));	//to determine a neighbor atom below
				Acceptor_List_Acceptor[N_Acceptor] = RealAtom_List[Acceptor_List_Acceptor[N_Acceptor]-1];
				if(fabs(netcharge_mol) > 0.2)	{	// charged species
					r_Min_QM_Acceptor[N_Acceptor] += SHIFT_QM_R_MIN_CHARGED;
				}
				else	{
					E_Min_QM_Acceptor[N_Acceptor] *= SCALE_QM_E_MIN;
					r_Min_QM_Acceptor[N_Acceptor] += SHIFT_QM_R_MIN;
				}
				ReadLine = fgets(szLine, 256, fIn);
				ReadItem = sscanf(szLine, "%s %d %lf %lf", szTmp, &(theta_Acceptor[N_Acceptor]), &(Para_Save_Acceptor[N_Acceptor][0]), &(Para_Save_Acceptor[N_Acceptor][1]));
				if(ReadItem != 4)	{
					sprintf(ErrorMsg, "Error in reading parameter from: %s\n%s\nQuit\n", szFile_MolWaterEnergy, szLine);
					Quit_With_Error_Msg(ErrorMsg);
				}
				Para_Save_Acceptor_Org[N_Acceptor][0] = Para_Save_Acceptor[N_Acceptor][0];	Para_Save_Acceptor_Org[N_Acceptor][1] = Para_Save_Acceptor[N_Acceptor][1];
				if(Acceptor_Set[Acceptor_List_Acceptor[N_Acceptor]] == 0)	{
					Acceptor_Set[Acceptor_List_Acceptor[N_Acceptor]] = 1;
					N_Acceptor++;
				}
			}
			else if(strncmp(szLine, "Donor", 5)==0)	{	// the entry for acceptor
				sscanf(szLine+5, "%d %lf %lf", &(Donor_List_Active[N_Donor]), &(E_Min_QM_Donor[N_Donor]), &(r_Min_QM_Donor[N_Donor]));	// to determine the atom connected with H below
				Donor_List_Active[N_Donor] = RealAtom_List[Donor_List_Active[N_Donor]-1];
				if(fabs(netcharge_mol) > 0.2)	{	// charged species
					r_Min_QM_Donor[N_Donor] += SHIFT_QM_R_MIN_CHARGED;
				}
				else	{
					E_Min_QM_Donor[N_Donor] *= SCALE_QM_E_MIN;
					r_Min_QM_Donor[N_Donor] += SHIFT_QM_R_MIN;
				}

				ReadLine = fgets(szLine, 256, fIn);
				ReadItem = sscanf(szLine, "%s %d %lf %lf", szTmp, &(theta_Donor[N_Donor]), &(Para_Save_Donor[N_Donor][0]), &(Para_Save_Donor[N_Donor][1]));
				if(ReadItem != 4)	{
					sprintf(ErrorMsg, "Error in reading parameter from: %s\n%s\nQuit\n", szFile_MolWaterEnergy, szLine);
					Quit_With_Error_Msg(ErrorMsg);
				}
				Para_Save_Donor_Org[N_Donor][0] = Para_Save_Donor[N_Donor][0];	Para_Save_Donor_Org[N_Donor][1] = Para_Save_Donor[N_Donor][1];	// saved

				if(Donor_Set[Donor_List_Active[N_Donor]] == 0)	{
					Donor_Set[Donor_List_Active[N_Donor]] = 1;
					N_Donor++;
				}
			}
			else	{
			}
		}
	}

	fclose(fIn);

	memcpy(theta_Acceptor_Save, theta_Acceptor, sizeof(int)*N_Acceptor);
	memcpy(theta_Donor_Save, theta_Donor, sizeof(int)*N_Donor);

	for(i=0; i<nAtom; i++)	{
		pMol_Wat->x[i] = pMol_ESP->x[i];
		pMol_Wat->y[i] = pMol_ESP->y[i];
		pMol_Wat->z[i] = pMol_ESP->z[i];
	}

	OptWater.Init(pMol_Wat);
}
#undef SCALE_QM_E_MIN
#undef SHIFT_QM_R_MIN
#undef SHIFT_QM_R_MIN_CHARGED

void CFitting::ReadParameters(char szName[])
{
	FILE *fIn;
	int i, nAtom;
	double *pCG, *pAlpha, *pThole;
	
//	fIn = fopen("para-opt-start.dat", "r");
	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		printf("Warning: %s does not exist. Using the charges in xpsf as the initial parameters.\n", szName);
		//		exit(1);
		// only for test
	}
	else	{
		for(i=0; i<n_Para; i++)	{
			fscanf(fIn, "%lf", &(Para_List[i]));
		}
		fclose(fIn);
		
		for(i=0; i<n_Para; i++)	{
			NLOpt_LowBound[i] = max(Para_List[i] - CG_Var_Range[i], LowBound[i]);
			NLOpt_UpBound[i] = min(Para_List[i] + CG_Var_Range[i], UpBound[i]);
		}


		DistributeData();	//update all parameters 
		fprintf(fFile_Run_Log, "Read the previous parameters as the initial parameters\n");
		fflush(fFile_Run_Log);
	}
	
	//start	to print parameters reading
	nAtom = pMol_ESP->nAtom - 1;
	pCG = pMol_ESP->CG;
	pAlpha = pMol_ESP->alpha;
	pThole = pMol_ESP->thole;
	for(i=0; i<nAtom; i++)	{
		if(pMol_ESP->IsDrude[i])	{	//skip drude
		}
		else	{
//			printf("%3d %4s ", i+1, pMol_ESP->AtomName[i]);
			fprintf(fFile_Run_Log, "%3d %4s ", i+1, pMol_ESP->AtomName[i]);
			fflush(fFile_Run_Log);
			if( ( (i+1) < nAtom) && (pMol_ESP->IsDrude[i+1]) )	{
//				printf("%10.7lf %10.7lf %10.7lf ", pCG[i]+pCG[i+1], pAlpha[i], pThole[i]);
				fprintf(fFile_Run_Log, "%10.7lf %10.7lf %10.7lf ", pCG[i]+pCG[i+1], pAlpha[i], pThole[i]);
				fflush(fFile_Run_Log);
			}
			else	{
//				printf("%10.7lf ", pCG[i]);
				fprintf(fFile_Run_Log, "%10.7lf ", pCG[i]);
				fflush(fFile_Run_Log);
			}
//			printf("\n");
			fprintf(fFile_Run_Log, "\n");
			fflush(fFile_Run_Log);
		}
	}
	//end	to print parameters reading
}

void CFitting::PrintInfo(void)
{
	int i, j, n_Dependent, AtomIdx;
	
	for(i=0; i<n_Type; i++)	{
		if(AtomType[i].IsRedundant == 0)	{	//parameters will be fitted
//			printf("%3d %5s  %3d  %4s ", 
			fprintf(fFile_Run_Log, "%3d %5s  %3d  %5s ", 
				i+1, AtomType[i].MolName, AtomType[i].AtomIdx+1, AtomType[i].pMol->AtomName[AtomType[i].AtomIdx]);
			if(AtomType[i].NotFitThole == 0)	{	//drude host, thole will be fitted
//				printf("   thole   ");
				fprintf(fFile_Run_Log, "   thole ");
			}
			else	{
//				printf("         ");
				fprintf(fFile_Run_Log, "         ");
			}
			
			if(AtomType[i].NotFitCharge == 0)	{	//charge will be fitted
//				printf(" CG-to-fit\n");
				fprintf(fFile_Run_Log, "   CG-to-fit\n");
			}
			else	{	//charge will not be fitted. 
//				printf(" CG = ");
				fprintf(fFile_Run_Log, "   CG = ");
				n_Dependent = AtomType[i].n_Atom_Neutral;
				for(j=0; j<n_Dependent; j++)	{
					AtomIdx = AtomType[AtomType[i].Cal_CG_List[j]].AtomIdx;
//					printf(" %2.0lf*(%s) ", AtomType[i].weight[j], AtomType[AtomType[i].Cal_CG_List[j]].pMol->AtomName[AtomIdx]);
					fprintf(fFile_Run_Log, " %4.1lf*(%s) ", AtomType[i].weight[j], AtomType[AtomType[i].Cal_CG_List[j]].pMol->AtomName[AtomIdx]);
				}
//				printf("\n");
				fprintf(fFile_Run_Log, "\n");
			}
		}
		else	{
//			printf("%3d %5s  %3d  %4s  redundant\n", 
			fprintf(fFile_Run_Log, "%3d %5s  %3d  %5s  redundant\n", 
				i+1, AtomType[i].MolName, AtomType[i].AtomIdx+1, AtomType[i].pMol->AtomName[AtomType[i].AtomIdx]);
		}
	}
	fflush(fFile_Run_Log);
}

int CFitting::FindTheRepresentative(int Idx)	//if Idx is an redundent atom type, find its representative; otherwise, just return Idx
{
	if(AtomType[Idx].IsRedundant)	{	// This is a redundent atom type
		return AtomType[Idx].CG_Represent;	//all equivalent atoms have same CG
	}
	else	{
		return Idx;
	}
}

int SplitString(char szBuff[], char ItemList[][32])
{
	int n_Item, ReadItem, Unfinished=1;
	int iPos, nLen;
	char szItem[256];
	
	n_Item = 0;
	iPos = 0;
	while(Unfinished)	{
		while( (szBuff[iPos] == ' ') || (szBuff[iPos] == '\t') )	{	//to find the first non-blank character
			if(szBuff[iPos] == 0)	{
				Unfinished = 0;
				break;
			}
			iPos++;
		}
		if(Unfinished == 0)	{
			break;
		}
		
		ReadItem = sscanf(szBuff+iPos, "%s", szItem);
		if(ReadItem != 1)	{
			break;
		}
		strcpy(ItemList[n_Item], szItem);
		nLen = strlen(szItem);
		iPos += nLen;
		n_Item++;
	}
	
	return n_Item;
}
#undef MAX_ITEM
#undef MAX_LEN_LINE


int CFitting::IndexAtomType(char szMolName[], int Idx)	// >=0 found this type; -1 - this is a new type
{
	int type;
	
	for(type=0; type<n_Type; type++)	{
		if( (strcmp(szMolName, AtomType[type].MolName)==0) && (Idx == AtomType[type].AtomIdx) )	{
			return type;
		}
	}
	
//	printf("Fail to find %d atom in molecule %s\nQuit.\n", Idx, szMolName);
	fprintf(fFile_Run_Log, "Fail to find %d atom in molecule %s\nQuit.\n", Idx, szMolName);
	fflush(fFile_Run_Log);
	exit(1);
	
	return (-1);	//not found. 
}

void CFitting::ReadDimerData(void)
{
	int i, j;
	FILE *fIn;
	double fTmp;
	char szTmp[256];

	fIn = fopen("./QM_ESP/ener-dimer.dat", "r");
	if(fIn == NULL)	{
//		printf("Fail to open file: ./QM_ESP/ener-dimer.dat\nQuit\n");
		fprintf(fFile_Run_Log, "Fail to open file: ./QM_ESP/ener-dimer.dat\nQuit\n");
		fflush(fFile_Run_Log);
		exit(1);
	}

	for(i=0; i<N_CONF_DIMER; i++)	{
		fscanf(fIn, "%lf %lf", &fTmp, &(E_Int_QM_Dimer[i]));
	}
	fclose(fIn);

	fIn = fopen("./QM_ESP/coor-dimer.dat", "r");
	if(fIn == NULL)	{
//		printf("Fail to open file: ./QM_ESP/coor-dimer.dat\nQuit\n");
		fprintf(fFile_Run_Log, "Fail to open file: ./QM_ESP/coor-dimer.dat\nQuit\n");
		fflush(fFile_Run_Log);
		exit(1);
	}
	for(i=0; i<N_CONF_DIMER; i++)	{
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][0 ]), &(y_Dimer[i][0 ]), &(z_Dimer[i][0 ]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][5 ]), &(y_Dimer[i][5 ]), &(z_Dimer[i][5 ]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][7 ]), &(y_Dimer[i][7 ]), &(z_Dimer[i][7 ]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][9 ]), &(y_Dimer[i][9 ]), &(z_Dimer[i][9 ]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][12]), &(y_Dimer[i][12]), &(z_Dimer[i][12]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][2 ]), &(y_Dimer[i][2 ]), &(z_Dimer[i][2 ]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][3 ]), &(y_Dimer[i][3 ]), &(z_Dimer[i][3 ]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][4 ]), &(y_Dimer[i][4 ]), &(z_Dimer[i][4 ]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][11]), &(y_Dimer[i][11]), &(z_Dimer[i][11]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][14]), &(y_Dimer[i][14]), &(z_Dimer[i][14]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][15]), &(y_Dimer[i][15]), &(z_Dimer[i][15]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][16]), &(y_Dimer[i][16]), &(z_Dimer[i][16]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][28]), &(y_Dimer[i][28]), &(z_Dimer[i][28]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][31]), &(y_Dimer[i][31]), &(z_Dimer[i][31]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][24]), &(y_Dimer[i][24]), &(z_Dimer[i][24]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][26]), &(y_Dimer[i][26]), &(z_Dimer[i][26]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][19]), &(y_Dimer[i][19]), &(z_Dimer[i][19]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][30]), &(y_Dimer[i][30]), &(z_Dimer[i][30]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][33]), &(y_Dimer[i][33]), &(z_Dimer[i][33]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][34]), &(y_Dimer[i][34]), &(z_Dimer[i][34]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][35]), &(y_Dimer[i][35]), &(z_Dimer[i][35]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][21]), &(y_Dimer[i][21]), &(z_Dimer[i][21]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][22]), &(y_Dimer[i][22]), &(z_Dimer[i][22]));
		fscanf(fIn, "%s %lf %lf %lf", szTmp, &(x_Dimer[i][23]), &(y_Dimer[i][23]), &(z_Dimer[i][23]));

		for(j=0; j<N_ATM_DIMER; j++)	{
			if(pNMA_Dimer->IsDrude[j])	{
				x_Dimer[i][j] = x_Dimer[i][j-1];
				y_Dimer[i][j] = y_Dimer[i][j-1];
				z_Dimer[i][j] = z_Dimer[i][j-1];
			}
		}
	}
	fclose(fIn);

}

int CFitting::To_Find_The_First_LP(int AtomIdx)
{
	int i, nLP;

	if(pMol_ESP == NULL)	{
		Quit_With_Error_Msg("pMol_ESP == NULL\nQuit.\n");
	}

	nLP = pMol_ESP->nLP;

	for(i=0; i<nLP; i++)	{
		if((pMol_ESP->LP_PosHost[i][0] == AtomIdx) && (fabs(pMol_ESP->LP_Dist[i]) > LOW_BOUND_LP_r) )	{	// not a dummy lone pair
			return i;
		}
	}

	fprintf(fFile_Run_Log, "Fail to find the associated lone pair for atom %d\nQuit\n", AtomIdx+1);
	exit(1);
	return -1;
}



CMol Mol_ESP, NMA_Single, NMA_Dimer, Mol_Ar, Mol_Water;
CForceField ForceField;

int main(int argc, char **argv)
{

  timebomb();

	int Time_1, Time_2;

	fFile_Run_Log = fopen(szName_Run_Log, "w");
	if(fFile_Run_Log==NULL)	{
		printf("Fail to create the log file.\n\n");
		exit(1);
	}

	if(argc == 2)	{
		strcpy(szName_Conf_File, argv[1]);
	}
	else if(argc == 3)	{
		if(strcmp(argv[2], "TEST_MODE") == 0)	{
			strcpy(szName_Conf_File, argv[1]);
			IsTestMode = 1;
		}
		else	{
			fprintf(fFile_Run_Log, "Unrecognize parameter: %s\nQuit\n", argv[2]);
			fflush(fFile_Run_Log);
			exit(1);
		}
	}
	else	{
		printf("Usage: fit fit.conf\n");
		fprintf(fFile_Run_Log, "Usage: fit fit.conf\n");
		fflush(fFile_Run_Log);
		exit(1);
	}

	ReadConfFile(szName_Conf_File);
/*
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProgID);
	MPI_Comm_size(MPI_COMM_WORLD, &nProc);
*/
	Time_1 = time(NULL);

	ForceField.ReadForceField(szName_Force_Field);
	fitting.myForceField = &ForceField;
	
	Mol_ESP.ReadPSF(szName_XPSF, 0);
	Get_Netcharge_From_Xpsf();
	Read_Rtf_File();

	for(int i=0; i<Mol_ESP.nDihedral; i++)	{
		Mol_ESP.Is_Phi_Constrained[i] = 1;
	}
	
	strcpy(Mol_ESP.MolName, Mol_ESP.ResName[0]);
	Mol_ESP.AssignForceFieldParameters(&ForceField);
//	Mol_ESP.ReadCRD("glyc.crd");
//	Mol_ESP.ReadCRD(szName_CRD);

	Mol_ESP.Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
	Mol_ESP.Is_Phi_Psi_Constrained = 0;
	Mol_ESP.E_CMap_On = 0;

//	Mol_ESP.Cal_E();
//	Mol_ESP.FullGeometryOptimization_LBFGS_drude();

	
	fitting.pMol_ESP = &Mol_ESP;

	fitting.RegisterAtomType(&Mol_ESP, Mol_ESP.nAtom-1);	//skip the test charge
	fitting.ProcessAtomTypeConstraint();

	if(Not_Fit_Elec_Para == 0)	{
		fitting.ReadESPData();	//read five structures of alad
	}

	if(Target_E_Int_Dimer)	{
		//start	to read psf file needed for dimer interactions
		NMA_Single.ReadPSF("nma-single.xpsf", 0);	// 37 atoms including the test charge
		NMA_Single.AssignForceFieldParameters(&ForceField);
		NMA_Single.Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
		NMA_Single.Is_Phi_Psi_Constrained = 0;
		NMA_Single.E_CMap_On = 0;
		fitting.pNMA_Single = &NMA_Single;

		NMA_Dimer.ReadPSF("nma-dimer.xpsf", 1);	// 37 atoms including the test charge
		NMA_Dimer.AssignForceFieldParameters(&ForceField);
		NMA_Dimer.Geo_Opt_Drude_Only = 1;	//only drude will be relaxed
		NMA_Dimer.Is_Phi_Psi_Constrained = 0;
		NMA_Dimer.E_CMap_On = 0;
		fitting.pNMA_Dimer = &NMA_Dimer;

		fitting.ReadDimerData();
		//end	to read psf file needed for dimer interactions
	}

	if(Target_E_Int_Water)	{
		fitting.pMol_Wat = &Mol_Water;
		Mol_Water.ReadPSF(szFile_MolWater_PSF, 1);
		Mol_Water.AssignForceFieldParameters(&ForceField);

		fitting.ReadMolWater_Emin_Rmin();
//		fitting.ReadMolWaterIntEnergyData();
		Mol_Water.Geo_Opt_Drude_Only = 1;	//only drude will be relaxed

		FILE *fIn;
		fIn = fopen("saved-para.dat", "r");
		if(fIn != NULL)	{
			fclose(fIn);
			Read_Soft_DihedralList();
			Read_Fitted_Torsion_Parameters();
		}
	}

	//start	to set up the psf of molecule plus Ar probe atom
	if(To_Fit_LJ)	{
		Mol_Ar.ReadPSF(szFile_MolAr_PSF, 1);
		Mol_Ar.AssignForceFieldParameters(&ForceField);
		fitting.pMol_Ar = &Mol_Ar;

//		fitting.ReadArPertForceData();
		fitting.ReadArPertEnergyData();
	}
	//end	to set up the psf of molecule plus Ar probe atom


	
	fitting.Optimization_NLOpt();
	
	
	Time_2 = time(NULL);
	
//	printf("Time = %d \n", Time_2-Time_1);
	fprintf(fFile_Run_Log, "Time = %d \n", Time_2-Time_1);
	fflush(fFile_Run_Log);
	
	fclose(fFile_Run_Log);

//	MPI_Finalize();


	return 0;
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

void Para_Aniso_K_A(double K11, double K22, double K33, double& A11, double& A22)
{
	K11 += (K_DRUDE + K33);
	K22 += (K_DRUDE + K33);
	K33 += (K_DRUDE);

	K11 /= K_DRUDE;
	K22 /= K_DRUDE;
	K33 /= K_DRUDE;

	A11 = 1.0/K11;
	A22 = 1.0/K22;
}

void Para_Aniso_A_K(double& K11, double& K22, double& K33, double A11, double A22)
{
	K11 = 1.0/A11;
	K22 = 1.0/A22;
	K33 = 1.0/(3.0-A11-A22);

	K11 *= K_DRUDE;
	K22 *= K_DRUDE;
	K33 *= K_DRUDE;

	K33 -= K_DRUDE;
	K11 -= (K_DRUDE+K33);
	K22 -= (K_DRUDE+K33);
}

void To_Upper_Case(char szBuff[], char szNew[])
{
	char gap;
	int i=0;

	gap = 'A'-'a';

	while(1)	{
		if(szBuff[i]==0)	{
			break;
		}
		if( (szBuff[i]>='a') && (szBuff[i]<='z') )	{
			szNew[i] = szBuff[i] + gap;
		}
		else	{
			szNew[i] = szBuff[i];
		}
		i++;
	}
	szNew[i]=0;
}

int CFitting::FindAtomsinMol(char szAtmName[], int IdxList[])
{
	int i, nLen, nMatch=0;
	char szAtmName_Upper[256];

	if(pMol_ESP == NULL)	{
//		printf("pMol_ESP == NULL\nQuit\n");
		fprintf(fFile_Run_Log, "pMol_ESP == NULL\nQuit\n");
		fflush(fFile_Run_Log);
		exit(1);
	}

	To_Upper_Case(szAtmName, szAtmName_Upper);
	nLen = strlen(szAtmName_Upper);

	if(szAtmName[nLen-1] == '*')	{	//wild character
		for(i=0; i<n_Type; i++)	{
			if(strncmp(pMol_ESP->AtomName[i], szAtmName_Upper, nLen-1)==0)	{
				IdxList[nMatch] = i;	// find one atom with matching name
				nMatch++;
			}
		}
	}
	else	{
		for(i=0; i<n_Type; i++)	{
			if(strcmp(pMol_ESP->AtomName[i], szAtmName_Upper)==0)	{
				IdxList[nMatch] = i;	// find one atom with matching name
				nMatch++;
			}
		}
	}

	if(nMatch == 0)	{
//		printf("nMatch == 0 for atom %s\nQuit\n", szAtmName);
		fprintf(fFile_Run_Log, "nMatch == 0 for atom %s\nQuit\n", szAtmName);
		fflush(fFile_Run_Log);
		exit(1);
	}

	return nMatch;
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
	int LineLen, i, ReadItem;
	FILE *fIn;
	char szValue[256], szBuff[256];

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
//		else if(strncmp(szReadConf_File[i], "FILE_XYZ", 8) == 0)	{
//			GetSecondItem(szReadConf_File[i], szName_CRD);
//		}
		else if(FindTags_As_First_String(szReadConf_File[i], "FILE_pot0"))	{
			GetSecondItem(szReadConf_File[i], szFile_pot0);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "FILE_pot"))	{
			GetSecondItem(szReadConf_File[i], szFile_pot);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "FILE_qpos"))	{
			GetSecondItem(szReadConf_File[i], szFile_qpos);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "FILE_MolAr_PSF"))	{
			GetSecondItem(szReadConf_File[i], szFile_MolAr_PSF);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "FILE_MolArForce"))	{
			GetSecondItem(szReadConf_File[i], szFile_MolArForce);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "FILE_MolArEnergy"))	{
			GetSecondItem(szReadConf_File[i], szFile_MolArEnergy);
		}

		else if(FindTags_As_First_String(szReadConf_File[i], "FILE_MolWater_PSF"))	{
			GetSecondItem(szReadConf_File[i], szFile_MolWater_PSF);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "FILE_MolWaterEnergy"))	{
			GetSecondItem(szReadConf_File[i], szFile_MolWaterEnergy);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "w_H_Donor_Acceptor"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			w_H_Donor_Acceptor = atof(szValue);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "w_charge"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			w_charge = atof(szValue);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "w_alpha"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			w_alpha = atof(szValue);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "w_thole"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			w_thole = atof(szValue);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "w_anisotropy"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			w_anisotropy = atof(szValue);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "Scale_Alpha"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			Scale_Alpha = atof(szValue);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "Scale_TargetAlpha"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			Scale_TargetAlpha = atof(szValue);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "FAR_DRUDE_R"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			FAR_DRUDE_R = atof(szValue);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "w_water_E_min"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			w_water_E_min = atof(szValue);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "w_water_R_min"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			w_water_R_min = atof(szValue);
		}

		else if(FindTags_As_First_String(szReadConf_File[i], "SCALE_QM_E_MIN_Drude"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			SCALE_QM_E_MIN = atof(szValue);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "SHIFT_QM_R_MIN_Drude"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			SHIFT_QM_R_MIN = atof(szValue);
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "SHIFT_QM_R_MIN_CHARGED_Drude"))	{
			GetSecondItem(szReadConf_File[i], szValue);
			SHIFT_QM_R_MIN_CHARGED = atof(szValue);
		}

		else if(FindTags_As_First_String(szReadConf_File[i], "To_Use_Miller_Alpha"))	{
			To_Use_Miller_Alpha = 1;
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "To_Fit_Aniso"))	{
			To_Fit_Aniso = 1;
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "To_Fit_LP_Para"))	{
			To_Fit_LP_Para= 1;
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "To_Fit_LJ"))	{
			To_Fit_LJ = 1;
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "Target_E_Int_Dimer"))	{
			Target_E_Int_Dimer = 1;
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "Target_E_Int_Water"))	{
			Target_E_Int_Water = 1;
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "Not_Fit_Elec_Para"))	{
			Not_Fit_Elec_Para = 1;
		}
		else if(FindTags_As_First_String(szReadConf_File[i], "DIPOLE_QM"))	{
			ReadItem = sscanf(szReadConf_File[i], "%s%lf%lf%lf%lf", szBuff, &QM_Dipole_x, &QM_Dipole_y, &QM_Dipole_z, &QM_Dipole);
			if(ReadItem != 5)	{	// something wrong
				QM_Dipole_x = QM_Dipole_y = QM_Dipole_z = QM_Dipole = 0.0;	// not valid dipole value
				printf("Warning> error in reading QM dipole.\n%s\n", szReadConf_File[i]);
			}
		}
		
		
	}
	//end	to assign some parameters 

	return n_Line_Conf_File;
}

int FindItemInConfFile(char szItem[], int LineStart)
{
	int Idx, Len;

	Len = strlen(szItem);

	for(Idx=LineStart; Idx<n_Line_Conf_File; Idx++)	{
		if(strncmp(szReadConf_File[Idx], szItem, Len)==0)	{
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
		if( (FindString(szTxt_Rtf[i], "RESI")>=0) && (FindString(szTxt_Rtf[i], "MOL")>=0) )	{
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

	nAnisotropy = 0;
	for(i=0; i<nLine_Rtf; i++)	{
		if(strncmp(szTxt_Rtf[i], "ANISOTROPY", 10)==0)	{
			Aniso_List[nAnisotropy] = i;
			nAnisotropy++;
		}
	}


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


void Quit_With_Error_Msg(char szMsg[])
{
	fprintf(fFile_Run_Log, "%s", szMsg);
	fflush(fFile_Run_Log);
	exit(1);
}


void Get_Netcharge_From_Xpsf(void)
{
	int i, nAtom;

	netcharge_mol = 0.0;
	nAtom = Mol_ESP.nAtom - 1;	// skip the test charge
	for(i=0; i<nAtom; i++)	{
		netcharge_mol += Mol_ESP.CG[i];
	}

	if(netcharge_mol > 0.3)	{
		netcharge_mol = 1.0*int(netcharge_mol + 0.2);	// assuming the netcharge is an integer !!!
	}
	else if(netcharge_mol < -0.3)	{
		netcharge_mol = 1.0*int(netcharge_mol - 0.2);	// assuming the netcharge is an integer !!!
	}
	else	{
		netcharge_mol = 0.0;
	}
}



void COptimizeWaterDimer::Gen_xyz(double x_A, double y_A, double z_A, 
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
	Normalize_drude(vz_x, vz_y, vz_z);	

	BA_x = x_A - x_B;
	BA_y = y_A - y_B;
	BA_z = z_A - z_B;
	Normalize_drude(BA_x, BA_y, BA_z);	

	vy_x = vz_y*BA_z - vz_z*BA_y;
	vy_y = vz_z*BA_x - vz_x*BA_z;
	vy_z = vz_x*BA_y - vz_y*BA_x;
	Normalize_drude(vy_x, vy_y, vy_z);

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

void COptimizeWaterDimer::GenerateWater(void)	// with six degree of freedom for water molecule
{
	double O_Wat_x, O_Wat_y, O_Wat_z, H1_Wat_x, H1_Wat_y, H1_Wat_z, H2_Wat_x, H2_Wat_y, H2_Wat_z, r;
	int Idx_O_Wat, Idx_H1_Wat, Idx_H2_Wat;

	Idx_O_Wat = fitting.pMol_ESP->nAtom-1;
	Idx_H1_Wat = Idx_O_Wat+3;
	Idx_H2_Wat = Idx_O_Wat+4;

	r = Para_List[0];

	if(IsAcceptor)	{
		H1_Wat_x = x[ActiveAtom] + r*v_x_X1_Acceptor;
		H1_Wat_y = y[ActiveAtom] + r*v_y_X1_Acceptor;
		H1_Wat_z = z[ActiveAtom] + r*v_z_X1_Acceptor;
		
		x_Dummy_3 = H1_Wat_x + (x_Dummy_2 - x[ActiveAtom]);
		y_Dummy_3 = H1_Wat_y + (y_Dummy_2 - y[ActiveAtom]);
		z_Dummy_3 = H1_Wat_z + (z_Dummy_2 - z[ActiveAtom]);
		
		Gen_xyz(x_Dummy_2, y_Dummy_2, z_Dummy_2, x_Dummy_3, y_Dummy_3, z_Dummy_3, 
			H1_Wat_x, H1_Wat_y, H1_Wat_z, b0_H_O_EXP, PI*0.5, PI, O_Wat_x, O_Wat_y, O_Wat_z);
		
		Gen_xyz(x_Dummy_3, y_Dummy_3, z_Dummy_3, H1_Wat_x, H1_Wat_y, H1_Wat_z, 
			O_Wat_x, O_Wat_y, O_Wat_z, b0_H_O_EXP, theta0_H_O_H_EXP*radianInv, Para_List[1], H2_Wat_x, H2_Wat_y, H2_Wat_z);
		
		pDimer->x[Idx_O_Wat] = O_Wat_x;		pDimer->y[Idx_O_Wat] = O_Wat_y;		pDimer->z[Idx_O_Wat] = O_Wat_z;
		pDimer->x[Idx_H1_Wat] = H1_Wat_x;	pDimer->y[Idx_H1_Wat] = H1_Wat_y;	pDimer->z[Idx_H1_Wat] = H1_Wat_z;
		pDimer->x[Idx_H2_Wat] = H2_Wat_x;	pDimer->y[Idx_H2_Wat] = H2_Wat_y;	pDimer->z[Idx_H2_Wat] = H2_Wat_z;
	}
	else	{
		O_Wat_x = x[ActiveHAtom] + r*v_x_X1_H;
		O_Wat_y = y[ActiveHAtom] + r*v_y_X1_H;
		O_Wat_z = z[ActiveHAtom] + r*v_z_X1_H;
		
		Gen_xyz(x_Dummy_2, y_Dummy_2, z_Dummy_2, x[ActiveHAtom], y[ActiveHAtom], z[ActiveHAtom], 
			O_Wat_x, O_Wat_y, O_Wat_z, b0_H_O_EXP, PI-0.5*theta0_H_O_H_EXP*radianInv, Para_List[1], H1_Wat_x, H1_Wat_y, H1_Wat_z);
		
		Gen_xyz(x_Dummy_3, y_Dummy_3, z_Dummy_3, H1_Wat_x, H1_Wat_y, H1_Wat_z, 
			O_Wat_x, O_Wat_y, O_Wat_z, b0_H_O_EXP, theta0_H_O_H_EXP*radianInv, 0.0, H2_Wat_x, H2_Wat_y, H2_Wat_z);
		
		pDimer->x[Idx_O_Wat] = O_Wat_x;		pDimer->y[Idx_O_Wat] = O_Wat_y;		pDimer->z[Idx_O_Wat] = O_Wat_z;
		pDimer->x[Idx_H1_Wat] = H1_Wat_x;	pDimer->y[Idx_H1_Wat] = H1_Wat_y;	pDimer->z[Idx_H1_Wat] = H1_Wat_z;
		pDimer->x[Idx_H2_Wat] = H2_Wat_x;	pDimer->y[Idx_H2_Wat] = H2_Wat_y;	pDimer->z[Idx_H2_Wat] = H2_Wat_z;
	}

	//start	to setup the coordinates for the drude
	pDimer->x[Idx_O_Wat+1] = pDimer->x[Idx_O_Wat];
	pDimer->y[Idx_O_Wat+1] = pDimer->y[Idx_O_Wat];
	pDimer->z[Idx_O_Wat+1] = pDimer->z[Idx_O_Wat];
	//end	to setup the coordinates for the drude

	pDimer->Position_LonePair();

//	pDimer->SavePdb("tmp.pdb");
}



double COptimizeWaterDimer::CalObjectiveFunction_Mol_Wat_Dimer_E_Int(void)	// only used for the optimization in MM
{
	GenerateWater();
	pDimer->Geo_Opt_Drude_Only = 1;	// only optimize the positions of drude
	pDimer->FullGeometryOptimization_LBFGS_drude(1);	// the total energy of molecule and water

	pDimer->Geo_Opt_Drude_Only = 0;	// only optimize the positions of drude
	E_Dimer_Int = pDimer->Cal_E(0);

//	pDimer->SavePdb("opt.pdb");

	return E_Dimer_Int;
}

#define DELTA_PARA_WAT	(2.0E-4)
double COptimizeWaterDimer::CalGradient_Wat_Opt(void)
{
	int i;
	double Para_Save[N_PARA_WAT];
	double f_Left, f_Right;
	
	for(i=0; i<N_PARA_WAT; i++)	{
		Para_Save[i] = Para_List[i];	//save current parameters
	}

	for(i=0; i<N_PARA_WAT; i++)	{
		Para_List[i] = Para_Save[i] - DELTA_PARA_WAT;
		f_Left = CalObjectiveFunction_Mol_Wat_Dimer_E_Int();
		
		Para_List[i] = Para_Save[i] + DELTA_PARA_WAT;
		f_Right = CalObjectiveFunction_Mol_Wat_Dimer_E_Int();
		
		Grad_List[i] = (f_Right-f_Left)/(2.0*DELTA_PARA_WAT);
//		printf("%10.5lf\n", Grad_List[i]);
		
		Para_List[i] = Para_Save[i];	//restore save (correct) parameter
	}
	
	CalObjectiveFunction_Mol_Wat_Dimer_E_Int();

	return E_Dimer_Int;
}
#undef DELTA_PARA_WAT

int FailCount_Wat_opt, Iteration_Wat_Opt;
double Callback_Eval_Gradient_Opt_Wat(unsigned n, const double *x, double *grad, void *my_func_data)
{
	int i, n_Local;

	n_Local = (int)n;
	for(i=0; i<n_Local; i++)	{
		OptWater.Para_List[i] = x[i];
	}

	if(grad)	{
		OptWater.CalGradient_Wat_Opt();	//gradient will be assigned into objGrad automatically
		Iteration_Wat_Opt++;
		for(i=0; i<n_Local; i++)	{
			grad[i] = OptWater.Grad_List[i];
		}
//		printf("Iteration %4d  E = %16.8lf\n", Iteration, E_Dimer_Int);

		if(OptWater.E_Dimer_Int < OptWater.E_Min)	{
			OptWater.E_Min = OptWater.E_Dimer_Int;
			FailCount_Wat_opt = 0;

			memcpy(OptWater.Para_Best, x, sizeof(double)*n);
		}
		else	{
//			if(FailCount_Wat_opt > 4)      {       // cannot converge due to limited accuracy of numerical gradient
//				printf("Too much failed tries.\nI have tried my best.\nQuit.\n");
//				memset(grad, 0, sizeof(double)*n);
//			}
			FailCount_Wat_opt++;
		}
	}
	else	{	// cal object function only
		OptWater.CalObjectiveFunction_Mol_Wat_Dimer_E_Int();
	}

    return OptWater.E_Dimer_Int;
}

int COptimizeWaterDimer::IsGradOptimized(double ftol)
{
	int i;

	for(i=0; i<N_PARA_WAT; i++)	{
		if(fabs(Grad_List[i]) > ftol)	{	// check the gradient
			if( (fabs(Low_Bound[i]-Para_List[i]) > 1.0E-3) && (fabs(Up_Bound[i]-Para_List[i]) > 1.0E-3)  )	{	// check touching boundaries or not
				return 0;
			}
		}
	}
	return 1;
}

double COptimizeWaterDimer::Optimize_Water_Pose(double Para_In[], double& E_Int_Min, double& r_Min)
{
	nlopt_opt opt;
	double E_Dimer_Min=1.0E100;
	int RetCode;

	memcpy(Para_List, Para_In, sizeof(double)*N_PARA_WAT);

	Low_Bound[0] = Lower_Bound_R_H_Bond;
	Up_Bound[0] = Upper_Bound_R_H_Bond;
	Low_Bound[1] = -10.0;
	Up_Bound[1] = 10.0;
// ELIOT
//	NLOPT_max_step_size = 0.2;
//	opt = nlopt_create(NLOPT_LD_LBFGS_GEO, N_PARA_WAT); /* algorithm and dimensionality */
	opt = nlopt_create(NLOPT_LD_LBFGS, N_PARA_WAT); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, Low_Bound);
	nlopt_set_upper_bounds(opt, Up_Bound);
	nlopt_set_min_objective(opt, Callback_Eval_Gradient_Opt_Wat, NULL);
	nlopt_set_ftol_rel(opt, 5.0E-10);

	E_Min=1.0E20;
	CalObjectiveFunction_Mol_Wat_Dimer_E_Int();

	FailCount_Wat_opt = 0;
	Iteration_Wat_Opt = 0;

	memcpy(Para_Best, Para_List, sizeof(double)*N_PARA_WAT);
	RetCode = nlopt_optimize(opt, Para_List, &E_Dimer_Min);
	if ( (RetCode != NLOPT_FTOL_REACHED) && (IsGradOptimized(1.0E-6) == 0) ) {
//	if ( (RetCode != NLOPT_FTOL_REACHED) && (RetCode != NLOPT_SUCCESS) ) {	// do optimization again
// ELIOT
//		NLOPT_max_step_size = 0.3;
		RetCode = nlopt_optimize(opt, Para_List, &E_Dimer_Min);
		if( (RetCode != NLOPT_FTOL_REACHED) && (RetCode != NLOPT_SUCCESS) )	{	// decrease the maximum step size and do optimization again
			printf("nlopt failed!\n");
		}
	}
	else {
//		printf("Obj_Min = %16.10lf\n", E_Dimer_Min);
	}

//	pDimer->SavePdb("opt.pdb");


	nlopt_destroy(opt);

	if(Para_Best[1] < -PI2)	{
		Para_Best[1] += PI2;
	}
	else if(Para_Best[1] > PI2)	{
		Para_Best[1] -= PI2;
	}

	memcpy(Para_List, Para_Best, sizeof(double)*N_PARA_WAT);
	memcpy(Para_In, Para_Best, sizeof(double)*N_PARA_WAT);
	E_Dimer_Min = CalObjectiveFunction_Mol_Wat_Dimer_E_Int();

	E_Int_Min = E_Dimer_Min;
	r_Min = Para_Best[0];	// r_min

//	pDimer->SavePdb("tmp.pdb");

	return E_Dimer_Min;
}

double Cal_Dist_Two_Atoms(int ia, int ib)
{
	double dx, dy, dz, r;
	CMol *Dimer;

	Dimer = fitting.pMol_Wat;
	dx = Dimer->x[ia] - Dimer->x[ib];
	dy = Dimer->y[ia] - Dimer->y[ib];
	dz = Dimer->z[ia] - Dimer->z[ib];
	r = sqrt(dx*dx + dy*dy + dz*dz);
	return r;
}


double Normalize_drude(double& vect_x, double& vect_y, double& vect_z)
{
	double dist, dist_Inv;

	dist = sqrt(vect_x*vect_x + vect_y*vect_y + vect_z*vect_z);
	if(dist < 1.0E-100)	{
		printf("Warning> Three atoms used in dihedral calculation are colinear!\n");
	}
	dist_Inv = 1.0/dist;
	vect_x *= dist_Inv;
	vect_y *= dist_Inv;
	vect_z *= dist_Inv;
	return dist;
}

void COptimizeWaterDimer::Init(CMol* pMol)
{
	int i;

	pDimer = pMol;
	x = pDimer->x;
	y = pDimer->y;
	z = pDimer->z;

	nAtom = fitting.pMol_ESP->nAtom;
	for(i=0; i<nAtom; i++)	{
		if( (pDimer->mass[i] > 0.99) && (pDimer->mass[i] < 1.01) )	{	// H atom
			Is_H[i] = 1;
		}
		else	{
			Is_H[i] = 0;
		}
	}

	for(i=0; i<nAtom; i++)	{
		GetBondList(i);	// setup Bond_Count[] and Bonded_Atom[]
	}

}

void COptimizeWaterDimer::GetBondList(int Idx)	// get the bond list for real atoms
{
	int i, nBond, iPos, *pBondList, Atom_2;
	double *mass;

	nBond = pDimer->nBond;
	pBondList = pDimer->BondList;
	mass = pDimer->mass;

	Bonded_Atom[Idx] = -1;
	Bond_Count[Idx] = 0;
	for(i=0; i<nBond; i++)	{
		iPos = 2*i;
		if(pBondList[iPos] == Idx)	{
			if( (mass[pBondList[iPos]] > 0.9) && (mass[pBondList[iPos+1]] > 0.9) )	{	// a bond connecting two real atoms
				Atom_2 = pBondList[iPos+1];
				BondList[Idx][Bond_Count[Idx]] = Atom_2;
				Bond_Count[Idx]++;

				if( (Bonded_Atom[Idx] < 0) && Is_H[Idx] )	{
					Bonded_Atom[Idx] = Atom_2;
				}
			}
		}
		if(pBondList[iPos+1] == Idx)	{
			if( (mass[pBondList[iPos]] > 0.9) && (mass[pBondList[iPos+1]] > 0.9) )	{	// a bond connecting two real atoms
				Atom_2 = pBondList[iPos];
				BondList[Idx][Bond_Count[Idx]] = Atom_2;
				Bond_Count[Idx]++;
				if( (Bonded_Atom[Idx] < 0) && Is_H[Idx] )	{
					Bonded_Atom[Idx] = Atom_2;
				}
			}
		}
	}
	return;
}

void COptimizeWaterDimer::Generate_Dummy_Atoms(int Idx, int bAcceptor, int theta)
{
	ActiveAtom = Idx;
	if(bAcceptor == 0)	{	// H donor
		IsAcceptor = 0;
		ActiveHAtom = Idx;
		ActiveAtom = Bonded_Atom[ActiveHAtom];
		Generate_Dummy_Atoms_Donor();
	}
	else	{
		IsAcceptor = 1;
		if(Bond_Count[ActiveAtom] == 1)	{
			Generate_Dummy_Atoms_Acceptor_1();
		}
		else if(Bond_Count[ActiveAtom] == 2)	{
			Generate_Dummy_Atoms_Acceptor_2();
		}
		else if(Bond_Count[ActiveAtom] == 3)	{
			Generate_Dummy_Atoms_Acceptor_3();
		}
		else {	// unexpected
		}
	}
	if(theta != 0)	{
		Rotate_Dummy_Atoms(theta);
	}
}

void COptimizeWaterDimer::Rotate_Dummy_Atoms(int theta)
{
	if(IsAcceptor)	{	// acceptor
		double v1_New_x, v1_New_y, v1_New_z;
		double v2_x, v2_y, v2_z, v2_New_x, v2_New_y, v2_New_z;
		double v3_x, v3_y, v3_z;
		double RotM[3][3], cos_theta, sin_theta;
		
		cos_theta = cos(theta*radianInv);
		sin_theta = sin(theta*radianInv);
		
		v_x_X1_Acceptor = x[ActiveAtom] - x_Dummy_1;	// X vector
		v_y_X1_Acceptor = y[ActiveAtom] - y_Dummy_1;
		v_z_X1_Acceptor = z[ActiveAtom] - z_Dummy_1;
		Normalize_drude(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor);
		
		v2_x = x[ActiveAtom] - x_Dummy_2;
		v2_y = y[ActiveAtom] - y_Dummy_2;
		v2_z = z[ActiveAtom] - z_Dummy_2;
		Normalize_drude(v2_x, v2_y, v2_z);
		
		Cross_Product(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z);
		Normalize_drude(v3_x, v3_y, v3_z);
		
		RotM[0][0] = v3_x*v3_x + (1-v3_x*v3_x)*cos_theta;
		RotM[0][1] = v3_x*v3_y*(1.0-cos_theta) - v3_z*sin_theta;
		RotM[0][2] = v3_x*v3_z*(1.0-cos_theta) + v3_y*sin_theta;
		
		RotM[1][0] = v3_x*v3_y*(1.0-cos_theta) + v3_z*sin_theta;
		RotM[1][1] = v3_y*v3_y + (1-v3_y*v3_y)*cos_theta;
		RotM[1][2] = v3_y*v3_z*(1.0-cos_theta) - v3_x*sin_theta;
		
		RotM[2][0] = v3_x*v3_z*(1.0-cos_theta) - v3_y*sin_theta;
		RotM[2][1] = v3_y*v3_z*(1.0-cos_theta) + v3_x*sin_theta;
		RotM[2][2] = v3_z*v3_z + (1-v3_z*v3_z)*cos_theta;
		
		
		v1_New_x = RotM[0][0] * v_x_X1_Acceptor + RotM[0][1] * v_y_X1_Acceptor + RotM[0][2] * v_z_X1_Acceptor;
		v1_New_y = RotM[1][0] * v_x_X1_Acceptor + RotM[1][1] * v_y_X1_Acceptor + RotM[1][2] * v_z_X1_Acceptor;
		v1_New_z = RotM[2][0] * v_x_X1_Acceptor + RotM[2][1] * v_y_X1_Acceptor + RotM[2][2] * v_z_X1_Acceptor;
		
		v2_New_x = RotM[0][0] * v2_x + RotM[0][1] * v2_y + RotM[0][2] * v2_z;
		v2_New_y = RotM[1][0] * v2_x + RotM[1][1] * v2_y + RotM[1][2] * v2_z;
		v2_New_z = RotM[2][0] * v2_x + RotM[2][1] * v2_y + RotM[2][2] * v2_z;
		
		x_Dummy_1 = x[ActiveAtom] - R_DUMMY_1*v1_New_x;
		y_Dummy_1 = y[ActiveAtom] - R_DUMMY_1*v1_New_y;
		z_Dummy_1 = z[ActiveAtom] - R_DUMMY_1*v1_New_z;
		
		x_Dummy_2 = x[ActiveAtom] - R_DUMMY_2*v2_New_x;
		y_Dummy_2 = y[ActiveAtom] - R_DUMMY_2*v2_New_y;
		z_Dummy_2 = z[ActiveAtom] - R_DUMMY_2*v2_New_z;
		
		v_x_X1_Acceptor = v1_New_x;
		v_y_X1_Acceptor = v1_New_y;
		v_z_X1_Acceptor = v1_New_z;
	}
	else	{	// donor
		double v1_New_x, v1_New_y, v1_New_z;
		double v2_x, v2_y, v2_z, v2_New_x, v2_New_y, v2_New_z;
		double v3_x, v3_y, v3_z;
		double RotM[3][3], cos_theta, sin_theta;
		
		cos_theta = cos(theta*radianInv);
		sin_theta = sin(theta*radianInv);
		
		v_x_X1_H = x[ActiveHAtom] - x_Dummy_1;	// X vector
		v_y_X1_H = y[ActiveHAtom] - y_Dummy_1;
		v_z_X1_H = z[ActiveHAtom] - z_Dummy_1;
		Normalize_drude(v_x_X1_H, v_y_X1_H, v_z_X1_H);
		
		v2_x = x[ActiveHAtom] - x_Dummy_2;
		v2_y = y[ActiveHAtom] - y_Dummy_2;
		v2_z = z[ActiveHAtom] - z_Dummy_2;
		Normalize_drude(v2_x, v2_y, v2_z);
		
		Cross_Product(v_x_X1_H, v_y_X1_H, v_z_X1_H, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z);
		Normalize_drude(v3_x, v3_y, v3_z);
		
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
		
		x_Dummy_1 = x[ActiveHAtom] - R_DUMMY_1*v1_New_x;
		y_Dummy_1 = y[ActiveHAtom] - R_DUMMY_1*v1_New_y;
		z_Dummy_1 = z[ActiveHAtom] - R_DUMMY_1*v1_New_z;
		
		x_Dummy_2 = x[ActiveHAtom] - R_DUMMY_2*v2_New_x;
		y_Dummy_2 = y[ActiveHAtom] - R_DUMMY_2*v2_New_y;
		z_Dummy_2 = z[ActiveHAtom] - R_DUMMY_2*v2_New_z;
		
		v_x_X1_H = v1_New_x;
		v_y_X1_H = v1_New_y;
		v_z_X1_H = v1_New_z;
		
		x_Dummy_3 = x[ActiveHAtom] + R_DUMMY_3 * v_x_X1_H;
		y_Dummy_3 = y[ActiveHAtom] + R_DUMMY_3 * v_y_X1_H;
		z_Dummy_3 = z[ActiveHAtom] + R_DUMMY_3 * v_z_X1_H;
	}
}


int COptimizeWaterDimer::Is_Three_Atoms_Colinear(int ia, int ib, int ic)
{
	double v1_x, v1_y, v1_z;
	double v2_x, v2_y, v2_z;
	double dot_v1_v2;

	v1_x = x[ib] - x[ia];
	v1_y = y[ib] - y[ia];
	v1_z = z[ib] - z[ia];

	v2_x = x[ic] - x[ib];
	v2_y = y[ic] - y[ib];
	v2_z = z[ic] - z[ib];

	Normalize_drude(v1_x, v1_y, v1_z);
	Normalize_drude(v2_x, v2_y, v2_z);

	dot_v1_v2 = fabs(Dot_Product(v1_x, v1_y, v1_z, v2_x, v2_y, v2_z));
	if(dot_v1_v2 > cos(15.0*radianInv))	{	// smaller than 10 degree
		return 1;
	}
	else	{
		return 0;
	}
}

void COptimizeWaterDimer::Generate_Dummy_Atoms_Acceptor_1(void)
{
	int i, i_O, i_A, v_Aux_Sel, Is_Colinear;
	double v_Aux[2][3]={{1, 1, 0}, {0, 1, 1}}, dot_1, dot_2;
	double v_x_Acceptor_X2, v_y_Acceptor_X2, v_z_Acceptor_X2;

	i_O = ActiveAtom;
	i_A = BondList[i_O][0];

	x_Dummy_1 = x[i_A];	// dummy atom 1
	y_Dummy_1 = y[i_A];	// dummy atom 1
	z_Dummy_1 = z[i_A];	// dummy atom 1

	v_x_X1_Acceptor = x[i_O] - x[i_A];
	v_y_X1_Acceptor = y[i_O] - y[i_A];
	v_z_X1_Acceptor = z[i_O] - z[i_A];
	Normalize_drude(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor);
	
	if(nAtom == 2)	{
		dot_1 = fabs(Dot_Product(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor, v_Aux[0][0], v_Aux[0][1], v_Aux[0][2]));
		dot_2 = fabs(Dot_Product(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor, v_Aux[1][0], v_Aux[1][1], v_Aux[1][2]));
		
		if(dot_1 <= dot_2)	{
			v_Aux_Sel = 0;
		}
		else	{
			v_Aux_Sel = 1;
		}
		
		Cross_Product(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor, v_Aux[v_Aux_Sel][0], v_Aux[v_Aux_Sel][1], v_Aux[v_Aux_Sel][2], v_x_Acceptor_X2, v_y_Acceptor_X2, v_z_Acceptor_X2);
	}
	else	{
		if(Bond_Count[i_A] == 0)	{
			GetBondList(i_A);
		}
		for(i=0; i<Bond_Count[i_A]; i++)	{
			if(BondList[i_A][i] != i_O)	{	// 2nd neighbor atom
				break;
			}
		}
		if(i >= Bond_Count[i_A])	{
			Quit_With_Error_Msg("Error in Generate_Dummy_Atoms_Acceptor_1().\nQuit\n");
		}
		i = BondList[i_A][i];

		//start	to check whether three atoms (Idx, i_A, i) are colinear or not
		Is_Colinear = Is_Three_Atoms_Colinear(i_O, i_A, i);
		if(Is_Colinear)	{
			dot_1 = fabs(Dot_Product(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor, v_Aux[0][0], v_Aux[0][1], v_Aux[0][2]));
			dot_2 = fabs(Dot_Product(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor, v_Aux[1][0], v_Aux[1][1], v_Aux[1][2]));
			
			if(dot_1 <= dot_2)	{
				v_Aux_Sel = 0;
			}
			else	{
				v_Aux_Sel = 1;
			}
			
			Cross_Product(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor, v_Aux[v_Aux_Sel][0], v_Aux[v_Aux_Sel][1], v_Aux[v_Aux_Sel][2], v_x_Acceptor_X2, v_y_Acceptor_X2, v_z_Acceptor_X2);
		}
		else	{
			Cross_Product(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor, x[i_A]-x[i], y[i_A]-y[i], z[i_A]-z[i], v_x_Acceptor_X2, v_y_Acceptor_X2, v_z_Acceptor_X2);
		}
		//end	to check whether three atoms (Idx, i_A, i) are colinear or not
	}

	Normalize_drude(v_x_Acceptor_X2, v_y_Acceptor_X2, v_z_Acceptor_X2);
	
	x_Dummy_2 = x[i_O] + R_DUMMY_2 * v_x_Acceptor_X2;
	y_Dummy_2 = y[i_O] + R_DUMMY_2 * v_y_Acceptor_X2;
	z_Dummy_2 = z[i_O] + R_DUMMY_2 * v_z_Acceptor_X2;
}

void COptimizeWaterDimer::Generate_Dummy_Atoms_Acceptor_2(void)
{
	int i_O, i_A, i_B;
	double vx_Acceptor_A, vy_Acceptor_A, vz_Acceptor_A, vx_Acceptor_B, vy_Acceptor_B, vz_Acceptor_B;
	double v_x_Acceptor_X2, v_y_Acceptor_X2, v_z_Acceptor_X2;

	i_O = ActiveAtom;
	i_A = BondList[i_O][0];
	i_B = BondList[i_O][1];

	vx_Acceptor_A = x[i_A] - x[i_O];
	vy_Acceptor_A = y[i_A] - y[i_O];
	vz_Acceptor_A = z[i_A] - z[i_O];
	Normalize_drude(vx_Acceptor_A, vy_Acceptor_A, vz_Acceptor_A);

	vx_Acceptor_B = x[i_B] - x[i_O];
	vy_Acceptor_B = y[i_B] - y[i_O];
	vz_Acceptor_B = z[i_B] - z[i_O];
	Normalize_drude(vx_Acceptor_B, vy_Acceptor_B, vz_Acceptor_B);

	v_x_X1_Acceptor = -(vx_Acceptor_A + vx_Acceptor_B);
	v_y_X1_Acceptor = -(vy_Acceptor_A + vy_Acceptor_B);
	v_z_X1_Acceptor = -(vz_Acceptor_A + vz_Acceptor_B);
	Normalize_drude(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor);

	x_Dummy_1 = x[i_O] - R_DUMMY_1 * v_x_X1_Acceptor;
	y_Dummy_1 = y[i_O] - R_DUMMY_1 * v_y_X1_Acceptor;
	z_Dummy_1 = z[i_O] - R_DUMMY_1 * v_z_X1_Acceptor;

	Cross_Product(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor, vx_Acceptor_A, vy_Acceptor_A, vz_Acceptor_A, v_x_Acceptor_X2, v_y_Acceptor_X2, v_z_Acceptor_X2);
	Normalize_drude(v_x_Acceptor_X2, v_y_Acceptor_X2, v_z_Acceptor_X2);

	x_Dummy_2 = x[i_O] + R_DUMMY_2 * v_x_Acceptor_X2;
	y_Dummy_2 = y[i_O] + R_DUMMY_2 * v_y_Acceptor_X2;
	z_Dummy_2 = z[i_O] + R_DUMMY_2 * v_z_Acceptor_X2;
}

void COptimizeWaterDimer::Generate_Dummy_Atoms_Acceptor_3(void)
{
	int i_O, i_A, i_B, i_C;
	double vx_Acceptor_A, vy_Acceptor_A, vz_Acceptor_A, vx_Acceptor_B, vy_Acceptor_B, vz_Acceptor_B, vx_Acceptor_C, vy_Acceptor_C, vz_Acceptor_C;
	double v_x_Acceptor_X2, v_y_Acceptor_X2, v_z_Acceptor_X2, r_Proj;

	i_O = ActiveAtom;
	i_A = BondList[i_O][0];
	i_B = BondList[i_O][1];
	i_C = BondList[i_O][2];

	vx_Acceptor_A = x[i_A] - x[i_O];
	vy_Acceptor_A = y[i_A] - y[i_O];
	vz_Acceptor_A = z[i_A] - z[i_O];
	Normalize_drude(vx_Acceptor_A, vy_Acceptor_A, vz_Acceptor_A);

	vx_Acceptor_B = x[i_B] - x[i_O];
	vy_Acceptor_B = y[i_B] - y[i_O];
	vz_Acceptor_B = z[i_B] - z[i_O];
	Normalize_drude(vx_Acceptor_B, vy_Acceptor_B, vz_Acceptor_B);

	vx_Acceptor_C = x[i_C] - x[i_O];
	vy_Acceptor_C = y[i_C] - y[i_O];
	vz_Acceptor_C = z[i_C] - z[i_O];
	Normalize_drude(vx_Acceptor_C, vy_Acceptor_C, vz_Acceptor_C);

	v_x_X1_Acceptor = -(vx_Acceptor_A + vx_Acceptor_B + vx_Acceptor_C);
	v_y_X1_Acceptor = -(vy_Acceptor_A + vy_Acceptor_B + vy_Acceptor_C);
	v_z_X1_Acceptor = -(vz_Acceptor_A + vz_Acceptor_B + vz_Acceptor_C);
	r_Proj = Normalize_drude(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor);

	if( r_Proj < 0.04 )	{	// four atoms in the same plane
		x_Dummy_2 = x[i_O] + R_DUMMY_2 * vx_Acceptor_A;
		y_Dummy_2 = y[i_O] + R_DUMMY_2 * vy_Acceptor_A;
		z_Dummy_2 = z[i_O] + R_DUMMY_2 * vz_Acceptor_A;

		Cross_Product(vx_Acceptor_A, vy_Acceptor_A, vz_Acceptor_A, vx_Acceptor_B, vy_Acceptor_B, vz_Acceptor_B, v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor);
		Normalize_drude(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor);
		x_Dummy_1 = x[i_O] - R_DUMMY_1 * v_x_X1_Acceptor;
		y_Dummy_1 = y[i_O] - R_DUMMY_1 * v_y_X1_Acceptor;
		z_Dummy_1 = z[i_O] - R_DUMMY_1 * v_z_X1_Acceptor;
	}
	else	{
		x_Dummy_1 = x[i_O] - R_DUMMY_1 * v_x_X1_Acceptor;
		y_Dummy_1 = y[i_O] - R_DUMMY_1 * v_y_X1_Acceptor;
		z_Dummy_1 = z[i_O] - R_DUMMY_1 * v_z_X1_Acceptor;

		Cross_Product(v_x_X1_Acceptor, v_y_X1_Acceptor, v_z_X1_Acceptor, vx_Acceptor_A, vy_Acceptor_A, vz_Acceptor_A, v_x_Acceptor_X2, v_y_Acceptor_X2, v_z_Acceptor_X2);
		Normalize_drude(v_x_Acceptor_X2, v_y_Acceptor_X2, v_z_Acceptor_X2);

		x_Dummy_2 = x[i_O] + R_DUMMY_2 * v_x_Acceptor_X2;
		y_Dummy_2 = y[i_O] + R_DUMMY_2 * v_y_Acceptor_X2;
		z_Dummy_2 = z[i_O] + R_DUMMY_2 * v_z_Acceptor_X2;
	}

}

void COptimizeWaterDimer::Generate_Dummy_Atoms_Donor(void)
{
	int i, v_Aux_Sel, Is_Colinear;
	double v_Aux[2][3]={{1, 1, 0}, {0, 1, 1}}, dot_1, dot_2;
	double v_x_H_X2, v_y_H_X2, v_z_H_X2;

	x_Dummy_1 = x[ActiveAtom];	// dummy atom 1
	y_Dummy_1 = y[ActiveAtom];	// dummy atom 1
	z_Dummy_1 = z[ActiveAtom];	// dummy atom 1

	v_x_X1_H = x[ActiveHAtom] - x_Dummy_1;
	v_y_X1_H = y[ActiveHAtom] - y_Dummy_1;
	v_z_X1_H = z[ActiveHAtom] - z_Dummy_1;
	Normalize_drude(v_x_X1_H, v_y_X1_H, v_z_X1_H);
	
	if(nAtom == 2)	{
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
			Cross_Product(v_x_X1_H, v_y_X1_H, v_z_X1_H, x[ActiveAtom]-x[i], y[ActiveAtom]-y[i], z[ActiveAtom]-z[i], v_x_H_X2, v_y_H_X2, v_z_H_X2);
		}
		//end	to check whether three atoms (Idx, i_A, i) are colinear or not
	}

	Normalize_drude(v_x_H_X2, v_y_H_X2, v_z_H_X2);
	
	x_Dummy_2 = x[ActiveHAtom] + R_DUMMY_2 * v_x_H_X2;
	y_Dummy_2 = y[ActiveHAtom] + R_DUMMY_2 * v_y_H_X2;
	z_Dummy_2 = z[ActiveHAtom] + R_DUMMY_2 * v_z_H_X2;

	x_Dummy_3 = x[ActiveHAtom] + R_DUMMY_3 * v_x_X1_H;
	y_Dummy_3 = y[ActiveHAtom] + R_DUMMY_3 * v_y_X1_H;
	z_Dummy_3 = z[ActiveHAtom] + R_DUMMY_3 * v_z_X1_H;
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

	return n_Phi;
}

void Read_Fitted_Torsion_Parameters(void)
{
	int i, j, ReadItem, Idx_Phi;
	FILE *fIn;
	double Para_Read[10];

	fIn = fopen("saved-para.dat", "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Error in opening file: saved-para.dat in BackupCurrent_Torsion_Parameters()\nQuit\n");
	}

	for(i=0; i<n_Phi; i++)	{
		Idx_Phi = IdxDihSelect[i];
		for(j=0; j<10; j++)	{
			ReadItem = fscanf(fIn, "%lf", &(Para_Read[j]));
			if(ReadItem != 1)	{
				fclose(fIn);
				Quit_With_Error_Msg("Error in reading saved-para.dat in BackupCurrent_Torsion_Parameters().\nQuit\n");
			}
		}
//		for(j=1; j<=6; j++)	{	// to assign the torsion parameters
		for(j=1; j<=4; j++)	{	// to assign the torsion parameters
			Mol_ESP.Para_k_Dih[Idx_Phi][j] = Para_Read[j*2-2];
			Mol_ESP.Para_phi[Idx_Phi][j] = Para_Read[j*2-1];

			Mol_Water.Para_k_Dih[Idx_Phi][j] = Para_Read[j*2-2];
			Mol_Water.Para_phi[Idx_Phi][j] = Para_Read[j*2-1];			
		}

		j=6;
		Mol_ESP.Para_k_Dih[Idx_Phi][j] = Para_Read[j*2-4];
		Mol_ESP.Para_phi[Idx_Phi][j] = Para_Read[j*2-3];
		Mol_Water.Para_k_Dih[Idx_Phi][j] = Para_Read[j*2-4];
		Mol_Water.Para_phi[Idx_Phi][j] = Para_Read[j*2-3];
	}

	fclose(fIn);
}



double CalRMSD_Drude(double *xa, double *ya, double *za, double *xb, double *yb, double *zb, int AtomNum_Drude, CMol *pMol)
{
	double RMSD, midxa, midya, midza, midxb, midyb, midzb;
	double xa_tmp[MAX_ATOM], ya_tmp[MAX_ATOM], za_tmp[MAX_ATOM];
	double xb_tmp[MAX_ATOM], yb_tmp[MAX_ATOM], zb_tmp[MAX_ATOM];
	int i, AtomNum=0;

	for(i=0; i<AtomNum_Drude; i++)	{
		if(pMol->mass[i] > 1.0)	{	// real atom
			AtomNum++;
			xa_tmp[AtomNum] = xa[i];
			ya_tmp[AtomNum] = ya[i];
			za_tmp[AtomNum] = za[i];

			xb_tmp[AtomNum] = xb[i];
			yb_tmp[AtomNum] = yb[i];
			zb_tmp[AtomNum] = zb[i];
		}
	}

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

//	for(i=1; i<=AtomNum; i++)	{	// generate the rotated coodinates
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

