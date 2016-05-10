//#define MD_COULOMB  (332.0636)	//(in NAMD)
#define MD_COULOMB  (332.0716)	//(in CHARMM)
#define K_DRUDE		(500.0)


#define N_LEN_CHEM_NAME	(12)
#define N_PARA_BOND		(2)
#define N_PARA_ANGLE	(4)	//k_theta, theta_0, k_UB, S_0
#define N_PARA_DIHEDRAL	(3)
#define N_PARA_IMPRDIHEDRAL	(3)
#define N_PARA_LJ	(6)	// non-1-4 and 1-4 parameters
#define N_PARA_NBFIX	(2)	// E_Min and r_Min


typedef struct	{
	char Chem[2][N_LEN_CHEM_NAME];
	double para[N_PARA_BOND];
}BOND_REC;

typedef struct	{
	char Chem[3][N_LEN_CHEM_NAME];
	double para[N_PARA_ANGLE];
}ANGLE_REC;

typedef struct	{
	char Chem[4][N_LEN_CHEM_NAME];
	double para[N_PARA_DIHEDRAL];
}DIHEDRAL_REC;

typedef struct	{
	char Chem[4][N_LEN_CHEM_NAME];
	double para[N_PARA_IMPRDIHEDRAL];
}IMPRODIHEDRAL_REC;

typedef struct	{
	char Chem[N_LEN_CHEM_NAME];
	double para[N_PARA_LJ];
}LJ_REC;

typedef struct	{
	char Chem_1[N_LEN_CHEM_NAME];
	char Chem_2[N_LEN_CHEM_NAME];
	double para[N_PARA_NBFIX];
}NBFix_REC;

#define N_REC_MAX	(2500)

class CForceField
{
public:
	int n_Rec_NBFix, n_Rec_LJ;
	NBFix_REC NBFix_Rec[N_REC_MAX];
	LJ_REC LJ_Rec[N_REC_MAX];
	double E14FAC;

private:
	int n_Rec_Bond, n_Rec_Angle, n_Rec_Dihedral, n_Rec_ImproDihedral;

	BOND_REC Bond_Rec[N_REC_MAX];
	ANGLE_REC Angle_Rec[N_REC_MAX];
	DIHEDRAL_REC Dihedral_Rec[N_REC_MAX];
	IMPRODIHEDRAL_REC ImproDihedral_Rec[N_REC_MAX];


public:

	CForceField();
	~CForceField();
	void ReadForceField(char szNamePara[]);
	void ReadUpdatedCMap(void);
	void GetPara_Bond(char szChemName[][N_LEN_CHEM_NAME], double Para[]);
	void GetPara_Angle(char szChemName[][N_LEN_CHEM_NAME], double Para[]);
	void GetPara_Dihedral(char szChemName[][N_LEN_CHEM_NAME], double Para[]);
	void GetPara_ImproDIhedral(char szChemName[][N_LEN_CHEM_NAME], double Para[]);
	int  GetPara_LJ(char szChemName[][N_LEN_CHEM_NAME], double Para[]);	// return the index in LJ record list
	int  QueryNBFix(char szChem_1[], char szChem_2[]);
};


#define MAX_ATOM	(256)

#define LEN_ATOM_NAME	(8)
#define MAX_FACTOR	(12)
#define MAX_DIH_ITEM	(7)	// 1, ..., 6
#define MAX_ATOM_BOND	(7)
#define MAX_LP_POS_HOST	(4)
#define MAX_DIHEDRAL	(MAX_ATOM*10)


class CMol
{
private:

public:
	CForceField *myForceField;
	int nAtom, nBond, nAngle, nDihedral, nImpro, nLP, nDrude, nAniso;
	int StartLast, EndLast;
	char MolName[8];
	double x[MAX_ATOM], y[MAX_ATOM], z[MAX_ATOM];	//coordinate
	double x_save[MAX_ATOM], y_save[MAX_ATOM], z_save[MAX_ATOM];	//coordinate
	double grad_x[MAX_ATOM], grad_y[MAX_ATOM], grad_z[MAX_ATOM];	//gradient
	double CG[MAX_ATOM];	//charge
	double alpha[MAX_ATOM];	//polarizability
	double thole[MAX_ATOM];	//thole
	double mass[MAX_ATOM];	//mass

	double E14FAC;

	int nRealAtom;
	int RealAtom_List[MAX_ATOM];

	int SegID[MAX_ATOM];
	int IsDrude[MAX_ATOM], DrudeList[MAX_ATOM], IsDrudeHost[MAX_ATOM];	//flag for drude

	//start	to data for lone pairs
	int IsLonePair[MAX_ATOM], LPList[MAX_ATOM], LP_Host[MAX_ATOM];	//flag for lone pair
	int LP_PosHost[MAX_ATOM][MAX_LP_POS_HOST];
	double LP_Dist[MAX_ATOM], LP_Theta[MAX_ATOM], LP_Phi[MAX_ATOM];
	double LP_Sin_Theta[MAX_ATOM], LP_Cos_Theta[MAX_ATOM], LP_Sin_Phi[MAX_ATOM], LP_Cos_Phi[MAX_ATOM];
	//end	to data for lone pairs

	int AnisoList[MAX_ATOM][4];
	double Para_Aniso[MAX_ATOM][3];

	int IsFixed[MAX_ATOM];	//flag for movablity
	int BondList[MAX_FACTOR*MAX_ATOM];	//the flag for bond list
	int AngleList[MAX_FACTOR*MAX_ATOM];	//the flag for angle list
	int DihedralList[MAX_FACTOR*MAX_DIHEDRAL];
	int ImprDihedralList[MAX_FACTOR*MAX_ATOM];
	char ResName[MAX_ATOM][LEN_ATOM_NAME];
	char AtomName[MAX_ATOM][LEN_ATOM_NAME];
	char ChemName[MAX_ATOM][N_LEN_CHEM_NAME];
	int Type[MAX_ATOM];		//atom type in psf file


	//start	to restrain on all dihedrals, only impose the forces not on energy
	int Soft_Restrain_Dihedrals;
	int Is_Phi_Constrained[MAX_DIHEDRAL];
	double All_Phi0[MAX_DIHEDRAL];
	//end	to restrain on all dihedrals, only impose the forces not on energy


	//start	to data used for geometry optimization when there are atoms fixed
	int Geo_Opt_Drude_Only;
	int Is_Phi_Psi_Constrained;
	int E_CMap_On;
	int nActive_Bond, Active_List_Bond[MAX_FACTOR*MAX_ATOM];
	int nActive_Angle, Active_List_Angle[MAX_FACTOR*MAX_ATOM];
	int nActive_Dihedral, Active_List_Dihedral[MAX_FACTOR*MAX_DIHEDRAL];
	int nActive_ImpDih, Active_List_ImpDih[MAX_FACTOR*MAX_ATOM];
	int nActive_NonBond, Active_List_NonBond[MAX_ATOM*MAX_ATOM];
	int nActive_CMap, Active_List_CMap[MAX_ATOM];
	int nActive_Aniso, Active_List_Aniso[MAX_ATOM];
	//end	to data used for geometry optimization when there are atoms fixed

	int AtomBond[MAX_ATOM], Bond_Array[MAX_ATOM][MAX_ATOM_BOND];	//the array to store all bonded atoms 
	char DistMatrix[MAX_ATOM][MAX_ATOM];	//distance matrix for all pair of atoms
	int n_NB_Pair, NB_List_i[MAX_ATOM*MAX_ATOM], NB_List_j[MAX_ATOM*MAX_ATOM];	//non-bond list
	double NB_Epsilon[MAX_ATOM*MAX_ATOM], NB_R_min[MAX_ATOM*MAX_ATOM];
	
	int NBFix_Rec[MAX_ATOM][MAX_ATOM];
	int LJ_Para_Rec[MAX_ATOM];

	int n_Thole_Pair, Thole_List_i[MAX_ATOM*MAX_ATOM], Thole_List_j[MAX_ATOM*MAX_ATOM];	//the list of thole
	double Para_Thole_aa[MAX_ATOM*MAX_ATOM], Para_Thole_qq[MAX_ATOM*MAX_ATOM];	//the parameters used for thole energy corrections

	double Para_k_b[MAX_FACTOR*MAX_ATOM], Para_b0[MAX_FACTOR*MAX_ATOM];	//parameters for bond
	double Para_k_Urey[MAX_FACTOR*MAX_ATOM], Para_b0_Urey[MAX_FACTOR*MAX_ATOM];	//parameters for Urey-Bradley
	double Para_k_a[MAX_FACTOR*MAX_ATOM], Para_theta0[MAX_FACTOR*MAX_ATOM];	//parameters for angle
	double Para_k_Dih[MAX_FACTOR*MAX_DIHEDRAL][MAX_DIH_ITEM], Para_phi[MAX_FACTOR*MAX_DIHEDRAL][MAX_DIH_ITEM];	//parameters for dihedral
	double Para_k_ImpDih[MAX_FACTOR*MAX_ATOM], Para_Imp_phi[MAX_FACTOR*MAX_ATOM];	//parameters for improper dihedral
	double Para_Type_ImpDih[MAX_FACTOR*MAX_ATOM];	// type of improper potential used. 0 - CHARMM; 2 - AMBER/GAFF
	double Para_LJ_Epsilon[MAX_ATOM], Para_LJ_Sigma[MAX_ATOM];
	double Para_LJ_Epsilon_14[MAX_ATOM], Para_LJ_Sigma_14[MAX_ATOM];
	double Para_LJ_Epsilon_IJ[MAX_ATOM*MAX_ATOM], Para_LJ_Sigma_IJ[MAX_ATOM*MAX_ATOM];
	double Para_LJ_Sigma_IJ_pow_6[MAX_ATOM*MAX_ATOM], Para_LJ_Sigma_IJ_pow_12[MAX_ATOM*MAX_ATOM];
	double Para_Elec_Pair[MAX_ATOM*MAX_ATOM];
	int nCMapTerm, CMapList[MAX_ATOM][8];

	double E_Total, E_Bond, E_UreyB, E_Angle, E_Dihedral, E_ImproDihedral, E_VDW, E_Elec, E_CMap, E_Aniso, E_Thole, E_DrudeHyper;
	double E_Constrain_Phi, E_Constrain_Psi;

	double Phi0_Constrain, Psi0_Constrain, Phi_Real, Psi_Real;

	void ReadPSF(char szNamePsf[], int MultiSegment);
	void ReadPDB(char szNamePdb[]);
	void ReadCRD(char szNameCrd[]);
	void AssignForceFieldParameters(CForceField* pForceField);
	double Cal_E(int PrintE);
	void Cal_E_Bond(void);
	void Cal_E_UreyB(void);
	void Cal_E_Angle(void);
	void Cal_E_Dihedral(int ToSave_Phi);
	void Cal_E_ImproperDihedral(void);
	void Cal_E_VDW_ELEC(void);
	void Cal_E_CMap(void);
	void Cal_E_DrudeHyper(void);
	void Cal_E_Constrain_Phi(void);	//design for alad only
	void Cal_E_Constrain_Psi(void);	//design for alad only

	void Cal_E_Anisotropy(void);
	void Cal_E_Thole(void);

	void TestFirstDerivative(void);

        void FullGeometryOptimization_LBFGS_drude(int);
        void FullGeometryOptimization_LBFGS(void);
	double FullGeometryOptimization_LBFGS(int Only_Opt_Drude);
	double Cal_Dipole(double& dipole_x, double& dipole_y, double& dipole_z, double& dipole);
	void Gen_Phi_Psi_Map(void);


	void BuildDistanceMatrix(void);
	void BuildActiveEnergyTermList(void);
	void Setup_NonBondParameters(void);
	void Setup_TholePairParameters(void);
	void FixAllRealAtoms(void);
	void TranslateLastSegment(void);
	void TranslateLastSegmentBack(void);
	void BackupCoordinates(void);
	void RestoreCoordinates(void);
	void WriteCRDFile(char szName[]);
	void Setup_NBFix(void);
	void SavePdb(char szName[]);

	void Position_LonePair(void);	//translated from charmm code
	void LonePairForceReDistribute(void);	//re-distribute force for lone pairs to hosts
	void CartCV_Charmm(int i, int j, int k, int l, int Idx_LP, double dist);	//translated from charmm

	int Query_Dihedral_Index(int ia, int ib, int ic, int id);
	void Restrain_All_Torsions(int Flag);
};

