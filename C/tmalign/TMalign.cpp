/*
===============================================================================
   This is a re-implementation of TM-align algorithm in C/C++. The code was 
   is written by Jianyi Yang and later updated by Jianjie Wu at The Yang Zhang 
   lab, Department of Computational Medicine and Bioinformatics, University of 
   Michigan, 100 Washtenaw Avenue, Ann Arbor, MI 48109-2218. Please report bugs 
   and questions to jianjiew@umich.edu or zhng@umich.edu

   DISCLAIMER:
     Permission to use, copy, modify, and distribute this program for 
     any purpose, with or without fee, is hereby granted, provided that
     the notices on the head, the reference information, and this
     copyright notice appear in all copies or substantial portions of 
     the Software. It is provided "as is" without express or implied 
     warranty.
   *************** updating history ********************************
   2012/01/24: A C/C++ code of TM-align was constructed by Jianyi Yang
   2016/05/21: Several updates of this program were made by Jianjie Wu, including
              (1) fixed several compiling bugs
	      (2) made I/O of C/C++ version consistent with the Fortran version
	      (3) added outputs including full-atom and ligand structures
	      (4) added options of '-i', '-I' and '-m'
   2016/05/25: fixed a bug on PDB file reading
===============================================================================
*/
#include "basic_define.h"

using namespace std;

#include "basic_fun.h"
#include "NW.h"
#include "Kabsch.h"
#include "TMalign.h"

// This define is necessary to unmangle C++ function names. Bullshit politics involved http://blog.vrplumber.com/b/2007/07/21/ctypes-for-c-is/
// http://wolfprojects.altervista.org/articles/dll-in-c-for-python/
#ifndef _WINDOWS
#define DLLEXPORT extern "C"
#endif

#ifdef _WINDOWS
#define DLLEXPORT extern "C" __declspec(dllexport)
#endif

#define Xf(atom, frame, nframes, nf3) atom*nf3+frame
#define Yf(atom, frame, nframes, nf3) atom*nf3+nframes+frame
#define Zf(atom, frame, nframes, nf3) atom*nf3+2*nframes+frame


void TMAlign::parameter_set4search(int xlen, int ylen)
{
	//parameter initilization for searching: D0_MIN, Lnorm, d0, d0_search, score_d8
	D0_MIN=0.5; 
	dcu0=4.25;                       //update 3.85-->4.25
 
	Lnorm=getmin(xlen, ylen);        //normaliz TMscore by this in searching
    if(Lnorm<=19)                    //update 15-->19
    {
        d0=0.168;                   //update 0.5-->0.168
    }
    else
    {
        d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    }
	D0_MIN=d0+0.8;              //this should be moved to above
    d0=D0_MIN;                  //update: best for search    


	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;


    score_d8=1.5*pow(Lnorm*1.0, 0.3)+3.5; //remove pairs with dis>d8 during search & final
}

void TMAlign::parameter_set4final(double len)
{
	D0_MIN=0.5; 
 
	Lnorm=len;            //normaliz TMscore by this in searching
    if(Lnorm<=21)         
    {
        d0=0.5;          
    }
    else
    {
        d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    }
    if(d0<D0_MIN) d0=D0_MIN;   

	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;  

}


void TMAlign::parameter_set4scale(int len, double d_s)
{
 
	d0=d_s;          
	Lnorm=len;            //normaliz TMscore by this in searching

	d0_search=d0;	
	if(d0_search>8) d0_search=8;
	if(d0_search<4.5) d0_search=4.5;  

}

int main(int argc, char *argv[]){
    /*********************************************************************************/
	/*                                get argument                                   */ 
    /*********************************************************************************/
    TMAlign* tma = new TMAlign();

    char fname_lign[MAXLEN], fname_matrix[MAXLEN];
    char Lnorm_ave[MAXLEN];
	int nameIdx = 0;
	for(int i = 1; i < argc; i++)
	{
		if ( !strcmp(argv[i],"-o") && i < (argc-1) ) 
		{
			strcpy(tma->out_reg, argv[i + 1]);      tma->o_opt = true; i++;
		}
		else if ( !strcmp(argv[i],"-u") && i < (argc-1) ) 
		{
			tma->Lnorm_ass = atof(argv[i + 1]); tma->u_opt = true; i++;
		}
		else if ( !strcmp(argv[i],"-a") && i < (argc-1) ) 
		{
			strcpy(Lnorm_ave, argv[i + 1]);     tma->a_opt = true; i++;
		}
		else if ( !strcmp(argv[i],"-d") && i < (argc-1) ) 
		{
			tma->d0_scale = atof(argv[i + 1]); tma->d_opt = true; i++;
		}
		else if ( !strcmp(argv[i],"-v") ) 
		{
			tma->v_opt = true; 
		}
		else if ( !strcmp(argv[i],"-i") && i < (argc-1) ) 
		{
			strcpy(fname_lign, argv[i + 1]);      tma->i_opt = true; i++;
		}
		else if (!strcmp(argv[i], "-m") && i < (argc-1) ) 
		{
			strcpy(fname_matrix, argv[i + 1]);      tma->m_opt = true; i++;
		}// get filename for rotation matrix
		else if (!strcmp(argv[i], "-I") && i < (argc-1) ) 
		{
			strcpy(fname_lign, argv[i + 1]);      tma->I_opt = true; i++;
		}
		else
		{
			if (nameIdx == 0)
			{
				strcpy(tma->xname, argv[i]);
			}
			else if (nameIdx == 1)
			{
				strcpy(tma->yname, argv[i]);
			}
			nameIdx++;
		}
	}

	if( tma->a_opt )
	{
		if(!strcmp(Lnorm_ave, "T"))
		{
		}
		else if(!strcmp(Lnorm_ave, "F"))
		{
			tma->a_opt=false;
		}
		else
		{
			cout << "Wrong value for option -a!  It should be T or F" << endl;
			exit(EXIT_FAILURE);
		}
	}
	if( tma->u_opt )
	{
		if(tma->Lnorm_ass<=0)
		{
			cout << "Wrong value for option -u!  It should be >0" << endl;
			exit(EXIT_FAILURE);
		}
	}
	if( tma->d_opt )
	{
		if(tma->d0_scale<=0)
		{
			cout << "Wrong value for option -d!  It should be >0" << endl;
			exit(EXIT_FAILURE);
		}
	}
    string basename = string(argv[0]);
	////// read initial alignment file from 'alignment.txt' //////
	int idx = basename.find_last_of("\\");
	basename = basename.substr(0, idx + 1);
	if (tma->i_opt || tma->I_opt)// Ask TM-align to start with an alignment,specified in fasta file 'align.txt'
	{
		if (fname_lign == "")
		{
			cout << "Please provide a file name for option -i!" << endl;
			exit(EXIT_FAILURE);
		}
		// open alignment file
		int n_p = 0;// number of structures in alignment file
		string line;

		string fullpath = basename + fname_lign;
		char path1[1000];
		strcpy(path1, fullpath.c_str());
		ifstream fileIn(path1);
		if (fileIn.is_open())
		{
			bool bContinue = true;
			while (fileIn.good() && bContinue)
			{
				getline(fileIn, line);
				if (line.compare(0, 1, ">") == 0)// Flag for a new structure
				{
					strcpy(tma->sequence[n_p], "");
					n_p++;
					if (n_p > 2)
						bContinue = false;
				}
				else// Read data
				{
					if (n_p > 0 && line!="")
					{
						strcat(tma->sequence[n_p-1], line.c_str());
					}
				}
			}
			fileIn.close();
		}
		else
		{
			cout << "\nAlignment file does not exist.\n";
			exit(EXIT_FAILURE);
		}
	
		if (n_p < 2)
		{
			cout << "\nERROR: Fasta format is wrong, two proteins should be included.\n";
			exit(EXIT_FAILURE);
		}
		else
		{
			if (strlen(tma->sequence[0]) != strlen(tma->sequence[1]))
			{
				cout << "\nWarning: FASTA format may be wrong, the length in alignment should be equal respectively to the aligned proteins.\n";
				exit(EXIT_FAILURE);
			}
		}
	}

    /*********************************************************************************/
	/*                                load data                                      */ 
    /*********************************************************************************/
    tma->load_PDB_allocate_memory(tma->xname, tma->yname); // this allocates nres1, nres2 (atoms per res), xa, ya (coords), seqx, seqy, xresno, yresno

    //for (int i=0; i<tma->ylen; i++){
    //    printf("P: %d %d %c %lf\n", tma->yresno[i], tma->nres2[tma->yresno[i]][32], tma->seqy[i], tma->ya[i][0]);
    //}
    double TM1, TM2, rmsd;
    tma->main(tma->xname, tma->yname, true, &TM1, &TM2, &rmsd);
    tma->free_all_memory();
    return 0;
}

DLLEXPORT void tmalign(int xlen, int ylen, int* xresno, int* yresno, char* seqx, char* seqy, float* xcoor, float* ycoor, int nframes, double *TM1, double *TM2, double *rmsd)
{
    double **xa, **ya;
    int nf3 = nframes * 3;
    
    TMAlign* tma = new TMAlign();
    NewArray(&(tma->xa), xlen, 3);
    NewArray(&(tma->ya), ylen, 3);
    tma->seqx = seqx;
    tma->seqy = seqy;
    tma->xresno = xresno;
    tma->yresno = yresno;
    tma->xlen = xlen;
    tma->ylen = ylen;
    tma->tempxlen = xlen;
    tma->tempylen = ylen;
    tma->minlen = min(tma->xlen, tma->ylen);

    for (int i=0; i<ylen; i++){
        tma->nres2[yresno[i]][32] = 1;
    }

    for (int i=0; i<xlen; i++){
        tma->xa[i][0] = (double)xcoor[Xf(i, 0, 1, 3)];
        tma->xa[i][1] = (double)xcoor[Yf(i, 0, 1, 3)];
        tma->xa[i][2] = (double)xcoor[Zf(i, 0, 1, 3)];
        tma->nres1[xresno[i]][32] = 1;
    }

    for (int f=0; f<nframes; f++){
        for (int i=0; i<ylen; i++){
            tma->ya[i][0] = (double)ycoor[Xf(i, f, nframes, nf3)];
            tma->ya[i][1] = (double)ycoor[Yf(i, f, nframes, nf3)];
            tma->ya[i][2] = (double)ycoor[Zf(i, f, nframes, nf3)];
            //fprintf(stderr, "ylen: %d X: %lf %lf %lf\n", ylen, tma->ya[i][0], tma->ya[i][1], tma->ya[i][2]);
        }
        tma->main(tma->xname, tma->yname, false, &TM1[f], &TM2[f], &rmsd[f]);
    }

    //fprintf(stderr, "E: %d %d\n", xlen, ylen);
    //for (int i=0; i<ylen; i++){
    //    fprintf(stderr, "X: %d %d %c %lf %lf %lf\n", yresno[i], tma->nres2[yresno[i]][32], seqy[i], tma->ya[i][0], tma->ya[i][1], tma->ya[i][2]);
    //}
    tma->free_rest_memory();
}


int TMAlign::main(char* xname, char* yname, bool outputres, double *r_TM1, double *r_TM2, double *r_rmsd)
{
	//clock_t t1, t2;
	//t1 = clock();

    temp_var_allocation();

    /*********************************************************************************/
	/*                                parameter set                                  */ 
    /*********************************************************************************/
	parameter_set4search(xlen, ylen);          //please set parameters in the function
    int simplify_step     = 40;               //for similified search engine
    int score_sum_method  = 8;                //for scoring method, whether only sum over pairs with dis<score_d8
        
	int i;
    int *invmap0          = new int[ylen+1]; 
    int *invmap           = new int[ylen+1]; 
    double TM, TMmax=-1;
	for(i=0; i<ylen; i++)
	{
		invmap0[i]=-1;
	}	

	double ddcc=0.4;
	if(Lnorm <= 40) ddcc=0.1;   //Lnorm was setted in parameter_set4search
	double local_d0_search = d0_search;

	//*********************************************************************************//
	//                  get initial alignment from user's input:                       //
	//                  Stick to the initial alignment                                 //
	//*********************************************************************************//
	char dest[1000];
	bool bAlignStick = false;
	if (I_opt)// if input has set parameter for "-I"
	{
		for (int j = 1; j < ylen; j++)// Set aligned position to be "-1"
			invmap[j] = -1;

		int i1 = -1;// in C version, index starts from zero, not from one
		int i2 = -1;
		int L1 = strlen(sequence[0]);
		int L2 = strlen(sequence[1]);
		int L = min(L1, L2);// Get positions for aligned residues
		for (int kk1 = 0; kk1 < L; kk1++)
		{
			if (sequence[0][kk1] != '-')
				i1++;
			if (sequence[1][kk1] != '-')
			{
				i2++;
				if (i2 >= ylen || i1 >= xlen)
					kk1 = L;
				else
				{
					if (sequence[0][kk1] != '-')
					{
						invmap[i2] = i1;
					}
				}
			}
		}
		
		//--------------- 2. Align proteins from original alignment
		double prevD0_MIN = D0_MIN;// stored for later use
		int prevLnorm = Lnorm;
		double prevd0 = d0;
		TM_ali = standard_TMscore(xa, ya, xlen, ylen, invmap, L_ali, rmsd_ali);
		D0_MIN = prevD0_MIN;
		Lnorm = prevLnorm;
		d0 = prevd0;
		TM = detailed_search_standard(xa, ya, xlen, ylen, invmap, t, u, 40, 8, local_d0_search, true);
		if (TM > TMmax)
		{
			TMmax = TM;
			for (i = 0; i<ylen; i++)
				invmap0[i] = invmap[i];
		}
		bAlignStick = true;
	}

    /*********************************************************************************/
	/*         get initial alignment with gapless threading                          */ 
    /*********************************************************************************/
	if (!bAlignStick)
	{
		get_initial(xa, ya, xlen, ylen, invmap0);
		//find the max TMscore for this initial alignment with the simplified search_engin
		TM = detailed_search(xa, ya, xlen, ylen, invmap0, t, u, simplify_step, score_sum_method, local_d0_search);
		if (TM>TMmax)
		{
			TMmax = TM;
		}
		//run dynamic programing iteratively to find the best alignment
		TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30, local_d0_search);
		if (TM>TMmax)
		{
			TMmax = TM;
			for (int i = 0; i<ylen; i++)
			{
				invmap0[i] = invmap[i];
			}
		}

		/*********************************************************************************/
		/*         get initial alignment based on secondary structure                    */
		/*********************************************************************************/
		get_initial_ss(xa, ya, xlen, ylen, invmap);
		TM = detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method, local_d0_search);
		if (TM>TMmax)
		{
			TMmax = TM;
			for (int i = 0; i<ylen; i++)
			{
				invmap0[i] = invmap[i];
			}
		}
		if (TM > TMmax*0.2)
		{
			TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30, local_d0_search);
			if (TM>TMmax)
			{
				TMmax = TM;
				for (int i = 0; i<ylen; i++)
				{
					invmap0[i] = invmap[i];
				}
			}
		}

		/*********************************************************************************/
		/*         get initial alignment based on local superposition                    */
		/*********************************************************************************/
		//=initial5 in original TM-align
		if (get_initial5(xa, ya, xlen, ylen, invmap))
		{
			TM = detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method, local_d0_search);
			if (TM>TMmax)
			{
				TMmax = TM;
				for (int i = 0; i<ylen; i++)
				{
					invmap0[i] = invmap[i];
				}
			}
			if (TM > TMmax*ddcc)
			{
				TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 2, local_d0_search);
				if (TM>TMmax)
				{
					TMmax = TM;
					for (int i = 0; i<ylen; i++)
					{
						invmap0[i] = invmap[i];
					}
				}
			}
		}
		else
		{
			cout << endl << endl << "Warning: initial alignment from local superposition fail!" << endl << endl << endl;
		}

		/*********************************************************************************/
		/*    get initial alignment based on previous alignment+secondary structure      */
		/*********************************************************************************/
		//=initial3 in original TM-align
		get_initial_ssplus(xa, ya, xlen, ylen, invmap0, invmap);
		TM = detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method, local_d0_search);
		if (TM>TMmax)
		{
			TMmax = TM;
			for (i = 0; i<ylen; i++)
			{
				invmap0[i] = invmap[i];
			}
		}
		if (TM > TMmax*ddcc)
		{
			TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30, local_d0_search);
			if (TM>TMmax)
			{
				TMmax = TM;
				for (i = 0; i<ylen; i++)
				{
					invmap0[i] = invmap[i];
				}
			}
		}

		/*********************************************************************************/
		/*        get initial alignment based on fragment gapless threading              */
		/*********************************************************************************/
		//=initial4 in original TM-align
		get_initial_fgt(xa, ya, xlen, ylen, xresno, yresno, invmap);
		TM = detailed_search(xa, ya, xlen, ylen, invmap, t, u, simplify_step, score_sum_method, local_d0_search);
		if (TM>TMmax)
		{
			TMmax = TM;
			for (i = 0; i<ylen; i++)
			{
				invmap0[i] = invmap[i];
			}
		}
		if (TM > TMmax*ddcc)
		{
			TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 1, 2, 2, local_d0_search);
			if (TM>TMmax)
			{
				TMmax = TM;
				for (i = 0; i<ylen; i++)
				{
					invmap0[i] = invmap[i];
				}
			}
		}

		//*********************************************************************************//
		//                  get initial alignment from user's input:                       //
		//*********************************************************************************//
		if (i_opt)// if input has set parameter for "-i"
		{
			for (int j = 0; j < ylen; j++)// Set aligned position to be "-1"
				invmap[j] = -1;
		
			int i1 = -1;// in C version, index starts from zero, not from one
			int i2 = -1;
			int L1 = strlen(sequence[0]);
			int L2 = strlen(sequence[1]);
			int L = min(L1, L2);// Get positions for aligned residues
			for (int kk1 = 0; kk1 < L; kk1++)
			{
				if (sequence[0][kk1] != '-')
					i1++;
				if (sequence[1][kk1] != '-')
				{
					i2++;
					if (i2 >= ylen || i1 >= xlen)
						kk1 = L;
					else
					{
						if (sequence[0][kk1] != '-')
						{
							invmap[i2] = i1;
						}
					}
				}
			}

			//--------------- 2. Align proteins from original alignment
			double prevD0_MIN = D0_MIN;// stored for later use
			int prevLnorm = Lnorm;
			double prevd0 = d0;
			TM_ali = standard_TMscore(xa, ya, xlen, ylen, invmap, L_ali, rmsd_ali);
			D0_MIN = prevD0_MIN;
			Lnorm = prevLnorm;
			d0 = prevd0;

			TM = detailed_search_standard(xa, ya, xlen, ylen, invmap, t, u, 40, 8, local_d0_search, true);
			if (TM > TMmax)
			{
				TMmax = TM;
				for (i = 0; i<ylen; i++)
					invmap0[i] = invmap[i];
			}
			TM = DP_iter(xa, ya, xlen, ylen, t, u, invmap, 0, 2, 30, local_d0_search);// Different from get_initial, get_initial_ss and get_initial_ssplus
			if (TM>TMmax)
			{
				TMmax = TM;
				for (i = 0; i<ylen; i++)
				{
					invmap0[i] = invmap[i];
				}
			}
		}
	}

    //*********************************************************************************//
    //     The alignment will not be changed any more in the following                 //
    //*********************************************************************************//
	//check if the initial alignment is generated approately	
	bool flag=false;
	for(i=0; i<ylen; i++)
	{
		if(invmap0[i]>=0)
		{
			flag=true;
			break;			
		}			
	}		
	if(!flag) 
	{
		cout << "There is no alignment between the two proteins!" << endl;
		cout << "Program stop with no result!" << endl;
		return 1;
	}

    //*********************************************************************************//
    //       Detailed TMscore search engine  --> prepare for final TMscore             //
    //*********************************************************************************//       
    //run detailed TMscore search engine for the best alignment, and 
	//extract the best rotation matrix (t, u) for the best alginment
	simplify_step=1;
    score_sum_method=8;
	TM = detailed_search_standard(xa, ya, xlen, ylen, invmap0, t, u, simplify_step, score_sum_method, local_d0_search, false);
		
	//select pairs with dis<d8 for final TMscore computation and output alignment
	int n_ali8, k=0;
	int n_ali=0;
	int *m1, *m2;
	double d;
	m1=new int[xlen]; //alignd index in x
	m2=new int[ylen]; //alignd index in y
	do_rotation(xa, xt, xlen, t, u);
	k=0;
    for(int j=0; j<ylen; j++)
    {
        i=invmap0[j];
        if(i>=0)//aligned
        {
			n_ali++;        
            d=sqrt(dist(&xt[i][0], &ya[j][0]));
			if (d <= score_d8 || (I_opt == true))
			{
				m1[k]=i;
				m2[k]=j;

				xtm[k][0]=xa[i][0];
                xtm[k][1]=xa[i][1];
                xtm[k][2]=xa[i][2];
                    
                ytm[k][0]=ya[j][0];
                ytm[k][1]=ya[j][1];
                ytm[k][2]=ya[j][2];			

				r1[k][0] = xt[i][0];
				r1[k][1] = xt[i][1];
				r1[k][2] = xt[i][2];
				r2[k][0] = ya[j][0];
				r2[k][1] = ya[j][1];
				r2[k][2] = ya[j][2];
				
				k++;
			}
		}
	}
	n_ali8=k;

	double rmsd0 = 0.0;
	Kabsch(r1, r2, n_ali8, 0, &rmsd0, t, u);// rmsd0 is used for final output, only recalculate rmsd0, not t & u
	rmsd0 = sqrt(rmsd0 / n_ali8);

    //*********************************************************************************//
    //                               Final TMscore                                     //
    //                     Please set parameters for output                            //
    //*********************************************************************************//
	double rmsd;
	double t0[3], u0[3][3];
	double TM1, TM2;
	double d0_out=5.0;  
    simplify_step=1;
    score_sum_method=0;

	double d0_0, TM_0;
	double Lnorm_0=ylen;
	
	
	//normalized by length of structure A
	parameter_set4final(Lnorm_0);
	d0A=d0;
	d0_0=d0A;
	local_d0_search = d0_search;
	TM1 = TMscore8_search(xtm, ytm, n_ali8, t0, u0, simplify_step, score_sum_method, &rmsd, local_d0_search);
	TM_0 = TM1;

	//normalized by length of structure B
	parameter_set4final(xlen+0.0);
	d0B=d0;
	local_d0_search = d0_search;
	TM2 = TMscore8_search(xtm, ytm, n_ali8, t, u, simplify_step, score_sum_method, &rmsd, local_d0_search);

	if(a_opt)
	{
		//normalized by average length of structures A, B
		Lnorm_0=(xlen+ylen)*0.5;
		parameter_set4final(Lnorm_0);
		d0a=d0;
		d0_0=d0a;
		local_d0_search = d0_search;

		TM3 = TMscore8_search(xtm, ytm, n_ali8, t0, u0, simplify_step, score_sum_method, &rmsd, local_d0_search);
		TM_0=TM3;
	}
	if(u_opt)
	{	
		//normalized by user assigned length		
		parameter_set4final(Lnorm_ass);		
		d0u=d0;		
		d0_0=d0u;
		Lnorm_0=Lnorm_ass;
		local_d0_search = d0_search;
		TM4 = TMscore8_search(xtm, ytm, n_ali8, t0, u0, simplify_step, score_sum_method, &rmsd, local_d0_search);
		TM_0=TM4;
	}
	if(d_opt)
	{	
		//scaled by user assigned d0
		parameter_set4scale(ylen, d0_scale);
		d0_out=d0_scale;
		d0_0=d0_scale;
		//Lnorm_0=ylen;
		Lnorm_d0=Lnorm_0;
		local_d0_search = d0_search;
		TM5 = TMscore8_search(xtm, ytm, n_ali8, t0, u0, simplify_step, score_sum_method, &rmsd, local_d0_search);
		TM_0=TM5;
	}

    if (outputres){
    	output_results(xname, yname, xlen, ylen, t0, u0, TM1, TM2, rmsd0, d0_out, m1, m2, n_ali8, n_ali, TM_0, Lnorm_0, d0_0);
    }
    *r_TM1 = TM1;
    *r_TM2 = TM2;
    *r_rmsd = rmsd0;

    free_temp_memory();

    delete [] invmap0;
    delete [] invmap;
	delete [] m1;
	delete [] m2;

    //t2 = clock();    
    //float diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    //printf("\nTotal running time is %5.2f seconds\n", diff);        
 	return 0;	
}
