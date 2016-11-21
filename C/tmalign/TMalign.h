/*
===============================================================================
   Implementation of TM-align in C/C++   

   This program is written by Jianyi Yang at
   Yang Zhang lab
   And it is updated by Jianjie Wu at
   Yang Zhang lab
   Department of Computational Medicine and Bioinformatics 
   University of Michigan 
   100 Washtenaw Avenue, Ann Arbor, MI 48109-2218 
           
   Please report bugs and questions to jianjiew@umich.edu or zhng@umich.edu
===============================================================================
*/
#include "basic_define.h"

class TMAlign {
    public:
        double D0_MIN;                             //for d0
        double Lnorm;                              //normalization length
        double score_d8, d0, d0_search, dcu0;      //for TMscore search
        double **score;            			       //Input score table for dynamic programming
        bool   **path;                             //for dynamic programming  
        double **val;                              //for dynamic programming  
        int    xlen, ylen, minlen;                 //length of proteins
        int tempxlen, tempylen, tempMinlen;
        double **xa, **ya;                         //for input vectors xa[0...xlen-1][0..2], ya[0...ylen-1][0..2]
                                                   //in general, ya is regarded as native structure --> superpose xa onto ya
        int    *xresno, *yresno;                   //residue numbers, used in fragment gapless threading 
        double **xtm, **ytm;                       //for TMscore search engine
        double **xt;                               //for saving the superposed version of r_1 or xtm
        char   *seqx, *seqy;                       //for the protein sequence 
        int    *secx, *secy;                       //for the secondary structure 
        double **r1, **r2;                         //for Kabsch rotation 
        double t[3], u[3][3];                      //Kabsch translation vector and rotation matrix

        int atomxlen, atomylen;        // length of atoms
        int *ia1, *ia2;                // ATOM indices, used for display
        char  **aa1, **aa2;           // "N", or "CA", or "C" 
        double **xyza1, **xyza2;      // for input vectors xa[0...xlen-1][0..2], ya[0...ylen-1][0..2], just for display
        char **ra1, **ra2;           // for the protein sequence 
        int  *ir1, *ir2;             // residue numbers, used in fragment gapless threading 

        char sequence[10][MAXLEN];// get value from alignment file
        double TM_ali, rmsd_ali;// TMscore and rmsd from standard_TMscore func, 
        int L_ali;// Aligned length from standard_TMscore func, 

        char *ins1, *ins2, *ains1, *ains2;// flag characters for data read from PDB file, which begins with "ATOM", and locates at s(27), usually are spaces
        int **nres1, **nres2;// number of atoms, nres1(i,j): the number of atoms for ith residue, j usually is 32 for a space


        //argument variables
        char out_reg[MAXLEN];
        double Lnorm_ass, Lnorm_d0, d0_scale, d0A, d0B, d0u, d0a;
        bool o_opt, a_opt, u_opt, d_opt, v_opt;
        bool i_opt;// flags for -i, with user given initial alignment file
        bool m_opt;// flags for -m, output rotation matrix
        bool I_opt;// flags for -I, stick to user given initial alignment file

        double TM3, TM4, TM5;

        char xname[MAXLEN], yname[MAXLEN];

        void parameter_set4search(int xlen, int ylen);
        void parameter_set4final(double len);
        void parameter_set4scale(int len, double d_s);
        int main(char* xname, char* yname, bool outputres, double *r_TM1, double *r_TM2, double *r_rmsd0);
        void print_help(char *arg);

        TMAlign()
        {
            // Set all input options to false
	        o_opt = a_opt = u_opt = d_opt = v_opt = i_opt = m_opt = I_opt = false;
            nres_allocation();
        }

        void nres_allocation(){
            //------load data------>  
	        NewArray(&nres1, NMAX, ASCIILimit);// Only data from nres1[0~ NMAX-1][32~122] is used
	        NewArray(&nres2, NMAX, ASCIILimit);
	        for (int i = 0; i < NMAX; i++)// Initialization
	        {
		        for (int j = 0; j < ASCIILimit; j++)
		        {
			        nres1[i][j] = 0;
			        nres2[i][j] = 0;
		        }
	        }	
        }

        void temp_var_allocation()
        {
            //------allocate memory for other temporary variables. Should be called after setting xlen, ylen, minlen------>
         	NewArray(&r1, minlen, 3);
	        NewArray(&r2, minlen, 3);
	        NewArray(&xtm, minlen, 3);
	        NewArray(&ytm, minlen, 3);
	        NewArray(&xt, xlen, 3);

	        NewArray(&score, xlen+1, ylen+1);
	        NewArray(&path, xlen+1, ylen+1);
	        NewArray(&val, xlen+1, ylen+1); 

            secx = new int[xlen];
            secy = new int[ylen];
        }

        void free_temp_memory()
        {
	        DeleteArray(&path, xlen+1);
	        DeleteArray(&val, xlen+1);
	        DeleteArray(&score, xlen+1);
	        DeleteArray(&xt, xlen);
	        DeleteArray(&r1, minlen);
	        DeleteArray(&r2, minlen);
	        DeleteArray(&xtm, minlen);
	        DeleteArray(&ytm, minlen);
            delete [] secx;
            delete [] secy;
        }

        void load_PDB_allocate_memory(char *xname, char *yname)
        {    
	        tempxlen = get_PDB_len(xname);// Get predicted length
	        tempylen = get_PDB_len(yname);

	        //------allocate memory for x and y------>
	        NewArray(&xa, tempxlen, 3);
	        seqx = new char[tempxlen + 1];
	        secx = new int[tempxlen];
	        xresno = new int[tempxlen];

	        NewArray(&ya, tempylen, 3);
	        seqy = new char[tempylen + 1];
	        yresno = new int[tempylen];
	        secy = new int[tempylen];

	        xlen = read_PDB(xname, xa, seqx, xresno, nres1);// Get exact length
	        ylen = read_PDB(yname, ya, seqy, yresno, nres2);
	        minlen = min(xlen, ylen);
        }

        void free_rest_memory()
        {
	        DeleteArray(&xa, tempxlen);
	        DeleteArray(&ya, tempylen);
	        DeleteArray(&nres1, NMAX);
	        DeleteArray(&nres2, NMAX);
        }

        // This function is only for commandline calls. We should not free python variables memory or it segfaults.
        void free_all_memory()
        {
            free_rest_memory();
	        delete [] seqx;
            delete [] seqy;
            delete [] xresno;
            delete [] yresno;

	        if (o_opt)
	        {
		        DeleteArray(&xyza1, tempxlen);
		        DeleteArray(&aa1, tempxlen);
		        DeleteArray(&ra1, tempxlen);
		        delete[] ia1;
		        delete[] ir1;

		        DeleteArray(&xyza2, tempylen);
		        DeleteArray(&aa2, tempylen);
		        DeleteArray(&ra2, tempylen);
		        delete[] ia2;
		        delete[] ir2;

		        delete[]ins1;
		        delete[]ins2;
		        delete[]ains1;
		        delete[]ains2;
	        }    
        }


        //     1, collect those residues with dis<d;
        //     2, calculate TMscore
        int score_fun8( double **xa, 
                        double **ya, 
                        int n_ali,
                        double d,
                        int i_ali[], 
                        double *score1,
                        int score_sum_method
                      )
        {
            double score_sum=0, di;
            double d_tmp=d*d;
	        double d02=d0*d0;
	        double score_d8_cut = score_d8*score_d8;
            
            int i, n_cut, inc=0;

            while(1)
            {
		        n_cut=0;
                score_sum=0;
                for(i=0; i<n_ali; i++)
                {
                    di = dist(xa[i], ya[i]);
                    if(di<d_tmp)
                    {
                        i_ali[n_cut]=i;
                        n_cut++;
                    }
                    if(score_sum_method==8)
                    {				
                        if(di<=score_d8_cut)
                        {					
                            score_sum += 1/(1+di/d02);
                        }                
                    }
                    else
                    {
				        score_sum += 1/(1+di/d02);
                    }
                }
                //there are not enough feasible pairs, reliefe the threshold 		
                if(n_cut<3 && n_ali>3)
                {
			        inc++;
			        double dinc=(d+inc*0.5);
			        d_tmp = dinc * dinc;
                }
                else
                {
                    break;
                }

            }  

            *score1=score_sum/Lnorm;
	
            return n_cut;
        }

        int score_fun8_standard(double **xa,
	        double **ya,
	        int n_ali,
	        double d,
	        int i_ali[],
	        double *score1,
	        int score_sum_method
	        )
        {
	        double score_sum = 0, di;
	        double d_tmp = d*d;
	        double d02 = d0*d0;
	        double score_d8_cut = score_d8*score_d8;

	        int i, n_cut, inc = 0;
	        while (1)
	        {
		        n_cut = 0;
		        score_sum = 0;
		        for (i = 0; i<n_ali; i++)
		        {
			        di = dist(xa[i], ya[i]);
			        if (di<d_tmp)
			        {
				        i_ali[n_cut] = i;
				        n_cut++;
			        }
			        if (score_sum_method == 8)
			        {
				        if (di <= score_d8_cut)
				        {
					        score_sum += 1 / (1 + di / d02);
				        }
			        }
			        else
			        {
				        score_sum += 1 / (1 + di / d02);
			        }
		        }
		        //there are not enough feasible pairs, reliefe the threshold 		
		        if (n_cut<3 && n_ali>3)
		        {
			        inc++;
			        double dinc = (d + inc*0.5);
			        d_tmp = dinc * dinc;
		        }
		        else
		        {
			        break;
		        }

	        }

	        *score1 = score_sum / n_ali;
	        return n_cut;
        }

        double TMscore8_search( double **xtm, 
                                double **ytm,
                                int Lali, 
                                double t0[3],
                                double u0[3][3],
                                int simplify_step,
                                int score_sum_method,
                                double *Rcomm,
						        double local_d0_search
                               )
        { 
            int i, m;
            double score_max, score, rmsd;    
            const int kmax=Lali;    
            int k_ali[kmax], ka, k;
            double t[3];
            double u[3][3];
	        double d;
	

	        //iterative parameters
	        int n_it=20;            //maximum number of iterations
            const int n_init_max=6; //maximum number of different fragment length 
            int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4   
            int L_ini_min=4;
            if(Lali<4) L_ini_min=Lali;   
            int n_init=0, i_init;      
            for(i=0; i<n_init_max-1; i++)
            {
                n_init++;
                L_ini[i]=(int) (Lali/pow(2.0, (double) i));
                if(L_ini[i]<=L_ini_min)
                {
                    L_ini[i]=L_ini_min;
                    break;
                }
            }
            if(i==n_init_max-1)
            {
                n_init++;
                L_ini[i]=L_ini_min;
            }
            
            score_max=-1;
            //find the maximum score starting from local structures superposition
            int i_ali[kmax], n_cut;
            int L_frag; //fragment length
            int iL_max; //maximum starting postion for the fragment
	
            for(i_init=0; i_init<n_init; i_init++)
            {
                L_frag=L_ini[i_init];
                iL_max=Lali-L_frag;
              
                i=0;   
                while(1)
                {
                    //extract the fragment starting from position i 
                    ka=0;
                    for(k=0; k<L_frag; k++)
                    {
				        int kk=k+i;
                        r1[k][0]=xtm[kk][0];  
                        r1[k][1]=xtm[kk][1]; 
                        r1[k][2]=xtm[kk][2];   
                        
                        r2[k][0]=ytm[kk][0];  
                        r2[k][1]=ytm[kk][1]; 
                        r2[k][2]=ytm[kk][2];
                        
                        k_ali[ka]=kk;
                        ka++;
                    }
                    
                    //extract rotation matrix based on the fragment
                    Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);
			        if (simplify_step != 1)
				        *Rcomm = 0;
			        do_rotation(xtm, xt, Lali, t, u);
                    
                    //get subsegment of this fragment
			        d = local_d0_search - 1;
                    n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score, score_sum_method);
                    if(score>score_max)
                    {
                        score_max=score;
                        
                        //save the rotation matrix
                        for(k=0; k<3; k++)
                        {
                            t0[k]=t[k];
                            u0[k][0]=u[k][0];
                            u0[k][1]=u[k][1];
                            u0[k][2]=u[k][2];
                        }
                    }
                    
                    //try to extend the alignment iteratively            
                    d = local_d0_search + 1;
                    for(int it=0; it<n_it; it++)            
                    {
                        ka=0;
                        for(k=0; k<n_cut; k++)
                        {
                            m=i_ali[k];
                            r1[k][0]=xtm[m][0];  
                            r1[k][1]=xtm[m][1]; 
                            r1[k][2]=xtm[m][2];
                            
                            r2[k][0]=ytm[m][0];  
                            r2[k][1]=ytm[m][1]; 
                            r2[k][2]=ytm[m][2];
                            
                            k_ali[ka]=m;
                            ka++;
                        } 
                        //extract rotation matrix based on the fragment                
                        Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);
                        do_rotation(xtm, xt, Lali, t, u);
                        n_cut=score_fun8(xt, ytm, Lali, d, i_ali, &score, score_sum_method);
                        if(score>score_max)
                        {
                            score_max=score;

                            //save the rotation matrix
                            for(k=0; k<3; k++)
                            {
                                t0[k]=t[k];
                                u0[k][0]=u[k][0];
                                u0[k][1]=u[k][1];
                                u0[k][2]=u[k][2];
                            }                     
                        }
                        
                        //check if it converges            
                        if(n_cut==ka)
                        {				
                            for(k=0; k<n_cut; k++)
                            {
                                if(i_ali[k]!=k_ali[k])
						        {
							        break;
						        }
                            }
                            if(k==n_cut)
                            {						
                                break; //stop iteration
                            }
                        }                                                               
                    } //for iteration            

			        if(i<iL_max)
			        {
				        i=i+simplify_step; //shift the fragment		
				        if(i>iL_max) i=iL_max;  //do this to use the last missed fragment
			        }
			        else if(i>=iL_max)
			        {
				        break;
			        }
                }//while(1)
                //end of one fragment
            }//for(i_init
            return score_max;
        }


        double TMscore8_search_standard(double **xtm,
	        double **ytm,
	        int Lali,
	        double t0[3],
	        double u0[3][3],
	        int simplify_step,
	        int score_sum_method,
	        double *Rcomm,
	        double local_d0_search
	        )
        {
	        int i, m;
	        double score_max, score, rmsd;
	        const int kmax = Lali;
	        int k_ali[kmax], ka, k;
	        double t[3];
	        double u[3][3];
	        double d;


	        //iterative parameters
	        int n_it = 20;            //maximum number of iterations
	        const int n_init_max = 6; //maximum number of different fragment length 
	        int L_ini[n_init_max];  //fragment lengths, Lali, Lali/2, Lali/4 ... 4   
	        int L_ini_min = 4;
	        if (Lali<4) L_ini_min = Lali;
	        int n_init = 0, i_init;
	        for (i = 0; i<n_init_max - 1; i++)
	        {
		        n_init++;
		        L_ini[i] = (int)(Lali / pow(2.0, (double)i));
		        if (L_ini[i] <= L_ini_min)
		        {
			        L_ini[i] = L_ini_min;
			        break;
		        }
	        }
	        if (i == n_init_max - 1)
	        {
		        n_init++;
		        L_ini[i] = L_ini_min;
	        }

	        score_max = -1;
	        //find the maximum score starting from local structures superposition
	        int i_ali[kmax], n_cut;
	        int L_frag; //fragment length
	        int iL_max; //maximum starting postion for the fragment

	        for (i_init = 0; i_init<n_init; i_init++)
	        {
		        L_frag = L_ini[i_init];
		        iL_max = Lali - L_frag;

		        i = 0;
		        while (1)
		        {
			        //extract the fragment starting from position i 
			        ka = 0;
			        for (k = 0; k<L_frag; k++)
			        {
				        int kk = k + i;
				        r1[k][0] = xtm[kk][0];
				        r1[k][1] = xtm[kk][1];
				        r1[k][2] = xtm[kk][2];

				        r2[k][0] = ytm[kk][0];
				        r2[k][1] = ytm[kk][1];
				        r2[k][2] = ytm[kk][2];

				        k_ali[ka] = kk;
				        ka++;
			        }
			        //extract rotation matrix based on the fragment
			        Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);
			        if (simplify_step != 1)
				        *Rcomm = 0;
			        do_rotation(xtm, xt, Lali, t, u);

			        //get subsegment of this fragment
			        d = local_d0_search - 1;
			        n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, &score, score_sum_method);

			        if (score>score_max)
			        {
				        score_max = score;

				        //save the rotation matrix
				        for (k = 0; k<3; k++)
				        {
					        t0[k] = t[k];
					        u0[k][0] = u[k][0];
					        u0[k][1] = u[k][1];
					        u0[k][2] = u[k][2];
				        }
			        }

			        //try to extend the alignment iteratively            
			        d = local_d0_search + 1;
			        for (int it = 0; it<n_it; it++)
			        {
				        ka = 0;
				        for (k = 0; k<n_cut; k++)
				        {
					        m = i_ali[k];
					        r1[k][0] = xtm[m][0];
					        r1[k][1] = xtm[m][1];
					        r1[k][2] = xtm[m][2];

					        r2[k][0] = ytm[m][0];
					        r2[k][1] = ytm[m][1];
					        r2[k][2] = ytm[m][2];

					        k_ali[ka] = m;
					        ka++;
				        }
				        //extract rotation matrix based on the fragment                
				        Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);
				        do_rotation(xtm, xt, Lali, t, u);
				        n_cut = score_fun8_standard(xt, ytm, Lali, d, i_ali, &score, score_sum_method);
				        if (score>score_max)
				        {
					        score_max = score;

					        //save the rotation matrix
					        for (k = 0; k<3; k++)
					        {
						        t0[k] = t[k];
						        u0[k][0] = u[k][0];
						        u0[k][1] = u[k][1];
						        u0[k][2] = u[k][2];
					        }
				        }

				        //check if it converges            
				        if (n_cut == ka)
				        {
					        for (k = 0; k<n_cut; k++)
					        {
						        if (i_ali[k] != k_ali[k])
						        {
							        break;
						        }
					        }
					        if (k == n_cut)
					        {
						        break; //stop iteration
					        }
				        }
			        } //for iteration            

			        if (i<iL_max)
			        {
				        i = i + simplify_step; //shift the fragment		
				        if (i>iL_max) i = iL_max;  //do this to use the last missed fragment
			        }
			        else if (i >= iL_max)
			        {
				        break;
			        }
		        }//while(1)
		        //end of one fragment
	        }//for(i_init
	        return score_max;
        }

        //Comprehensive TMscore search engine
        // input:   two vector sets: x, y
        //          an alignment invmap0[] between x and y
        //          simplify_step: 1 or 40 or other integers
        //          score_sum_method: 0 for score over all pairs
        //                            8 for socre over the pairs with dist<score_d8          
        // output:  the best rotaion matrix t, u that results in highest TMscore
        double detailed_search( double **x,
                                double **y, 
                                int x_len, 
                                int y_len, 
                                int invmap0[],
                                double t[3],
                                double u[3][3],
                                int simplify_step,
                                int score_sum_method,
						        double local_d0_search
                               )
        {
            //x is model, y is template, try to superpose onto y
            int i, j, k;     
            double tmscore;
            double rmsd;



            k=0;
            for(i=0; i<y_len; i++) 
            {
                j=invmap0[i];
                if(j>=0) //aligned
                {
                    xtm[k][0]=x[j][0];
                    xtm[k][1]=x[j][1];
                    xtm[k][2]=x[j][2];
                        
                    ytm[k][0]=y[i][0];
                    ytm[k][1]=y[i][1];
                    ytm[k][2]=y[i][2];
                    k++;
                }
            }

            //detailed search 40-->1
	        tmscore = TMscore8_search(xtm, ytm, k, t, u, simplify_step, score_sum_method, &rmsd, local_d0_search);
            return tmscore;
        }

        double detailed_search_standard(double **x,
                                double **y, 
                                int x_len, 
                                int y_len, 
                                int invmap0[],
                                double t[3],
                                double u[3][3],
                                int simplify_step,
                                int score_sum_method,
						        double local_d0_search,
						        const bool& bNormalize
                               )
        {
            //x is model, y is template, try to superpose onto y
            int i, j, k;     
            double tmscore;
            double rmsd;



            k=0;
            for(i=0; i<y_len; i++) 
            {
                j=invmap0[i];
                if(j>=0) //aligned
                {
                    xtm[k][0]=x[j][0];
                    xtm[k][1]=x[j][1];
                    xtm[k][2]=x[j][2];
                        
                    ytm[k][0]=y[i][0];
                    ytm[k][1]=y[i][1];
                    ytm[k][2]=y[i][2];
                    k++;
                }
            }

            //detailed search 40-->1
	        tmscore = TMscore8_search_standard(xtm, ytm, k, t, u, simplify_step, score_sum_method, &rmsd, local_d0_search);
	        if (bNormalize)// "-i", to use standard_TMscore, then bNormalize=true, else bNormalize=false; 
		        tmscore = tmscore * k / Lnorm;

            return tmscore;
        }

        //compute the score quickly in three iterations
        double get_score_fast(double **x, double **y, int x_len, int y_len, int invmap[])
        {
            double rms, tmscore, tmscore1, tmscore2;
            int i, j, k;

            k=0;
            for(j=0; j<y_len; j++)
            {
                i=invmap[j];
                if(i>=0)
                {
                    r1[k][0]=x[i][0];
                    r1[k][1]=x[i][1];
                    r1[k][2]=x[i][2];

                    r2[k][0]=y[j][0];
                    r2[k][1]=y[j][1];
                    r2[k][2]=y[j][2];
                    
                    xtm[k][0]=x[i][0];
                    xtm[k][1]=x[i][1];
                    xtm[k][2]=x[i][2];
                    
                    ytm[k][0]=y[j][0];
                    ytm[k][1]=y[j][1];
                    ytm[k][2]=y[j][2];                  
                    
                    k++;
                }
                else if(i!=-1)
                {
                    PrintErrorAndQuit("Wrong map!\n");
                }       
            }
            Kabsch(r1, r2, k, 1, &rms, t, u);
            
            //evaluate score   
            double di;
	        const int len=k;
            double dis[len];    
	        double d00=d0_search;
	        double d002=d00*d00;
	        double d02=d0*d0;
	
            int n_ali=k;
	        double xrot[3];
	        tmscore=0;
	        for(k=0; k<n_ali; k++)
	        {
                transform(t, u, &xtm[k][0], xrot);        
                di=dist(xrot, &ytm[k][0]);
                dis[k]=di;
                tmscore += 1/(1+di/d02);
            }
	
           
           
           //second iteration 
            double d002t=d002;
            while(1)
            {
		        j=0;
                for(k=0; k<n_ali; k++)
                {            
                    if(dis[k]<=d002t)
                    {
                        r1[j][0]=xtm[k][0];
                        r1[j][1]=xtm[k][1];
                        r1[j][2]=xtm[k][2];
                        
                        r2[j][0]=ytm[k][0];
                        r2[j][1]=ytm[k][1];
                        r2[j][2]=ytm[k][2];
                        
                        j++;
                    }
                }
                //there are not enough feasible pairs, relieve the threshold 
                if(j<3 && n_ali>3)
                {
                    d002t += 0.5;
                }
                else
                {
                    break;
                }
            }
            
            if(n_ali!=j)
            {
                Kabsch(r1, r2, j, 1, &rms, t, u);
            	tmscore1=0;
            	for(k=0; k<n_ali; k++)
            	{
                    transform(t, u, &xtm[k][0], xrot);        
                    di=dist(xrot, &ytm[k][0]);
                    dis[k]=di;
                    tmscore1 += 1/(1+di/d02);
                }
                
                //third iteration
                d002t=d002+1;
               
                while(1)
                {
			        j=0;
                    for(k=0; k<n_ali; k++)
                    {            
                        if(dis[k]<=d002t)
                        {
                            r1[j][0]=xtm[k][0];
                            r1[j][1]=xtm[k][1];
                            r1[j][2]=xtm[k][2];
                            
                            r2[j][0]=ytm[k][0];
                            r2[j][1]=ytm[k][1];
                            r2[j][2]=ytm[k][2];
                                                
                            j++;
                        }
                    }
                    //there are not enough feasible pairs, relieve the threshold 
                    if(j<3 && n_ali>3)
                    {
                        d002t += 0.5;
                    }
                    else
                    {
                        break;
                    }
                }

                //evaluate the score
                Kabsch(r1, r2, j, 1, &rms, t, u);
                tmscore2=0;
                for(k=0; k<n_ali; k++)
                {
                    transform(t, u, &xtm[k][0], xrot);
                    di=dist(xrot, &ytm[k][0]);
                    tmscore2 += 1/(1+di/d02);
                }    
            }
            else
            {
                tmscore1=tmscore;
                tmscore2=tmscore;
            }
            
              
            if(tmscore1>=tmscore) tmscore=tmscore1;
            if(tmscore2>=tmscore) tmscore=tmscore2;

            
            return tmscore; // no need to normalize this score because it will not be used for latter scoring
        }


        //perform gapless threading to find the best initial alignment
        //input: x, y, x_len, y_len
        //output: y2x0 stores the best alignment: e.g., 
        //y2x0[j]=i means:
        //the jth element in y is aligned to the ith element in x if i>=0 
        //the jth element in y is aligned to a gap in x if i==-1
        double get_initial( double **x, 
                            double **y, 
                            int x_len,
                            int y_len, 
                            int *y2x
                           )
        {
            int min_len=getmin(x_len, y_len);
            if(min_len<=5) PrintErrorAndQuit("Sequence is too short <=5!\n");
            
            int min_ali= min_len/2;              //minimum size of considered fragment 
            if(min_ali<=5)  min_ali=5;    
            int n1, n2;
            n1 = -y_len+min_ali; 
            n2 = x_len-min_ali;

            int i, j, k, k_best;
            double tmscore, tmscore_max=-1;

            k_best=n1;
            for(k=n1; k<=n2; k++)
            {
                //get the map
                for(j=0; j<y_len; j++)
                {
                    i=j+k;
                    if(i>=0 && i<x_len)
                    {
                        y2x[j]=i;
                    }
                    else
                    {
                        y2x[j]=-1;
                    }
                }
                
                //evaluate the map quickly in three iterations
		        //this is not real tmscore, it is used to evaluate the goodness of the initial alignment
                tmscore=get_score_fast(x, y, x_len, y_len, y2x); 
                if(tmscore>=tmscore_max)
                {
                    tmscore_max=tmscore;
                    k_best=k;
                }
            }
            
            //extract the best map
            k=k_best;
            for(j=0; j<y_len; j++)
            {
                i=j+k;
                if(i>=0 && i<x_len)
                {
                    y2x[j]=i;
                }
                else
                {
                    y2x[j]=-1;
                }
            }    

            return tmscore_max;
        }

        void smooth(int *sec, int len)
        {
	        int i, j;
	        //smooth single  --x-- => -----
	        for(i=2; i<len-2; i++)
	        {
		        if(sec[i]==2 || sec[i]==4)
		        {
			        j=sec[i];
			        if(sec[i-2] != j)
			        {
				        if(sec[i-1] != j)
				        {
					        if(sec[i+1] != j)
					        {
						        if(sec[i+2] != j) 
						        {
							        sec[i]=1;
						        }
					        }
				        }
			        }
		        }
	        }

	        //   smooth double 
	        //   --xx-- => ------

	        for(i=0; i<len-5; i++)
	        {
		        //helix
		        if(sec[i] != 2)
		        {
			        if(sec[i+1] != 2)
			        {
				        if(sec[i+2] == 2)
				        {
					        if(sec[i+3] == 2)
					        {
						        if(sec[i+4] != 2)
						        {
							        if(sec[i+5] != 2)
							        {
								        sec[i+2]=1;
								        sec[i+3]=1;
							        }
						        }
					        }
				        }
			        }
		        }

		        //beta
		        if(sec[i] != 4)
		        {
			        if(sec[i+1] != 4)
			        {
				        if(sec[i+2] ==4)
				        {
					        if(sec[i+3] == 4)
					        {
						        if(sec[i+4] != 4)
						        {
							        if(sec[i+5] != 4)
							        {
								        sec[i+2]=1;
								        sec[i+3]=1;
							        }
						        }
					        }
				        }
			        }
		        }
	        }

	        //smooth connect
	        for(i=0; i<len-2; i++)
	        {		
		        if(sec[i] == 2)
		        {
			        if(sec[i+1] != 2)
			        {
				        if(sec[i+2] == 2)
				        {
					        sec[i+1]=2;
				        }
			        }
		        }
		        else if(sec[i] == 4)
		        {
			        if(sec[i+1] != 4)
			        {
				        if(sec[i+2] == 4)
				        {
					        sec[i+1]=4;
				        }
			        }
		        }
	        }

        }

        int sec_str(double dis13, double dis14, double dis15, double dis24, double dis25, double dis35)
        {
	        int s=1;
	
	        double delta=2.1;
	        if(fabs(dis15-6.37)<delta)
	        {
		        if(fabs(dis14-5.18)<delta)
		        {
			        if(fabs(dis25-5.18)<delta)
			        {
				        if(fabs(dis13-5.45)<delta)
				        {
					        if(fabs(dis24-5.45)<delta)
					        {
						        if(fabs(dis35-5.45)<delta)
						        {
							        s=2; //helix						
							        return s;
						        }
					        }
				        }
			        }
		        }
	        }

	        delta=1.42;
	        if(fabs(dis15-13)<delta)
	        {
		        if(fabs(dis14-10.4)<delta)
		        {
			        if(fabs(dis25-10.4)<delta)
			        {
				        if(fabs(dis13-6.1)<delta)
				        {
					        if(fabs(dis24-6.1)<delta)
					        {
						        if(fabs(dis35-6.1)<delta)
						        {
							        s=4; //strand
							        return s;
						        }
					        }
				        }
			        }
		        }
	        }

	        if(dis15 < 8)
	        {
		        s=3; //turn
	        }	  


	        return s;
        }


        //1->coil, 2->helix, 3->turn, 4->strand
        void make_sec(double **x, int len, int *sec)
        {
            int j1, j2, j3, j4, j5;
	        double d13, d14, d15, d24, d25, d35;
            for(int i=0; i<len; i++)
            { 	
		        sec[i]=1;
                j1=i-2;
                j2=i-1;
                j3=i;
                j4=i+1;
                j5=i+2;		
                
                if(j1>=0 && j5<len)
		        {
			        d13=sqrt(dist(x[j1], x[j3]));
			        d14=sqrt(dist(x[j1], x[j4]));
			        d15=sqrt(dist(x[j1], x[j5]));
			        d24=sqrt(dist(x[j2], x[j4]));
			        d25=sqrt(dist(x[j2], x[j5]));
			        d35=sqrt(dist(x[j3], x[j5]));
			        sec[i]=sec_str(d13, d14, d15, d24, d25, d35);			
		        }    
            } 
        }




        //get initial alignment from secondary structure alignment
        //input: x, y, x_len, y_len
        //output: y2x stores the best alignment: e.g., 
        //y2x[j]=i means:
        //the jth element in y is aligned to the ith element in x if i>=0 
        //the jth element in y is aligned to a gap in x if i==-1
        void get_initial_ss(  double **x, 
					          double **y, 
					          int x_len,
					          int y_len, 
					          int *y2x
                              )
        {
            //assign secondary structures
	        make_sec(x, x_len, secx);
	        make_sec(y, y_len, secy);

	        double gap_open=-1.0;
	        NWDP_TM(secx, secy, x_len, y_len, gap_open, y2x, val, path);    
        }


        // get_initial5 in TMalign
        //get initial alignment of local structure superposition
        //input: x, y, x_len, y_len
        //output: y2x stores the best alignment: e.g., 
        //y2x[j]=i means:
        //the jth element in y is aligned to the ith element in x if i>=0 
        //the jth element in y is aligned to a gap in x if i==-1
        bool get_initial5(double **x,
	        double **y,
	        int x_len,
	        int y_len,
	        int *y2x
	        )
        {
	        double GL, rmsd;
	        double t[3];
	        double u[3][3];

	        double d01 = d0 + 1.5;
	        if (d01 < D0_MIN) d01 = D0_MIN;
	        double d02 = d01*d01;

	        double GLmax = 0;
	        int aL = getmin(x_len, y_len);
	        int *invmap = new int[y_len + 1];

	        // jump on sequence1-------------->
	        int n_jump1 = 0;
	        if (x_len > 250)
		        n_jump1 = 45;
	        else if (x_len > 200)
		        n_jump1 = 35;
	        else if (x_len > 150)
		        n_jump1 = 25;
	        else
		        n_jump1 = 15;
	        if (n_jump1 > (x_len / 3))
		        n_jump1 = x_len / 3;

	        // jump on sequence2-------------->
	        int n_jump2 = 0;
	        if (y_len > 250)
		        n_jump2 = 45;
	        else if (y_len > 200)
		        n_jump2 = 35;
	        else if (y_len > 150)
		        n_jump2 = 25;
	        else
		        n_jump2 = 15;
	        if (n_jump2 > (y_len / 3))
		        n_jump2 = y_len / 3;

	        // fragment to superimpose-------------->
	        int n_frag[2] = { 20, 100 };
	        if (n_frag[0] > (aL / 3))
		        n_frag[0] = aL / 3;
	        if (n_frag[1] > (aL / 2))
		        n_frag[1] = aL / 2;

	        // start superimpose search-------------->
	        bool flag = false;
	        for (int i_frag = 0; i_frag < 2; i_frag++)
	        {
		        int m1 = x_len - n_frag[i_frag] + 1;
		        int m2 = y_len - n_frag[i_frag] + 1;

		        //for (int i = 1; i<m1; i = i + n_jump1) //index starts from 0, different from FORTRAN
		        // for debug
		        for (int i = 0; i<m1; i = i + n_jump1) //index starts from 0, different from FORTRAN
		        {
			        //for (int j = 1; j<m2; j = j + n_jump2)
			        for (int j = 0; j<m2; j = j + n_jump2)
			        {
				        for (int k = 0; k<n_frag[i_frag]; k++) //fragment in y
				        {
					        r1[k][0] = x[k + i][0];
					        r1[k][1] = x[k + i][1];
					        r1[k][2] = x[k + i][2];

					        r2[k][0] = y[k + j][0];
					        r2[k][1] = y[k + j][1];
					        r2[k][2] = y[k + j][2];
				        }

				        // superpose the two structures and rotate it
				        Kabsch(r1, r2, n_frag[i_frag], 1, &rmsd, t, u);

				        double gap_open = 0.0;
				        NWDP_TM(x, y, x_len, y_len, t, u, d02, gap_open, invmap, val, path);
				        GL = get_score_fast(x, y, x_len, y_len, invmap);
				        if (GL>GLmax)
				        {
					        GLmax = GL;
					        for (int ii = 0; ii<y_len; ii++)
					        {
						        y2x[ii] = invmap[ii];
					        }
					        flag = true;
				        }
			        }
		        }
	        }
	        delete[] invmap;
	        return flag;
        }

        //with invmap(i) calculate score(i,j) using RMSD rotation
        void score_matrix_rmsd(  double **x, 
						         double **y, 
						         int x_len,
						         int y_len,
						         int *y2x
						         )
        {
	        double t[3], u[3][3];
	        double rmsd, dij;
	        double d01=d0+1.5;
	        if(d01 < D0_MIN) d01=D0_MIN;
	        double d02=d01*d01;

	        double xx[3];
	        int i, k=0;
	        for(int j=0; j<y_len; j++)
	        {
		        i=y2x[j];
		        if(i>=0)
		        {
			        r1[k][0]=x[i][0];  
			        r1[k][1]=x[i][1]; 
			        r1[k][2]=x[i][2];   
                    
			        r2[k][0]=y[j][0];  
			        r2[k][1]=y[j][1]; 
			        r2[k][2]=y[j][2];
			
			        k++;
		        }
	        }
	        Kabsch(r1, r2, k, 1, &rmsd, t, u);
	        //do_rotation(x, xt, x_len, t, u);
	
	
	        for(int ii=0; ii<x_len; ii++)
	        {		
		        transform(t, u, &x[ii][0], xx);
		        for(int jj=0; jj<y_len; jj++)
		        {
			        //dij=dist(&xt[ii][0], &y[jj][0]);   
			        dij=dist(xx, &y[jj][0]); 
			        score[ii+1][jj+1] = 1.0/(1+dij/d02);
			        //	cout << ii+1 << " " << jj+1 << " " << score[ii+1][jj+1]<< endl;
		        }
	        }		
        }


        void score_matrix_rmsd_sec(  double **x, 
							         double **y, 
							         int x_len,
							         int y_len,
							         int *y2x
							         )
        {
	        double t[3], u[3][3];
	        double rmsd, dij;
	        double d01=d0+1.5;
	        if(d01 < D0_MIN) d01=D0_MIN;
	        double d02=d01*d01;

	        double xx[3];
	        int i, k=0;
	        for(int j=0; j<y_len; j++)
	        {
		        i=y2x[j];
		        if(i>=0)
		        {
			        r1[k][0]=x[i][0];  
			        r1[k][1]=x[i][1]; 
			        r1[k][2]=x[i][2];   
                    
			        r2[k][0]=y[j][0];  
			        r2[k][1]=y[j][1]; 
			        r2[k][2]=y[j][2];
			
			        k++;
		        }
	        }
	        Kabsch(r1, r2, k, 1, &rmsd, t, u);

	
	        for(int ii=0; ii<x_len; ii++)
	        {		
		        transform(t, u, &x[ii][0], xx);
		        for(int jj=0; jj<y_len; jj++)
		        {
			        dij=dist(xx, &y[jj][0]); 
			        if(secx[ii]==secy[jj])
			        {
				        score[ii+1][jj+1] = 1.0/(1+dij/d02) + 0.5;
			        }
			        else
			        {
				        score[ii+1][jj+1] = 1.0/(1+dij/d02);
			        }		
		        }
	        }		
        }


        //get initial alignment from secondary structure and previous alignments
        //input: x, y, x_len, y_len
        //output: y2x stores the best alignment: e.g., 
        //y2x[j]=i means:
        //the jth element in y is aligned to the ith element in x if i>=0 
        //the jth element in y is aligned to a gap in x if i==-1
        void get_initial_ssplus( double **x, 
						         double **y, 
						         int x_len,
						         int y_len,
						         int *y2x0,
						         int *y2x						
						         )
        {

	        //create score matrix for DP
	        score_matrix_rmsd_sec(x, y, x_len, y_len, y2x0);
	
	        double gap_open=-1.0;
	        NWDP_TM(x_len, y_len, gap_open, y2x, val, path, score);
        }


        void find_max_frag(double **x, int *resno, int len, int *start_max, int *end_max)
        {
	        int r_min, fra_min=4;           //minimum fragment for search
	        double d;
	        int start;
	        int Lfr_max=0, flag;

	        r_min= (int) (len*1.0/3.0); //minimum fragment, in case too small protein
	        if(r_min > fra_min) r_min=fra_min;
	
	        int inc=0;
	        double dcu0_cut=dcu0*dcu0;;
	        double dcu_cut=dcu0_cut;

	        while(Lfr_max < r_min)
	        {		
		        Lfr_max=0;			
		        int j=1;    //number of residues at nf-fragment
		        start=0;
		        for(int i=1; i<len; i++)
		        {			
			        d = dist(x[i-1], x[i]);
			        flag=0;
			        if(dcu_cut>dcu0_cut)
			        {
				        if(d<dcu_cut)
				        {
					        flag=1;
				        }
			        }
			        else if(resno[i] == (resno[i-1]+1)) //necessary??
			        {
				        if(d<dcu_cut)
				        {
					        flag=1;
				        }
			        }

			        if(flag==1)
			        {
				        j++;

				        if(i==(len-1))
				        {
					        if(j > Lfr_max) 
					        {
						        Lfr_max=j;
						        *start_max=start;
						        *end_max=i;						
					        }
					        j=1;
				        }
			        }
			        else
			        {
				        if(j>Lfr_max) 
				        {
					        Lfr_max=j;
					        *start_max=start;
					        *end_max=i-1;										
				        }

				        j=1;
				        start=i;
			        }
		        }// for i;
		
		        if(Lfr_max < r_min)
		        {
			        inc++;
			        double dinc=pow(1.1, (double) inc) * dcu0;
			        dcu_cut= dinc*dinc;
		        }
	        }//while <;	
        }

        //perform fragment gapless threading to find the best initial alignment
        //input: x, y, x_len, y_len
        //output: y2x0 stores the best alignment: e.g., 
        //y2x0[j]=i means:
        //the jth element in y is aligned to the ith element in x if i>=0 
        //the jth element in y is aligned to a gap in x if i==-1
        double get_initial_fgt( double **x, 
						        double **y, 
						        int x_len,
						        int y_len, 
						        int *xresno,
						        int *yresno,
						        int *y2x
						        )
        {
	        int fra_min=4;           //minimum fragment for search
	        int fra_min1=fra_min-1;  //cutoff for shift, save time

	        int xstart=0, ystart=0, xend=0, yend=0;

	        find_max_frag(x, xresno, x_len,  &xstart, &xend);
	        find_max_frag(y, yresno, y_len, &ystart, &yend);


	        int Lx = xend-xstart+1;
	        int Ly = yend-ystart+1;
	        int *ifr, *y2x_;
	        int L_fr=getmin(Lx, Ly);
	        ifr= new int[L_fr];
	        y2x_= new int[y_len+1];

	        //select what piece will be used (this may araise ansysmetry, but
	        //only when L1=L2 and Lfr1=Lfr2 and L1 ne Lfr1
	        //if L1=Lfr1 and L2=Lfr2 (normal proteins), it will be the same as initial1

	        if(Lx<Ly || (Lx==Ly && x_len<=y_len))
	        {		
		        for(int i=0; i<L_fr; i++)
		        {
			        ifr[i]=xstart+i;
		        }
	        }
	        else if(Lx>Ly || (Lx==Ly && x_len>y_len))
	        {		
		        for(int i=0; i<L_fr; i++)
		        {
			        ifr[i]=ystart+i;
		        }	
	        }

	
	        int L0=getmin(x_len, y_len); //non-redundant to get_initial1
	        if(L_fr==L0)
	        {
		        int n1= (int)(L0*0.1); //my index starts from 0
		        int n2= (int)(L0*0.89);

		        int j=0;
		        for(int i=n1; i<= n2; i++)
		        {
			        ifr[j]=ifr[i];
			        j++;
		        }
		        L_fr=j;
	        }


	        //gapless threading for the extracted fragment
	        double tmscore, tmscore_max=-1;

	        if(Lx<Ly || (Lx==Ly && x_len<=y_len))
	        {
		        int L1=L_fr;
	            int min_len=getmin(L1, y_len);    
		        int min_ali= (int) (min_len/2.5);              //minimum size of considered fragment 
		        if(min_ali<=fra_min1)  min_ali=fra_min1;    
		        int n1, n2;
		        n1 = -y_len+min_ali; 
		        n2 = L1-min_ali;

		        int i, j, k;
		        for(k=n1; k<=n2; k++)
		        {
			        //get the map
			        for(j=0; j<y_len; j++)
			        {
				        i=j+k;
				        if(i>=0 && i<L1)
				        {				
					        y2x_[j]=ifr[i];
				        }
				        else
				        {
					        y2x_[j]=-1;
				        }
			        }

			        //evaluate the map quickly in three iterations
			        tmscore=get_score_fast(x, y, x_len, y_len, y2x_);

			        if(tmscore>=tmscore_max)
			        {
				        tmscore_max=tmscore;
				        for(j=0; j<y_len; j++)
				        {
					        y2x[j]=y2x_[j];
				        }
			        }
		        }
	        }
	        else
	        {
		        int L2=L_fr;
	            int min_len=getmin(x_len, L2);    
		        int min_ali= (int) (min_len/2.5);              //minimum size of considered fragment 
		        if(min_ali<=fra_min1)  min_ali=fra_min1;    
		        int n1, n2;
		        n1 = -L2+min_ali; 
		        n2 = x_len-min_ali;

		        int i, j, k;	

		        for(k=n1; k<=n2; k++)
		        {
			        //get the map
			        for(j=0; j<y_len; j++)
			        {
				        y2x_[j]=-1;
			        }

			        for(j=0; j<L2; j++)
			        {
				        i=j+k;
				        if(i>=0 && i<x_len)
				        {
					        y2x_[ifr[j]]=i;
				        }
			        }
                
			        //evaluate the map quickly in three iterations
			        tmscore=get_score_fast(x, y, x_len, y_len, y2x_);
			        if(tmscore>=tmscore_max)
			        {
				        tmscore_max=tmscore;
				        for(j=0; j<y_len; j++)
				        {
					        y2x[j]=y2x_[j];
				        }
			        }
		        }
	        }    


	        delete [] ifr;
	        delete [] y2x_;
            return tmscore_max;
        }





        //heuristic run of dynamic programing iteratively to find the best alignment
        //input: initial rotation matrix t, u
        //       vectors x and y, d0
        //output: best alignment that maximizes the TMscore, will be stored in invmap
        double DP_iter( double **x,
                        double **y, 
                        int x_len, 
                        int y_len, 
                        double t[3],
                        double u[3][3],
				        int invmap0[],
				        int g1,
				        int g2,
				        int iteration_max,
				        double local_d0_search
				        )
        {
            double gap_open[2]={-0.6, 0};
            double rmsd; 
            int *invmap=new int[y_len+1];
            
            int iteration, i, j, k;
            double tmscore, tmscore_max, tmscore_old=0;    
            int score_sum_method=8, simplify_step=40;
            tmscore_max=-1;

	        //double d01=d0+1.5;
            double d02=d0*d0;
            for(int g=g1; g<g2; g++)
            {
                for(iteration=0; iteration<iteration_max; iteration++)
                {           
			        NWDP_TM(x, y, x_len, y_len, t, u, d02, gap_open[g], invmap, val, path);
                    
                    k=0;
                    for(j=0; j<y_len; j++) 
                    {
                        i=invmap[j];

                        if(i>=0) //aligned
                        {
                            xtm[k][0]=x[i][0];
                            xtm[k][1]=x[i][1];
                            xtm[k][2]=x[i][2];
                            
                            ytm[k][0]=y[j][0];
                            ytm[k][1]=y[j][1];
                            ytm[k][2]=y[j][2];
                            k++;
                        }
                    }

                    //tmscore=TMscore8_search(xtm, ytm, k, t, u, simplify_step, score_sum_method, &rmsd);
			        tmscore = TMscore8_search(xtm, ytm, k, t, u, simplify_step, score_sum_method, &rmsd, local_d0_search);

                   
                    if(tmscore>tmscore_max)
                    {
                        tmscore_max=tmscore;
                        for(i=0; i<y_len; i++) 
                        {                
                            invmap0[i]=invmap[i];                                      
                        }				
                    }
            
                    if(iteration>0)
                    {
                        if(fabs(tmscore_old-tmscore)<0.000001)
                        {     
                            break;       
                        }
                    }
                    tmscore_old=tmscore;
                }// for iteration           
                
            }//for gapopen
            
            
            delete []invmap;
            return tmscore_max;
        }


        void output_superpose(char *xname,
					          char *yname,
					          int x_len,
					          int y_len,
					          double t[3],
					          double u[3][3],
					          double rmsd,
					          double d0_out,
					          int m1[], 
					          int m2[],
					          int n_ali8,
					          double seq_id,
					          double TM_0,
					          double Lnorm_0,
					          double d0_0
					         )
        {
	        int i, j, j1;
	        double dis2;

	        int max=5000;


	        //aligned region
	        FILE *fp = fopen(out_reg, "w");
	        fprintf(fp, "load inline\n");
	        fprintf(fp, "select *A\n");
	        fprintf(fp, "wireframe .45\n");
	        fprintf(fp, "select *B\n");
	        fprintf(fp, "wireframe .20\n");
	        fprintf(fp, "select all\n");
	        fprintf(fp, "color white\n");


	        do_rotation(xa, xt, x_len, t, u);
	        for(i=0; i<n_ali8; i++)
	        {
		        j=m1[i];
		        j1=m2[i];
		        dis2=sqrt(dist(&xt[j][0], &ya[j1][0])) ;
		        if(dis2<=d0_out)
		        {
			        fprintf(fp, "select %4d:A,%4d:B\n", xresno[j], yresno[j1]);
			        fprintf(fp, "color red\n");
		        }
	        }

	        fprintf(fp, "select all\n");
	        fprintf(fp, "exit\n");
	        string basename(xname);
	        int idx = basename.find_last_of("\\");
	        string xtempbase = basename.substr(idx + 1, basename.length());
	        fprintf(fp, "REMARK Structure A:%s   Size=%4d\n", xtempbase.c_str(), x_len);
	        basename = string(yname);
	        idx = basename.find_last_of("\\");
	        string ytempbase = basename.substr(idx + 1, basename.length());
	        fprintf(fp, "REMARK Structure B:%s   Size=%4d\n", ytempbase.c_str(), y_len);
	        fprintf(fp, "(TM-score is normalized by %4d, d0=%6.2f)\n", int(Lnorm_0), d0_0);
	        fprintf(fp, "REMARK Aligned length=%4d, RMSD=%6.2f, TM-score=%7.5f, ID=%5.3f\n", n_ali8, rmsd, TM_0, seq_id);


	        char AA[4];
	        //superposed structure B
	        for(i=0; i<n_ali8; i++)
	        {
		        j=m1[i];
		        AAmap3(seqx[j], AA);
		        fprintf(fp, "ATOM  %5d  CA  %3s A%4d%c   %8.3f%8.3f%8.3f\n", j + 1, AA, xresno[j], ins1[j], xt[j][0], xt[j][1], xt[j][2]);
	        }
	        fprintf(fp, "TER\n");

	        for(i=1; i<n_ali8; i++)
	        {
		        j=m1[i-1]+1;
		        j1=m1[i]+1;
		        fprintf(fp, "CONECT%5d%5d\n", j, j1);
	        }
	        //structure A
	        for(i=0; i<n_ali8; i++)
	        {
		        j=m2[i];
		        AAmap3(seqy[j], AA);
		        fprintf(fp, "ATOM  %5d  CA  %3s B%4d%c   %8.3f%8.3f%8.3f\n", j + max+ 1, AA, yresno[j], ins2[j], ya[j][0], ya[j][1], ya[j][2]);
	        }
	        fprintf(fp, "TER\n");
	        for(i=1; i<n_ali8; i++)
	        {
		        j=max+m2[i-1]+1;
		        j1=max+m2[i]+1;
		        fprintf(fp, "CONECT%5d%5d\n", j, j1);
	        }

	        fclose(fp);







	        //  output CA - trace of whole chain in 'TM.sup_all' -------->
	        //all regions
	        char str[3000];
	        sprintf(str, "%s_all", out_reg);
	        fp=fopen(str, "w");
	        fprintf(fp, "load inline\n");
	        fprintf(fp, "select *A\n");
	        fprintf(fp, "wireframe .45\n");
	        fprintf(fp, "select none\n");
	        fprintf(fp, "select *B\n");
	        fprintf(fp, "wireframe .20\n");
	        fprintf(fp, "color white\n");

	        for(i=0; i<n_ali8; i++)
	        {
		        j=m1[i];
		        j1=m2[i];
		        dis2=sqrt(dist(&xt[j][0], &ya[j1][0])) ;
		        if(dis2<=d0_out)
		        {
			        fprintf(fp, "select %4d:A,%4d:B\n", xresno[j], yresno[j1]);
			        fprintf(fp, "color red\n");
		        }
	        }


	        fprintf(fp, "select all\n");
	        fprintf(fp, "exit\n");
	        fprintf(fp, "REMARK Structure A:%s   Size=%4d\n", xtempbase.c_str(), x_len);
	        fprintf(fp, "REMARK Structure B:%s   Size=%4d ", ytempbase.c_str(), y_len);
	        fprintf(fp, "(TM-score is normalized by %4d, d0=%6.2f)\n", int(Lnorm_0), d0_0);
	        fprintf(fp, "REMARK Aligned length=%4d, RMSD=%6.2f, TM-score=%7.5f, ID=%5.3f\n", n_ali8, rmsd, TM_0, seq_id);



	        //superposed structure B
	        for(i=0; i<x_len; i++)
	        {
		        j=i;
		        AAmap3(seqx[j], AA);
		        fprintf(fp, "ATOM  %5d  CA  %3s A%4d%c   %8.3f%8.3f%8.3f\n", j + 1, AA, xresno[j], ins1[j], xt[j][0], xt[j][1], xt[j][2]);
	        }
	        fprintf(fp, "TER\n");

	        for(i=1; i<x_len; i++)
	        {
		        fprintf(fp, "CONECT%5d%5d\n", i, i+1);
	        }
	        //structure A
	        for(i=0; i<y_len; i++)
	        {
		        j=i;
		        AAmap3(seqy[j], AA);
		        fprintf(fp, "ATOM  %5d  CA  %3s B%4d%c   %8.3f%8.3f%8.3f\n", j + max + 1, AA, yresno[j], ins2[j], ya[j][0], ya[j][1], ya[j][2]);
	        }
	        fprintf(fp, "TER\n");
	        for(i=1; i<y_len; i++)
	        {
		        fprintf(fp, "CONECT%5d%5d\n", max+i, max+i+1);
	        }
	        fclose(fp);


	        ///////  output full - atomic structure of whole chain in 'TM.sup_atm' -------->
	        sprintf(str, "%s_atm", out_reg);
	        fp = fopen(str, "w");
	        fprintf(fp, "load inline\n");
	        fprintf(fp, "select *A\n");
	        fprintf(fp, "color blue\n");
	        fprintf(fp, "select *B\n");
	        fprintf(fp, "color red\n");
	        fprintf(fp, "select all\n");
	        fprintf(fp, "cartoon\n");
	        fprintf(fp, "exit\n");

	        fprintf(fp, "REMARK Structure A:%s   Size=%4d\n", xtempbase.c_str(), x_len);
	        fprintf(fp, "REMARK Structure B:%s   Size=%4d ", ytempbase.c_str(), y_len);
	        fprintf(fp, "(TM-score is normalized by %4d, d0=%6.2f)\n", int(Lnorm_0), d0_0);
	        fprintf(fp, "REMARK Aligned length=%4d, RMSD=%6.2f, TM-score=%7.5f, ID=%5.3f\n", n_ali8, rmsd, TM_0, seq_id);

	        // chain_1: structure A
	        double **coord;
	        NewArray(&coord, atomxlen, 3);
	        do_rotation(xyza1, coord, atomxlen, t, u);
	        for (i = 0; i < atomxlen; i++)
	        {
		        for (j = 0; j < n_ali8; j++)
		        {
			        if (ir1[i] == xresno[m1[j]])
			        {
				        if (ains1[i] == ins1[m1[j]])
					        fprintf(fp, "ATOM  %5d  %-4s%3s A%4d%c   %8.3f%8.3f%8.3f\n", ia1[i], aa1[i], ra1[i], ir1[i], ains1[i], coord[i][0], coord[i][1], coord[i][2]);
			        }
		        }
	        }
	        fprintf(fp, "TER\n");

	        // chain_2: structure B, coordinates don't change at all
	        for (i = 0; i < atomylen; i++)
	        {
		        for (j = 0; j < n_ali8; j++)
		        {
			        if (ir2[i] == yresno[m2[j]])
			        {
				        if (ains2[i] == ins2[m2[j]])
					        fprintf(fp, "ATOM  %5d  %-4s%3s B%4d%c   %8.3f%8.3f%8.3f\n", ia2[i], aa2[i], ra2[i], ir2[i], ains2[i], xyza2[i][0], xyza2[i][1], xyza2[i][2]);
			        }
		        }
	        }
	        fprintf(fp, "TER\n");

	        fclose(fp);


	        ///////  output full - atomic structure of whole chain in 'TM.sup_all_atm' -------->
	        sprintf(str, "%s_all_atm", out_reg);
	        fp = fopen(str, "w");
	        fprintf(fp, "load inline\n");
	        fprintf(fp, "select *A\n");
	        fprintf(fp, "color blue\n");
	        fprintf(fp, "select *B\n");
	        fprintf(fp, "color red\n");
	        fprintf(fp, "select all\n");
	        fprintf(fp, "cartoon\n");
	        fprintf(fp, "exit\n");

	        fprintf(fp, "REMARK Structure A:%s   Size=%4d\n", xtempbase.c_str(), x_len);
	        fprintf(fp, "REMARK Structure B:%s   Size=%4d ", ytempbase.c_str(), y_len);
	        fprintf(fp, "(TM-score is normalized by %4d, d0=%.2f)\n", int(Lnorm_0), d0_0);
	        fprintf(fp, "REMARK Aligned length=%4d, RMSD=%6.2f, TM-score=%7.5f, ID=%5.3f\n", n_ali8, rmsd, TM_0, seq_id);

	        // chain_1: structure A
	        for (i = 0; i < atomxlen; i++)
		        fprintf(fp, "ATOM  %5d  %-4s%3s A%4d%c   %8.3f%8.3f%8.3f\n", ia1[i], aa1[i], ra1[i], ir1[i], ains1[i], coord[i][0], coord[i][1], coord[i][2]);
	        fprintf(fp, "TER\n");
	        DeleteArray(&coord, atomxlen);

	        // chain_2: structure B, coordinates don't change at all
	        for (i = 0; i < atomylen; i++)
		        fprintf(fp, "ATOM  %5d  %-4s%3s B%4d%c   %8.3f%8.3f%8.3f\n", ia2[i], aa2[i], ra2[i], ir2[i], ains2[i], xyza2[i][0], xyza2[i][1], xyza2[i][2]);
	        fprintf(fp, "TER\n");
	        fclose(fp);


	        ///////  output full - atomic structure of whole chain in 'TM.sup_all_atm_lig' -------->
	        sprintf(str, "%s_all_atm_lig", out_reg);
	        fp = fopen(str, "w");
	        fprintf(fp, "load inline\n");
	        fprintf(fp, "select all\n");
	        fprintf(fp, "cartoon\n");
	        fprintf(fp, "select *A\n");
	        fprintf(fp, "color blue\n");
	        fprintf(fp, "select *B\n");
	        fprintf(fp, "color red\n");
	        fprintf(fp, "select ligand\n");
	        fprintf(fp, "wireframe 0.25\n");
	        fprintf(fp, "select solvent\n");
	        fprintf(fp, "spacefill 0.25\n");
	        fprintf(fp, "select all\n");
	        fprintf(fp, "exit\n");

	        fprintf(fp, "REMARK Structure A:%s   Size=%4d\n", xtempbase.c_str(), x_len);
	        fprintf(fp, "REMARK Structure B:%s   Size=%4d ", ytempbase.c_str(), y_len);
	        fprintf(fp, "(TM-score is normalized by %4d, d0=%.2f)\n", int(Lnorm_0), d0_0);
	        fprintf(fp, "REMARK Aligned length=%4d, RMSD=%6.2f, TM-score=%7.5f, ID=%5.3f\n", n_ali8, rmsd, TM_0, seq_id);

	        // chain_1: structure A
	        int ligLen = get_ligand_len(xname);
	        NewArray(&coord, ligLen, 3);
	        string* strSeq = new string[ligLen];
	        double** coord2;
	        NewArray(&coord2, ligLen, 3);
	        read_ligand(xname, coord, strSeq);
	        do_rotation(coord, coord2, ligLen, t, u);
	        for (i = 0; i < ligLen; i++)
	        {
		        fprintf(fp, "%21sA%8s%8.3f%8.3f%8.3f\n", strSeq[i].substr(0, 21).c_str(), strSeq[i].substr(22, 8).c_str(), coord2[i][0], coord2[i][1], coord2[i][2]);
	        }
	        fprintf(fp, "TER\n");
	        DeleteArray(&coord, ligLen);
	        DeleteArray(&coord2, ligLen);
	        delete []strSeq;

	        // chain_2: structure B, coordinates don't change at all
	        ligLen = get_ligand_len(yname);
	        NewArray(&coord, ligLen, 3);
	        strSeq = new string[ligLen];
	        read_ligand(yname, coord, strSeq);
	        for (i = 0; i < ligLen; i++)
	        {
		        fprintf(fp, "%21sB%8s%8.3f%8.3f%8.3f\n", strSeq[i].substr(0, 21).c_str(), strSeq[i].substr(22, 8).c_str(), coord[i][0], coord[i][1], coord[i][2]);
	        }
	        fprintf(fp, "TER\n");
	        DeleteArray(&coord, ligLen);
	        delete []strSeq;

	        fclose(fp);
        }


        //output the final results
        void output_results(char *xname,
					         char *yname,
					         int x_len,
					         int y_len,
					         double t[3],
					         double u[3][3],
					         double TM1,
					         double TM2,
					         double rmsd,
					         double d0_out,
					         int m1[], 
					         int m2[],
					         int n_ali8,
					         int n_ali,
					         double TM_0,
					         double Lnorm_0,
					         double d0_0
					         )
        {
            double seq_id;          
            int i, j, k;
            double d;
            int ali_len=x_len+y_len; //maximum length of alignment
	        char *seqM, *seqxA, *seqyA;
	        seqM=new char[ali_len];
	        seqxA=new char[ali_len];
	        seqyA=new char[ali_len];
	

	        do_rotation(xa, xt, x_len, t, u);

	        seq_id=0;
	        int kk=0, i_old=0, j_old=0;
	        for(k=0; k<n_ali8; k++)
	        {
		        for(i=i_old; i<m1[k]; i++)
		        {
			        //align x to gap
			        seqxA[kk]=seqx[i];
			        seqyA[kk]='-';
			        seqM[kk]=' ';					
			        kk++;
		        }

		        for(j=j_old; j<m2[k]; j++)
		        {
			        //align y to gap
			        seqxA[kk]='-';
                    seqyA[kk]=seqy[j];
                    seqM[kk]=' ';
                    kk++;
		        }

		        seqxA[kk]=seqx[m1[k]];
		        seqyA[kk]=seqy[m2[k]];
		        if(seqxA[kk]==seqyA[kk])
		        {
			        seq_id++;
		        }
		        d=sqrt(dist(&xt[m1[k]][0], &ya[m2[k]][0]));
		        if(d<d0_out)
		        {
			        seqM[kk]=':';
		        }
		        else
		        {
			        seqM[kk]='.';
		        } 
		        kk++;  
		        i_old=m1[k]+1;
		        j_old=m2[k]+1;
	        }

	        //tail
	        for(i=i_old; i<x_len; i++)
	        {
		        //align x to gap
		        seqxA[kk]=seqx[i];
		        seqyA[kk]='-';
		        seqM[kk]=' ';					
		        kk++;
	        }    
	        for(j=j_old; j<y_len; j++)
	        {
		        //align y to gap
		        seqxA[kk]='-';
		        seqyA[kk]=seqy[j];
		        seqM[kk]=' ';
		        kk++;
	        }
         
            seqxA[kk]='\0';
            seqyA[kk]='\0';
            seqM[kk]='\0';

	        seq_id=seq_id/( n_ali8+0.00000001); //what did by TMalign, but not reasonable, it should be n_ali8    

	        cout <<endl;	
	
	        printf("\nName of Chain_1: %s (to be superimposed onto Chain_2)\n", xname); 
	        printf("Name of Chain_2: %s\n", yname);
	        printf("Length of Chain_1: %d residues\n", x_len);
	        printf("Length of Chain_2: %d residues\n\n", y_len);

	        if (i_opt || I_opt)
	        {
		        printf("User-specified initial alignment: TM/Lali/rmsd = %7.5lf, %4d, %6.3lf\n", TM_ali, L_ali, rmsd_ali);
	        }

	        printf("Aligned length= %d, RMSD= %6.2f, Seq_ID=n_identical/n_aligned= %4.3f\n", n_ali8, rmsd, seq_id); 
	        printf("TM-score= %6.5f (if normalized by length of Chain_1, i.e., LN=%d, d0=%.2f)\n", TM2, x_len, d0B);
	        printf("TM-score= %6.5f (if normalized by length of Chain_2, i.e., LN=%d, d0=%.2f)\n", TM1, y_len, d0A);
	        if(a_opt)
	        {
	          double L_ave=(x_len+y_len)*0.5;
	          printf("TM-score= %6.5f (if normalized by average length of two structures, i.e., LN= %.2f, d0= %.2f)\n", TM3, L_ave, d0a);
	        }
	        if(u_opt)
	        {		
	          printf("TM-score= %6.5f (if normalized by user-specified LN=%.2f and d0=%.2f)\n", TM4, Lnorm_ass, d0u);
	        }
	        if(d_opt)
	          {		
	            printf("TM-score= %6.5f (if scaled by user-specified d0= %.2f, and LN= %.2f)\n", TM5, d0_scale, Lnorm_0);
	        }
	        printf("(You should use TM-score normalized by length of the reference protein)\n");

            
            //output structure alignment
            printf("\n(\":\" denotes residue pairs of d < %4.1f Angstrom, ", d0_out);
            printf("\".\" denotes other aligned residues)\n");
            printf("%s\n", seqxA);
            printf("%s\n", seqM);
            printf("%s\n", seqyA);

	        cout << endl;

	        if(o_opt)
	        {
		        output_superpose(xname, yname, x_len, y_len, t, u, rmsd, d0_out, m1, m2, n_ali8, seq_id, TM_0, Lnorm_0, d0_0);
	        }


	        delete [] seqM;
	        delete [] seqxA;
	        delete [] seqyA;
            
        }

        double standard_TMscore(double **x, double **y, int x_len, int y_len, int invmap[], int& L_ali, double& RMSD )
        {
	        D0_MIN = 0.5;
	        Lnorm = y_len;
	        if (Lnorm > 21)
		        d0 = (1.24*pow((Lnorm*1.0 - 15), 1.0 / 3) - 1.8);
	        else
		        d0 = D0_MIN;
	        if (d0 < D0_MIN)
		        d0 = D0_MIN;
	        double d0_input = d0;// Scaled by seq_min

	        double tmscore;// collected alined residues from invmap
	        int n_al = 0;
	        int i;
	        for (int j = 0; j<y_len; j++)
	        {
		        i = invmap[j];
		        if (i >= 0)
		        {
			        xtm[n_al][0] = x[i][0];
			        xtm[n_al][1] = x[i][1];
			        xtm[n_al][2] = x[i][2];

			        ytm[n_al][0] = y[j][0];
			        ytm[n_al][1] = y[j][1];
			        ytm[n_al][2] = y[j][2];

			        r1[n_al][0] = x[i][0];
			        r1[n_al][1] = x[i][1];
			        r1[n_al][2] = x[i][2];

			        r2[n_al][0] = y[j][0];
			        r2[n_al][1] = y[j][1];
			        r2[n_al][2] = y[j][2];

			        n_al++;
		        }
		        else if (i != -1)
		        {
			        PrintErrorAndQuit("Wrong map!\n");
		        }
	        }
	        L_ali = n_al;

	        Kabsch(r1, r2, n_al, 0, &RMSD, t, u);
	        RMSD = sqrt( RMSD/(1.0*n_al) );
	
	        int temp_simplify_step = 1;
	        int temp_score_sum_method = 0;
	        d0_search = d0_input;
	        double rms = 0.0;
	        tmscore = TMscore8_search_standard(xtm, ytm, n_al, t, u, temp_simplify_step, temp_score_sum_method, &rms, d0_input);
	        tmscore = tmscore * n_al / (1.0*Lnorm);

	        return tmscore;
        }
};
