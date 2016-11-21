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

#define ATOMNMAX 90000//upper limit for number of ATOM, 
#define CANMAX 5000//upper limit for number of C-alpha, 


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <malloc.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>

#include "basic_define.h"



void PrintErrorAndQuit(const char* sErrorString)
{
	cout << sErrorString << endl;
	exit(1);
}


template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
  *array=new A* [Narray1];
  for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
};

template <class A> void DeleteArray(A *** array, int Narray)
{
  for(int i=0; i<Narray; i++)
    if(*(*array+i)) delete [] *(*array+i);
  if(Narray) delete [] (*array);
  (*array)=NULL;
};


char AAmap(string AA)
{
    char A=' ';
    if(     AA.compare("BCK")==0)   A='X';
    else if(AA.compare("GLY")==0)   A='G';
    else if(AA.compare("ALA")==0)   A='A';
    else if(AA.compare("SER")==0)   A='S';
    else if(AA.compare("CYS")==0)   A='C';
    else if(AA.compare("VAL")==0)   A='V';     
    else if(AA.compare("THR")==0)   A='T';
    else if(AA.compare("ILE")==0)   A='I';
    else if(AA.compare("PRO")==0)   A='P';
    else if(AA.compare("MET")==0)   A='M';
    else if(AA.compare("ASP")==0)   A='D';
    else if(AA.compare("ASN")==0)   A='N';
    else if(AA.compare("LEU")==0)   A='L';
    else if(AA.compare("LYS")==0)   A='K';
    else if(AA.compare("GLU")==0)   A='E';
    else if(AA.compare("GLN")==0)   A='Q';
    else if(AA.compare("ARG")==0)   A='R';
    else if(AA.compare("HIS")==0)   A='H';
    else if(AA.compare("PHE")==0)   A='F';
    else if(AA.compare("TYR")==0)   A='Y';
    else if(AA.compare("TRP")==0)   A='W';    
    else if(AA.compare("CYX")==0)   A='C';
    else
        A='Z'; //ligand
        
    return A;
}

void AAmap3(char A, char AA[3])
{
    if     ( A=='X')   strcpy(AA, "BCK");
	else if( A=='G')   strcpy(AA, "GLY");
	else if( A=='A')   strcpy(AA, "ALA");
	else if( A=='S')   strcpy(AA, "SER");
	else if( A=='C')   strcpy(AA, "CYS");
	else if( A=='V')   strcpy(AA, "VAL");
	else if( A=='T')   strcpy(AA, "THR");
	else if( A=='I')   strcpy(AA, "ILE");
	else if( A=='P')   strcpy(AA, "PRO");
	else if( A=='M')   strcpy(AA, "MET");
	else if( A=='D')   strcpy(AA, "ASP");
	else if( A=='N')   strcpy(AA, "ASN");
	else if( A=='L')   strcpy(AA, "LEU");
	else if( A=='K')   strcpy(AA, "LYS");
	else if( A=='E')   strcpy(AA, "GLU");
	else if( A=='Q')   strcpy(AA, "GLN");
	else if( A=='R')   strcpy(AA, "ARG");
	else if( A=='H')   strcpy(AA, "HIS");
	else if( A=='F')   strcpy(AA, "PHE");
	else if( A=='Y')   strcpy(AA, "TYR");
	else if( A=='W')   strcpy(AA, "TRP");
	else if( A=='C')   strcpy(AA, "CYX");
    else
        strcpy(AA, "UNK");           
}


void get_xyz(string line, double *x, double *y, double *z, char *resname, int *no)
{
    char cstr[50];    
    
    strcpy(cstr, (line.substr(30, 8)).c_str());
    sscanf(cstr, "%lf", x);
    
    strcpy(cstr, (line.substr(38, 8)).c_str());
	sscanf(cstr, "%lf", y);
    
    strcpy(cstr, (line.substr(46, 8)).c_str());
	sscanf(cstr, "%lf", z);
    
    strcpy(cstr, (line.substr(17, 3)).c_str());
	*resname = AAmap(cstr);

    strcpy(cstr, (line.substr(22, 4)).c_str());
	sscanf(cstr, "%d", no);
}

string Trim(string inputString)
{
	string result = inputString;
	int idxBegin = inputString.find_first_not_of(" ");
	int idxEnd = inputString.find_last_not_of(" ");
	if (idxBegin >= 0 && idxEnd >= 0)
		result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
	return result;
}

// This function is used to get data for output, the data will not be used for alignment
void get_xyz(string line, int *ia1, char *aa, char *ra, int *ir, double *x, double *y, double *z)
{
	string tempLine = line.substr(6, 5);
	tempLine = Trim(tempLine);
	sscanf(tempLine.c_str(), "%d", ia1);

	tempLine = line.substr(12, 4);
	tempLine = Trim(tempLine);
	strcpy(aa, tempLine.c_str());

	strcpy(ra, (line.substr(17, 3)).c_str());

	char cstr[50];
	strcpy(cstr, (line.substr(22, 4)).c_str());
	sscanf(cstr, "%d", ir);

	strcpy(cstr, (line.substr(30, 8)).c_str());
	sscanf(cstr, "%lf", x);

	strcpy(cstr, (line.substr(38, 8)).c_str());
	sscanf(cstr, "%lf", y);

	strcpy(cstr, (line.substr(46, 8)).c_str());
	sscanf(cstr, "%lf", z);
}

void get_xyz(string line, double *x, double *y, double *z, string *resname)
{
	char cstr[50];

	strcpy(cstr, (line.substr(30, 8)).c_str());
	sscanf(cstr, "%lf", x);

	strcpy(cstr, (line.substr(38, 8)).c_str());
	sscanf(cstr, "%lf", y);

	strcpy(cstr, (line.substr(46, 8)).c_str());
	sscanf(cstr, "%lf", z);

	(*resname) = line.substr(0, 30);
}

int get_PDB_len(char *filename)
{
	int i = 0;
	string line;
	string atom("ATOM ");

	ifstream fin(filename);
	if (fin.is_open())
	{
		while (fin.good())
		{
			getline(fin, line);
			if (line.compare(0, atom.length(), atom) == 0)
				i++;
		}
		fin.close();
	}
	else
	{
		char message[5000];
		sprintf(message, "Can not open file: %s\n", filename);
		PrintErrorAndQuit(message);
	}

	return i;
}

int read_PDB(char *filename, double **a, char *seq, int *resno, int **nres)
{
    int i=0;
    string line, str;    
    string atom ("ATOM "); 
	string ter("TER");
	string du1, i8;
    
	int mk = 1;
    ifstream fin (filename);
    if (fin.is_open())
    {
		bool bGoRead = true;
		while (fin.good() && bGoRead)
        {
            getline(fin, line);
			if (i > 0 && line.compare(0, 3, ter) == 0)
				bGoRead = false;
			else
			{
				if(line.compare(0, atom.length(), atom)==0)
				{
					if( line.compare(12, 4, "CA  ")==0 ||\
						line.compare(12, 4, " CA ")==0 ||\
						line.compare(12, 4, "  CA")==0 )
					{
						du1 = line.substr(26, 1);
						int nDu1 = *(du1.c_str());// the ASCII code of du1

						mk = 1;
						if (line.compare(16, 1, "") != 0)
						{
							i8 = line.substr(22, 4);// get index of residue
							const char* pi8 = i8.c_str();
							int numi8 = atoi(pi8);
							if (nres[numi8][nDu1] >= 1)// atom0[i][j]: i,j index begin from one
								mk = -1;
						}
						
						if (mk == 1)
						{  
							get_xyz(line, &a[i][0], &a[i][1], &a[i][2], &seq[i], &resno[i]);
							
							nres[resno[i]][nDu1] ++;
							i++;
						}
					}                  
				}    			
			}        
        }
        fin.close();
    }
    else
    {
		char message[5000];
		sprintf(message, "Can not open file: %s\n", filename);
        PrintErrorAndQuit(message);
    } 
    seq[i]='\0';   
    
    return i;
}

int read_PDB_fullatom(char *filename, double **a, char *seq, int *resno, int *ia, char **aa, char **ra, int *ir, double **xyza,
	int **nres, string **atom0, char* ains, char* ins, int& atomlen)
{
	char dest[1000];


	int i = 0;
	string line, str;
	string atom("ATOM ");
	string ter("TER");

	int na = 0;
	int ntt = 0;
	int mk = 1;
	string du1, du2, i8;
	ifstream fin(filename);
	if (fin.is_open())
	{
		bool bGoRead = true;
		while (fin.good() && bGoRead)
		{
			getline(fin, line);
			if (i > 0 && line.compare(0, 3, ter) == 0)
			{
				bGoRead = false;
			}
			else
			{			
				if (line.compare(0, atom.length(), atom) == 0)
				{
					du1 = line.substr(26, 1);
					int nDu1 = *(du1.c_str());// the ASCII code of du1

					mk = 1;// Get flag for determination of nres, and atom0
					if (line.compare(16, 1, " ") != 0)//if(s(17:17).ne.' '), find alternate atom
					{
						du2 = line.substr(12, 4);//get name of atom0
						int idxBegin = du2.find_first_not_of(" ");
						int idxEnd = du2.find_last_not_of(" ");
						if (idxBegin>=0 && idxEnd >=0)
							du2 = du2.substr(idxBegin, idxEnd + 1 - idxBegin);

						i8 = line.substr(22, 4);// get index of residue
						const char* pi8 = i8.c_str();
						int numi8 = atoi(pi8);
						int tempNo = nres[numi8][nDu1];// nres[i][j]: i index begins from one, j index begins from 32 to 122
						for (int i1 = 1; i1 <= tempNo; i1++)
						{
							if (du2 == atom0[numi8][i1])// atom0[i][j]: i,j index begin from one
								mk = -1;
						}
					}

					if (mk == 1)
					{
						ntt++;
						if (ntt < NMAX2)
						{
							get_xyz(line, &ia[na], &aa[na][0], &ra[na][0], &ir[na], &xyza[na][0], &xyza[na][1], &xyza[na][2]);// Get data for output	
						
							int numi8 = ir[na];
							int tempNo = nres[numi8][nDu1];
							if (tempNo < AtomLenMax - 1)
							{
								tempNo++;
								nres[numi8][nDu1] = tempNo;
								atom0[numi8][tempNo] = aa[na];
							}
							strcpy(&ains[na], du1.c_str());
							na++;

							if (line.compare(12, 4, "CA  ") == 0 || \
								line.compare(12, 4, " CA ") == 0 || \
								line.compare(12, 4, "  CA") == 0)
							{
								get_xyz(line, &a[i][0], &a[i][1], &a[i][2], &seq[i], &resno[i]);
								strcpy(&ins[i], du1.c_str());
								i++;
							}// end if("CA")
						}// end if (ntt < NMAX2)
						else
							bGoRead = false;// Stop reading data from PDB file
					}				
				}// end if (line.compare(0, atom.length(), atom) == 0), "ATOM"
			}
		}
		fin.close();
	}
	else
	{
		char message[5000];
		sprintf(message, "Can not open file: %s\n", filename);
		PrintErrorAndQuit(message);
	}
	seq[i] = '\0';
	
	ins[i] = '\0';
	ains[na] = '\0';

	atomlen = na;
	return i;
}

int get_ligand_len(char *filename)
{
    int i=0;
    string line;    
	string atom1("ATOM  ");
	string atom2("HETATM");
    string finish ("END"); 
    
    ifstream fin (filename);
    if (fin.is_open())
    {
        while ( fin.good() )
        {
            getline(fin, line);
			if (line.compare(0, atom1.length(), atom1) == 0 || line.compare(0, atom2.length(), atom2) == 0)
            {
               i++;
            }
            else if(line.compare(0, finish.length(), finish)==0) 
            {
                break;
            }          
        }
        fin.close();
    }
    else
    {
		char message[5000];
		sprintf(message, "Can not open file: %s\n", filename);
        PrintErrorAndQuit(message);
    } 
    
    return i;
}

int read_ligand(char *filename, double **a, string *seq)
{
    int i=0;
    char cstr[100];    
    string line, du;    
	string atom1("ATOM  ");
	string atom2 ("HETATM"); 
    string finish ("END"); 
    
    
    ifstream fin (filename);
    if (fin.is_open())
    {
        while ( fin.good() )
        {
            getline(fin, line);
            if(line.compare(0, atom1.length(), atom1)==0 || line.compare(0, atom2.length(), atom2)==0)
            {
				get_xyz(line, &a[i][0], &a[i][1], &a[i][2], &seq[i]);
                i++;
            }
            else if(line.compare(0, finish.length(), finish)==0) 
            {
                break;
            }          
        }
        fin.close();
    }
    else
    {
		char message[5000];
		sprintf(message, "Can not open file: %s\n", filename);
        PrintErrorAndQuit(message);     
    } 
    
    return i;
}



double dist(double x[3], double y[3])
{
	double d1=x[0]-y[0];
	double d2=x[1]-y[1];
	double d3=x[2]-y[2];	
 
    return (d1*d1 + d2*d2 + d3*d3);
}

double dot(double *a, double *b)
{
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

void transform(double t[3], double u[3][3], double *x, double *x1)
{
    x1[0]=t[0]+dot(&u[0][0], x);
    x1[1]=t[1]+dot(&u[1][0], x);
    x1[2]=t[2]+dot(&u[2][0], x);
}

void do_rotation(double **x, double **x1, int len, double t[3], double u[3][3])
{
    for(int i=0; i<len; i++)
    {
        transform(t, u, &x[i][0], &x1[i][0]);
    }    
}



void output_align1(int *invmap0, int len)
{
	for(int i=0; i<len; i++)
	{
		if(invmap0[i]>=0)
		{
			cout << invmap0[i]+1 << " ";
		}			
		else
			cout << invmap0[i] << " ";

	}	
	cout << endl << endl;	
}


int output_align(int *invmap0, int len)
{
	int n_ali=0;
	for(int i=0; i<len; i++)
	{		
		cout <<  invmap0[i] << " ";
		n_ali++;		
	}	
	cout << endl << endl;	

	return n_ali;
}
// output aligned residues to string
string output_align_to_string(int *invmap, int len)
{
	string result = "";
	int step = 10;
	int rest = len % step;
	int loop = len / step;
	int i,k;
	char temp[1000];
	for ( i = 0; i<loop; i++)
	{
		k = i * step;
		for (int j = 0; j < step; j++)
		{
			sprintf(temp, "%4d ", invmap[k + j]);
			result += string(temp);
		}
		result += "\n";
	}
	k = i * step;
	for (int j = 0; j < rest; j++)
	{
		sprintf(temp, "%4d ", invmap[k + j]);
		result += string(temp);
	}
	result += "\n";

	return result;
}

// output aligned coords to string
string output_coord_to_string(double **x, double **y, int len)
{
	string result = "";
	int i, k;
	char temp[1000];
	result += "xtm coords:\n";
	for (i = 0; i<len; i++)
	{
		sprintf(temp, "%8.3f %8.3f %8.3f\n", x[i][0], x[i][1], x[i][2]);
		result += string(temp);
	}
	result += "ytm coords:\n";
	for (i = 0; i<len; i++)
	{
		sprintf(temp, "%8.3f %8.3f %8.3f\n", y[i][0], y[i][1], y[i][2]);
		result += string(temp);
	}
	return result;
}
