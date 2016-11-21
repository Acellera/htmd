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
#ifndef BASICDEFINE_H
#define BASICDEFINE_H

#define getmax(a,b) a>b?a:b
#define getmin(a,b) a>b?b:a

#define AtomLenMax 30

#define dimmax 30000
#define charmax 10

#define NMAX 5000
#define NMAX2 90000
#define ASCIILimit 123

#define MAXLEN 10000 //maximum length of filenames
#include <string>// TO use string variables

using namespace std;

string **atom1, **atom2;// atom name, atom1(i,j): the name of jth atom for ith residue

#endif
