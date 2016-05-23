/* (c) Benoit Roux, Lei Huang    */
/* Licensed under GPL version 2  */
/* Modifications by Acellera Ltd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ff.h"

#ifndef FALSE
#define FALSE	0
#define TRUE	1
#endif

#ifndef MINAB
#define MINAB 
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif

#define N_MAX_DIH	(32)
#define MAX_RAND_INI	(50)
#define MAX_POINT	(6000)
#define MAX_CLUSTER	(8)
double Centers[MAX_CLUSTER], withinss[MAX_CLUSTER];
int n_Point_List[MAX_CLUSTER];

int DihList[N_MAX_DIH][4], n_Phi;

int n_Point=0, iseed=-1275;
double *pData=NULL, *pSilhouette;
int *GroupID=NULL;
int n_Cluster;

char szPhiToScan[]="soft-dih-list.txt";

FILE *fFile_Run_Log;     // will be shared by other source code just for compiling

int Read_1D_Data(char szName[]);
int kmeans_Lloyd_1D(double *x, int n, double *cen, int k, int *cl, int maxiter, int *nc, double *wss, double& wss_Sum);
double rand(int &iseed);
void Output_Data(int Idx);
double Cal_Silhouette(void);
double Cal_Silhouette_a_b(int i, int Idx);
int To_Find_Closest_Cluster(int i);
int Is_A_Single_Cluster(void);
double Do_Clustering(int k);
void Setup_Centers(int k);
void Shuffle_Data(int nStep);
int Read_Soft_DihedralList(void);
void Quit_With_Error_Msg(char szMsg[]);
void Split_Phi_Traj(void);
void SortCenters(void);

int main(void)
{
	int n_Correct=-1, i, Idx_Phi;
	double Silhouette, Silhouette_Max=-1.0E100;
	char szName[256];
	FILE *fOut;

  timebomb();

	Read_Soft_DihedralList();
	Split_Phi_Traj();


	fOut = fopen("soft-dih-list-new.txt", "w");
	for(Idx_Phi=0; Idx_Phi<n_Phi; Idx_Phi++)	{
		sprintf(szName, "phi-traj-%d.dat", Idx_Phi+1);
		Read_1D_Data(szName);
		
		if(Is_A_Single_Cluster())	{	// a single cluster
			fprintf(fOut, "%3d %3d %3d %3d  %.0lf\n", DihList[Idx_Phi][0], DihList[Idx_Phi][1], DihList[Idx_Phi][2], DihList[Idx_Phi][3], Centers[0]);
			continue;
		}
		
		Silhouette_Max=-1.0E100;
		for(int k=2; k<=4; k++)	{
			Silhouette = Do_Clustering(k);
			printf("k = %2d ,  %.3lf\n", k, Silhouette);
			if(Silhouette > (Silhouette_Max+0.01) )	{	// we prefer smaller cluster
				Silhouette_Max = Silhouette;
				n_Correct = k;
			}
		}
		
		printf("Phi %2d   n_Correct = %d\n\n", Idx_Phi+1, n_Correct);
		
		Silhouette = Do_Clustering(n_Correct);

		fprintf(fOut, "%3d %3d %3d %3d ", DihList[Idx_Phi][0], DihList[Idx_Phi][1], DihList[Idx_Phi][2], DihList[Idx_Phi][3]);

		SortCenters();
		for(i=0; i<n_Cluster; i++)	{
			fprintf(fOut, "  %.0lf", Centers[i]);
		}
		fprintf(fOut, "\n");


//		for(i=0; i<n_Cluster; i++)	{
//			Output_Data(i+1);
//		}
		
		
		if(pData!=NULL)	{
			delete []pData;
			delete []GroupID;
			delete []pSilhouette;
		}
	}
	fclose(fOut);


	return 0;
}

void SortCenters(void)
{
	int i, j;
	double fTmp;

	for(i=0; i<n_Cluster; i++)	{
		for(j=i+1; j<n_Cluster; j++)	{
			if(Centers[i] > Centers[j])	{
				fTmp = Centers[i];
				Centers[i] = Centers[j];
				Centers[j] = fTmp;
			}
		}
	}
}

void Setup_Centers(int k)
{
	int i;
	double d_Phi;

	d_Phi = 360.0/k;

	Centers[0] = -180.0 + d_Phi * rand(iseed);

	for(i=1; i<k; i++)	{
		Centers[i] = Centers[0] + i*d_Phi + d_Phi * 0.4 * (rand(iseed)-0.5);
	}
}

double Do_Clustering(int k)
{
	double Centers_Save[MAX_CLUSTER], wss_Sum, wss_Sum_Min=1.0E100, Avg_Silhouette;
	int i_Run;

	n_Cluster = k;
	for(i_Run=0; i_Run<MAX_RAND_INI; i_Run++)	{
		Setup_Centers(n_Cluster);
		kmeans_Lloyd_1D(pData, n_Point, Centers, n_Cluster, GroupID, 100, n_Point_List, withinss, wss_Sum);
		if(wss_Sum_Min > wss_Sum)	{
			wss_Sum_Min = wss_Sum;
			memcpy(Centers_Save, Centers, sizeof(double)*k);
		}
	}
	memcpy(Centers, Centers_Save, sizeof(double)*k);
	kmeans_Lloyd_1D(pData, n_Point, Centers, n_Cluster, GroupID, 100, n_Point_List, withinss, wss_Sum);

	Avg_Silhouette = Cal_Silhouette();

	return Avg_Silhouette;
}



int Read_1D_Data(char szName[])
{
	FILE *fIn;
	double Phi=0.0;
	int ReadItem, i;


	fIn = fopen(szName, "r");
	n_Point = 0;
	while(1)	{
		ReadItem = fscanf(fIn, "%lf", &Phi);
		if(ReadItem == 1)	{
			n_Point++;
		}
		else	{
			break;
		}
	}

	pData = new double[n_Point];
	GroupID = new int[n_Point];
	pSilhouette = new double[n_Point];

	fseek(fIn, 0, SEEK_SET);
	for(i=0; i<n_Point; i++)	{
		fscanf(fIn, "%lf", &(pData[i]));
		GroupID[i] = -1;
	}

	fclose(fIn);

	if(n_Point >= MAX_POINT)	{	// to control the number of data points
		Shuffle_Data(n_Point*8);
		n_Point = MAX_POINT;
	}
	

	return n_Point;
}

void Shuffle_Data(int nStep)
{
	int Step, i, j;
	double Phi_Tmp;

	for(Step=0; Step<nStep; Step++)	{
		i = (int)(n_Point*rand(iseed));
		j=i;
		while(j==i)	{
			j = (int)(n_Point*rand(iseed));
		}
		Phi_Tmp = pData[i];
		pData[i] = pData[j];
		pData[j] = Phi_Tmp;
	}
}

int kmeans_Lloyd_1D(double *x, int n, double *cen, int k, int *cl, int maxiter, int *nc, double *wss, double& wss_Sum)	// periodical data clustering
{
    int iter, i, j, it, inew = 0;
    double best, dd, tmp, cen_tmp, dx;
    int updated;

    for(i = 0; i < n; i++) cl[i] = -1;

    for(iter = 0; iter < maxiter; iter++) {
		updated = FALSE;
		for(i = 0; i < n; i++) {
			/* find nearest centre for each point */
			best = 1.0E100;
			for(j = 0; j < k; j++) {
				dd = 0.0;

				tmp = x[i] - cen[j];
				if(tmp > 180.0)	{
					tmp -= 360.0;
				}
				else if(tmp < -180.0)	{
					tmp += 360.0;
				}

				dd += tmp * tmp;

				if(dd < best) {
					best = dd;
					inew = j+1;
				}
			}
			if(cl[i] != inew) {
				updated = TRUE;
				cl[i] = inew;
			}
		}
		if(!updated) break;
		/* update each centre */
		for(j = 0; j < k; j++) cen[j] = 0.0;
		for(j = 0; j < k; j++) nc[j] = 0;
		for(i = 0; i < n; i++) {
			it = cl[i] - 1;
			nc[it]++;
			cen_tmp = cen[it]/nc[it];
			dx = cen_tmp - x[i];

			if(dx > 180.0)	{
				cen[it] += (x[i]+360.0);
			}
			else if(dx < -180.0)	{
				cen[it] += (x[i]-360.0);
			}
			else	{
				cen[it] += x[i];
			}
		}
		for(j = 0; j < k; j++)	{
			if(nc[j] > 0)	{
				cen[j] /= nc[j];
			}
			else	{
				cen[j] = 0.0;
			}

			if(cen[j] > 180.0)	{
				cen[j] -= 360.0;
			}
			else if(cen[j] < -180.0)	{
				cen[j] += 360.0;
			}
		}
    }
	
    for(j = 0; j < k; j++) wss[j] = 0.0;
    for(i = 0; i < n; i++) {
		it = cl[i] - 1;
		tmp = x[i] - cen[it];
		if(tmp > 180.0)	{
			tmp -= 360.0;
		}
		else if(tmp < -180.0)	{
			tmp += 360.0;
		}
		wss[it] += tmp * tmp;
    }
	wss_Sum = 0.0;
    for(j = 0; j < k; j++) wss_Sum +=wss[j];

	if(iter >= maxiter)	{
		return 0;
	}
	else	{
		return 1;
	}
}

/*

int kmeans_Lloyd_1D(double *x, int n, double *cen, int k, int *cl, int maxiter, int *nc, double *wss, double& wss_Sum)
{
    int iter, i, j, it, inew = 0;
    double best, dd, tmp;
    int updated;

    for(i = 0; i < n; i++) cl[i] = -1;

    for(iter = 0; iter < maxiter; iter++) {
		updated = FALSE;
		for(i = 0; i < n; i++) {
			//find nearest centre for each point
			best = 1.0E100;
			for(j = 0; j < k; j++) {
				dd = 0.0;

				tmp = x[i] - cen[j];
				dd += tmp * tmp;

				if(dd < best) {
					best = dd;
					inew = j+1;
				}
			}
			if(cl[i] != inew) {
				updated = TRUE;
				cl[i] = inew;
			}
		}
		if(!updated) break;
		//update each centre
		for(j = 0; j < k; j++) cen[j] = 0.0;
		for(j = 0; j < k; j++) nc[j] = 0;
		for(i = 0; i < n; i++) {
			it = cl[i] - 1; nc[it]++;
			cen[it] += x[i];
		}
		for(j = 0; j < k; j++) cen[j] /= nc[j % k];
    }
	
    for(j = 0; j < k; j++) wss[j] = 0.0;
    for(i = 0; i < n; i++) {
		it = cl[i] - 1;
		tmp = x[i] - cen[it];
		wss[it] += tmp * tmp;
    }
	wss_Sum = 0.0;
    for(j = 0; j < k; j++) wss_Sum +=wss[j];

	if(iter >= maxiter)	{
		return 0;
	}
	else	{
		return 1;
	}
}
*/


#define IADD   453806245
#define IMUL   314159269
#define MASK   2147483647
#define SCALE  0.4656612873e-9

double rand(int &iseed)
{
  iseed = (iseed * IMUL + IADD) & MASK;
  return (iseed * SCALE);
}

void Output_Data(int Idx)
{
	FILE *fOut;
	char szName[256];
	int i;

	sprintf(szName, "cluster-%d.dat", Idx);
	fOut = fopen(szName, "w");
	for(i=0; i<n_Point; i++)	{
		if(GroupID[i] == Idx)	{
			fprintf(fOut, "%.0lf %.3lf\n", pData[i], pSilhouette[i]);
//			fprintf(fOut, "%.0lf\n", pData[i]);
		}
	}
	fclose(fOut);
}

inline int To_Find_Closest_Cluster(int i)
{
	int ActiveCluster, j, Neighbor=-1;
	double dist, dist_Min=1.0E100, Phi;

	ActiveCluster = GroupID[i] - 1;
	Phi = pData[i];

	for(j=0; j<n_Cluster; j++)	{
		if(j == ActiveCluster)	{
			continue;
		}

		dist =  fabs(Phi - Centers[j]);
		if(dist > 180.0)	{
			dist = 360.0 - dist;
		}
		if(dist < dist_Min)	{
			dist_Min = dist;
			Neighbor = j;
		}
	}

	return (Neighbor+1);
}

inline double Cal_Silhouette_a_b(int i, int Idx)
{
	int Count=0, j;
	double Phi_0, dist, dist_Sum=0.0;

	Phi_0 = pData[i];

	for(j=0; j<n_Point; j++)	{
		if(GroupID[j] == Idx)	{
			dist = fabs(Phi_0 - pData[j]);
			if(dist > 180.0)	{
				dist = 360.0 - dist;
			}
			dist_Sum += dist;
			Count++;
		}
	}
	return (dist_Sum/Count);
}

double Cal_Silhouette(void)
{
	int i, ActiveCluster, Neighbor_Cluster;
	double a, b, max_a_b, Mean;

	memset(pSilhouette, 0, sizeof(double)*n_Point);
	for(ActiveCluster=1; ActiveCluster<=n_Cluster; ActiveCluster++)	{
		for(i=0; i<n_Point; i++)	{
			if(GroupID[i] == ActiveCluster)	{	//
				a = Cal_Silhouette_a_b(i, ActiveCluster);
				Neighbor_Cluster = To_Find_Closest_Cluster(i);
				b = Cal_Silhouette_a_b(i, Neighbor_Cluster);
				max_a_b = max(a, b);
				pSilhouette[i] = (b-a)/max_a_b;
			}
		}
	}

	Mean = 0.0;
	for(i=0; i<n_Point; i++)	{
		Mean += pSilhouette[i];
	}
	return (Mean/n_Point);
}

int Is_A_Single_Cluster(void)
{
	double cen=0.0, cen_Local=0.0, dx, dist_SQ_Sum, sigma;
	int i;

	for(i = 0; i < n_Point; i++) {
		cen_Local = cen/(i+1);
		dx = cen_Local - pData[i];
		
		if(dx > 180.0)	{
			cen += (pData[i]+360.0);
		}
		else if(dx < -180.0)	{
			cen += (pData[i]-360.0);
		}
		else	{
			cen += pData[i];
		}
	}
	cen /= n_Point;
	
	if(cen > 180.0)	{
		cen -= 360.0;
	}
	else if(cen < -180.0)	{
		cen += 360.0;
	}
	
	dist_SQ_Sum = 0.0;
	for(i = 0; i < n_Point; i++) {
		dx = fabs(cen - pData[i]);
		if(dx > 180.0)	{
			dx -= 360.0;
		}
		dist_SQ_Sum += (dx*dx);
	}

	sigma = sqrt(dist_SQ_Sum/n_Point);

	if(sigma < 30.0)	{	// the criteria for a single cluster
		n_Cluster = 1;
		Centers[0] = cen;
		return 1;
	}
	else	{
		return 0;
	}
}


int Read_Soft_DihedralList(void)
{
	FILE *fIn;
	int ReadItem;

	n_Phi = 0;

	fIn = fopen(szPhiToScan, "r");

	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open file in Read_Soft_DihedralList().\n");
	}

	while(1)	{
		ReadItem = fscanf(fIn, "%d %d %d %d", &(DihList[n_Phi][0]), &(DihList[n_Phi][1]), &(DihList[n_Phi][2]), &(DihList[n_Phi][3]));
		if(ReadItem == 4)	{
			n_Phi++;
		}
		else	{
			break;
		}
	}

	fclose(fIn);

	return n_Phi;
}

void Quit_With_Error_Msg(char szMsg[])
{
	FILE *fOut;
	fOut = fopen("error.txt", "a+");
	fseek(fOut, 0, SEEK_END);
	fprintf(fOut, "Error in clustering.cpp\n");
	fprintf(fOut, "%s\n", szMsg);
	fclose(fOut);

	exit(1);
}

#define E_MAX	(150.0)
void Split_Phi_Traj(void)
{
	FILE *fIn, *fOut[N_MAX_DIH];
	double Phi_List[N_MAX_DIH], fTmp, Emin=1.0E200, E_Mol;
	char szLine[256], *ReadLine, szName[256];
	int i, ReadItem, n_Conf=0, n_Conf_Used=0;

	fIn = fopen("E-phi-mm-pes.txt", "r");
	if(fIn == NULL)	{
		Quit_With_Error_Msg("Fail to open file: E-phi-mm-pes.txt in Split_Phi_Traj()\n");
	}

	for(i=0; i<n_Phi; i++)	{
		sprintf(szName, "phi-traj-%d.dat", i+1);
		fOut[i] = fopen(szName, "w");
	}

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		if(strlen(ReadLine) < 22)	{
			break;
		}
		ReadItem = sscanf(szLine, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
			&fTmp, &E_Mol, &(Phi_List[0]), &(Phi_List[1]), &(Phi_List[2]), &(Phi_List[3]), &(Phi_List[4]), 
			&(Phi_List[5]), &(Phi_List[6]), &(Phi_List[7]), &(Phi_List[8]), &(Phi_List[9]), &(Phi_List[10]), &(Phi_List[11]), 
			&(Phi_List[12]), &(Phi_List[13]), &(Phi_List[14]), &(Phi_List[15]), &(Phi_List[16]), &(Phi_List[17]), &(Phi_List[18]), &(Phi_List[19]), 
			&(Phi_List[20]), &(Phi_List[21]), &(Phi_List[22]), &(Phi_List[23]), &(Phi_List[24]), &(Phi_List[25]), &(Phi_List[26]), &(Phi_List[27]), 
			&(Phi_List[28]), &(Phi_List[29]), &(Phi_List[30]), &(Phi_List[31]));
		if(ReadItem == (n_Phi+2))	{
			if(E_Mol < Emin)	{
				Emin = E_Mol;
			}
			n_Conf++;
		}
		else	{
			printf("Warning> Reading: %s\nin E-phi-mm-pes.txt.\n", szLine);
		}
	}

	printf("Read data for %d configurations.\n", n_Conf);
	fseek(fIn, 0, SEEK_SET);


	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		if(strlen(ReadLine) < 22)	{
			break;
		}
		ReadItem = sscanf(szLine, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
			&fTmp, &E_Mol, &(Phi_List[0]), &(Phi_List[1]), &(Phi_List[2]), &(Phi_List[3]), &(Phi_List[4]), 
			&(Phi_List[5]), &(Phi_List[6]), &(Phi_List[7]), &(Phi_List[8]), &(Phi_List[9]), &(Phi_List[10]), &(Phi_List[11]), 
			&(Phi_List[12]), &(Phi_List[13]), &(Phi_List[14]), &(Phi_List[15]), &(Phi_List[16]), &(Phi_List[17]), &(Phi_List[18]), &(Phi_List[19]), 
			&(Phi_List[20]), &(Phi_List[21]), &(Phi_List[22]), &(Phi_List[23]), &(Phi_List[24]), &(Phi_List[25]), &(Phi_List[26]), &(Phi_List[27]), 
			&(Phi_List[28]), &(Phi_List[29]), &(Phi_List[30]), &(Phi_List[31]));
		if(ReadItem == (n_Phi+2))	{
			if( (E_Mol-Emin) < (E_MAX) )	{
				for(i=0; i<n_Phi; i++)	{
					fprintf(fOut[i], "%.0lf\n", Phi_List[i]);
				}
				n_Conf_Used++;
			}
		}
		else	{
			printf("Warning> Reading: %s\nin E-phi-mm-pes.txt.\n", szLine);
		}
	}
	printf("Using data for %d configurations.\n", n_Conf_Used);

	for(i=0; i<n_Phi; i++)	{
		fclose(fOut[i]);
	}

	fclose(fIn);
}
#undef E_MAX

