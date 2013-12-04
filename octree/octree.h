// version 1.0 11/11/06
// Author: Thierry Sousbie.
// For any information contact tsousbie@obs.univ-lyon1.fr

#ifndef __NEI_FIND__
#define __NEI_FIND__

#include "oct_struct.h"

#define FLAG_MIN(x) (1<<(2*x))
#define FLAG_MAX(x) (1<<(2*x+1))

#define FLAG_CROSS(x) (1<<x)

int freeINT(int *p);
int FreeGroupList(group_list *gl);
int FreeGroupList_Y(group_list *gl, long *index, int N);
int FreeOctData(OcTree_data *data);
int FreeOctDataKeepId(OcTree_data *data);

OcTree_data *BuildNeiTree(float *Pos,int n,int inplace,float *min,float *max,int ndims,int nmaxincell);
int *BuildFoF (OcTree_data *data, double p_dist, int **fof_id,int periodic);

void *ComputeGroups_Y(OcTree_data *data,int *Id,float *P,float *V,int N,int Nmin,float ll,snapshot_data *snap,float *delta,int ndims_p,int cp);
group_list *ComputeGroupsAttributes(int *Id,int Nmin,float ll, snapshot_data *snap,int cp);
group_list *ComputeGroupsFromTree(OcTree_data *data, float *P, float *V, int Nmin,int cp);
group_list *ComputeGroupsFromArrays(int *Id,int Nmin,float ll, float *P,float *V,int N,float *delta,int ndims,int cp);

int GetIndexClosest(OcTree_data *data,int index,double *dist,int periodic);
int GetPosClosest(OcTree_data *data,float *p_RefPos,double *dist,int periodic);

int GetIndexNeighbours(OcTree_data *,int,int,Nei_List **,int);
int GetPosNeighbours(OcTree_data *,float *,int,Nei_List **,int);

int GetIndexNeighbours_Y(OcTree_data *data,int index,int N,long *id,double *d,int periodic);
int GetPosNeighbours_Y(OcTree_data *data,float *RefPos,int N,long *id,double *d,int periodic);

void GetNeighbours_Y(OcTree_data *octdata,void *data,int ndata, int dataisint,int N,long *id,double *d,int periodic);

int GetPointsInSphere(OcTree_data *data,float *Center,float radius,int **id,int periodic);

#endif

