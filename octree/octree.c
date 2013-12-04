// version 1.0 11/11/06
// Author: Thierry Sousbie.
// For any information contact tsousbie@obs.univ-lyon1.fr

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "octree.h"

//int TREE_MAXINCELL;
#define MAX_TREE_LEVEL 50

#ifndef PI
#define PI 3.1415926535897932384
#endif

#define SQRT3INV 0.57735027

#define DSQUARE(_vect,_index,_ndims,_result) \
	for (_index=0,_result=0;_index<_ndims;_index++) \
_result += (_vect)[_index]*(_vect)[_index ]

int freeINT(int *p)
{
	//printf ("pointer : %d\n",(int)p);
	free (p);
	return 1;
}

int FreeGroupList(group_list *gl)
{
	int i;

	for (i=0;i<gl->NGroups;i++)
	{
		free(gl->group[i].index);
		free(gl->group[i].VMean);
		free(gl->group[i].Center);
		free(gl->group[i].Spin);
	}

	free(gl->group);

	return 1;
}

int FreeGroupList_Y(group_list *gl, long *index, int N)
{
	int i;

	for (i=0;i<N;i++)
	{
		free((int *) index[i]);
	}

	free(gl->group);

	return 1;
}


int FreeOctData(OcTree_data *data)
{
	free(data->Pos);
	free(data->Id);
	free(data->Inv_Id);
	free(data->Tree_ind);
	free(data->limits);
	free(data->Tree);
	free(data->Fof_Id);
	free(data->min);
	free(data->max);

	return 0;
}

int FreeOctDataKeepId(OcTree_data *data)
{
	free(data->Pos);
	free(data->Id);
	free(data->Inv_Id);
	free(data->Tree_ind);
	free(data->limits);
	free(data->Tree);
	free(data->min);
	free(data->max);

	return 0;
}


int Diagonalise3x3_oct(double *M,double *L,double *vec)
{
	double a,b,c,u,v,w;
	double B,C,D;
	double p,q,r;
	double phi;
	double tmp;
	int i;

	a=M[0];u=M[1];v=M[2];
	b=M[3];c=M[5];w=M[4];

	//if ((a==b)&&(b==c)&&(c==0)) return -1;

	//Compute characteristic polynome
	B=-(a+b+c);
	C=a*b+a*c+b*c-u*u-v*v-w*w;
	D=u*u*c+v*v*b+w*w*a-a*b*c-2*u*v*w;

	//solve 3rd degree equation
	p=(3.*C-B*B)/9.;
	q=(2.*B*B*B/27.-B*C/3.+D)/2.;

	if (q*q+p*p*p>=0)
	{
		L[0]=L[1]=L[2]=0;
		return -1;
	}


	r=((q<0)?-1:1)*sqrt(fabs(p));
	phi = acos(q/(r*r*r))/3.;

	r*=2;//B/=3;
	L[0]=-r*cos(phi)-B/3.;
	L[1]=r*cos(PI/3-phi)-B/3.;
	L[2]=r*cos(PI/3+phi)-B/3.;



	if ((L[2])>(L[1])) {tmp=L[1];L[1]=L[2];L[2]=tmp;}
	if ((L[1])>(L[0])) {tmp=L[1];L[1]=L[0];L[0]=tmp;}
	if ((L[2])>(L[1])) {tmp=L[1];L[1]=L[2];L[2]=tmp;}

	//Some special case ...

	//matrix is already diagonal
	if ((u==0)&&(v==0)&&(w==0))
	{
		vec[0]=0;vec[1]=0;vec[2]=0;
		//associate each axis (x,y and z) with it s eigenvalues
		for (i=0;i<3;i++)
		{
			if (fabs(L[i])>1.E-8)
			{
				vec[3*i+0]=0;vec[3*i+1]=0;vec[3*i+2]=0;
				if (fabs((a-L[i])/L[i])<1.E-4) vec[3*i+0]=1;
				else if (fabs((b-L[i])/L[i])<1.E-4) vec[3*i+1]=1;
				else if (fabs((c-L[i])/L[i])<1.E-4) vec[3*i+2]=1;
				else printf ("OOPSSSSSS\n");
			}
			else {L[i]=0;vec[3*i+0]=0;vec[3*i+1]=0;vec[3*i+2]=0;}
		}
	}

	//We already have 1 eigenvalue (matrix is block diagonal)
	else if (((u==0)&&(v==0))||((v==0)&&(w==0))||((u==0)&&(w==0)))
	{
		double tmpa=0,tmpb=0,tmpu=0;
		int tmpd=0;
		//Finds which axis are already eigenvector (x,y or z)
		vec[0]=0;vec[1]=0;vec[2]=0;
		for (i=0;i<3;i++)
		{
			if (fabs(L[i])>1.E-8)
			{
				vec[3*i+0]=0;vec[3*i+1]=0;vec[3*i+2]=0;
				//Check if we already have this eigenvector
				if (fabs((a-L[i])/L[i])<1.E-4) vec[3*i+0]=1;
				else if (fabs((b-L[i])/L[i])<1.E-4) vec[3*i+1]=1;
				else if (fabs((c-L[i])/L[i])<1.E-4) vec[3*i+2]=1;
				//no ...
				//diagonalise 2x2 block ...
				else
				{
					if ((u==0)&&(v==0)) {tmpa=b;tmpb=c;tmpu=w;tmpd=1;}
					if ((v==0)&&(w==0)) {tmpa=a;tmpb=b;tmpu=u;tmpd=0;}
					if ((u==0)&&(w==0)) {tmpa=a;tmpb=c;tmpu=v;tmpd=2;}
					if ((tmpa==tmpb)&&(tmpb==0))
						B=C=0;
					else
					{
						B=(L[i]+tmpu-tmpa)/(L[i]+tmpu-tmpb);
						C=1./sqrt(1+B*B);
					}
					if (tmpd==2)
					{
						vec[3*i]=C;
						vec[3*i+2]=B*C;
					}
					else
					{
						vec[3*i+tmpd]=C;
						vec[3*i+tmpd+1]=B*C;
					}
				}
			}
			else {L[i]=0;vec[3*i+0]=0;vec[3*i+1]=0;vec[3*i+2]=0;}
		}
	}

	//General case ...
	else for (i=0;i<3;i++)
	{

		B=L[i]-b-u*u/(L[i]-a);
		B=(w+u*v/(L[i]-a))/B;

		C=u/(L[i]-a)*B+v/(L[i]-a);
		D=1./sqrt(B*B+C*C+1.);

		vec[3*i+2]=D;
		vec[3*i+1]=B*D;
		vec[3*i]=C*D;
	}

	//Now rearrange vectors to have a direct eigenbasis
	/*
	   B=vec[0]*(vec[3+1]*vec[6+2]-vec[3+2]*vec[6+1])-
	   vec[3]*(vec[0+1]*vec[6+2]-vec[0+2]*vec[6+1])+
	   vec[6]*(vec[0+1]*vec[3+2]-vec[0+2]*vec[3+1]);

	   if (B<0) {vec[0]=-vec[0];vec[1]=-vec[1];vec[2]=-vec[2];}
	   */


	return 0;
}


//dim is which coordinate to look at
//tab should be (kind of) sorted in ascending order
//return the index of the first value >=val
int Get_geq_index(float *pos,int n,int dim,float val,int ndims)
{
	long imin,imax;
	float Pos;

	imin=0;
	imax=n-1;

	if (n<=0) return 0;

	while (abs(imax-imin)>1)
	{
		Pos = pos[(size_t)ndims*(int)((imin+imax)/2)+dim];
		if (Pos>=val)
			imax=(imin+imax)/2;
		else
			imin=(imin+imax)/2;
	}

	Pos = pos[(size_t)ndims*imin+dim];

	if (pos[(size_t)ndims*imax+dim]<val) return imax+1;

	if (pos[(size_t)ndims*imin+dim] >= val)
		return imin;
	else
		return imin+1;
}

//Sorts Id and Pos so that all values <val are on left and returns
//sort is followind dim componant of pos (dim=0->x, 1->y and 2->z)
int PivotSort(int *Id,float *Pos,int n,int dim,float val,int ndims)
{
	int i=0;
	int j=n-1;
	int k;
	float tmp[ndims];

	if (n<=0) return 0;

	for (;;)
	{

		//while ((Pos[(size_t)ndims*i+dim]<val)&&(i<n)) i++;
		//while ((Pos[(size_t)ndims*j+dim]>=val)&&(j>=0)) j--;
		while ((i<n)&&(Pos[(size_t)ndims*i+dim]<val)) i++;
		while ((j>=0)&&(Pos[(size_t)ndims*j+dim]>=val)) j--;
		//printf ("%d,%d (dim=%d) ",i,j,dim);

		if (i>j)
			return i;
		else
		{
			k=Id[i];
			Id[i]=Id[j];
			Id[j]=k;
			memcpy(tmp,&Pos[(size_t)ndims*i],(size_t)ndims*sizeof(float));
			memcpy(&Pos[(size_t)ndims*i],&Pos[(size_t)ndims*j],(size_t)ndims*sizeof(float));
			memcpy(&Pos[(size_t)ndims*j],tmp,(size_t)ndims*sizeof(float));
			i++;j--;
		}
	}

}

// Sorts data depending on Position
// This is a nested sort so that data in two nodes close to each other are close to each other in Tree
// prev =-1 and PosIndex=TreeIndex=0 to initialize
int SortTree(int *Id,float *Pos,int n,int level,OcTree_data *data,int TreeIndex, int PosIndex,int prev,
		float *min,float *max,int ndims)
{
	long index[1<<ndims];
	/* Passer en static ... mais attention a l alloc !!! */
	int sortindex[1<<ndims];
	float tmpmin[ndims];
	float tmpmax[ndims];
	/**************************/
	int nextindex[1<<ndims];
	int n_nodes=0;
	double delta[ndims];

	long i,j,k,k0,l;

	static int Cur_Tree_index;

	for (i=0;i<ndims;i++)
	{
		delta[i] = ((double)max[i]-(double)min[i])/2;
		//printf ("delta[%d]=%f\n",i,delta[i]);
	}

	//This branch is empty ...
	if (n==0)
	{
		Cur_Tree_index--;
		return 0;
	}

	if (level<0)
	{
		//Tree has been allocated, fill it up
		Nei_OcTree *T = data->Tree;

		level = -2;
		T[TreeIndex].Prev = &T[prev];
		T[TreeIndex].Cur = TreeIndex;
		T[TreeIndex].N = n;
		T[TreeIndex].index = PosIndex;
		T[TreeIndex].NNext = 0;

		for (i=0;i<ndims;i++)
		{
			data->limits[2*ndims*TreeIndex+2*i]=min[i];
			data->limits[2*ndims*TreeIndex+2*i+1]=max[i];
		}
	}

	//this branch does not need refinement
	if (n<=data->maxincell) return 1;

	if (level>MAX_TREE_LEVEL)
	{
		printf ("Warning: Maximum level (%d) reached.\n",MAX_TREE_LEVEL);
		return 1;
	}

	//index = malloc (sizeof(long) * (1<<ndims));
	index[0]=0;

	sortindex[0]=0;
	for (j=1;j<=ndims;j++)
		for (i=0;i<(1<<(j-1));i++)
			sortindex[(1<<(ndims-j))*(2*i+1)]=(1<<(j-1))+i;

	for (i=0;i<(1<<ndims);i++)
		nextindex[sortindex[i]]=i;

	if (level >= 0)
	{
		int u;

		//Sort for the (1<<ndims) sub cells
		for (i=0,k=1;i<ndims;i++)
		{
			k0=k;
			for (j=0,l=0;j<(1<<i);j++,k++,l++)
			{
				if (l+1==k0)
				{
					index[k] = index[l] + PivotSort(&Id[index[l]],&Pos[(size_t)ndims*index[l]],n-index[l],i,min[i]+delta[i],ndims);
					//printf ("%d-%d\n",n,index[l]);
				}
				else
				{
					u=sortindex[l*(1<<(ndims-i))];

					//printf ("%d-%d(%d-%d)\n",index[sortindex[nextindex[u]+(1<<(ndims-i))]],index[u],sortindex[nextindex[u]+(1<<(ndims-i))],u);
					index[k] = index[u] + PivotSort(&Id[index[u]],&Pos[(size_t)ndims*index[u]],index[sortindex[nextindex[u]+(1<<(ndims-i))]]-index[u],i,min[i]+delta[i],ndims);

				}
			}
		}
	}
	else
	{
		int u;

		for (i=0,k=1;i<ndims;i++)
		{
			k0=k;
			for (j=0,l=0;j<(1<<i);j++,k++,l++)
			{
				u=sortindex[l*(1<<(ndims-i))];
				if (l+1==k0)
					index[k] = index[u] +  Get_geq_index(&Pos[(size_t)ndims*index[u]],n-index[u],i,min[i]+delta[i],ndims);
				else
					index[k] = index[u] +  Get_geq_index(&Pos[(size_t)ndims*index[u]],index[sortindex[nextindex[u]+(1<<(ndims-i))]]-index[u],i,min[i]+delta[i],ndims);
			}
		}
	}

	for (i=0;i<(1<<ndims);i++)
	{
		if (nextindex[i]+1 == (1<<ndims))
			nextindex[i] = -1;
		else
			nextindex[i] = sortindex[nextindex[i]+1];
	}
	//for (k=0;k<(1<<ndims);k++) printf ("%d ",index[k]);
	//printf ("\n");

	for (k=0;k<(1<<ndims);k++)
	{
		for (l=0;l<ndims;l++)
		{
			if (k&(1<<(ndims-l-1)))
			{
				tmpmin[l]=min[l]+ delta[l];
				tmpmax[l]=max[l];
				//printf ("%d(%f): [1,2],[%f,%f] ",l,delta[l], tmpmin[l],tmpmax[l]);
			}
			else
			{
				tmpmin[l]=min[l];
				tmpmax[l]=min[l]+ delta[l];
				//printf ("%d(%f): [0,1],[%f,%f] ",l,delta[l], tmpmin[l],tmpmax[l]);
			}
		}

		i = sortindex[k];

		Cur_Tree_index++;
		if (nextindex[i]<0)
		{
			//printf ("%d-%d\n",n,index[i]);
			n_nodes+=
				SortTree(&Id[index[i]],&Pos[(size_t)ndims*index[i]],
						n-index[i],level+1,data,Cur_Tree_index,PosIndex+index[i],
						TreeIndex,tmpmin,tmpmax,ndims);
		}
		else
		{
			//printf ("%d-%d\n",index[nextindex[i]],index[i]);
			n_nodes+=
				SortTree(&Id[index[i]],&Pos[(size_t)ndims*index[i]],
						index[nextindex[i]]-index[i],level+1,data,Cur_Tree_index,
						PosIndex+index[i],TreeIndex,tmpmin,tmpmax,ndims);
		}

	}

	//for (i=0;i<ndims;i++) free(index[i]);

	if (level == 0)
	{
		Nei_OcTree *T;
		Nei_OcTree *Tprev;

		printf ("Filling up...");fflush(0);


		//Now we can fill the tree up
		data->NLeaves = n_nodes+1;
		T = data->Tree = (Nei_OcTree *) calloc (data->NLeaves,sizeof(Nei_OcTree));
		if (T==NULL)
		{
			fprintf (stderr,"Not enough memory to allocate data->Tree.\n");
			exit(0);
		}

		data->limits = (float *) malloc ((long)2*data->NLeaves*ndims*sizeof(float));
		if (data->limits==NULL)
		{
			fprintf (stderr,"Not enough memory to allocate data->limits.\n");
			exit(0);
		}

		Cur_Tree_index=0;

		//do the fill up
		SortTree(Id,Pos,n,-1,data,0,0,-1,min,max,ndims);

		//associate the Next with the right node
		//printf ("Associating next\n");
		for (i=1;i<data->NLeaves;i++)
		{
			Tprev = T[i].Prev;
			if ((Tprev->Next==NULL)||((Tprev->NNext>>3)!=((Tprev->NNext-1)>>3)))
			{
				//printf ("%d-%d\n",Tprev->NNext,(1<<3));
				Tprev->Next = realloc(Tprev->Next,(Tprev->NNext+(1<<3))*sizeof(Nei_OcTree*));
			}

			Tprev->Next[Tprev->NNext]=&T[i];
			T[i].Prev_NextIndex = Tprev->NNext;
			Tprev->NNext++;
		}

		printf ("local sort...");fflush(0);

		//And finally sort the pos in deepest nodes following R (for TINN method)
		for (i=1;i<data->NLeaves;i++)
		{
			//Check if it is a deepest node
			if (T[i].NNext==0)
			{
				int j,k,l;
				int tmpint;
				float min;
				float tmp[ndims];

				//And dirty sort it (should be fine for at most TREE_MAXINCELL points)
				for (j=T[i].index;j<T[i].index+T[i].N;j++)
				{
					min = data->Pos[(size_t)ndims*j];l=j;
					for (k=j;k<T[i].index+T[i].N;k++)
					{
						if (data->Pos[(size_t)ndims*k]<min) {min=data->Pos[(size_t)ndims*k];l=k;}
					}
					if (l!=j)
					{
						tmpint=data->Id[l];data->Id[l]=data->Id[j];data->Id[j]=tmpint;
						memcpy(tmp,&data->Pos[(size_t)ndims*l],(size_t)ndims*sizeof(float));
						memcpy(&data->Pos[(size_t)ndims*l],&data->Pos[(size_t)ndims*j],(size_t)ndims*sizeof(float));
						memcpy(&data->Pos[(size_t)ndims*j],tmp,(size_t)ndims*sizeof(float));
					}
				}
			}
		}

		//Compute the Inv_ID
		data->Inv_Id = (int *) malloc (data->Npart*sizeof(int));
		if (data->Inv_Id==NULL)
		{
			fprintf (stderr,"Not enough memory to allocate data->Inv_Id.\n");
			exit(0);
		}

		for (i=0;i<data->Npart;i++)
		{
			data->Inv_Id[data->Id[i]]=i;
		}

		//Compute the Tree_ind
		data->Tree_ind = (int *) malloc (data->Npart*sizeof(int));
		if (data->Tree_ind==NULL)
		{
			fprintf (stderr,"Not enough memory to allocate data->Tree_ind.\n");
			exit(0);
		}

		for (i=0;i<data->NLeaves;i++)
		{
			int j;
			if (data->Tree[i].NNext==0)
			{
				for (j=0;j<data->Tree[i].N;j++)
					data->Tree_ind[data->Tree[i].index+j]=i;
			}
		}


	}

	return n_nodes+1;
}

// Build the Octree with n particles at positiuon  Pos[x1,y1,z1,x2,y2,...],
// with particle identity Id.
// if inplace is true, data in Pos is used and modified but no mem is allocated for it
OcTree_data *BuildNeiTree(float *Pos,int n,int inplace,float *min,float *max,int ndims,int nmaxincell)
{
	//Nei_OcTree *Tree;
	OcTree_data *data;
	int *Id;
	long i,j;
	clock_t time1;
	clock_t time2;

	time1=clock();
	printf ("Building %dD tree for %d particles ...",ndims,n);fflush(0);

	Id = malloc (n*sizeof(int));

	//Fill up data structure
	for (i=0;i<n;i++)
		Id[i] = i;

	data = calloc (1,sizeof(OcTree_data));


	if (inplace)
		data->Pos = Pos;
	else
	{
		data->Pos = (float *) malloc ((size_t)ndims*n*sizeof(float));
		memcpy(data->Pos,Pos,(size_t)ndims*n*sizeof(float));
	}

	data->Id = Id;
	data->Npart = n;
	data->maxincell = nmaxincell;
	data->ndims=ndims;

	data->min = malloc(ndims*sizeof(float));
	data->max = malloc(ndims*sizeof(float));

	memcpy(data->min,min,ndims*sizeof(float));
	memcpy(data->max,max,ndims*sizeof(float));

	printf ("Sorting...");fflush(0);

	SortTree(data->Id,data->Pos,n,0,data,0,0,-1,min,max,ndims);

	for (i=0,j=0;i<data->NLeaves;i++)
	{
		if (data->Tree[i].NNext!=0)
			j+=(1<<3)+((data->Tree[i].NNext-1)&(~((1<<3)-1)));
	}

	//printf ("En moyenne: %f\n",(double)j/(double)data->NLeaves);

	i=(long)data->NLeaves*((long)sizeof(Nei_OcTree) + (long)sizeof(Nei_OcTree)*(double)j/(double)data->NLeaves);
	i+= n*ndims*sizeof(float) + n*2*sizeof(int) + 2*ndims*sizeof(float)*data->NLeaves;

	time2=clock();
	printf (" done. (%d nodes (~%.2f Mo) in %.2fs)\n",data->NLeaves,((double)(i))/(1024.*1024.),(double)(time2-time1)/CLOCKS_PER_SEC);

	//printf ("long(data)=%d (%d)\n",(long)data,sizeof(long));
	return data;
}

//periodic is true or false
int RecursiveSetFofPart(OcTree_data *data,Nei_OcTree *Tree,float *RefPos,int Prev_NextIndex,int Id,double dist,int *NInNode,int ndims, int periodic)
{

	float min[ndims];
	long i,j,k;
	long index;
	double d2;
	float x[ndims];
	double u[ndims];
	double a;
	Nei_OcTree *TmpTree;

	int NInGroup = 0;
	int per_flags=0;

	memcpy (x,RefPos,ndims*sizeof(float));
	d2=dist*dist;

	//Checks if all the particles in that tree branch are not already in the fof
	if (Tree->N>NInNode[Tree->Cur])
	{
		if (Tree->NNext == 0)
		{
			index = Get_geq_index(&data->Pos[(size_t)ndims*Tree->index],Tree->N,0,x[0]-dist,ndims);
			j=index+Tree->index;

			while ((j<Tree->index+Tree->N)&&(data->Pos[(size_t)ndims*j] <= x[0]+dist))
			{

				for (k=0;k<ndims;k++)
					u[k] = data->Pos[(size_t)ndims*j+k] - x[k];
				DSQUARE(u,k,ndims,a);

				//Can set the particle group
				if ((a<=d2)&&(data->Fof_Id[j] == 0))
				{

					//Add the new point
					NInGroup++;
					data->Fof_Id[j]=Id;

					TmpTree=Tree;
					while (TmpTree->Cur!=0)
					{
						NInNode[TmpTree->Cur]++;
						TmpTree=TmpTree->Prev;
					}

					NInGroup += RecursiveSetFofPart(data,Tree,&data->Pos[(size_t)ndims*j],0,Id,dist,NInNode,ndims,periodic);


					/* "periodic boundary condition" */

					//Check periodic boundary conditions ...
					if (periodic)
					{
						per_flags=0;

						for (k=0;k<ndims;k++)
						{
							if (data->Pos[(size_t)ndims*j+k]-dist<data->min[k]) per_flags|=(1<<(2*k));
							else if (data->Pos[(size_t)ndims*j+k]+dist>data->max[k]) per_flags|=(1<<(2*k+1));
						}

						//sphere crosses a boundary
						if (per_flags)
						{
							int per_index;
							int inter;
							float New_RefPos[ndims];

							for (per_index=1;per_index<=per_flags;per_index++)
							{
								inter = per_index&per_flags;

								if ( (per_index&(~per_flags)) != 0) continue;

								if (inter)
								{
									memcpy(New_RefPos,&data->Pos[(size_t)ndims*j],(size_t)ndims*sizeof(float));

									for (k=0;k<ndims;k++)
									{
										if (inter&(1<<(2*k))) New_RefPos[k]+=data->max[k]-data->min[k];
										else if (inter&(1<<(2*k+1))) New_RefPos[k]-=data->max[k]-data->min[k];
									}

									NInGroup+=RecursiveSetFofPart(data,data->Tree,New_RefPos,-1,Id,dist,NInNode,ndims,periodic);
								}
							}
						}
					}

					/* end of "periodic boundary condition" */
				}
				j++;
			}

		}
		else
		{
			double dmin2[1<<ndims];

			//Compute the min distances for the next ...
			for (i=0;i<Tree->NNext;i++)
			{
				if (i!=Prev_NextIndex)
				{

					index = 2*ndims*Tree->Next[i]->Cur;

					for (k=0;k<ndims;k++)
						if (x[k]<data->limits[index+2*k])
							min[k] = data->limits[index+2*k] - x[k];
						else if (x[k]<data->limits[index+2*k+1])
							min[k] = 0;
						else
							min[k] = x[k] - data->limits[index+2*k+1];

					DSQUARE(min,k,ndims,dmin2[i]);
				}
				else
				{
					dmin2[i] = 1.E20;
				}
			}

			for (i=0;i<Tree->NNext;i++)
			{
				if (dmin2[i]<d2)
					NInGroup += RecursiveSetFofPart(data,Tree->Next[i],RefPos,-1,Id,dist,NInNode,ndims,periodic);
			}
		}
	}


	if ((Prev_NextIndex >= 0)&&(Tree->Cur!=0))
	{
		index = ndims*Tree->Cur;

		for (i=0,k=0;i<2*ndims;i++)
			if (fabs(x[i]-data->limits[index+i])>dist)
				return (NInGroup + RecursiveSetFofPart(data,Tree->Prev,RefPos,Tree->Prev_NextIndex,Id,dist,NInNode,ndims,periodic));
	}

	return NInGroup;
}

int *BuildFoF (OcTree_data *data, double p_dist, int **fof_id,int periodic)
{
	long i,n;
	long N;
	int Cur_Id=1;
	int NGroups=0;
	double dist = p_dist;
	Nei_OcTree *Tree;
	Nei_OcTree *TmpTree;
	int *NInNode;
	clock_t time1;
	clock_t time2;
	int *tmp_id;
	int ndims = data->ndims;

	printf ("Building FOF with l=%.2e ...",p_dist);fflush(0);

	time1 = clock();

	data->Fof_ll = dist;

	if (fof_id!=NULL)
	{
		data->Fof_Id = *fof_id;
	}

	data->Fof_Id = (int *) realloc (data->Fof_Id,(size_t)data->Npart*sizeof(int));
	memset (data->Fof_Id,0,(size_t)data->Npart*sizeof(int));

	NInNode = calloc(data->NLeaves,sizeof(int));
	N=0;

	for (i=0;i<data->Npart;i++)
	{
		if (data->Fof_Id[i]==0)
		{
			n=1;
			data->Fof_Id[i]=Cur_Id;
			Cur_Id++;

			TmpTree = Tree = &data->Tree[data->Tree_ind[i]];
			while (TmpTree->Cur!=0)
			{
				NInNode[TmpTree->Cur]++;
				TmpTree=TmpTree->Prev;
			}


			n+=RecursiveSetFofPart(data,Tree,&data->Pos[(size_t)ndims*i],0,data->Fof_Id[i],dist,NInNode,ndims,periodic);

			if (n>=20) N++;

			NGroups ++;
		}
	}

	time2 = clock();
	printf ("\rBuilding FOF with l=%.2e ... %d groups found (>=20) in %.5fs.\n",p_dist,(int)N,(double)(time2-time1)/CLOCKS_PER_SEC);

	free(NInNode);

	tmp_id = (int *) malloc((size_t)data->Npart*sizeof(int));
	memcpy(tmp_id,data->Fof_Id,(size_t)data->Npart*sizeof(int));

	for (i=0;i<data->Npart;i++)
	{
		data->Fof_Id[data->Id[i]]=tmp_id[i];
	}

	free(tmp_id);

	if (fof_id!=NULL)
		*fof_id = data->Fof_Id;

	return data->Fof_Id;

}

/*
   snap_p is expected to be of type (snapshot_data *)
   it is casted to void for yorick wrapper compatibility
   */
/*
   group_list *ComputeGroupsAttributes_Y(int *Id,int Nmin,float ll, void *snap_p,int cp)
   {
   return ComputeGroupsAttributes(Id,Nmin,ll, (snapshot_data *) snap_p,cp);
   }
   */

#define ALLOCANDRETGRP(__gr,__grp,__grpl,__index,__dim)\
{\
	Group_list##__dim *  __grpl;\
	Group_type##__dim *  __grp;\
	__grpl=calloc(sizeof(Group_list##__dim),1);\
	memcpy( __grpl, __gr,sizeof(group_list));\
	__grp= __grpl->group = calloc(sizeof(Group_type##__dim), __gr->NGroups);\
	for ( __index=0; __index< __gr->NGroups; __index++)\
	{\
		__grp[ __index].N= __gr->group[ __index].N;\
		__grp[ __index].cut= __gr->group[ __index].cut;\
		__grp[ __index].index= __gr->group[ __index].index;\
		if ( __gr->group[ __index].Center!=NULL)\
		memcpy( __grp[ __index].Center, __gr->group[ __index].Center, __gr->ndims*sizeof(float));\
		if ( __gr->group[ __index].Spin!=NULL)\
		memcpy( __grp[ __index].Spin, __gr->group[ __index].Spin, __gr->ndims*sizeof(float));\
		if ( __gr->group[ __index].VMean!=NULL)\
		memcpy( __grp[ __index].VMean, __gr->group[ __index].VMean, __gr->ndims*sizeof(float));\
		if ( __gr->ndims==3)\
		{\
			memcpy( __grp[ __index].sig_evec, __gr->group[ __index].sig_evec, __gr->ndims*__gr->ndims*sizeof(float));\
			memcpy( __grp[ __index].sig_eval, __gr->group[ __index].sig_eval, __gr->ndims*sizeof(float));\
		}\
	}\
	for ( __index=0; __index < __gr->NGroups; __index++)\
	{\
		free( __gr->group[ __index].VMean);\
		free( __gr->group[ __index].Center);\
		free( __gr->group[ __index].Spin);\
	}\
	free( __gr->group);free( __gr);\
	return (void *) __grpl;\
}\

void *ComputeGroups_Y(OcTree_data *data,int *Id,float *P,float *V,int N,int Nmin,float ll,snapshot_data *snap,float *delta,int ndims_p,int cp)
{
	group_list *gr;
	int i,ndims;
	if (data!=NULL)
	{gr = ComputeGroupsFromTree(data, P,V,Nmin,cp);ndims=data->ndims;}
	else if (snap!=NULL)
	{gr = ComputeGroupsAttributes(Id,Nmin,ll,snap,cp);ndims=3;}
	else
	{gr = ComputeGroupsFromArrays(Id,Nmin,ll,P,V,N,delta,ndims_p,cp);ndims=ndims_p;}

	switch (ndims)
	{
		case 1:ALLOCANDRETGRP(gr,grp,grpl,i,1);break;
		case 2:ALLOCANDRETGRP(gr,grp,grpl,i,2);break;
		case 3:ALLOCANDRETGRP(gr,grp,grpl,i,3);break;
		case 4:ALLOCANDRETGRP(gr,grp,grpl,i,4);break;
		case 5:ALLOCANDRETGRP(gr,grp,grpl,i,5);break;
		case 6:ALLOCANDRETGRP(gr,grp,grpl,i,6);break;
		case 7:ALLOCANDRETGRP(gr,grp,grpl,i,7);break;
		case 8:ALLOCANDRETGRP(gr,grp,grpl,i,8);break;
		case 9:ALLOCANDRETGRP(gr,grp,grpl,i,9);break;
		case 10:ALLOCANDRETGRP(gr,grp,grpl,i,10);break;
		case 11:ALLOCANDRETGRP(gr,grp,grpl,i,11);break;
		case 12:ALLOCANDRETGRP(gr,grp,grpl,i,12);break;
		case 13:ALLOCANDRETGRP(gr,grp,grpl,i,13);break;
		case 14:ALLOCANDRETGRP(gr,grp,grpl,i,14);break;
		case 15:ALLOCANDRETGRP(gr,grp,grpl,i,15);break;
		case 16:ALLOCANDRETGRP(gr,grp,grpl,i,16);break;
	}

	return NULL;
}

group_list *ComputeGroupsAttributes(int *Id,int Nmin,float ll, snapshot_data *snap,int cp)
{
	float delta[3];
	int k;

	for (k=0;k<3;k++)
		delta[k] = snap->header.BoxSize;

	return ComputeGroupsFromArrays(Id,Nmin,ll, snap->Pos,snap->Vel,snap->N,delta,3,cp);
}

group_list *ComputeGroupsFromTree(OcTree_data *data, float *P,float *V, int Nmin,int cp)
{
	int k;
	float delta[data->ndims];

	for (k=0;k<data->ndims;k++)
		delta[k] = data->max[k]-data->min[k];

	return ComputeGroupsFromArrays(data->Fof_Id,Nmin,data->Fof_ll,P,V,data->Npart,delta,data->ndims,cp);
}

group_list *ComputeGroupsFromArrays(int *Id,int Nmin,float ll, float *P,float *V,int N,float *delta,int ndims,int cp)
{
	long i,j,k;
	int *NInGroup;
	int NGroupsTot;
	int NGroups;
	group_list *list;
	group_type *group;

	float *tmp;
	double tmpmat[2*ndims];
	float pos[ndims];
	float *vel;
	int *GroupIndex;

	// So far only works for ndims=3
	double evec[9];
	double eval[3];

	float min[ndims];
	float max[ndims];
	//snapshot_data *snap;

	list = calloc(1,sizeof(group_list));

	if (cp)
	{
		if (P!=NULL)
		{
			list->Pos=malloc((size_t)sizeof(float)*ndims*N);
			memcpy(list->Pos,P,(size_t)sizeof(float)*ndims*N);
		}
		if (V!=NULL)
		{
			list->Vel=malloc((size_t)sizeof(float)*ndims*N);
			memcpy(list->Vel,V,(size_t)sizeof(float)*ndims*N);
		}
	}
	else
	{
		list->Pos = P;
		list->Vel = V;
	}

	list->NMin=Nmin;
	list->NPart=N;
	list->ll = ll;
	list->ndims = ndims;

	NGroupsTot=0;
	for (i=0;i<N;i++)
		if (Id[i]>NGroupsTot) NGroupsTot = Id[i];

	NGroupsTot++;

	NInGroup=calloc (NGroupsTot,sizeof(int));
	for (i=0;i<N;i++)
		NInGroup[Id[i]]++;


	GroupIndex = calloc (NGroupsTot,sizeof(int));

	NGroups=0;

	for (i=0;i<NGroupsTot;i++)
		if (NInGroup[i] >= Nmin)
			GroupIndex[i] = NGroups++;
		else
			GroupIndex[i] = -1;

	list->NGroups = NGroups;
	list->group   = calloc (NGroups,sizeof(group_type));

	for (i=0;i<NGroupsTot;i++)
	{
		if (GroupIndex[i]>=0)
		{
			list->group[GroupIndex[i]].index = calloc (NInGroup[i],sizeof(int));
			list->group[GroupIndex[i]].N=0;
		}
	}

	free(NInGroup);

	for (i=0;i<N;i++)
	{
		if (GroupIndex[Id[i]]>=0)
		{
			group = &list->group[GroupIndex[Id[i]]];
			group->index[group->N] = i;
			group->N++;
		}
	}

	/*
	   group->Center = calloc(ndims,sizeof(float));
	   group->Spin = calloc (ndims,sizeof(float));
	   group->VMean = calloc(ndims,sizeof(float));
	   */
	for (i=0;i<list->NGroups;i++)
	{
		group=&list->group[i];

		group->Spin = NULL;
		group->Center = NULL;
		group->VMean = NULL;
	}

	if (P!=NULL)
	{

		//Check boundaries
		for (i=0;i<list->NGroups;i++)
		{
			group=&list->group[i];

			group->Center = calloc(ndims,sizeof(float));


			for (j=0;j<ndims;j++)
			{
				min[j]=1.E15;
				max[j]=-1.E15;
			}

			for (j=0;j<group->N;j++)
			{
				for (k=0;k<ndims;k++)
				{
					if (P[(size_t)ndims*group->index[j]+k]<min[k])
						min[k] = P[(size_t)ndims*group->index[j]+k];

					if (P[(size_t)ndims*group->index[j]+k]>max[k])
						max[k] = P[(size_t)ndims*group->index[j]+k];
				}
			}

			for (k=0;k<ndims;k++)
				if (max[k]-min[k] >= delta[k] - 2*ll)
					group->cut|=FLAG_MIN(k);

			//Compute center
			for (j=0;j<group->N;j++)
			{
				for (k=0;k<ndims;k++)
					if ((group->cut&FLAG_MIN(k))&&
							(P[(size_t)ndims*group->index[j]+k]>max[k]-delta[k]/2))
						group->Center[k] += P[(size_t)ndims*group->index[j]+k]-delta[k];
					else
						group->Center[k] +=P[(size_t)ndims*group->index[j]+k];
			}

			for (j=0;j<ndims;j++)
				group->Center[j]/=group->N;

			//Check boundary conditions
			for (k=0;k<ndims;k++)
				if (group->Center[k]<min[k])
				{
					group->cut^=(FLAG_MIN(k)|FLAG_MAX(k));
					group->Center[k]+=delta[k];
				}

		}

		//Compute moments
		if ((V!=NULL)&&(ndims==3))
		{
			//   printf ("Pos and vel \n");
			for (i=0;i<list->NGroups;i++)
			{
				group=&list->group[i];
				group->Spin = calloc (ndims,sizeof(float));
				tmp = group->Spin;
				for (j=0;j<group->N;j++)
				{
					for (k=0;k<ndims;k++)
					{
						pos[k] = P[(size_t)ndims*group->index[j]+k] - group->Center[k];

						if (pos[k]>delta[k]/2)
							pos[k]-=delta[k];
						else if (pos[k]<-delta[k]/2)
							pos[k]+=delta[k];
					}

					vel = &V[(size_t)ndims*group->index[j]];

					tmp[0] += + pos[1]*vel[2] - pos[2]*vel[1];
					tmp[1] += - pos[0]*vel[2] + pos[2]*vel[0];
					tmp[2] += + pos[0]*vel[1] - pos[1]*vel[0];
				}
			}

		}
	}

	if (V!=NULL)
	{
		//Compute VMean and VSig;
		for (i=0;i<list->NGroups;i++)
		{
			group=&list->group[i];
			group->VMean = calloc(ndims,sizeof(float));
			tmp = group->VMean;
			for (j=0;j<group->N;j++)
			{
				for (k=0;k<ndims;k++)
					tmp[k] += V[(size_t)ndims*group->index[j]+k];
			}
			for (k=0;k<ndims;k++)
				tmp[k]/=group->N;
		}

		if (ndims==3)
		{
			for (i=0;i<list->NGroups;i++)
			{
				group=&list->group[i];
				memset(tmpmat,0,sizeof(double)*6);
				for (j=0;j<group->N;j++)
				{

					tmpmat[0]+=(
							(V[(size_t)3*group->index[j]]-group->VMean[0])*
							(V[(size_t)3*group->index[j]]-group->VMean[0]));
					tmpmat[1]+=(
							(V[(size_t)3*group->index[j]]-group->VMean[0])*
							(V[(size_t)3*group->index[j]+1]-group->VMean[1]));
					tmpmat[2]+=(
							(V[(size_t)3*group->index[j]]-group->VMean[0])*
							(V[(size_t)3*group->index[j]+2]-group->VMean[2]));
					tmpmat[3]+=(
							(V[(size_t)3*group->index[j]+1]-group->VMean[1])*
							(V[(size_t)3*group->index[j]+1]-group->VMean[1]));
					tmpmat[4]+=(
							(V[(size_t)3*group->index[j]+1]-group->VMean[1])*
							(V[(size_t)3*group->index[j]+2]-group->VMean[2]));
					tmpmat[5]+=(
							(V[(size_t)3*group->index[j]+2]-group->VMean[2])*
							(V[(size_t)3*group->index[j]+2]-group->VMean[2]));
				}

				for (j=0;j<6;j++)
					tmpmat[j]=(tmpmat[j]/group->N);

				Diagonalise3x3_oct(tmpmat,eval,evec);
				for (j=0;j<9;j++)
					group->sig_evec[j]=(float)evec[j];
				for (j=0;j<3;j++)
				{
					if (eval[j]<0)
					{
						group->sig_eval[j]=0;
						printf ("Problem computing eigenvalues. Too few particles in group ?\n");
					}
					else group->sig_eval[j]=(float)sqrt(eval[j]);
				}
			}
		}
	}

	free(GroupIndex);

	return list;
}


// give the particle closest to coordinates given in point ([x,y,z])
// 2 modes:
//  1- [x,y,z] is arbitrary (not in the list of Pos)
//    Pivot should be initialised to any particle close to point and
//    p_dmax is the exact distance^2 between the two points.
//  2- [x,y,z] belongs to Pos (it s a point in the list)
//    Pivot should be the corresponding index in Pos i.
//    In the last case, p_dmax should be initialised to 0.

//Pos is sorting in ascending order depending on X (=Pos[(size_t)3*i]) value
//the index of the closest point is returned and it s distance is in p_dmax
int CASSort (float *Pos,int n,float *point,int Pivot,double *p_dmax,int ndims)
{
	double dmax;
	int ic=Pivot;//index of the currently closest
	long i,j,k,l;
	int delta_max;
	double a[ndims],d;
	float x[ndims];

	memcpy(x,point,ndims*sizeof(float));

	//printf ("N=%d\n",n);
	if (Pivot<0)
		ic=-1;
	else if (Pivot>n-1)
		ic=n;

	if (*p_dmax<=0)
		dmax=1.E30;
	else
	{
		dmax=*p_dmax;
	}

	delta_max = (n-1-ic);
	if (ic>delta_max) delta_max=ic;

	//Starting from pivot goes from closest to closest in R
	for (i=ic-1,j=ic+1;(i>=0)||(j<n);i--,j++)
	{
		//printf ("%d,%d\n",i,j);
		if (i>=0)
		{
			k=(long)ndims*i;
			a[0]=x[0]-Pos[k];//point_R-R[i];
			//printf ("i:%e (%e)\n",a,dmax);
			//closest than dmax, have to check
			if (a[0]<=dmax)
			{
				for (l=1;l<ndims;l++)
					a[l] = x[l]-Pos[k+l];

				DSQUARE(a,l,ndims,d);

				if (d<dmax*dmax)
				{
					dmax=sqrt(d);
					ic=i;
				}
			}
			//triangle inequality ==> can stop directly
			else {i=-1;}//printf ("istop\n");}
	}

	if (j<n)
	{
		k=(long)ndims*j;
		a[0]=Pos[k]-x[0];//R[j]-point_R;
		//printf ("j:%e (%e)\n",a,dmax);
		if (a[0]<dmax)
		{
			for (l=1;l<ndims;l++)
				a[l] = x[l]-Pos[k+l];

			DSQUARE(a,l,ndims,d);

			if (d<dmax*dmax)
			{
				dmax=sqrt(d);
				ic=j;
			}
		}
		else {j=n;}//printf ("jstop\n");}
}

//printf ("dmax = %e\n",dmax);

}
//printf ("i=%d, j=%d\n",i,j);

if ((*p_dmax<=0)||(dmax<*p_dmax))
{
	*p_dmax = dmax;
	return ic;
}
else
return Pivot;

return ic;
}

//Pos is sorting in ascending order depending on X (=Pos[(size_t)3*i]) value
//the index of the closest point is returned and it s distance is in p_dmax
int CASNeighbours (float *Pos,int n,float *point,int Pivot,Nei_List *list,Nei_List **lastinlist,int tab_start,int Nnb,int ndims)
{
	double dmax;
	long ic=Pivot;//index of the currently closest
	long i,j,k,l;
	long delta_max;
	double a[ndims];
	double d;

	float x[ndims];
	Nei_List *TmpList=list;
	Nei_List *NewLast;

	memcpy(x,point,ndims*sizeof(float));

	if (Pivot<0)
		ic=-1;
	else if (Pivot>n-1)
		ic=n;

	if ((*lastinlist)->d<=0)
		dmax=1.E30;
	else
	{
		dmax=(*lastinlist)->d;
	}

	delta_max = (n-1-ic);
	if (ic>delta_max) delta_max=ic;

	//Starting from pivot goes from closest to closest in R
	for (i=ic-1,j=ic+1;(i>=0)||(j<n);i--,j++)
	{
		if (i>=0)
		{
			k=(long)ndims*i;
			a[0]=x[0]-Pos[k];//point_R-R[i];

			//closest than dmax, have to check
			if (a[0]<dmax)
			{
				for (l=1;l<ndims;l++)
					a[l] = x[l]-Pos[k+l];

				DSQUARE(a,l,ndims,d);

				if (d<dmax*dmax)
				{
					d=sqrt(d);

					if (TmpList->d<d)
						while ((TmpList->d<d))
							TmpList=TmpList->Next;
					else
						while ((TmpList->Prev!=NULL)&&(TmpList->Prev->d>d))
							TmpList=TmpList->Prev;

					(*lastinlist)->index = tab_start+i;
					(*lastinlist)->d = d;

					if (TmpList->Next!=NULL)
					{
						NewLast=(*lastinlist)->Prev;
						(*lastinlist)->Next = TmpList;

						NewLast->Next = NULL;

						if (TmpList->Prev!=NULL)
						{
							TmpList->Prev->Next=(*lastinlist);
							(*lastinlist)->Prev = TmpList->Prev;

						}
						else (*lastinlist)->Prev = NULL;

						TmpList->Prev = (*lastinlist);
						*lastinlist = NewLast;

					}
					dmax=(*lastinlist)->d;
				}

			}
			//triangle inequality ==> can stop directly
			else {i=-1;}//printf ("istop\n");}
	}

	if (j<n)
	{

		k=(long)ndims*j;
		a[0]=Pos[k]-x[0];//R[j]-point_R;

		if (a[0]<dmax)
		{
			for (l=1;l<ndims;l++)
				a[l] = x[l]-Pos[k+l];
			DSQUARE(a,l,ndims,d);

			if (d<dmax*dmax)
			{
				d=sqrt(d);

				if (TmpList->d<d)
					while ((TmpList->d<d))
						TmpList=TmpList->Next;
				else
					while ((TmpList->Prev!=NULL)&&(TmpList->Prev->d>d))
						TmpList=TmpList->Prev;

				(*lastinlist)->index = tab_start+j;
				(*lastinlist)->d = d;
				if (TmpList->Next!=NULL)
				{
					NewLast=(*lastinlist)->Prev;
					(*lastinlist)->Next = TmpList;
					NewLast->Next = NULL;

					if (TmpList->Prev!=NULL)
					{
						TmpList->Prev->Next=(*lastinlist);
						(*lastinlist)->Prev = TmpList->Prev;
					}
					else
						(*lastinlist)->Prev = NULL;

					TmpList->Prev = (*lastinlist);
					*lastinlist = NewLast;

				}
				dmax=(*lastinlist)->d;
			}

		}
		else {j=n;}//printf ("jstop\n");}
}

}

return 1;
}



int RecursiveGetClosest(OcTree_data *data,Nei_OcTree *Tree,float *RefPos,int Prev_NextIndex,double *dmin,int *min_index)
{
	long i,j;
	long index;
	double min[data->ndims];
	double dmin2[1<<data->ndims];
	double d2;
	long tab_start;
	float x[data->ndims];

	//Computes the minimal distance to every node
	for (i=0;i<Tree->NNext;i++)
	{
		if (i!=Prev_NextIndex)
		{
			index = (long)2*data->ndims*Tree->Next[i]->Cur;

			for (j=0;j<data->ndims;j++)
			{
				if (RefPos[j]<data->limits[index+2*j]) min[j] = data->limits[index+2*j]-RefPos[j];
				else if (RefPos[j]<data->limits[index+2*j+1]) min[j] = 0;
				else min[j] = RefPos[j] - data->limits[index+2*j+1];
			}
			DSQUARE(min,j,data->ndims,dmin2[i]);

			//dmin2[i] = xmin*xmin+ymin*ymin+zmin*zmin;
		}
		else
		{
			dmin2[i] = 1.E30;
		}
	}

	d2=(*dmin)*(*dmin);

	for (i=0;i<Tree->NNext;i++)
	{
		if (dmin2[i]<d2)
		{

			if (Tree->Next[i]->NNext==0)
			{
				tab_start = Tree->Next[i]->index;
				*min_index = tab_start + CASSort (&data->Pos[(size_t)data->ndims*tab_start],Tree->Next[i]->N,RefPos,(*min_index)-tab_start,dmin,data->ndims);
			}
			else
			{
				RecursiveGetClosest(data,Tree->Next[i],RefPos,-1,dmin,min_index);
			}
		}
	}

	if ((Prev_NextIndex >= 0)&&(Tree->Cur!=0))
	{
		//Peut etre une optimisation ici ... a voir

		index = (long)2*data->ndims*Tree->Cur;

		memcpy(x,RefPos,data->ndims*sizeof(float));

		for (j=0,i=1;j<data->ndims;j++)
		{
			if (((RefPos[j]-data->limits[index+2*j])<*dmin)||((data->limits[index+2*j+1]-RefPos[j])<*dmin))
			{
				i=0;
				j=data->ndims+1;
			}
		}
		if (i) return 1;

		return RecursiveGetClosest(data,Tree->Prev,RefPos,Tree->Prev_NextIndex,dmin,min_index);
	}
	else
		return 1;
}

//Returns the index of the point closest from the point of index "index"
//"index" is the index in the original array not the one modified when
//building the tree.
int GetIndexClosest(OcTree_data *data,int index,double *dist,int periodic)
{
	int my_index=data->Inv_Id[index];
	int tree_index=data->Tree_ind[my_index];
	Nei_OcTree *Tree = &data->Tree[tree_index];
	int tab_start = Tree->index;

	double d;
	int closest;
	long i,j;
	double x[data->ndims];

	d=0;

	if (Tree->N==1)
	{

		closest=my_index+1;

		for (i=0;i<data->ndims;i++)
			x[i]=data->Pos[(size_t)data->ndims*closest+i]-data->Pos[(size_t)data->ndims*my_index+i];

		DSQUARE(x,i,data->ndims,d);
		d=sqrt(d);
		//RecursiveGetClosest(data,Tree->Prev,&data->Pos[(size_t)3*my_index],Tree->Prev_NextIndex,&d,&closest);
		//*dist = d*data->ScaleFactor;

		//return data->Id[closest];
	}
	else
	{

		closest = CASSort (&data->Pos[(size_t)data->ndims*tab_start],Tree->N,&data->Pos[(size_t)data->ndims*my_index],my_index-tab_start,&d,data->ndims);
		closest +=tab_start;
	}

	for (i=0;i<data->ndims;i++)
		x[i]=data->Pos[(size_t)data->ndims*my_index+i];


	for (i=0;i<data->ndims;i++)
		if (((x[i]-data->limits[(long)2*data->ndims*tree_index+2*i])<d)||
				((data->limits[(long)2*data->ndims*tree_index+2*i+1]-x[i])<d))
		{
			RecursiveGetClosest(data,Tree->Prev,&data->Pos[(size_t)data->ndims*my_index],Tree->Prev_NextIndex,&d,&closest);
			i=data->ndims+1;
		}

	if (periodic)
	{
		long per_flags=0;
		float *RefPos = &data->Pos[(size_t)data->ndims*my_index];

		for (i=0;i<data->ndims;i++)
		{
			if (RefPos[i]-d<data->min[i])
				per_flags|=(1<<(2*i));
			else if (RefPos[i]+d>data->max[i])
				per_flags|=(1<<(2*i+1));
		}

		//sphere crosses a boundary
		if (per_flags)
		{
			int per_index;
			int inter;
			float New_RefPos[data->ndims];
			int found;
			long locindex;
			//printf ("per:%d\n",per_flags);
			for (per_index=1;per_index<=per_flags;per_index++)
			{
				inter = per_index&per_flags;
				if ( (per_index&(~per_flags)) != 0) continue;
				if (inter)
				{
					memcpy(New_RefPos,RefPos,(size_t)data->ndims*sizeof(float));

					for (i=0;i<data->ndims;i++)
					{
						if (inter&(1<<(2*i)))
							New_RefPos[i]+=(data->max[i]-data->min[i]);
						else if (inter&(1<<(2*i+1)))
							New_RefPos[i]-=(data->max[i]-data->min[i]);
					}

					for (i=0;i<data->ndims;i++)
					{
						if (New_RefPos[i]<data->limits[2*i])
							x[i]=data->limits[2*i];
						else if (New_RefPos[i]>data->limits[2*i+1])
							x[i]=data->limits[2*i+1];
						else x[i]=New_RefPos[i];
					}

					Tree = &data->Tree[0];

					found=1;
					while (found)
					{
						for (i=0;(i<Tree->NNext)&&(found);i++)
						{
							locindex = (long)2*data->ndims*Tree->Next[i]->Cur;
							for (j=0;j<data->ndims;j++)
								if ((x[j]<data->limits[locindex+2*j])||(x[j]>data->limits[locindex+2*j+1]))
								{
									found=0;
								}

							if (found) Tree=Tree->Next[i];
						}

						if (Tree->NNext==0) found=0;
					}

					if (Tree->NNext!=0)
						Tree=Tree->Next[0];

					RecursiveGetClosest(data,Tree->Prev,New_RefPos,1<<30,&d,&closest);
				}
			}
		}
	}

	*dist = d;

	return data->Id[closest];
}

//Returns the index of the point closest from the point of index "index"
//"index" is the index in the original array not the one modified when
//building the tree.
int GetPosClosest(OcTree_data *data,float *p_RefPos,double *dist,int periodic)
{
	float RefPos[data->ndims];
	Nei_OcTree *Tree = &data->Tree[0];
	Nei_OcTree *TmpTree = Tree;
	double x[data->ndims];
	int i,j;
	int index;
	double d;
	int closest;
	int found;
	int nfound=0;

	for (i=0;i<data->ndims;i++)
		x[i] = RefPos[i] = p_RefPos[i];

	if (periodic)
	{
		for (i=0;i<data->ndims;i++)
		{
			if (x[i]<data->min[i])
				x[i]+=data->max[i]-data->min[i];
			else if (x[i]>data->max[i])
				x[i]-=(data->max[i]-data->min[i]);
		}
	}
	else
	{
		//This is if the point is out of the tree
		//the closest point should not be far from this one ...
		//good enought for a first guess
		for (i=0;i<data->ndims;i++)
		{
			if (x[i]<data->min[i])
				x[i]=data->min[i];
			else if (x[i]>data->max[i])
				x[i]=data->max[i];
		}
	}
	//printf("P=%f %f %f\n",x[0],x[1],x[2]);
	//Find the node it belongs to
	found=1;nfound=0;

	while (found)
	{
		for (i=0;i<Tree->NNext;i++)
		{
			index = (long)2*data->ndims*Tree->Next[i]->Cur;
			for (j=0,found=1;(j<data->ndims)&&(found);j++)
				if ((x[j]<data->limits[index+2*j])||
						(x[j]>data->limits[index+2*j+1]))
				{
					found=0;
				}

			if (found)
			{
				Tree=Tree->Next[i];
				nfound=1;
				break;
			}
		}

		if (Tree->NNext==0) found=0;
	}

	while (Tree->NNext!=0) Tree=Tree->Next[0];

	//Finds the point that has closest X value
	index = Get_geq_index(&data->Pos[(size_t)data->ndims*Tree->index],Tree->N,0,RefPos[0],data->ndims);

	if (index>=Tree->N) index = Tree->index+Tree->N-1;
	else index+=Tree->index;

	for (i=0;i<data->ndims;i++)
		x[i] = RefPos[i]-data->Pos[(size_t)data->ndims*index+i];

	DSQUARE(x,i,data->ndims,d);
	d=sqrt(d);

	//finds the closest in current node
	closest = Tree->index + CASSort (&data->Pos[(size_t)data->ndims*Tree->index],Tree->N,RefPos,index - Tree->index,&d,data->ndims);

	for (i=0;i<data->ndims;i++)
		x[i]=RefPos[i];

	//and go on like in point closest

	if (nfound)
	{
		i=1;
		for (i=0;i<data->ndims;i++)
			if (((x[i]-data->limits[(long)2*data->ndims*Tree->Cur+2*i])<=d)||((data->limits[(long)2*data->ndims*Tree->Cur+2*i+1]-x[i])<=d))
			{
				i=1<<30;
			}

		if (i<(1<<30))
		{
			*dist = d;
			return data->Id[closest];
		}

		RecursiveGetClosest(data,Tree->Prev,RefPos,Tree->Prev_NextIndex,&d,&closest);
	}
	else
		RecursiveGetClosest(data,TmpTree,RefPos,TmpTree->Prev_NextIndex,&d,&closest);

	if (periodic)
	{
		long per_flags=0;

		for (i=0;i<data->ndims;i++)
		{
			if (RefPos[i]-d<data->min[i])
				per_flags|=(1<<(2*i));
			else if (RefPos[i]+d>data->max[i])
				per_flags|=(1<<(2*i+1));
		}

		//sphere crosses a boundary
		if (per_flags)
		{
			int per_index;
			int inter;
			float New_RefPos[3];
			int found;

			for (per_index=1;per_index<=per_flags;per_index++)
			{
				inter = per_index&per_flags;
				if ( (per_index&(~per_flags)) != 0) continue;
				if (inter)
				{
					memcpy(New_RefPos,RefPos,(size_t)data->ndims*sizeof(float));

					for (i=0;i<data->ndims;i++)
					{
						if (inter&(1<<(2*i)))
							New_RefPos[i]+=(data->max[i]-data->min[i]);
						else if (inter&(1<<(2*i+1)))
							New_RefPos[i]-=(data->max[i]-data->min[i]);
					}

					for (i=0;i<data->ndims;i++)
					{
						if (New_RefPos[i]<data->limits[2*i])
							x[i]=data->limits[2*i];
						else if (New_RefPos[i]>data->limits[2*i+1])
							x[i]=data->limits[2*i+1];
						else x[i]=New_RefPos[i];
					}

					Tree = &data->Tree[0];

					found=1;
					while (found)
					{
						for (i=0;(i<Tree->NNext)&&(found);i++)
						{
							index = (long)2*data->ndims*Tree->Next[i]->Cur;
							for (j=0;j<data->ndims;j++)
								if ((x[j]<data->limits[index+2*j])||(x[j]>data->limits[index+2*j+1]))
								{
									found=0;
								}

							if (found) Tree=Tree->Next[i];
						}
						if (Tree->NNext==0) found=0;
					}

					if (Tree->NNext!=0)
						Tree=Tree->Next[0];

					RecursiveGetClosest(data,Tree->Prev,New_RefPos,1<<30,&d,&closest);
				}
			}
		}
	}

	*dist = d;
	return data->Id[closest];

}

int RecursiveGetNeighbours(OcTree_data *data,Nei_OcTree *Tree,float *RefPos,int RefIndex,int Prev_NextIndex,Nei_List *list,Nei_List **p_lastinlist,int Nnb)
{

	long i,j;
	long index;
	double min[data->ndims];
	double dmin2[1<<data->ndims];
	double d2;
	long tab_start;

	//return 0;
	index = (long)2*data->ndims*Tree->Cur;

	//Computes the minimal distance to every node
	for (i=0;i<Tree->NNext;i++)
	{
		if (i!=Prev_NextIndex)
		{
			index = (long)2*data->ndims*Tree->Next[i]->Cur;

			for (j=0;j<data->ndims;j++)
			{
				if (RefPos[j]<data->limits[index+2*j]) min[j] = data->limits[index+2*j]-RefPos[j];
				else if (RefPos[j]<data->limits[index+2*j+1]) min[j] = 0;
				else min[j] = RefPos[j] - data->limits[index+2*j+1];
			}
			DSQUARE(min,j,data->ndims,dmin2[i]);

		}
		else
		{
			dmin2[i] = 1.E30;
		}
	}

	if (Tree->NNext==0)
	{
		tab_start = Tree->index;
		CASNeighbours(&data->Pos[(size_t)data->ndims*tab_start],Tree->N,RefPos,RefIndex-tab_start,list,p_lastinlist,tab_start,Nnb,data->ndims);
	}
	else
	{
		d2=(*p_lastinlist)->d*(*p_lastinlist)->d;
		for (i=0;i<Tree->NNext;i++)
			if (dmin2[i]<d2)
				RecursiveGetNeighbours(data,Tree->Next[i],RefPos,RefIndex,-1,list,p_lastinlist,Nnb);
	}

	if ((Prev_NextIndex >= 0)&&(Tree->Cur!=0))
	{
		double a,c;

		index = (long)2*data->ndims*Tree->Cur;

		d2 = (*p_lastinlist)->d*(*p_lastinlist)->d;

		for (j=0,c=0,i=1;j<data->ndims;j++)
		{
			if (RefPos[j]<data->limits[index+2*j]) a = data->limits[index+2*j]-RefPos[j];
			else if (RefPos[j]<data->limits[index+2*j+1]) a = 0;
			else a = RefPos[j] - data->limits[index+2*j+1];

			c+=a*a;
		}
		if (c>d2) return 1;

		return RecursiveGetNeighbours(data,Tree->Prev,RefPos,RefIndex,Tree->Prev_NextIndex,list,p_lastinlist,Nnb);
	}
	else
		return 1;

}

int GetIndexNeighbours_Y(OcTree_data *data,int index,int N,long *id,double *d,int periodic)
{
	Nei_List *plist=NULL;
	int i;

	GetIndexNeighbours(data,index-1,N,&plist,periodic);

	if (id!=NULL)
		for (i=0;i<N;i++)
			id[i]=plist[i].index+1;

	if (d!=NULL)
		for (i=0;i<N;i++)
			d[i]=plist[i].d;

	free(plist);

	return 0;
}

int CompList(const void *pa,const void *pb)
{
	Nei_List *a=(Nei_List *)pa;
	Nei_List *b=(Nei_List *)pb;

	if (a->d>b->d) return 1;
	else if (a->d<b->d) return -1;
	else return 0;
}


// N is the number of neighbours
// list is the list of the neighbours
int GetIndexNeighbours(OcTree_data *data,int index,int N,Nei_List **plist,int periodic)
{
	int my_index=data->Inv_Id[index];
	int tree_index=data->Tree_ind[my_index];
	Nei_OcTree *Tree = &data->Tree[tree_index];
	Nei_List *list=NULL;
	Nei_List *lastinlist=NULL;

	long i,j,k;
	int NInList=0;
	//double dmax;
	Nei_OcTree *TmpTree;
	int LastNext=0;
	double x[data->ndims];

	if (Tree->N>N)
		list = calloc (Tree->N,sizeof(Nei_List));
	else
		list = calloc (N,sizeof(Nei_List));

	//Add the current leave to the list
	for (i=Tree->index;i<Tree->index+Tree->N;i++)
	{
		//Uncomment this to exclude the point itself from the list
		if (i!=my_index)
		{
			list[NInList].index = i;
			for (j=0;j<data->ndims;j++)
				x[j]=data->Pos[(size_t)data->ndims*my_index+j]-data->Pos[(size_t)data->ndims*i+j];

			DSQUARE(x,j,data->ndims,list[NInList].d);
			NInList++;
		}
	}

	//need to get at least N neighbours ...
	while (NInList<N)
	{

		for (i=0;(i<Tree->Prev->NNext)&&(NInList<N);i++)
		{
			if (i!=Tree->Prev_NextIndex)
			{
				TmpTree=Tree->Prev->Next[i];
				if (NInList+TmpTree->N>N)
					list = realloc (list,(NInList+TmpTree->N)*sizeof(Nei_List));

				for (j=TmpTree->index;j<TmpTree->index+TmpTree->N;j++)
				{

					list[NInList].index = j;
					for (k=0;k<data->ndims;k++)
						x[k]=data->Pos[(size_t)data->ndims*my_index+k]-data->Pos[(size_t)data->ndims*j+k];

					DSQUARE(x,k,data->ndims,list[NInList].d);
					NInList++;
				}

			}
		}

		if (NInList<N) Tree=Tree->Prev;
		else LastNext = i;
	}

	qsort (list,NInList,sizeof(Nei_List),CompList);

	NInList = N;
	list = realloc (list,N*sizeof(Nei_List));

	list[0].Next=&list[1];
	list[0].Prev=NULL;
	list[N-1].Next=NULL;
	list[N-1].Prev=&list[N-2];

	lastinlist = &list[N-1];

	for (i=1;i<N-1;i++)
	{
		list[i].Prev = &list[i-1];
		list[i].Next = &list[i+1];
	}

	for (i=0;i<N;i++) list[i].d=sqrt(list[i].d);

	//Finish inspecting the top reached cell
	for (i=LastNext;i<Tree->Prev->NNext;i++)
	{
		if (i!=Tree->Prev_NextIndex)
		{
			//lance la recursion ...
			RecursiveGetNeighbours(data,Tree->Prev->Next[i],&data->Pos[(size_t)data->ndims*my_index],my_index,-1,list,&lastinlist,N);
		}
	}

	Tree=Tree->Prev;

	//On the whole tree
	if (Tree->Cur!=0)
		RecursiveGetNeighbours(data,Tree->Prev,&data->Pos[(size_t)data->ndims*my_index],my_index,Tree->Prev_NextIndex,list,&lastinlist,N);


	if (periodic)
	{
		int per_flags=0;
		float *RefPos = &data->Pos[(size_t)data->ndims*my_index];
		double d=list[N-1].d;
		//float x,y,z;

		for (i=0;i<data->ndims;i++)
		{
			if (RefPos[i]-d<data->min[i])
				per_flags|=(1<<(2*i));
			else if (RefPos[i]+d>data->max[i])
				per_flags|=(1<<(2*i+1));
		}

		//sphere crosses a boundary
		if (per_flags)
		{
			int per_index;
			int inter;
			float New_RefPos[data->ndims];
			int found;
			long locindex;

			//printf ("per:%d\n",per_flags);
			for (per_index=1;per_index<=per_flags;per_index++)
			{
				inter = per_index&per_flags;
				if ( (per_index&(~per_flags)) != 0) continue;
				if (inter)
				{
					memcpy(New_RefPos,RefPos,(size_t)data->ndims*sizeof(float));

					for (i=0;i<data->ndims;i++)
					{
						if (inter&(1<<(2*i)))
							New_RefPos[i]+=(data->max[i]-data->min[i]);
						else if (inter&(1<<(2*i+1)))
							New_RefPos[i]-=(data->max[i]-data->min[i]);
					}

					for (i=0;i<data->ndims;i++)
					{
						if (New_RefPos[i]<data->limits[2*i])
							x[i]=data->limits[2*i];
						else if (New_RefPos[i]>data->limits[2*i+1])
							x[i]=data->limits[2*i+1];
						else x[i]=New_RefPos[i];
					}

					Tree = &data->Tree[0];

					found=1;
					while (found)
					{
						for (i=0;(i<Tree->NNext)&&(found);i++)
						{
							locindex = (long)2*data->ndims*Tree->Next[i]->Cur;
							for (j=0;j<data->ndims;j++)
								if ((x[j]<data->limits[locindex+2*j])||(x[j]>data->limits[locindex+2*j+1]))
								{
									found=0;
								}

							if (found) Tree=Tree->Next[i];
						}

						if (Tree->NNext==0) found=0;
					}

					if (Tree->NNext!=0)
						Tree=Tree->Next[0];

					RecursiveGetNeighbours(data,Tree->Prev,New_RefPos,my_index,Tree->Prev_NextIndex,list,&lastinlist,N);
				}
			}
		}
	}

	qsort (list,N,sizeof(Nei_List),CompList);


	for (i=0;i<N;i++)
	{
		list[i].index=data->Id[list[i].index];
	}

	if (*plist != NULL) free(*plist);
	*plist=list;

	return 0;
}

void GetNeighbours_Y(OcTree_data *octdata,void *data,int ndata, int dataisint,int N,long *id,double *d,int periodic)
{
	float *pos;
	int *index;
	int i,j;
	Nei_List *plist=NULL;
	double dist;

	if (N>1)
	{
		if (dataisint)
		{
			index = data;
			for (j=0;j<ndata;j++)
			{
				GetIndexNeighbours(octdata,index[j]-1,N,&plist,periodic);
				if (id!=NULL)
					for (i=0;i<N;i++)
						id[j*N+i]=plist[i].index+1;

				if (d!=NULL)
					for (i=0;i<N;i++)
						d[j*N+i]=plist[i].d;
			}
		}
		else
		{
			pos = data;
			for (j=0;j<ndata;j++)
			{
				GetPosNeighbours(octdata,&pos[(long)octdata->ndims*j],N,&plist,periodic);

				if (id!=NULL)
					for (i=0;i<N;i++)
						id[j*N+i]=plist[i].index+1;

				if (d!=NULL)
					for (i=0;i<N;i++)
						d[j*N+i]=plist[i].d;
			}
		}
	}
	else
	{
		if (dataisint)
		{
			index = data;
			for (j=0;j<ndata;j++)
			{
				i=GetIndexClosest(octdata,index[j]-1,&dist,periodic);
				if (id!=NULL) id[j]=i+1;
				if (d!=NULL) d[j]=dist;
			}
		}
		else
		{
			pos = data;
			for (j=0;j<ndata;j++)
			{
				i=GetPosClosest(octdata,&pos[(long)octdata->ndims*j],&dist,periodic);
				//printf("pos =(%f %f %f) -> %d\n",pos[(long)octdata->ndims*j],pos[(long)octdata->ndims*j+1],pos[(long)octdata->ndims*j+2],i);fflush(0);
				if (id!=NULL) id[j]=i+1;
				if (d!=NULL) d[j]=dist;
			}
		}
	}

	free(plist);
}

int GetPosNeighbours_Y(OcTree_data *data,float *RefPos,int N,long *id,double *d,int periodic)
{
	Nei_List *plist=NULL;
	int i;

	GetPosNeighbours(data,RefPos,N,&plist,periodic);

	if (id!=NULL)
		for (i=0;i<N;i++)
			id[i]=plist[i].index+1;

	if (d!=NULL)
		for (i=0;i<N;i++)
			d[i]=plist[i].d;

	free(plist);

	return 1;
}

// N is the number of neighbours
// list is the list of the neighbours
int GetPosNeighbours(OcTree_data *data,float *p_RefPos,int N,Nei_List **plist,int periodic)
{
	int my_index;
	Nei_OcTree *Tree = &data->Tree[0];
	Nei_List *list=NULL;
	Nei_List *lastinlist=NULL;

	long i,j,k;
	int NInList=0;
	Nei_OcTree *TmpTree;
	int LastNext=0;
	float RefPos[data->ndims];
	double x[data->ndims];
	double u[data->ndims];
	int found;
	long index;

	for (i=0;i<data->ndims;i++)
		x[i] = RefPos[i] = p_RefPos[i];

	if (periodic)
	{
		for (i=0;i<data->ndims;i++)
		{
			if (x[i]<data->min[i])
				x[i]+=data->max[i]-data->min[i];
			else if (x[i]>data->max[i])
				x[i]-=(data->max[i]-data->min[i]);
		}
	}
	else
	{
		//This is if the point is out of the tree
		//the closest point should not be far from this one ...
		//good enough for a first guess
		for (i=0;i<data->ndims;i++)
		{
			if (x[i]<data->min[i])
				x[i]=data->min[i];
			else if (x[i]>data->max[i])
				x[i]=data->max[i];
		}
	}

	//Find the node it belongs to
	found=1;
	while (found)
	{
		for (i=0,found=0;(i<Tree->NNext)&&(!found);i++)
		{
			index = (long)2*data->ndims*Tree->Next[i]->Cur;
			for (j=0,found=1;(j<data->ndims)&&(found);j++)
				if ((x[j]<data->limits[index+2*j])||
						(x[j]>data->limits[index+2*j+1]))
				{
					found=0;
					break;
				}

			if (found)
			{
				Tree=Tree->Next[i];
				break;
			}
		}
		if ((Tree->NNext==0)||(Tree->N<N)) found=0;
	}

	my_index = Get_geq_index(&data->Pos[(size_t)data->ndims*Tree->index],Tree->N,0,RefPos[0],data->ndims);

	if (my_index>=Tree->N) my_index =  Tree->index+Tree->N-1;
	else my_index+=Tree->index;

	//Tree = data->Tree[data->Tree_ind[my_index]];

	if (Tree->N>N)
		list = calloc (Tree->N,sizeof(Nei_List));
	else
		list = calloc (N,sizeof(Nei_List));

	//Add the current leave to the list
	for (i=Tree->index;i<Tree->index+Tree->N;i++)
	{
		list[NInList].index = i;
		for (j=0;j<data->ndims;j++)
			u[j]=x[j]-data->Pos[(size_t)data->ndims*i+j];

		DSQUARE(u,j,data->ndims,list[NInList].d);

		NInList++;
	}
	//printf ("N=%d\n",Tree->N);
	//printf("ninlist %d\n",NInList);
	//need to get at least N neighbours ...
	while (NInList<N)
	{

		for (i=0;(i<Tree->Prev->NNext)&&(NInList<N);i++)
		{
			if (i!=Tree->Prev_NextIndex)
			{
				TmpTree=Tree->Prev->Next[i];
				if (NInList+TmpTree->N>N)
					list = realloc (list,(NInList+TmpTree->N)*sizeof(Nei_List));

				for (j=TmpTree->index;j<TmpTree->index+TmpTree->N;j++)
				{
					list[NInList].index = j;
					for (k=0;k<data->ndims;k++)
						u[k]=x[k]-data->Pos[(size_t)data->ndims*j+k];

					DSQUARE(u,k,data->ndims,list[NInList].d);
					NInList++;
				}
			}
		}

		if (NInList<N) Tree=Tree->Prev;
		else LastNext = i;
	}
	//printf("ninlist %d\n",NInList);
	qsort (list,NInList,sizeof(Nei_List),CompList);

	NInList = N;

	list = realloc (list,N*sizeof(Nei_List));

	list[0].Next=&list[1];
	list[0].Prev=NULL;

	list[N-1].Next=NULL;
	list[N-1].Prev=&list[N-2];

	lastinlist = &list[N-1];

	for (i=1;i<N-1;i++)
	{
		list[i].Prev = &list[i-1];
		list[i].Next = &list[i+1];
	}

	for (i=0;i<N;i++) list[i].d=sqrt(list[i].d);

	//Finish inspecting the top reached cell
	for (i=LastNext;i<Tree->Prev->NNext;i++)
	{
		if (i!=Tree->Prev_NextIndex)
		{
			//lance la recursion ...
			RecursiveGetNeighbours(data,Tree->Prev->Next[i],RefPos,my_index,-1,list,&lastinlist,N);
		}
	}

	Tree=Tree->Prev;

	//On the whole tree
	if (Tree->Cur!=0)
		RecursiveGetNeighbours(data,Tree->Prev,RefPos,my_index,Tree->Prev_NextIndex,list,&lastinlist,N);

	if (periodic)
	{
		int per_flags=0;
		double d=list[N-1].d;

		for (i=0;i<data->ndims;i++)
		{
			if (RefPos[i]-d<data->min[i])
				per_flags|=(1<<(2*i));
			else if (RefPos[i]+d>data->max[i])
				per_flags|=(1<<(2*i+1));
		}

		//sphere crosses a boundary
		if (per_flags)
		{
			int per_index;
			int inter;
			float New_RefPos[data->ndims];
			int found;

			for (per_index=1;per_index<=per_flags;per_index++)
			{
				inter = per_index&per_flags;

				if ( (per_index&(~per_flags)) != 0) continue;

				if (inter)
				{
					memcpy(New_RefPos,RefPos,(size_t)data->ndims*sizeof(float));

					for (i=0;i<data->ndims;i++)
					{
						if (inter&(1<<(2*i)))
							New_RefPos[i]+=(data->max[i]-data->min[i]);
						else if (inter&(1<<(2*i+1)))
							New_RefPos[i]-=(data->max[i]-data->min[i]);
					}

					for (i=0;i<data->ndims;i++)
					{
						if (New_RefPos[i]<data->limits[2*i])
							x[i]=data->limits[2*i];
						else if (New_RefPos[i]>data->limits[2*i+1])
							x[i]=data->limits[2*i+1];
						else x[i]=New_RefPos[i];
					}

					Tree = &data->Tree[0];

					found=1;
					while (found)
					{
						for (i=0;i<Tree->NNext;i++)
						{
							index = (long)2*data->ndims*Tree->Next[i]->Cur;
							for (j=0;j<data->ndims;j++)
								if ((x[j]<data->limits[index+2*j])||(x[j]>data->limits[index+2*j+1]))
								{
									found=0;
									break;
								}

							if (found)
							{
								Tree=Tree->Next[i];
								break;
							}
						}

						if (Tree->NNext==0) found=0;
					}

					if (Tree->NNext!=0)
						Tree=Tree->Next[0];

					RecursiveGetNeighbours(data,Tree->Prev,New_RefPos,Tree->index,Tree->Prev_NextIndex,list,&lastinlist,N);
				}
			}
		}
	}


	qsort (list,N,sizeof(Nei_List),CompList);


	for (i=0;i<N;i++)
	{
		list[i].index=data->Id[list[i].index];
	}

	if (*plist != NULL) free(*plist);
	*plist=list;

	return 0;
}

int RecursiveGetPointsInSphere(OcTree_data *data,Nei_OcTree *Tree,float *Center,float radius,int **id,int *curid,int *nalloc)
{
	int i,j;
	double dmax;
	double dmin;
	double d;
	float *text;

	text = &data->limits[(long)2*data->ndims*Tree->Cur];

	dmax=dmin=0;
	for (i=0;i<data->ndims;i++)
		if (fabs(Center[i]-text[2*i])<fabs(Center[i]-text[2*i+1]))
		{
			dmax += (Center[i]-text[2*i+1])*(Center[i]-text[2*i+1]);
			if (Center[i]<text[2*i])
				dmin += (Center[i]-text[2*i])*(Center[i]-text[2*i]);
		}
		else
		{
			if (Center[i]>text[2*i+1])
				dmin += (Center[i]-text[2*i+1])*(Center[i]-text[2*i+1]);
			dmax += (Center[i]-text[2*i])*(Center[i]-text[2*i]);
		}

	if (dmax<radius*radius)
	{
		//add all
		if ((*curid)+Tree->N > *nalloc)
		{
			if (Tree->N >(*nalloc))
				(*nalloc) += Tree->N/2;

			(*nalloc)=(*nalloc)*2;
			*id =realloc (*id,sizeof(int)*(*nalloc));
		}

		for (i=0;i<Tree->N;i++)
		{
			(*id)[(*curid)+i] = Tree->index + i;
		}
		(*curid)+=Tree->N;

		return (*curid);
	}


	if (dmin>radius*radius)
	{
		for (i=0,j=1;i<data->ndims;i++)
			if (((Center[i]-text[(long)2*i])*(text[(long)2*i+1] - Center[i]))<0)
				return (*curid);
	}


	if (Tree->NNext==0)
	{
		double x[data->ndims];

		dmax = radius*radius;
		for (i=0;i<Tree->N;i++)
		{
			text=&data->Pos[(long)data->ndims*(Tree->index+i)];
			for (j=0;j<data->ndims;j++)
				x[j] = text[j]-Center[j];

			DSQUARE(x,j,data->ndims,d);

			if (d<dmax)
			{
				if ((*curid)+1 > *nalloc)
				{
					(*nalloc)=(*nalloc)*2;
					*id =realloc (*id,sizeof(int)*(*nalloc));
				}
				(*id)[(*curid)++] = Tree->index + i;
			}
		}
		return (*curid);
	}
	else
	{
		for (i=0;i<Tree->NNext;i++)
		{
			RecursiveGetPointsInSphere(data,Tree->Next[i],Center,radius,id,curid,nalloc);
		}
	}

	return (*curid);
}

int *GetPointsInSphere_Y(OcTree_data *data,float *Center,float radius,int periodic, int *N)
{
	int *id=NULL;
	int i;

	*N = GetPointsInSphere(data,Center,radius,&id,periodic);

	for (i=0;i<*N;i++)
		id[i]++;

	//printf ("allocated pointer : %d\n",(int)id);

	return id;
}

//returns the number of points
int GetPointsInSphere(OcTree_data *data,float *Center,float radius,int **id,int periodic)
{
	Nei_OcTree *Tree = &data->Tree[0];
	long i,j,k,l;
	int inside=0;
	int intersect =0;
	float min[data->ndims];
	float max[data->ndims];
	float *text;
	float ce[data->ndims];
	float oldce[data->ndims];

	int *nid;
	int np;
	int napprox;
	int curid=0;

	for (i=0;i<data->ndims;i++)
		ce[i]=data->max[i]-data->min[i];

	DSQUARE(ce,i,data->ndims,ce[0]);
	//this is the approximate number of particles we expect (for a homogeneous distribution)
	napprox = (int)((4.1887902*radius*radius*radius)/ce[0]);

	if (napprox<=10) napprox = 10;
	nid = malloc (napprox*sizeof(int));

	for (i=0;i<data->ndims;i++)
		ce[i]=Center[i];

	if (periodic)
	{
		//Make sure the center is in the box
		for (i=0;i<data->ndims;i++)
		{
			if (ce[i]<data->min[i])
				ce[i]+=data->max[i]-data->min[i];
			else if (ce[i]>data->max[i])
				ce[i]-=(data->max[i]-data->min[i]);
		}
	}

	for (i=0;i<data->ndims;i++)
	{
		min[i]=ce[i]-radius;
		max[i]=ce[i]+radius;
	}

	text = &data->limits[(long)2*data->ndims*Tree->Cur];

	for (i=0;i<data->ndims;i++)
		if ((min[i]<text[2*i])||(max[i]>text[2*i+1]))
		{
			//totally outside
			for (j=0;j<data->ndims;j++)
				if ((max[j]<text[2*j])||(min[j]>text[2*j+1]))
					return 0;

			for (j=0;j<data->ndims;j++)
			{
				if (min[j]<text[2*j])   {min[j]=text[2*j];intersect|=(1<<(2*j));}
				if (max[j]>text[2*j+1]) {max[j]=text[2*j+1];intersect|=(1<<(2*j+1));}
			}
			break;
		}

	//Finds the widest Tree cell containing the whole sphere
	do {
		inside=0;i=0;
		while ((i<Tree->NNext)&&(!inside))
		{
			text = &data->limits[(long)2*data->ndims*Tree->Next[i]->Cur];
			for (j=0,k=1;j<data->ndims;j++)
				if ((min[j]<text[2*j])||(max[j]>text[2*j+1]))
				{k=0;break;}

			if (k) inside = i+1;
			i++;
		}
		if (inside) Tree=Tree->Next[inside-1];
	} while (inside);


	//Now Tree contains the whole sphere, use recursion ...

	np=RecursiveGetPointsInSphere(data,Tree,ce,radius,&nid,&curid,&napprox);

	//Periodic boundary conditions
	if ((periodic)&&(intersect))
	{
		int inter[data->ndims];
		int ninter=0;
		float delta;
		int m;

		memcpy(oldce,ce,data->ndims*sizeof(float));
		for (k=0;k<2*data->ndims;k++)
		{
			if (intersect&(1<<k))
			{

				inter[ninter]=k;
				ninter++;
			}
		}

		for (i=1;i<(1<<ninter);i++)
		{
			memcpy(ce,oldce,data->ndims*sizeof(float));
			for (j=0;j<ninter;j++)
				if (i&(1<<j))
				{
					delta = data->max[inter[j]>>1]-data->min[inter[j]>>1];
					if (inter[j]&1) ce[inter[j]>>1]-=delta;
					else ce[inter[j]>>1]+=delta;
				}

			Tree = &data->Tree[0];
			//Finds the widest Tree cell containing the whole sphere
			do {
				inside=0;k=0;
				while ((k<Tree->NNext)&&(!inside))
				{
					text = &data->limits[(long)2*data->ndims*Tree->Next[i]->Cur];
					for (j=0,m=1;j<data->ndims;j++)
						if ((min[j]<text[2*j])||(max[j]>text[2*j+1]))
						{m=0;break;}

					if (m) inside = k+1;
					k++;
				}
				if (inside) Tree=Tree->Next[inside-1];
			} while (inside);

			np=RecursiveGetPointsInSphere(data,Tree,ce,radius,&nid,&curid,&napprox);
		}
	}
	if (np!=0)
		nid =realloc (nid,sizeof(int)*np);

	for (i=0;i<np;i++)
		nid[i] = data->Id[nid[i]];

	if (*id != NULL) free(*id);
	*id=nid;

	return np;
}
