#include "types.h"

Part  Part_add(Part a, Part b)
{
	b.x  += a.x;
	b.y  += a.y;
	b.z  += a.z;
	b.vx += a.vx;
	b.vy += a.vy;
	b.vz += a.vz;
	b.v  += a.v;
	b.r  += a.r;

	return b;
}

Part* Part1d(const int n)
{
	Part * ptf=NULL;
	if (n >0) {
		ptf = (Part *) calloc(n, sizeof(Part)) ;
	}
	if (ptf==NULL) {
		fprintf(stderr, "erreur allocation Part1d\n");
		exit(EXIT_FAILURE);
	}
	return ptf;
}

void Part1d_libere(Part *ptf)
{
	/* liberation d'un tableau 1D de Parts */
	free(ptf);
	return ;
}

int qsort_partstr(const void *a, const void *b)
{
	const Part *p1 = (const Part*)a,
		   *p2 = (const Part*)b;

	if( p1->r > p2->r )//sqrt(p1->x*p1->x + p1->y*p1->y) > sqrt(p2->x*p2->x + p2->y*p2->y))
		return 1;
	else if( p1->r < p2->r )
		return -1;
	else
		return 0;
}

int qsort_partaxe(const void *a, const void *b)
{
	const Part *p1 = (const Part*)a,
		   *p2 = (const Part*)b;

	if( p1->x > p2->x )
		return 1;
	else if( p1->x < p2->x )
		return -1;
	else
		return 0;
}

void Echange(Part *a, Part *b)
{
	Part tmp;

	memcpy(&tmp, a, sizeof(Part));
	memcpy(a,    b, sizeof(Part));
	memcpy(b, &tmp, sizeof(Part));
/*	Part tmp;

	tmp.x  = a->x;
	tmp.y  = a->y;
	tmp.z  = a->z;
	tmp.r  = a->r;
	tmp.id = a->id;

	a->x   = b->x;
	a->y   = b->y;
	a->z   = b->z;
	a->r   = b->r;
	a->id  = b->id;

	b->x   = tmp.x;
	b->y   = tmp.y;
	b->z   = tmp.z;
	b->r   = tmp.r;
	b->id  = tmp.id;
// */
}


