#ifdef USE_TIMER
#include <time.h>
#endif

#include "tree.h"

int Level_Max = 50;
int NB_bro = 8;
int Zc = 0;
double G = 6.67e-11;
//#define NEAREST(x, boxhalf, boxsize) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))

//#include "tree_create.c"
//#include "tree_voisin.c"

void Tree_SetG(double nG)
{
	G = nG;
}
double Tree_GetG(void)
{
	return G;
}

int NotIn(Part cherch, Part *Tab, const int NbVois)
{
	for(int i = 0; i < NbVois; i++)
		if( Tab[i].id == cherch.id )
			return 0;
	return 1;
}

void Tree_var(const int LMax)
{
	Level_Max = LMax;

	NB_bro    = 8;
}

TNoeud Tree_Init(int NbPart, double xc, double yc, double zc, double cote)
{
	TNoeud root  = NULL;

	root         = malloc(sizeof(*root));
	if( root == NULL )
	{
		perror("malloc a rencontré un probléme ");
		fprintf(stderr, "\033[31m%s::%s::%d ==> Échec de l'allocation du noeud\033[00m\n",
			__FILE__, __func__, __LINE__
		       );
		return NULL;
	}
	root->N      = NbPart;
	root->level  = 0;
	root->x      = xc;
	root->y      = yc;
	root->z      = zc;
	root->cote   = cote;

	root->first  = NULL;
	root->parent = NULL;
	root->frere  = NULL;
	root->fils   = NULL;

	return root;
}

void Tree_Save(TNoeud root, FILE *fich)
{
//	if( root->frere == NULL )
//	{
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x - root->cote / 2.0), (root->y + root->cote / 2.0), (root->z + root->cote / 2.0), root->level, root->N);
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x + root->cote / 2.0), (root->y + root->cote / 2.0), (root->z + root->cote / 2.0), root->level, root->N);
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x + root->cote / 2.0), (root->y - root->cote / 2.0), (root->z + root->cote / 2.0), root->level, root->N);
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x - root->cote / 2.0), (root->y - root->cote / 2.0), (root->z + root->cote / 2.0), root->level, root->N);
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x - root->cote / 2.0), (root->y + root->cote / 2.0), (root->z + root->cote / 2.0), root->level, root->N);
		fprintf(fich, "\n");

		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x - root->cote / 2.0), (root->y + root->cote / 2.0), (root->z - root->cote / 2.0), root->level, root->N);
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x + root->cote / 2.0), (root->y + root->cote / 2.0), (root->z - root->cote / 2.0), root->level, root->N);
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x + root->cote / 2.0), (root->y - root->cote / 2.0), (root->z - root->cote / 2.0), root->level, root->N);
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x - root->cote / 2.0), (root->y - root->cote / 2.0), (root->z - root->cote / 2.0), root->level, root->N);
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x - root->cote / 2.0), (root->y + root->cote / 2.0), (root->z - root->cote / 2.0), root->level, root->N);
		fprintf(fich, "\n");

		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x - root->cote / 2.0), (root->y - root->cote / 2.0), (root->z - root->cote / 2.0), root->level, root->N);
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x - root->cote / 2.0), (root->y - root->cote / 2.0), (root->z + root->cote / 2.0), root->level, root->N);
		fprintf(fich, "\n");

		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x + root->cote / 2.0), (root->y - root->cote / 2.0), (root->z - root->cote / 2.0), root->level, root->N);
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x + root->cote / 2.0), (root->y - root->cote / 2.0), (root->z + root->cote / 2.0), root->level, root->N);
		fprintf(fich, "\n");

		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x + root->cote / 2.0), (root->y + root->cote / 2.0), (root->z - root->cote / 2.0), root->level, root->N);
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x + root->cote / 2.0), (root->y + root->cote / 2.0), (root->z + root->cote / 2.0), root->level, root->N);
		fprintf(fich, "\n");

		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x - root->cote / 2.0), (root->y + root->cote / 2.0), (root->z - root->cote / 2.0), root->level, root->N);
		fprintf(fich, "%.16g %.16g %.16g %d %d\n", (root->x - root->cote / 2.0), (root->y + root->cote / 2.0), (root->z + root->cote / 2.0), root->level, root->N);
		fprintf(fich, "\n");
		fprintf(fich, "\n");
//	}

//	if( root->level > Level_Max / 3 )
//		fprintf (stderr, "\033[32m|_ %s:: Level = %d ; x in ]%.16g, %.16g] ; y in ]%.16g, %.16g]\033[00m\n",
//				__func__,
//				root->level,
//				root->x - root->cote / 2.0, root->x + root->cote / 2.0,
//				root->y - root->cote / 2.0, root->y + root->cote / 2.0
//			);

	if( root->frere != NULL )
		Tree_Save(root->frere, fich);
	if( root->fils != NULL )
		Tree_Save(root->fils, fich);
}

int Tree_Read(TNoeud root, FILE *fich)
{
	int res = 0, err = 0;
	if( fscanf(fich, "%lf %lf %lf %d %lf %d", &root->x, &root->y, &root->z, &root->level, &root->cote, &res) != 6 )
		return -6;

	if( root->frere != NULL )
	{
		err = Tree_Read(root->frere, fich);
		if( err != 0 )
			return err;
	}
	if( res == 1 )
	{
		err = Tree_Read(root->fils, fich);
		if( err != 0 )
			return err;
	}

	return 0;
}

void Tree_Free(TNoeud root)
{
	if( root->fils != NULL )
		Tree_Free(root->fils);
	if( root->frere != NULL )
		Tree_Free(root->frere);
	free(root), root = NULL;
}

void Tree_Calc(TNoeud t1, const int NbPart)
{
	t1->CM    = 0.0;
	t1->cm.x  = 0.0;
	t1->cm.y  = 0.0;
	t1->cm.z  = 0.0;

	for(int i = 0; i < NbPart; i++)
	{
#ifdef P_DBG_TREECODE_P_CALC
		//for(int j = 0; j<=t1->level; j++) fprintf(stderr, " ");
		fprintf (stderr, "\033[33m|_ %s:: Level = %d ; x in ]%.16g, %.16g] ; y in ]%.16g, %.16g] ; z in ]%.16g, %.16g] :: P = (%.16g, %.16g, %.16g)::%d\033[00m\n",
				__func__,
				t1->level,
				t1->x - t1->cote / 2.0, t1->x + t1->cote / 2.0,
				t1->y - t1->cote / 2.0, t1->y + t1->cote / 2.0,
				t1->z - t1->cote / 2.0, t1->z + t1->cote / 2.0,
				t1->first[i].x, t1->first[i].y, t1->first[i].z, t1->first[i].id
			);
#endif
		if( (t1->first[i].x > t1->x - t1->cote / 2.0 && t1->first[i].x <= t1->x + t1->cote / 2.0) &&
		    (t1->first[i].y > t1->y - t1->cote / 2.0 && t1->first[i].y <= t1->y + t1->cote / 2.0) &&
		    (t1->first[i].z > t1->z - t1->cote / 2.0 && t1->first[i].z <= t1->z + t1->cote / 2.0)
		  )
		{
			Echange(&t1->first[t1->N], &t1->first[i]);
			t1->CM   += t1->first[t1->N].m;
			t1->cm.x += t1->first[t1->N].x;
			t1->cm.y += t1->first[t1->N].y;
			t1->cm.z += t1->first[t1->N].z;
			t1->N++;
#ifdef P_DBG_TREECODE_P_CALC
			//for(int j = 0; j<=t1->level; j++) fprintf(stderr, " ");
			fprintf(stderr, "\033[36m|_ %s:: Level = %d ==> Prise : %d\033[00m\n", __func__, t1->level, t1->N);
#endif
		}
	}

	if( t1->N != 0 )
	{
		t1->cm.x /= t1->N;
		t1->cm.y /= t1->N;
		t1->cm.z /= t1->N;
#ifdef TREE_CM_BUILD_DEBUG_
		fprintf(stderr, "\033[35mCM :: (%g, %g, %g) ; %g (N == %d)\033[00m\n", t1->x, t1->y, t1->z, t1->CM, t1->N);
#endif
	}
}

TNoeud tmp_Build2(TNoeud root, int NbPart, int bro)
{
	TNoeud t1  = NULL;

	double xc, yc, zc;

	switch(bro)
	{
		case 1:
			xc = root->x - root->cote/4.0;
			yc = root->y + root->cote/4.0;
			zc = root->z + root->cote/4.0;
			break;
		case 2:
			xc = root->x + root->cote/4.0;
			yc = root->y + root->cote/4.0;
			zc = root->z + root->cote/4.0;
			break;
		case 3:
			xc = root->x + root->cote/4.0;
			yc = root->y - root->cote/4.0;
			zc = root->z + root->cote/4.0;
			break;
		case 4:
			xc = root->x - root->cote/4.0;
			yc = root->y - root->cote/4.0;
			zc = root->z + root->cote/4.0;
			break;
		case 5:
			xc = root->x - root->cote/4.0;
			yc = root->y - root->cote/4.0;
			zc = root->z - root->cote/4.0;
			break;
		case 6:
			xc = root->x + root->cote/4.0;
			yc = root->y - root->cote/4.0;
			zc = root->z - root->cote/4.0;
			break;
		case 7:
			xc = root->x + root->cote/4.0;
			yc = root->y + root->cote/4.0;
			zc = root->z - root->cote/4.0;
			break;
		case 8:
			xc = root->x - root->cote/4.0;
			yc = root->y + root->cote/4.0;
			zc = root->z - root->cote/4.0;
			break;
		default:
			fprintf(stderr, "\033[31m%s::Erreur : %d n'est pas prévu.\033[00m\n", __func__, bro);
			exit(EXIT_FAILURE);
			break;
	}

	t1         = Tree_Init( NbPart,
				xc,
				yc,
				zc,
				root->cote / 2.0
			      );
	if( t1 == NULL )
	{
		perror("malloc a rencontré un probléme ");
		fprintf(stderr, "\033[31m%s::%s::%d => Échec de l'allocation du noeud t1 de niveau %d\033[00m\n",
			__FILE__, __func__, __LINE__, root->level + 1);
		return NULL;
	}
	// Question : On garde un indice de la premiére case ou carément l'adresse de
	// la premiére case.
	t1->level  = root->level + 1;
	t1->cote   = root->cote / 2.0; //root->cote / (root->level + 1);
	t1->parent = root;

	return t1;
}

int Tree_Build2(TNoeud root, int NbPart, int NbMin)
{
	if( root->level >= Level_Max )
		return 0;

	int Nb_use = 0,
	    bro    = 0;

	TNoeud t1   = NULL;
	t1          = tmp_Build2(root, 0, bro + 1);
	root->fils  = t1;

	if( Nb_use < NbPart )
	{
		t1->first  = &(root->first[Nb_use]);
		Tree_Calc(t1, NbPart);
	}

#ifdef P_DBG_TREECODE_P_CALC2
	fprintf(stderr, "\033[32m|-%s:: t1->N = %d ; deb = %d ; NbPart = %d, NbMin = %d\033[00m\n", __func__, t1->N, 0, NbPart, NbMin);
#endif

	if( t1->N == 0 )
		t1->first = NULL;

	Nb_use += t1->N;
	bro++;

	while( bro < NB_bro )
	{
//		t1	    = NULL;
//		t1          = tmp_Build2(root, 0, bro + 1);
//		root->frere = t1;
		t1->frere   = tmp_Build2(root, 0, bro + 1);
		t1          = t1->frere;

		if( Nb_use < NbPart )
		{
			t1->first  = &(root->first[Nb_use]);

			Tree_Calc(t1, NbPart - Nb_use);
		}

#ifdef P_DBG_TREECODE_P_CALC2
		fprintf(stderr, "\033[32m|-%s:: t%d->N = %d ; deb = %d ; NbPart = %d, NbMin = %d\033[00m\n", __func__, bro + 1, t1->N, Nb_use, NbPart, NbMin);
#endif

		if( t1->N == 0 )
			t1->first = NULL;

		Nb_use += t1->N;
		bro++;
	}

	if( Nb_use != NbPart )
	{
		fprintf(stderr, "\033[31m%s::Erreur :: Toute les particules n'ont pas été prise au niveau %d (%d au lieu de %d)!!!\033[00m\n", __func__, root->level, Nb_use, NbPart);
		exit(EXIT_FAILURE);
	}

	t1 = root->fils;
	do
	{
		if( t1->N > NbMin )
			Tree_Build2(t1, t1->N, NbMin);
		t1 = t1->frere;
	}
	while(t1 != NULL);

	return 1;
}

#ifdef PERIODIC
void CalcVois(Part *insert, const int N, Part *Tab, const int NbVois, const Part *part, const double BS)
#else
void CalcVois(Part *insert, const int N, Part *Tab, const int NbVois, const Part *part)
#endif
{
#ifdef DOUBLE_BOUCLE
#	ifdef __DEBUG_CALCVOIS_TREECODE_P__
	fprintf(stderr, "\033[35m%s:: Vérification :: (%d, %d)\033[00m\n", __func__, N, NbVois);
#	endif
	Part *di = NULL;
	di       = Part1d(N);
	//Calcul des distances :
	for(int i=0; i<N; i++)
	{
#	ifdef PERIODIC
		di[i].r  = sqrt(  pow( (NEAREST( (insert[i].x - part->x), (BS/2.0), BS )), 2.0 ) +
				  pow( (NEAREST( (insert[i].y - part->y), (BS/2.0), BS )), 2.0 ) +
				  pow( (NEAREST( (insert[i].z - part->z), (BS/2.0), BS )), 2.0 )
			    );
#	else
		di[i].r  = sqrt(  pow( (insert[i].x - part->x), 2.0 ) +
				  pow( (insert[i].y - part->y), 2.0 ) +
				  pow( (insert[i].z - part->z), 2.0 )
			    );
#	endif
		di[i].id = insert[i].id;
#	ifdef __DEBUG_CALCVOIS_TREECODE_P__
		fprintf(stderr, "\033[36m%s::di :: %.16g (%.16g)\033[00m\n", __func__, di[i].r, Tab[NbVois - 1].r);
#	endif
		// Il faut conserver le tableau ordonné, ou on fait un qsort après la boucle :
#	ifndef USE_VOIS_QSORT
		for (int j = i-1/*N-2*/; j >= 0; j--)
		{
			if( di[j].r > di[j+1].r )
			{
				Echange(&di[j], &di[j+1]);
			}
			else
				break;
		}
#	endif
	}
#	ifdef USE_VOIS_QSORT
#		warning "Use of qsort in neighbourhood research : performance will be reduced."
	qsort(Tab, (size_t)NbVois, sizeof(Part), qsort_partstr);
#	endif

	for(int i=0; i<N; i++)
	{
		if( di[i].r < Tab[NbVois - 1].r && di[i].id != part->id && NotIn(di[i], Tab, NbVois) )
		{
			Tab[NbVois - 1].id = di[i].id;
			Tab[NbVois - 1].r  = di[i].r;
#	ifdef __DEBUG_CALCVOIS_TREECODE_P__
			fprintf(stderr, "\033[38m%s::selection :: %.16g\033[00m\n",
					__func__,
					Tab[NbVois - 1].r);
#	endif
			//On garde le tableau des voisins ordonné :
#	ifdef USE_VOIS_QSORT
#		warning "Use of qsort in neighbourhood research : performance will be reduced."
			qsort(Tab, (size_t)NbVois, sizeof(Part), qsort_partstr);
#	else
			for (int j = NbVois-2; j >= 0; j--)
			{
				if( Tab[j].r > Tab[j+1].r )
				{
					Echange(&Tab[j], &Tab[j+1]);
				}
				else
					break;
			}
#	endif
		}
		else if( di[i].r > Tab[NbVois - 1].r )
			break;
	}
	free(di);
#else
	double di = 0.0;
	for(int i = 0; i < N; i++)
	{
#	ifdef PERIODIC
		di  = sqrt(  pow( (NEAREST( (insert[i].x - part->x), (BS/2.0), BS )), 2.0 ) +
				  pow( (NEAREST( (insert[i].y - part->y), (BS/2.0), BS )), 2.0 ) +
				  pow( (NEAREST( (insert[i].z - part->z), (BS/2.0), BS )), 2.0 )
			    );
#	else
		di  = sqrt(  pow( (insert[i].x - part->x), 2.0 ) +
				  pow( (insert[i].y - part->y), 2.0 ) +
				  pow( (insert[i].z - part->z), 2.0 )
			    );
#	endif
		if( di < Tab[NbVois - 1].r && insert[i].id != part->id && NotIn(insert[i], Tab, NbVois) )
		{
			Tab[NbVois - 1].id = insert[i].id;
			Tab[NbVois - 1].r  = di;
			for (int j = NbVois-2; j >= 0; j--)
			{
				if( Tab[j].r > Tab[j+1].r )
				{
					Echange(&Tab[j], &Tab[j+1]);
				}
				else
					break;
			}
		}
	}
#endif
}

#ifdef PERIODIC
/*inline*/ double Tree_Dist(const TNoeud root, const Part *part, const double BS)
{
	// Calcul des distances particules--cube :
	double dx = 0.0,
	       dy = 0.0,
	       dz = 0.0,
	       d  = 0.0;

	//Faire par rapport au centre. max d-cote/2, 0
//	d  = part->x - root->x;
//	d  = NEAREST(d, (BS/2.0), BS);
//	dx = fmax( fmin(d, d - root->cote/2.0), 0.0);
	d  = fabs( root->x - part->x );
	dx = fmax( 0., fmin( d, BS - d ) - root->cote/2.0 );

//	d  = part->x - root->x;
//	d  = NEAREST(d, (BS/2.0), BS);
//	dy = fmax( fmin(d, d - root->cote/2.0), 0.0);
	d  = fabs( root->y - part->y );
	dy = fmax( 0., fmin( d, BS - d ) - root->cote/2.0 );

//	d  = part->x - root->x;
//	d  = NEAREST(d, (BS/2.0), BS);
//	dz = fmax( fmin(d, d - root->cote/2.0), 0.0);
	d  = fabs( root->z - part->z );
	dz = fmax( 0., fmin( d, BS - d ) - root->cote/2.0 );

	return sqrt(dx*dx + dy*dy + dz*dz);
}
#else
/*inline*/ double Tree_Dist(const TNoeud root, const Part *part)
{
	double dx = 0.0,
	       dy = 0.0,
	       dz = 0.0,
	       d1 = 0.0,
	       d2 = 0.0;

//	d1 =
	d1 = part->x - (root->x - root->cote/2.0);
	d2 = part->x - (root->x + root->cote/2.0);

	if( d1 > 0.0 && d2 <= 0.0 )
		dx = 0.0;
	else
		dx = fmin(fabs(d1), fabs(d2));

	d1 = part->y - (root->y - root->cote/2.0);
	d2 = part->y - (root->y + root->cote/2.0);

	if( d1 > 0.0 && d2 <= 0.0 )
		dy = 0.0;
	else
		dy = fmin(fabs(d1), fabs(d2));

	d1 = part->z - (root->z - root->cote/2.0);
	d2 = part->z - (root->z + root->cote/2.0);

	if( d1 > 0.0 && d2 <= 0.0 )
		dz = 0.0;
	else
		dz = fmin(fabs(d1), fabs(d2));

	return sqrt( dx*dx + dy*dy + dz*dz );
}
#endif

#ifdef PERIODIC
void Tree_Voisin(TNoeud root, Part *Tab, int NbVois, const Part *part, const double BS)
#else
void Tree_Voisin(TNoeud root, Part *Tab, int NbVois, const Part *part)
#endif
{
#ifdef __DEBUG_CALCVOIS_TREECODE_P1__
	fprintf(stderr, "%s:: level -> %d ; indice -> %d\n", __func__, root->level, 0);
#endif
	/************************************************************************\
	 *	1ére étape : Vérifier qu'au moins une particule est candidate	*
	\************************************************************************/
#ifdef PERIODIC
	if( Tree_Dist(root, part, BS) > Tab[NbVois-1].r )
#else
	if( Tree_Dist(root, part) > Tab[NbVois-1].r )
#endif
	{
#ifdef __DEBUG_CALCVOIS_TREECODE_P2__
		fprintf(stderr, "\033[37m%d::Inutile de descendre plus bas (level : %d) !!!\033[00m\n", part->id, root->level);
#endif
		return ;
	}
#ifdef __DEBUG_CALCVOIS_TREECODE_P3__
	fprintf(stderr, "\033[37mDescendre plus bas (level : %d) !!!\033[00m\n", root->level);
#endif


	if(root->fils != NULL)
	{
		TNoeud t1 = root->fils;
		do
		{
#ifdef PERIODIC
			Tree_Voisin(t1, Tab, NbVois, part, BS);
#else
			Tree_Voisin(t1, Tab, NbVois, part);
#endif
			t1 = t1->frere;
		}
		while(t1 != NULL);
	}
	else if( root->N !=0 )
	{
#ifdef __DEBUG_CALCVOIS_TREECODE_P4__
		fprintf(stderr, "\033[35mParcours (level : %d) (%p, %d, %d)!!!\033[00m\n", root->level, root->first, root->N, NbVois);
#endif
#ifdef PERIODIC
		CalcVois(root->first, root->N, Tab, NbVois, part, BS);
#else
		CalcVois(root->first, root->N, Tab, NbVois, part);
#endif
	}

	return ;
}


double Tree_ExactPot(const TNoeud root, const Part *part, const double soft)
{
	double pot = 0.0;

	for(int i=0; i<root->N; i++)
	{
		if( root->first[i].id != part->id )
		{
			pot += -G * root->first[i].m / (
				sqrt( pow(root->first[i].x - part->x, 2.0) + pow(root->first[i].y - part->y, 2.0) + pow(root->first[i].z - part->z, 2.0) )
				       	+ soft);
		}
	}
	return pot;
}

double Tree_ApproxPot(const TNoeud root, const Part *part, const double soft)
{
	double x = 0.0,
	       y = 0.0,
	       z = 0.0;
	x = root->cm.x;
	y = root->cm.y;
	z = root->cm.z;

	if( root->N != 0 )
		return -G * root->CM / (sqrt( pow(x - part->x, 2.0) + pow(y - part->y, 2.0) + pow(z - part->z, 2.0) ) + soft);
	else
		return 0.0;
}

#ifdef PERIODIC
#	ifdef __bool_true_false_are_defined
/*inline*/ bool Tree_Accept(const TNoeud root, const Part const * part, const double accept, const double BS)
#	else
/*inline*/ int Tree_Accept(const TNoeud root, const Part const * part, const double accept, const double BS)
#	endif
#else
#	ifdef __bool_true_false_are_defined
/*inline*/ bool Tree_Accept(const TNoeud root, const Part const * part, const double accept)
#	else
/*inline*/ int Tree_Accept(const TNoeud root, const Part const * part, const double accept)
#	endif
#endif
{
#ifdef PERIODIC
	if( Tree_Dist(root, part, BS) == 0.0)
#else
	if( Tree_Dist(root, part) == 0.0)
#endif
#ifdef __bool_true_false_are_defined
		return true;
#else
		return 1;
#endif
#ifdef PERIODIC
	return ( root->cote / Tree_Dist(root, part, BS) ) > accept;
#else
	return ( root->cote / Tree_Dist(root, part) ) > accept;
#endif
}

#ifdef PERIODIC
double Tree_CalcPot(TNoeud root, const Part *part, const double accept, const double soft, const double BS)
#else
double Tree_CalcPot(TNoeud root, const Part *part, const double accept, const double soft)
#endif
{
	double pot = 0.0;

#ifdef PERIODIC
	if( Tree_Accept(root, part, accept, BS) //( ( root->cote / Tree_Dist(root, part) ) > accept )
		&& root->fils != NULL
	  )
#else
	if( Tree_Accept(root, part, accept) //( ( root->cote / Tree_Dist(root, part) ) > accept )
		&& root->fils != NULL
	  )
#endif
	{
		TNoeud t1 = root->fils;
		do
		{
			//Tree_Voisin(t1, Tab, NbVois, part);
#ifdef PERIODIC
			pot += Tree_CalcPot(t1, part, accept, soft, BS);
#else
			pot += Tree_CalcPot(t1, part, accept, soft);
#endif
			t1   = t1->frere;
		}
		while(t1 != NULL);
		return pot;
	}

	if( root->fils == NULL )
		return pot + Tree_ExactPot(root, part, soft);
	else
	{
#ifdef TREE_CALCPOT_DEBUG_
	fprintf(stderr, "\033[31mUtilisation de ApproxPot :: %d (val : %g ; %g) (critere : %g, soft : %g)\033[00m\n",
#ifdef PERIODIC
			Tree_Accept(root, part, accept, BS),
#else
			Tree_Accept(root, part, accept),
#endif
			root->cote,
#ifdef PERIODIC
			Tree_Dist(root, part, BS),
#else
			Tree_Dist(root, part),
#endif
			accept,
			soft);
#endif
		return pot + Tree_ApproxPot(root, part, soft);
	}
}

