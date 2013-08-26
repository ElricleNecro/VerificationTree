#include "Verif_tools.h"

#include <assert.h>

#ifdef __DEBUG_VOIS_LOG
Part MoreDenseParticule(const TNoeud root, const int NbVois, const double BS, const char * fname)
#else
#ifdef PERIODIC
Part MoreDenseParticule(const TNoeud root, const int NbVois, const double BS)
#else
Part MoreDenseParticule(const TNoeud root, const int NbVois)
#endif
#endif
{
	Part *Vois  = NULL;

	Vois        = Part1d(NbVois);

#ifdef __DEBUG_DENSCENTER_P
	fprintf(stderr, "\033[35mRoot de l'arbre : %d particules, nombre de voisin voulu : %d\033[00m\n", root->N, NbVois);
#endif

	double rholoc    = 0.0,
	       rholoctot = 0.0;
	int    NbPart    = root->N,
	       ind       = 0;

#ifdef __DEBUG_VOIS_LOG
	FILE *fich = NULL;
	if( (fich = fopen(fname, "w") ) == NULL )
	{
		perror("Impossible d'ouvrir le fichier : ");
		exit(EXIT_FAILURE);
	}
#endif

	for (int i = 0; i < NbPart; i++)
	{
		//On initialise le tableau des voisins :
		for(int j = 0, k = 0; j < NbVois && k < NbPart; j++,k++)
		{
			if( j == i )
				k++;
			Vois[j].x  = root->first[k].x;
			Vois[j].y  = root->first[k].y;
			Vois[j].z  = root->first[k].z;
			Vois[j].r  = sqrt( pow(root->first[i].x - root->first[k].x, 2.0) + pow(root->first[i].y - root->first[k].y, 2.0) + pow(root->first[i].z - root->first[k].z, 2.0) );
			Vois[j].id = root->first[k].id;
			Vois[j].vx = root->first[k].vx;
			Vois[j].vy = root->first[k].vy;
			Vois[j].vz = root->first[k].vz;
			Vois[j].v  = sqrt( pow(root->first[i].vx - root->first[k].vx, 2.0) + pow(root->first[i].vy - root->first[k].vy, 2.0) + pow(root->first[i].vz - root->first[k].vz, 2.0) );
			if( k == NbPart -1)
				k = 0;
		}
//		for(int j = 0; j < NbVois; j++)
//			Vois[j].r  = 2.0*BS;


		if( i < 1 )
		{
			fprintf(stderr, "Part : %d\n", root->first[i].id);
			for (int j = 0; j < NbVois; j++)
			{
				fprintf(stderr, "%d ", Vois[j].id);
			}
			fprintf(stderr, "\n");
		}

		qsort(Vois, (size_t)NbVois, sizeof(Part), qsort_partstr);

		if( i < 1 )
		{
			for (int j = 0; j < NbVois; j++)
			{
				fprintf(stderr, "%d ", Vois[j].id);
			}
			fprintf(stderr, "\n");
		}
		//Calcul des voisins :
#ifdef PERIODIC
		Tree_Voisin(root, Vois, NbVois, &root->first[i], BS);
#else
		Tree_Voisin(root, Vois, NbVois, &root->first[i]);
#endif

		//Calcul de la densité locale :
		rholoc     = (NbVois) * root->first[0].m / ( (4.0/3.0) * acos(-1.0) * Vois[NbVois -1].r * Vois[NbVois -1].r * Vois[NbVois -1].r );

#ifdef __DEBUG_VOIS_LOG
		fprintf(fich, "%g %g %g %g\n", root->first[i].x, root->first[i].y, root->first[i].z, rholoc);
#endif

		if( rholoc > rholoctot )
		{
			rholoctot = rholoc;
			ind = i;
		}

		if( i < 1 )
		{
			for (int j = 0; j < NbVois; j++)
			{
				fprintf(stderr, "%d ", Vois[j].id);
			}
			fprintf(stderr, "\n");
		}

#ifndef __DENSITYCENTER_NOPROGRESS_P
		fprintf(stderr, "\r\033[31m%s:: %03.3f%%\033[00m", __func__, ( (double)i + 1.0)/( (double)NbPart ) * 100.0);
#endif
	}
#ifndef __DENSITYCENTER_NOPROGRESS_P
	fputs("\n", stderr);
#endif
	Part1d_libere(Vois);

#ifdef __DEBUG_VOIS_LOG
	fclose(fich);
#endif
	fprintf(stderr, "%d %d %g %g %g %g\n", ind, root->first[ind].id, root->first[ind].x, root->first[ind].y, root->first[ind].z, rholoctot);

	return root->first[ind];
}

Part ReCentre(TNoeud root, Part *posvits, const int NbPart, const int NbVois, const int NbMin, const double BoxSize)
{
	/****************************************************************\
	 *		Correction des conditions périodiques		*
	\****************************************************************/
	double Bord    = 0.5 * BoxSize,
	       taille  = BoxSize;
	root           = Tree_Init(NbPart, BoxSize /2.0, BoxSize /2.0, BoxSize /2.0, taille);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);

	root->first = posvits;

	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partaxe);
	fprintf(stderr, "(%g, %g)\t<=>\t(%g, %g, %g)\n", Bord, BoxSize, posvits[NbPart-1].x, posvits[NbPart-1].y, posvits[NbPart-1].z);

	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);
	Tree_Build2(root, NbPart, NbMin);

	Part   TotMove = {.x  = 0.0,
		          .y  = 0.0,
			  .z  = 0.0,
			  .r  = 0.0,
			  .vx = 0.0,
			  .vy = 0.0,
			  .vz = 0.0,
			  .v  = 0.0},
	       Center;

#ifdef USE_TIMER2
	time_t t1, t2;
	clock_t start, finish;
	start = clock();
	time(&t1);
#endif
#ifdef __DEBUG_VOIS_LOG
	Center         = MoreDenseParticule(root, NbVois, BoxSize, "debug_vois_recentre.log");
#else
#ifdef PERIODIC
	Center         = MoreDenseParticule(root, NbVois, BoxSize);
#else
	Center         = MoreDenseParticule(root, NbVois);
#endif
#endif

#ifdef TEST_VOISIN_LOG_SHUFFLE
#warning "Macro TEST_VOISIN_LOG_SHUFFLE activated, performance will be reduced."
	for (int i = 0; i < NbPart; i++)
	{
		int j;
		do {
			j = rand()%NbPart;
		}while( j == i );
		if( i == 0 )
		{
			fprintf(stderr, "%g, %g, %g, %g, %g, %g, %d\n", posvits[i].x, posvits[i].y, posvits[i].z, posvits[i].vx, posvits[i].vy, posvits[i].vz, posvits[i].id );
			fprintf(stderr, "%g, %g, %g, %g, %g, %g, %d\n", posvits[j].x, posvits[j].y, posvits[j].z, posvits[j].vx, posvits[j].vy, posvits[j].vz, posvits[j].id );
		}
		Part  tmp  = posvits[i];
		posvits[i] = posvits[j];
		posvits[j] = tmp;
		if( i == 0 )
		{
			fprintf(stderr, "%g, %g, %g, %g, %g, %g, %d\n", posvits[i].x, posvits[i].y, posvits[i].z, posvits[i].vx, posvits[i].vy, posvits[i].vz, posvits[i].id );
			fprintf(stderr, "%g, %g, %g, %g, %g, %g, %d\n", posvits[j].x, posvits[j].y, posvits[j].z, posvits[j].vx, posvits[j].vy, posvits[j].vz, posvits[j].id );
		}
	}
	Tree_Free(root), root = NULL;
	root           = Tree_Init(NbPart, BoxSize /2.0, BoxSize /2.0, BoxSize /2.0, taille);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);

	root->first = posvits;
	Tree_Build2(root, NbPart, NbMin);

#	ifdef __DEBUG_VOIS_LOG
	Center         = MoreDenseParticule(root, NbVois, BoxSize, "debug_vois_recentre-Manuel.log");
#	else
	Center         = MoreDenseParticule(root, NbVois, BoxSize);
#	endif
#endif
#ifdef USE_TIMER2
	finish = clock();
	time(&t2);
	fprintf(stderr, "\033[32mTemps d'exécution de la fonction %s :: \033[33m%f (%.3f) secondes\033[00m\n", "MoreDenseParticule", difftime(t2,t1), (double)(finish - start) / (double)CLOCKS_PER_SEC);
#endif
	TotMove = Part_add(Center, TotMove);

	for (int i = 0; i < NbPart; i++) {
		posvits[i].x -= Center.x;
		posvits[i].y -= Center.y;
		posvits[i].z -= Center.z;
		posvits[i].r  = sqrt( pow(posvits[i].x, 2.0) + pow(posvits[i].y, 2.0) + pow(posvits[i].z, 2.0) );
	}

#ifdef TEST_INFLUENCE_MODIF_MARCHE_ARBRE
#warning "Macro TEST_INFLUENCE_MODIF_MARCHE_ARBRE activated, performance will be reduced."
	Tree_Free(root), root = NULL;
	root           = Tree_Init(NbPart, 0., 0., 0., taille*20.0);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);
	root->first = posvits;
	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);
	Tree_Build2(root, NbPart, NbMin);

#	ifdef __DEBUG_VOIS_LOG
	Center         = MoreDenseParticule(root, NbVois, BoxSize, "debug_vois_recentre2.log");
#	else
	Center         = MoreDenseParticule(root, NbVois, BoxSize);
#	endif
#endif
	for (int i = 0; i < NbPart; i++) {
		if(	 posvits[i].x > Bord )
			 posvits[i].x -= BoxSize;
		else if( posvits[i].x <= -Bord )
			 posvits[i].x += BoxSize;

		if(	 posvits[i].y > Bord )
			 posvits[i].y -= BoxSize;
		else if( posvits[i].y <= -Bord )
			 posvits[i].y += BoxSize;

		if(	 posvits[i].z > Bord )
			 posvits[i].z -= BoxSize;
		else if( posvits[i].z <= -Bord)
			 posvits[i].z += BoxSize;

		posvits[i].r  = sqrt( pow(posvits[i].x, 2.0) + pow(posvits[i].y, 2.0) + pow(posvits[i].z, 2.0) );
	}
#ifdef TEST_INFLUENCE_MODIF_MARCHE_ARBRE
#warning "Macro TEST_INFLUENCE_MODIF_MARCHE_ARBRE activated, performance will be reduced."
	Tree_Free(root), root = NULL;
	root           = Tree_Init(NbPart, 0., 0., 0., taille*2.0);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);
	root->first = posvits;
	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);
	Tree_Build2(root, NbPart, NbMin);

#	ifdef __DEBUG_VOIS_LOG
	Center         = MoreDenseParticule(root, NbVois, BoxSize, "debug_vois_recentre3.log"); // * 3.086e16 );
#	else
	Center         = MoreDenseParticule(root, NbVois, BoxSize); // * 3.086e16 );
#	endif
#endif

//	root           = Tree_Init(NbPart, 0., 0., 0., taille*2.0);
//	if( root == NULL )
//		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);
//
//	root->first = posvits;
//	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);
//	Tree_Build2(root, NbPart, NbMin);
//
//
//	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partaxe);
//	fprintf(stderr, "(%g, %g)\t<=>\t(%g, %g, %g)\n", Bord, BoxSize, posvits[NbPart-1].x, posvits[NbPart-1].y, posvits[NbPart-1].z);
//
	Tree_Free(root), root = NULL;

	return Center;
}

Part     GravityCenter(const TNoeud root, const int N)
{
	Part center = {.x = 0.0, .y = 0.0, .z = 0.0, .r = 0.0, .vx = 0.0, .vy = 0.0, .vz = 0.0, .v = 0.0};

	for(int i = 0; i< root->N; i++)
	{
		center.x  += root->first[i].x;
		center.y  += root->first[i].y;
		center.z  += root->first[i].z;
		center.vx += root->first[i].vx;
		center.vy += root->first[i].vy;
		center.vz += root->first[i].vz;
	}

	center.x  /= root->N;
	center.y  /= root->N;
	center.z  /= root->N;
	center.vx /= root->N;
	center.vy /= root->N;
	center.vz /= root->N;

	printf("\033[31m%g\t%g\t%g\n\033[00m", center.x, center.y, center.z);
#ifdef __DEBUG_GRAVCENTER_P
	fprintf(stderr, "%g %g %g %g %g %g\n", center.x, center.y, center.z, center.vx, center.vy, center.vz);
#endif

	(void)N;

	return center;
}

#ifdef PERIODIC
Part     DensityCenter(const TNoeud root, const int NbVois, const double BS)
#else
Part     DensityCenter(const TNoeud root, const int NbVois)
#endif
{
	Part center = {.x = 0.0, .y = 0.0, .z = 0.0, .r = 0.0, .vx = 0.0, .vy = 0.0, .vz = 0.0, .v = 0.0},
	     *Vois  = NULL;

	Vois        = Part1d(NbVois);

#ifdef __DEBUG_DENSCENTER_P
	fprintf(stderr, "\033[35mRoot de l'arbre : %d particules, nombre de voisin voulu : %d\033[00m\n", root->N, NbVois);
#endif

	/***************************************************************************\
	 *			Calcul du centre de densité			   *
	\***************************************************************************/
	double rholoc    = 0.0,
	       rholoctot = 0.0;
	int    NbPart    = root->N;

#ifdef __DEBUG_VOIS_LOG
	FILE *fich = NULL;
	if( (fich = fopen("debug_vois_density.log", "w") ) == NULL )
	{
		perror("Impossible d'ouvrir le fichier : ");
		exit(EXIT_FAILURE);
	}
#endif

	for (int i = 0; i < NbPart; i++)
	{
		//On initialise le tableau des voisins :
		for (int j = 0, k = 0; j < NbVois && k < NbPart; j++,k++)
		{
			if( j == i )
				k++;
			//memcpy(&Vois[j], &root->first[k], sizeof(Part));
			Vois[j].x = root->first[k].x;
			Vois[j].y = root->first[k].y;
			Vois[j].z = root->first[k].z;
			Vois[j].r = sqrt( pow(root->first[i].x - root->first[k].x, 2.0) + pow(root->first[i].y - root->first[k].y, 2.0) + pow(root->first[i].z - root->first[k].z, 2.0) );
			Vois[j].vx = root->first[k].vx;
			Vois[j].vy = root->first[k].vy;
			Vois[j].vz = root->first[k].vz;
			Vois[j].v = sqrt( pow(root->first[i].vx - root->first[k].vx, 2.0) + pow(root->first[i].vy - root->first[k].vy, 2.0) + pow(root->first[i].vz - root->first[k].vz, 2.0) );
			if( k == NbPart -1)
				k = 0;
		}

		qsort(Vois, (size_t)NbVois, sizeof(Part), qsort_partstr);
		//Calcul des voisins :
#ifdef PERIODIC
		Tree_Voisin(root, Vois, NbVois, &root->first[i], BS);
#else
		Tree_Voisin(root, Vois, NbVois, &root->first[i]);
#endif

		qsort(Vois, (size_t)NbVois, sizeof(Part), qsort_partstr);
		//Calcul de la densité locale :
		rholoc     = (NbVois - 1) * root->first[0].m / ( (4.0/3.0) * acos(-1.0) * Vois[NbVois -1].r * Vois[NbVois -1].r * Vois[NbVois -1].r );
		rholoctot += rholoc;

#ifdef __DEBUG_VOIS_LOG
		fprintf(fich, "%d %g %g %g %g\n", root->first[i].id, root->first[i].x, root->first[i].y, root->first[i].z, rholoc);
#endif

		center.x  += root->first[i].x  * rholoc;
		center.y  += root->first[i].y  * rholoc;
		center.z  += root->first[i].z  * rholoc;
		center.vx += root->first[i].vx * rholoc;
		center.vy += root->first[i].vy * rholoc;
		center.vz += root->first[i].vz * rholoc;
#ifndef __DENSITYCENTER_NOPROGRESS_P
		fprintf(stderr, "\r\033[31m%s:: %03.3f%%\033[00m", __func__, ( (double)i + 1.0)/( (double)NbPart ) * 100.0);
#endif
	}
#ifndef __DENSITYCENTER_NOPROGRESS_P
	fputs("\n", stderr);
#endif

	center.x  /= rholoctot;
	center.y  /= rholoctot;
	center.z  /= rholoctot;
	center.vx /= rholoctot;
	center.vy /= rholoctot;
	center.vz /= rholoctot;

	printf("\033[31m%g\t%g\t%g\n\033[00m", center.x, center.y, center.z);

	Part1d_libere(Vois);

#ifdef __DEBUG_VOIS_LOG
	fclose(fich);
#endif

	return center;
}

void     Axial_ratio(const TNoeud root, double R_ori, double *grand, double *petit)
{
	//Axial Ratio :
	double  a1	   = 0.0,
		a2	   = 0.0,
		a3	   = 0.0,
		a4	   = 0.0,
		xy	   = 0.0,
		yz	   = 0.0,
		xz	   = 0.0,
		xx	   = 0.0,
		yy	   = 0.0,
		zz	   = 0.0;

	/********************************************************************************************************************************\
	 *					Création de l'histogramme de masse :							*
	\********************************************************************************************************************************/
	for(int i = 0; i < root->N; i++)
	{
		if( root->first[i].r < R_ori )
		{
			//********************************************
			// Calcul des termes de la matrice d'inertie :
			//********************************************
			xx += root->first[i].y * root->first[i].y + root->first[i].z * root->first[i].z;
			yy += root->first[i].x * root->first[i].x + root->first[i].z * root->first[i].z;
			zz += root->first[i].x * root->first[i].x + root->first[i].y * root->first[i].y;
			xy -= root->first[i].x * root->first[i].y;
			xz -= root->first[i].x * root->first[i].z;
			yz -= root->first[i].y * root->first[i].z;
		}
	}

	a1 = -1.0;
	a2 = xx + yy + zz;
	a3 = xy * xy + yz * yz + xz * xz - xx * yy - xx * zz - yy * zz;
	a4 = xx * yy * zz + 2.0 * xy * yz * xz - xx * yz * yz - yy * xz * xz - zz * xy * xy;

	//***************************************
	// Utilisation de la Méthode de Cardan pour résoudre le polynôme d'ordre 3 :
	//**************************************************************************
	double p       = a3/a1- a2/a1 * a2/a1 / 3.0,
	       q       = /*a2/a1 / 27.0 * ( 2.0 * a2/a1 * a2/a1 - 9.0 * a3/a1) + a4/a1,			// */ 2.0 * a2/a1 * a2/a1 * a2/a1 / 27.0 - a2/a1 * a3/a1 / 3.0 + a4/a1,
	       Rac     = /*q*q + 4.0 * p * p * p / 27.0,						// */ (q / 2.0) * (q / 2.0) + (p / 3.0) * (p / 3.0) * (p / 3.0),
	       lambda1 = 0.0,
	       lambda2 = 0.0,
	       lambda3 = 0.0,
	       g_ratio = 0.0,
	       p_ratio = 0.0;

	p        = p/3.0;
	q        = q/2.0;

	if( Rac == 0 )
	{
		double aux1 = -3.0 * q / (2.0 * p)/*,							// */ - a2/a1 / 3.0,
		       aux2 = -3.0 * q / (2.0 * p)/*,							// */ - a2/a1 / 3.0,
		       aux3 = 3.0 * q / p;
		lambda1 = (aux1 > aux3)?aux1:aux3;

		if(lambda1 == aux1)
		{
			lambda2 = aux2;
			lambda3 = aux3;
		}
		else
		{
			lambda2 = aux2;
			lambda3 = aux3;
		}
	}
	else if ( Rac < 0 )
	{
		double phi  = /*acos(-q * sqrt(27.0 / (- p*p*p)) / 2.0)/3.0,				// */ acos(q / (-p * sqrt(-p))),
		       aux1 = /*fabs(2.0*sqrt(-p/3.0) * cos(phi)),					// */ -2.0 * sqrt(-p) * cos(phi / 3.0) - a2/a1 / 3.0,
		       aux2 = /*fabs(2.0*sqrt(-p/3.0) * cos(phi + 2.0 * PI / 3.0)),			// */ 2.0 * sqrt(-p) * cos((PI - phi) / 3.0) - a2/a1 / 3.0,
		       aux3 = /*fabs(2.0*sqrt(-p/3.0) * cos(phi + 4.0 * PI / 3.0));			// */ 2.0 * sqrt(-p) * cos((PI + phi) / 3.0) - a2/a1 / 3.0;
		lambda1     = max(aux1, aux2, aux3);							//(aux1 > aux2)? ( (aux1 > aux3)?aux1:aux3 ):( (aux2 > aux3)?aux2:aux3);
		lambda3     = min(aux1, aux2, aux3);							//(aux1 < aux2)? ( (aux1 < aux3)?aux1:aux3 ):( (aux2 < aux3)?aux2:aux3);
		lambda2     = aux1 + aux2 + aux3 - lambda1 - lambda3;

#ifdef __DEBUG_AXIALRATIO_P
		printf("Lambda1 : %g\tLambda2 : %g\tLambda3 : %g\naux1 : %g\taux2 : %g\taux3 : %g\n", lambda1, lambda2, lambda3, aux1, aux2, aux3);
#endif

		g_ratio     = lambda1 / lambda2;
		p_ratio     = lambda3 / lambda2;
	}
	else
		lambda1     = pow((-q / 2.0 + sqrt(Rac)/2.0), 1.0/3.0) + pow((-q / 2.0 - sqrt(Rac)/2.0), 1.0/3.0)/*; // */ - a2/a1 / 3.0;

	*grand              = g_ratio;
	*petit              = p_ratio;
}

#ifdef PERIODIC
double** CalcPotentiel(const TNoeud root, const double theta, const double rsoft, const double BS)
#else
double** CalcPotentiel(const TNoeud root, const double theta, const double rsoft)
#endif
{
	double **potentiel;
	potentiel  = double2d(root->N, 2);
	for (int i = 0; i < root->N/*1000*/; i++)
	{
#ifdef PERIODIC
		potentiel[i][1] = Tree_CalcPot(root, &root->first[i], theta, rsoft, BS);
#else
		potentiel[i][1] = Tree_CalcPot(root, &root->first[i], theta, rsoft);
#endif
		potentiel[i][0] = root->first[i].r;
#ifndef __POTENTIALCENTER_NOPROGRESS_P
		fprintf(stderr, "\r\033[31mCalcul du potentiel :: %03.3f%%\033[00m", ( (double)i + 1.0)/( (double)root->N ) * 100.0);
#endif
	}
#ifndef __POTENTIALCENTER_NOPROGRESS_P
	fputs("\n", stderr);
#endif

//	qsort(potentiel, (size_t)root->N, sizeof(double*), qsort_comp); //trie_Pot);

	return potentiel;
}

double*  CalcMasse(const TNoeud root)
{
	double *masse = NULL;
	masse = double1d(root->N);

	for(int i = 0; i<root->N; i++)
		masse[i]    = ((i==0)?(0.0):(masse[i-1])) + root->first[i].m; //posvits[i].m;

	return masse;
}

double*  CalcDensite(const TNoeud root, const int NbBin, const double dr, const double rmax)
{
	int     ind     = 0;
	double *densite = NULL;
	densite         = double1d(NbBin);

	for(int i = 0; i < root->N; i++)
	{
		ind = (int)(root->first[i].r/dr);
		if( ind >= NbBin )
		{
			ind = NbBin - 1;
			densite[NbBin-1]   += root->first[i].m;
		}
		else
			densite[ind]       += root->first[i].m;
		assert( ind < NbBin );
	}

	for(int i = 0; i < NbBin; i++)
		densite[i]   = 3.0* densite[i] / (4.0 * PI * dr*dr*dr * (3*i*i + 3*i + 1));

	(void)rmax;

	return densite;
}

double** CalcLogDensite(const TNoeud root, const int NbBin, const double rmin, const double rmax, const double rnorm)
{
	double *densite = NULL,
	       *lr	= NULL,
	        cte     = 0.0;

	densite         = double1d(NbBin);
	lr              = double1d(NbBin + 1);

	lr[0]           = log10((rmin)/rnorm);
	lr[NbBin]       = log10(rmax/rnorm);
	cte             = ( lr[NbBin] - lr[0] ) / NbBin;

	printf("(%g -> %g) ; (%g -> %g)\n", (rmin)/rnorm, lr[0], rmax/rnorm, lr[NbBin]);

	for (int i = NbBin; i > 1; i--)
		lr[i-1] = lr[i] - cte;

	for(int i = 0; i < root->N; i++)
		for(int j = 0; j < NbBin; j++)
			if( root->first[i].r/rnorm > pow(10.0, lr[j]) && root->first[i].r/rnorm <= pow(10.0, lr[j+1]) )
			{
				densite[j]++;
				break;
			}

//	for(int i = 0; i < NbBin; i++)
//		densite[i] *= 3.0 * root->first[0].m / ( 4.0 * PI * (		//pow(10.0, 3.0*lr[i+1]) - pow(10.0, 3.0*lr[i]) )*rnorm*rnorm*rnorm ) ; // root->N;
//					  pow(10.0, 2.0*lr[i+1])		/*lr[i+1]*lr[i+1]*lr[i+1]*/
//					- pow(10.0, 2.0*lr[i])			/*lr[i]*lr[i]*lr[i]*/
//					)*rnorm*rnorm*(pow(10.0, lr[i+1]) - pow(10.0, lr[i]))*rnorm  );

	double **res = double2d(NbBin, 2);

	for (int i = 0; i < NbBin; i++)
	{
		res[i][0] = pow(10., lr[i+1]);
		res[i][1] = densite[i] * root->first[0].m / cte / (4.0 * PI * rnorm*rnorm*rnorm *log(10.0)*pow(10.0, 3.0*lr[i+1]));
	}

	double1d_libere(densite);
	double1d_libere(lr);

	return res;
}

void     CalcEnergie(const TNoeud root, double *energie_c, double *energie_t, double **potentiel, double *Ectot, double *Eptot)
{
	FILE *fich = NULL;
	fich = fopen("log_ener.log", "w");

	*Ectot = 0.0;
	*Eptot = 0.0;

	for (int i = 0; i < root->N; i++)
	{
		energie_c[i] = 0.5 * root->first[i].m * root->first[i].v * root->first[i].v;
		energie_t[i] = energie_c[i] + root->first[i].m * potentiel[i][1];
		if(energie_t[i] == 0.0)
			fprintf(stderr, "%d::v = %g, psi = %g, m = %g\n", i, root->first[i].v, potentiel[i][1], root->first[i].m);

		fprintf(fich, "%g %g %g %g %g %g\n", root->first[i].m, root->first[i].vx, root->first[i].vy, root->first[i].vz, root->first[i].v, energie_c[i]);

		*Ectot += 0.5 * root->first[i].m * root->first[i].v * root->first[i].v;
		*Eptot += root->first[i].m * potentiel[i][1];
	}
	*Eptot = *Eptot / 2.0;
	fclose(fich);
#ifdef __DEBUG_VIRIEL_P
	fprintf(stderr, "%s:: Viriel = %g (%g, %g)\n", __func__, 2.0 * *Ectot / *Eptot, *Ectot, 0.5 * *Eptot);
#endif
}

double*  CalcTemperature(const TNoeud root, const double *densite, const int nb_bin, const double dr, double *Tmoy)
{
	double *temperature = NULL,
	       //*Deltatemp   = NULL,
	        d_all       = 0.0;
	int    *compteur    = NULL,
	       ind          = 0;

	*Tmoy               = 0.0;

	temperature = double1d(nb_bin);
	//Deltatemp   = double1d(nb_bin);
	compteur    = int1d(nb_bin);

//	for(int i = 0; i < nb_bin; i++)
//		temperature[i] = 0.0;

	for(int i = 0; i < root->N; i++)
	{
		ind = (int)(root->first[i].r/dr);
		if( /*root->first[i].r >= rmax || (int)(root->first[i].r/dr)*/ind >= nb_bin )
		{
			temperature[nb_bin-1] += root->first[i].v*root->first[i].v;
			compteur[nb_bin-1]++;
			ind = nb_bin-1;
		}
		else
		{
			temperature[ind /*(int)(root->first[i].r/dr)*/] += root->first[i].v*root->first[i].v;
			compteur[ind/*(int)(root->first[i].r/dr)*/]++;
		}
		assert( ind < nb_bin );
		*Tmoy += root->first[i].v*root->first[i].v;
	}

	for(int i = 0; i < nb_bin; i++)
	{
		temperature[i] *= densite[i] / ( (double)( (compteur[i] != 0)?compteur[i]:1.0 ) );
		//temperature[i] /= (double)( (compteur[i] != 0)?compteur[i]:1.0 );
		d_all += densite[i];
	}

	*Tmoy /= d_all;//root->N;

	for (int i = 0; i < nb_bin; i++) {
		//Deltatemp[i] = temperature[i];
		temperature[i] /= d_all;
	}

//	for(int i = 0; i < root->N; i++)
//	{
//		// Calcul du profil de température :
//		if( root->first[i].r == rmax || (int)(root->first[i].r/dr) == nb_bin )
//			Deltatemp[nb_bin-1] += pow(root->first[i].v*root->first[i].v - temperature[nb_bin - 1], 2.0);
//		else
//			Deltatemp[(int)(root->first[i].r/dr)] += pow(root->first[i].v*root->first[i].v - temperature[(int)(root->first[i].r/dr)], 2.0);
//
//		// Calcul de la température moyenne :
//		DeltaTmoy    += pow(root->first[i].v*root->first[i].v - *Tmoy, 2.0);
//		//DeltaTmoy    += pow(root->first[i].v * masse / 2.0 - *Tmoy, 2.0);
//	}
//
//	// Normalisation de l'écart-type :
//	for(int i = 0; i < nb_bin; i++)
//		Deltatemp[i] /= (double)( (compteur[i] > 1)?(compteur[i]-1):1 ); //(root->N - 1.0);
//	DeltaTmoy /= (root->N - 1.0);

//	*Tmoy      = DeltaTmoy;

	//free(temperature);
	int1d_libere(compteur);

	//return Deltatemp;
	return temperature;
}

double*  CalcJacobien(const TNoeud root, const int NbBin, const double *energie_t, double **potentiel, const double Emin, const double Emax, const double dE, double *distrib)
{
	int NbPart  = root->N,
	    *nb     = NULL;
	double /*Jdr  = 0.0,*/
	       *Jac = NULL,
	       rmax = maxdouble2d(potentiel, NbPart, 0),
	       dr   = rmax / NbBin;

	double *Pot = NULL;

//	fprintf(stderr, "%s::rmax = %g ; dr = %g\n", __func__, rmax, dr);

	Jac = double1d(NbBin);
	Pot = (double *)calloc(NbBin, sizeof(double));
	nb  = int1d(NbBin);

	// Calcul de la distribution en énergie ainsi que de la température :
	for (int i = 0; i < NbPart; i++)
	{
		if( (energie_t[i]-Emin) == (Emax-Emin) || abs((int)( (energie_t[i] - Emin)/dE)) == NbBin )
			distrib[NbBin-1]++;
		else
			distrib[abs((int)( (energie_t[i] - Emin)/dE))]++;

		if( (potentiel[i][0] == rmax) || abs( (int)(potentiel[i][0]/dr) ) == NbBin )
		{
			Pot[NbBin-1] += potentiel[i][1];
			nb[NbBin-1]++;
		}
		else
		{
			Pot[abs( (int)(potentiel[i][0]/dr) )] += potentiel[i][1];
			nb[abs( (int)(potentiel[i][0]/dr) )]++;
		}
	}

	for (int i = 0; i < NbBin; i++) {
		if( nb[i] != 0 )
			Pot[i] /= nb[i];
	}

	lissage(Pot, NbBin);

	FILE *fich = NULL;
	fich = fopen("Potentiel-bin.log", "w");
	for (int i = 0; i < NbBin; i++) {
		fprintf(fich, "%g %g\n", i*dr, Pot[i]);
	}
	fclose(fich);

	// Calcul du jacobien nécessaire à la comparaison avec la théorie plus calcul de la densité :
	for(int i = 0; i < NbBin; i++)
	{
		Jac[i] = 0.0;
#ifdef PREV_NB_NEG_P
		fprintf(stderr, "%s::Négatif[%d] (%d/%d)\r", __func__, i, 0, NbBin);
#endif
		// Intégation du terme en Jacobien :
		for(int z = 0; z < NbBin - 1; z ++)
		{
			if( Emin + i*dE - root->first[0].m*Pot[z] >= 0.0 )
				Jac[i] += 16.0 * PI * PI * root->first[0].m * (z*dr) * (z*dr) * sqrt((Emin + i*dE - root->first[0].m*Pot[z])*2.0*root->first[0].m) * dr;
			else
			{
#ifdef PREV_NB_NEG_P
				fprintf(stderr, "%s::Négatif[%d] (%d/%d)\r", __func__, i, z, NbBin);
#endif
				Jac[i] += 0.0;
			}
		}
#ifdef PREV_NB_NEG_P
		fprintf(stderr, "\n");
#endif

		// On vérifie que le term dont on veut prendre la racine soit bien positif, pour la derniére étape de l'intégration :
		if( (Emin + i*dE - root->first[0].m*Pot[NbBin-1]) > 0 )
			Jac[i]         += 16.0 * PI * PI * root->first[0].m * (NbBin)*dr*(NbBin)*dr * sqrt((Emin + i*dE - root->first[0].m*Pot[NbBin-1])*2.0*root->first[0].m) * dr/2.0;
		// Sinon on met zéro :
		else
			Jac[i]		+= 0.0;
	}
	free(Pot);
	free(nb);

	return Jac;
}

double*  CalcAnisotropie(const TNoeud root, const int NbBin, const double dr, double *Coeff)
{
	double *aniso = NULL,
	       *vr    = NULL,
	       *vt    = NULL,
	       *vp    = NULL,
	       *vrn   = NULL,
	       *vtn   = NULL,
	       *vpn   = NULL,
	       cvr    = 0.0,
	       cvt    = 0.0,
	       cvp    = 0.0,
	       cvrm   = 0.0,
	       cvpm   = 0.0,
	       cvtm   = 0.0;
	int    *nb    = NULL,
	       ind    = 0;

	vr    = double1d(NbBin);
	vt    = double1d(NbBin);
	vp    = double1d(NbBin);
	vrn   = double1d(NbBin);
	vtn   = double1d(NbBin);
	vpn   = double1d(NbBin);
	aniso = double1d(NbBin);

	nb    = int1d(NbBin);

	//Calcul de la moyenne :
	for (int i = 0; i < root->N; i++)
	{
		double svr, vtheta, vphi, theta, phi;

		theta = (root->first[i].z>0.0)?(asin(sqrt(root->first[i].x*root->first[i].x + root->first[i].y*root->first[i].y)/root->first[i].r)):(acos(-1.0) - asin(sqrt(root->first[i].x*root->first[i].x + root->first[i].y*root->first[i].y)/root->first[i].r));
		phi   = atan2(root->first[i].x, root->first[i].y);

		svr = root->first[i].vx * cos(theta)*cos(phi)
		       + root->first[i].vy * cos(theta)*sin(phi)
			- root->first[i].vz * sin(theta);
		/* \f$ \frac{d\phi}{dt} = \frac{1}{y} \left( x \frac{x \frac{dx}{dt} + y \frac{dy}{dt}}{x^2 + y^2} - \frac{dx}{dt} \right) \f$ */
		vphi = -root->first[i].vx * sin(phi) + root->first[i].vy * cos(phi);
		/* \f$ \frac{d\theta}{dt} = \frac{1}{\sqrt{x^2 + y^2}} \left( \frac{dr}{dt}\frac{z}{r} - \frac{dz}{dt} \right) \f$ */
		vtheta = root->first[i].vx * cos(phi)*sin(theta)
				+ root->first[i].vy * sin(phi)*sin(theta)
				+ root->first[i].vz * cos(theta);

		cvr += svr;
		cvp += vphi;
		cvt += vtheta;

		root->first[i].r = sqrt(root->first[i].x*root->first[i].x + root->first[i].y*root->first[i].y + root->first[i].z*root->first[i].z);

		ind = (int)(root->first[i].r/dr);
		if( /*root->first[i].r == rmax || (int)(root->first[i].r/dr)*/ ind >= NbBin )
		{
			vr[NbBin - 1] += svr;
			vp[NbBin - 1] += vphi;
			vt[NbBin - 1] += vtheta;
			nb[NbBin - 1] ++;
			ind = NbBin - 1;
		}
		else
		{
			vr[ind]                               += svr;		// --> r
			vp[ind/*(int)(root->first[i].r/dr)*/] += vphi; 	// --> phi (xOy)
			vt[ind/*(int)(root->first[i].r/dr)*/] += vtheta;	// --> theta (rOz)
			nb[ind/*(int)(root->first[i].r/dr)*/] ++;
		}
		assert( ind < NbBin );
	}
	for (int i = 0; i < NbBin; i++) {
		vr[i] /= ( (nb[i]!=0)?nb[i]:1 );//root->N;
		vt[i] /= ( (nb[i]!=0)?nb[i]:1 );//root->N;
		vp[i] /= ( (nb[i]!=0)?nb[i]:1 );//root->N;
	}
	cvr /= root->N;
	cvt /= root->N;
	cvp /= root->N;

#ifdef __DEBUG_ANISO_P
	fprintf(stderr, "%s:: %g\t%g\t%g\n", __func__, cvt, cvp, cvr);
#endif

#ifdef __LOG_NB_ANISO_P
	FILE *fich = NULL;
	fich = fopen("aniso_nb.log", "w");
	for(int i = 0; i < NbBin; i++)
		fprintf(fich, "%d %d\n", i, nb[i]);
	fclose(fich);
#endif

	ind = 0;
	//Calcul de l'écart-type :
	for (int i = 0; i < root->N; i++)
	{
		double svr, vtheta, vphi, theta, phi;

		theta = (root->first[i].z>0.0)?(asin(sqrt(root->first[i].x*root->first[i].x + root->first[i].y*root->first[i].y)/root->first[i].r)):(acos(-1.0) - asin(sqrt(root->first[i].x*root->first[i].x + root->first[i].y*root->first[i].y)/root->first[i].r));
		phi   = atan2(root->first[i].x, root->first[i].y);

		svr = root->first[i].vx * cos(theta)*cos(phi)
		       + root->first[i].vy * cos(theta)*sin(phi)
			- root->first[i].vz * sin(theta);
		vphi = -root->first[i].vx * sin(phi) + root->first[i].vy * cos(phi);
		vtheta = root->first[i].vx * cos(phi)*sin(theta)
				+ root->first[i].vy * sin(phi)*sin(theta)
				+ root->first[i].vz * cos(theta);

		cvrm += pow(svr - cvr, 2.0);
		cvpm += pow(vphi - cvp, 2.0);
		cvtm += pow(vtheta - cvt, 2.0);

		ind = (int)(root->first[i].r/dr);
		if( /*root->first[i].r == rmax || (int)(root->first[i].r/dr)*/ ind >= NbBin )
		{
			vrn[NbBin - 1] += pow(svr - vr[NbBin - 1], 2.0);
			vpn[NbBin - 1] += pow(vphi - vp[NbBin - 1], 2.0);
			vtn[NbBin - 1] += pow(vtheta - vt[NbBin - 1], 2.0);
			ind = NbBin - 1;
		}
		else
		{
			vrn[ind/*(int)(root->first[i].r/dr)*/] += pow(svr - vr[ind/*(int)(root->first[i].r/dr)*/], 2.0);
			vpn[ind/*(int)(root->first[i].r/dr)*/] += pow(vphi - vp[ind/*(int)(root->first[i].r/dr)*/], 2.0);
			vtn[ind/*(int)(root->first[i].r/dr)*/] += pow(vtheta - vt[ind/*(int)(root->first[i].r/dr)*/], 2.0);
		}
		assert( ind < NbBin );
	}
	for (int i = 0; i < NbBin; i++) {
		vrn[i] /= ( (nb[i]-1>0)?(nb[i]-1):1 );
		vtn[i] /= ( (nb[i]-1>0)?(nb[i]-1):1 );
		vpn[i] /= ( (nb[i]-1>0)?(nb[i]-1):1 );
#ifdef PREV_NB_WEAK_P
		if( nb[i] <= 1 )
			fprintf(stderr, "%s:: Nombre de particules trop faible :: nb[%d] = %d\n", __func__, i, nb[i]);
#endif
	}
	cvrm /= (root->N-1.0);
	cvtm /= (root->N-1.0);
	cvpm /= (root->N-1.0);


	//Calcul de l'anisotropie 1 - 2*sig_r^2 / (sig_t^2 + sig_p^2) :
	for (int i = 0; i < NbBin; i++) {
		if( vrn[i] != 0.0 )
			aniso[i] = 1.0 - ((vtn[i]) + (vpn[i])) / (2.0 * (vrn[i]));
		else
			aniso[i] = 0.0;
	}

	*Coeff = 1.0 - ( (cvtm) + (cvpm) ) / (2.0 * (cvrm));

#ifdef __DEBUG_ANISO_P
	fprintf(stderr, "%s:: %g\t%g\t%g\n\t\t%g\n", __func__, cvtm, cvpm, cvrm, *Coeff);
#endif

	int1d_libere(nb);
	double1d_libere(vr);
	double1d_libere(vt);
	double1d_libere(vp);
	double1d_libere(vrn);
	double1d_libere(vtn);
	double1d_libere(vpn);

	return aniso;
}

/**
 * Retourne le maximum entre 3 paramètres
 * @param a 1éres valeurs à comparer
 * @param b 2nde valeurs à comparer
 * @param c 3émes valeurs à comparer
 * @return Le maximum
 */
double max(double a, double b, double c)
{
	double _max = a;
	if( _max < b )
		_max = b;
	if( _max < c )
		_max = c;
	return _max;
}

/**
 * Retourne le minimum entre 3 paramètres
 * @param a 1éres valeurs à comparer
 * @param b 2nde valeurs à comparer
 * @param c 3émes valeurs à comparer
 * @return Le minimum
 */
double min(double a, double b, double c)
{
	double _min = a;
	if( _min > b )
		_min = b;
	if( _min > c )
		_min = c;
	return _min;
}

