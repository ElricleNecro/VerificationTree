/****************************************************************************************\
 *	But    : Obtenir les densité et fonction de distribution à partir de position   *
 *	         et vitesse de particule.						*
 *	Entrée : Fichier snapshot Gadget2						*
 *	Sortie : 3 fichiers : densité(r), potntiel(r) et distribution(E)		*
 * Auteur      : Guillaume Plum								*
 * Date        : Vendredi 13 Mai 2011							*
 * Cadre       : Stage de M1 sur la stabilité des systèmes auto-gravitants		*
 * Modification: Thése									*
 * 											*
 * Dépendances :									*
 * 		-- Librairies standards stdio.h et stdlib.h				*
 * 		-- Librairies mathématiques (fonctions sqrt, exp) tgmath.h		*
 *		-- Librairies associées au Tree Code (livrées avec les sources		*
 * 											*
 * TODO : Calculer Nombre de particule 	au-delà du rayon d'origine. 			*
 * 500000 part -- horizon8 -- eps ~= 0.05 t0 = 0.0 et tmax = 1.0 et 20 snappshot	*
 * 100000 part -- horizon  -- pour verif, bcp snapshot			Ok !!!		*
 * 5000	       -- t = 0.1 et bcp bcp de snapshot (~200)			Ok !!!		*
 *											*
 * TODO2 : * Facteur d'anisotropie (1 +/- 2*sig^2_r / sig^2_t				*
 *	   * Température contre rayon							*
 *	   * Fonction de distribution f(E)						*
 *	   * Rapport du Viriel								*
\****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

#ifdef USE_TIMER
#include <time.h>
#endif

#ifdef USE_MCHECK
#include <mcheck.h>
#endif

#include "tree.h"
#include "utils.h"
#include "types.h"
#include "cte_phys.h"

#ifdef IO_PERSO
#include "iogadget2.h"
#else
#include "io_snapshot.h"
#endif

#include "Verif_tools.h"

// Flags de debug : 1 activé, 0 sinon
#define __DBG_SORT_MYSORT_P 0
#define __DBG_SORT__QSORT_P 0

#ifdef __DEBUG_VOIS_LOG
Part MoreDenseParticule(const TNoeud root, const int NbVois, const double BS, const char * fname)
#else
Part MoreDenseParticule(const TNoeud root, const int NbVois, const double BS)
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
			Vois[j].vx = root->first[k].vx;
			Vois[j].vy = root->first[k].vy;
			Vois[j].vz = root->first[k].vz;
			Vois[j].v  = sqrt( pow(root->first[i].vx - root->first[k].vx, 2.0) + pow(root->first[i].vy - root->first[k].vy, 2.0) + pow(root->first[i].vz - root->first[k].vz, 2.0) );
			if( k == NbPart -1)
				k = 0;
		}
//		for(int j = 0; j < NbVois; j++)
//			Vois[j].r  = 2.0*BS;


		qsort(Vois, (size_t)NbVois, sizeof(Part), qsort_partstr);
		//Calcul des voisins :
		Tree_Voisin(root, Vois, NbVois, &root->first[i], BS);

		//Calcul de la densité locale :
		rholoc     = (NbVois - 1) * root->first[0].m / ( (4.0/3.0) * acos(-1.0) * Vois[NbVois -1].r * Vois[NbVois -1].r * Vois[NbVois -1].r );

#ifdef __DEBUG_VOIS_LOG
		fprintf(fich, "%g %g %g %g\n", root->first[i].x, root->first[i].y, root->first[i].z, rholoc);
#endif

		if( rholoc > rholoctot )
		{
			rholoctot = rholoc;
			ind = i;
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

	//Center = GravityCenter(root, NbVois);
	//Center = PeriodicGravityCenter(posvits, NbPart, header.BoxSize*3.086e16);
#ifdef USE_TIMER2
	time_t t1, t2;
	clock_t start, finish;
	start = clock();
	time(&t1);
#endif
#ifdef __DEBUG_VOIS_LOG
	Center         = MoreDenseParticule(root, NbVois, BoxSize, "debug_vois_recentre.log");
#else
	Center         = MoreDenseParticule(root, NbVois, BoxSize);
#endif
//	Tree_Free(root), root = NULL;
//	root           = Tree_Init(NbPart, BoxSize /2.0, BoxSize /2.0, BoxSize /2.0, taille);
//	if( root == NULL )
//		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);
//
//	root->first = posvits;
//	Tree_Build2(root, NbPart, NbMin);
//
//#ifdef __DEBUG_VOIS_LOG
//	Center         = MoreDenseParticule(root, NbVois, BoxSize, "debug_vois_recentre.log");
//#else
//	Center         = MoreDenseParticule(root, NbVois, BoxSize);
//#endif
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

	Tree_Free(root), root = NULL;
	root           = Tree_Init(NbPart, 0., 0., 0., taille*20.0);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);
	root->first = posvits;
	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);
	Tree_Build2(root, NbPart, NbMin);

#ifdef __DEBUG_VOIS_LOG
	Center         = MoreDenseParticule(root, NbVois, BoxSize, "debug_vois_recentre2.log");
#else
	Center         = MoreDenseParticule(root, NbVois, BoxSize);
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
	Tree_Free(root), root = NULL;
	root           = Tree_Init(NbPart, 0., 0., 0., taille*2.0);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);
	root->first = posvits;
	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);
	Tree_Build2(root, NbPart, NbMin);

#ifdef __DEBUG_VOIS_LOG
	Center         = MoreDenseParticule(root, NbVois, BoxSize, "debug_vois_recentre3.log"); // * 3.086e16 );
#else
	Center         = MoreDenseParticule(root, NbVois, BoxSize); // * 3.086e16 );
#endif

	root           = Tree_Init(NbPart, 0., 0., 0., taille*2.0);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);

	root->first = posvits;
	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);
	Tree_Build2(root, NbPart, NbMin);


	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partaxe);
	fprintf(stderr, "(%g, %g)\t<=>\t(%g, %g, %g)\n", Bord, BoxSize, posvits[NbPart-1].x, posvits[NbPart-1].y, posvits[NbPart-1].z);

	Tree_Free(root), root = NULL;

	return Center;
}

int trie_rayon(const void *a, const void *b)
{
	const double *p1 = (const double*)a,
	             *p2 = (const double*)b;

	if( *p1 > *p2 )//sqrt(p1->x*p1->x + p1->y*p1->y) > sqrt(p2->x*p2->x + p2->y*p2->y))
		return 1;
	else if( *p1 < *p2 )
		return -1;
	else
		return 0;
}

void Save_Part(const char const *fn, const Part const *posvits, const int NbPart)
{
	FILE *fich     = NULL;
	if( (fich      = fopen(fn, "w")) == NULL )
	{
		perror("fopen a échouer ");
		fprintf(stderr, "\033[31m%s::%s::%d ==> Impossible d'ouvrir le fichier : %s\n", __FILE__, __func__, __LINE__, fn);
		exit(EXIT_FAILURE);
	}
	for(int i      = 0; i<NbPart; i++)
		fprintf(fich, "%g %g %g %g %g %g %g %g %d\n", posvits[i].x, posvits[i].y, posvits[i].z, posvits[i].r, posvits[i].vx, posvits[i].vy, posvits[i].vz, posvits[i].v, posvits[i].id);
	fclose(fich);
}

void usage(const char *exec)
{
	printf( "\033[04mUtilisation\033[00m :\n"
		"\t\033[34m%s [-g] [-p] [-h] [-n Entier] [-s Réel] [-R Réel] [-G Réel] \033[04mFichier\033[00m\n"
	        "\tFichier : Fichier de données au format Gadget 1\n"
		"\n"
		"\033[04mOption du programme\033[00m :\n"
		"\t\033[35m-G \033[04mRéel\033[00m   : Valeur de G dans les unités du fichier (\033[36;03mdéfaut : 6.67e-11\033[00m). \033[31mInutile pour le moment.\033[00m\n"
		"\t\033[35m-R \033[04mRéel\033[00m   : Rayon maximum souhaité, toutes les particules au-delà seront ignorées (\033[36;03mdéfaut : -1.0 => désactivé\033[00m).\n"
		"\t\033[35m-h \033[00m       : Affiche cette aide.\n"
		"\t\033[35m-s \033[04mRéel\033[00m   : Paramètre de lissage à utiliser (\033[36;03mdéfaut : 0.0\033[00m).\n"
		"\t\033[35m-o \033[04mRéel\033[00m   : Angle d'ouverture de Tree Code (\033[36;03mdéfaut : 0.0 si Nombre de Particules < 1000, 0.5 sinon\033[00m, \033[01mvoir Barnes & Hut, 1986\033[00m).\n"
		"\t\033[35m-n \033[04mEntier\033[00m : Nombre de bin à utiliser (\033[36;03mdéfaut : 100\033[00m).\n"
		"\t\033[35m-t \033[04mEntier\033[00m : Type Gadget à charger (\033[36;03mdéfaut : 4\033[00m). \033[01mNe peut travailler que sur un type à la fois.\033[00m\n"
		"\t\033[35m-v \033[04mEntier\033[00m : Nombre de Voisin à utiliser pour le centre de gravité (\033[36;03mdéfaut : 20\033[00m)\n"
		"\t\033[35m-m \033[04mEntier\033[00m : Nombre minimum de particule à mettre dans un nœud de l'arbre (\033[36;03mdéfaut : 15\033[00m).\n"
		"\t\033[35m-p       \033[00m : La simulation a t elle été faîtes avec des conditions périodiques.\n"
		"\t\033[35m-g \033[00m       : Active l'utilisation du centre de gravité au lieu du centre de densité (\033[36;03mdéfaut : désactivé\033[00m).\n"
		"\n"
		"\033[04mOption longue du programme\033[00m :\n"
		"\t\033[35m--help      \033[00m: Affiche cette aide.\n"
		"\t\033[34m--timeparam \033[00m: Nom du fichier de paramètre.\n"
		"\n"
		"\033[04mÀ venir\033[00m :\n"
		"\t- Support du type de fichier Gadget 2.\n"
		"\n",
		exec
	      );
}

int main(int argc, char **argv)
{
#ifdef USE_MCHECK
	mtrace();
#endif
	if( argc < 2 )
	{
		fprintf(stderr, "\033[31mUsage : (%s) (fichier de données) (Valeur de G) (Rayon d'origine) (softening en pc) [(Nombre de bin)]\n%d arguments au lieu de %d\033[00m\n", argv[0], argc - 1, 4);
		return EXIT_FAILURE;
	}

	int     NbPart     = 0,
		NbPartOri  = 0,
		NbMin      = 15,
		NbVois     = 20,
	        nb_bin     = 100,
		nb_hors_Ro = 0,
		nbfiles    = 1,
		type       = 4;
	bool    gc         = false,
		periodic   = false;
	char   *filename   = NULL,
	       *timeparam  = "TimeParam.dat";
	double  G          = 6.67384e-11,
		R_ori      = -1.0,
		taille     = 0.0,
		theta      = -1.0,
		simu_time  = 0.0,
		rmax       = 0.0,
		rsoft      = 0.0,
		PosFact    = 3.086e16,
		VitFact    = 1.0;

	for (int i = 1; i < argc; i++)
	{
		if( argv[i][0] == '-' )
		{
			switch(argv[i][1])
			{
				case 'G':
					G = atof(argv[i+1]);
					i++;
					break;
				case 'R':
					R_ori = atof(argv[i+1]);
					i++;
					break;
				case 'h':
					usage(argv[0]);
					exit(EXIT_SUCCESS);
					break;
				case 's':
					rsoft = atof(argv[i+1]);
					i++;
					break;
				case 'n':
					nb_bin = atoi(argv[i+1]);
					i++;
					break;
				case 'm':
					NbMin = atoi(argv[i+1]);
					i++;
					break;
				case 'g':
					gc = true;
					break;
				case 'p':
					periodic = true;
					break;
				case 'o':
					theta = atof(argv[i+1]);
					i++;
					break;
				case 't':
					type = atoi(argv[i+1]);
					i++;
					break;
				case 'v':
					NbVois = atoi(argv[i+1]);
					i++;
					break;
				case '-':
					if( !strcmp("--help", argv[i]) )
					{
						usage(argv[0]);
						exit(EXIT_SUCCESS);
					}
					else if( !strcmp("--timeparam", argv[i]) )
					{
						timeparam = argv[i+1];
						i++;
					}
					else if( !strcmp("--posfact", argv[i]) )
					{
						PosFact = atof(argv[i+1]);
						i++;
					}
					else if( !strcmp("--vitfact", argv[i]) )
					{
						VitFact = atof(argv[i+1]);
						i++;
					}
					else
					{
						fprintf(stderr, "\033[00mArgument '%s' invalide\033[00m\n", argv[i]);
						usage(argv[0]);
						exit(EXIT_FAILURE);
					}
					break;
				default:
					fprintf(stderr, "\033[00mArgument '%s' invalide\033[00m\n", argv[i]);
					usage(argv[0]);
					exit(EXIT_FAILURE);
					break;
			}
		}
		else
		{
			filename = argv[i];
		}
	}
	if( filename == NULL )
	{
		fprintf(stderr, "\033[31mLe nom du fichier à lire est obligatoire !!!\033[00m\n");
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}
	if( theta < 0.0 )
	{
		if( NbPart <= 1000 )
			theta = 0.0;
		else
			theta = 0.5;
	}

	printf("Option du programme :\n");
	for (int i = 0; i < argc; i++) {
		printf("\t%d :: %s\n", i, argv[i]);
	}

	/****************************************************************************\
	 *		Chargement des données et allocation :			    *
	\****************************************************************************/

	IO_Header header;
	FILE   *fich   = NULL;
	TNoeud  root   = NULL;
	Part  *posvits = NULL,
	       Center  = {.x  = 0.0,
		          .y  = 0.0,
			  .z  = 0.0,
			  .r  = 0.0,
			  .vx = 0.0,
			  .vy = 0.0,
			  .vz = 0.0,
			  .v  = 0.0},
	       TotMove = {.x  = 0.0,
		          .y  = 0.0,
			  .z  = 0.0,
			  .r  = 0.0,
			  .vx = 0.0,
			  .vy = 0.0,
			  .vz = 0.0,
			  .v  = 0.0};

	Tree_var(50);

	posvits         = read_snapshot(filename, nbfiles, type, PosFact, VitFact, &NbPart, &simu_time, &header);
	header.BoxSize *= PosFact;

	if( posvits == NULL )
		fprintf(stderr, "Erreur avec le tableau de particules !!!\n"),exit(EXIT_FAILURE);
	NbPartOri = NbPart;
	printf("Chargement de %d particule de type %d effectué.\n", NbPart, type);

	Save_Part("Particule-BeCorr.dat", posvits, NbPart);

	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);
#ifdef USE_TIMER
	time_t t1, t2;
	clock_t start, finish;
	start = clock();
	time(&t1);
#endif
	TotMove = ReCentre(root, posvits, NbPart, NbVois, NbMin, header.BoxSize); // * PosFact);
#ifdef USE_TIMER
	finish = clock();
	time(&t2);
	fprintf(stderr, "\033[32mTemps d'exécution de la fonction %s :: \033[33m%f (%.3f) secondes\033[00m\n", "ReCentre", difftime(t2,t1), (double)(finish - start) / (double)CLOCKS_PER_SEC);
#endif
	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);

	rmax        = posvits[NbPart-1].r;
	taille      = /*header.BoxSize; / */ 2.0 * posvits[NbPart-1].r;
	root        = Tree_Init(NbPart, 0.0, 0.0, 0.0, taille);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);

	Save_Part("Particule-Classe.dat", posvits, NbPart);

	root->first = posvits;
	Tree_Build2(root, NbPart, NbMin);

	printf("\033[36mNombre de particule : %d de masse : %g, La plus lointaine : %g (rayon maximum = %g, taille de la boîte = %g)\nsoftening : %g, critère : %g\nNombre de particule par feuille de Tree Code : %d, Nombre de voisin à chercher : %d\033[00m\n",
					NbPart,		posvits[0].m,		  rmax,		R_ori,			taille,		rsoft,		theta,						NbMin,				NbVois);
	if( R_ori > 0.0 )
		printf("\033[36mTaille du Cut-Off : %g\033[00m\n", R_ori);

	// Soustraire Centre gravité.
	// Retirer particule lointaine
	/****************************************************************************\
	 * 			Calcul du centre de gravité :			    *
	\****************************************************************************/
#ifdef USE_TIMER
	start = clock();
	time(&t1);
#endif
	if( gc )
		Center = GravityCenter(root, NbVois);
	else
		Center = DensityCenter(root, NbVois, header.BoxSize); // * PosFact); //3.086e16 );
	printf("\033[31m%g\t%g\t%g\n\033[00m", Center.x, Center.y, Center.z);
#ifdef USE_TIMER
	finish = clock();
	time(&t2);
	fprintf(stderr, "\033[32mTemps d'exécution de la fonction %sCenter :: \033[33m%f (%.3f) secondes\033[00m\n", (gc)?"Gravity":"Density", difftime(t2,t1), (double)(finish - start) / (double)CLOCKS_PER_SEC);
#endif
	for (int i = 0; i < NbPart; i++)
	{
		root->first[i].x  -= Center.x;
		root->first[i].y  -= Center.y;
		root->first[i].z  -= Center.z;
		root->first[i].vx -= Center.vx;
		root->first[i].vy -= Center.vy;
		root->first[i].vz -= Center.vz;
		root->first[i].r   = sqrt( pow(root->first[i].x, 2.0) + pow(root->first[i].y, 2.0) + pow(root->first[i].z, 2.0) );
		root->first[i].v   = sqrt( pow(root->first[i].vx, 2.0) + pow(root->first[i].vy, 2.0) + pow(root->first[i].vz, 2.0) );
		if( R_ori > 0.0 && root->first[i].r < R_ori )
			NbPart--;
	}
	TotMove = Part_add(Center, TotMove);

	Tree_Free(root), root = NULL;

	qsort(posvits, (size_t)NbPartOri, sizeof(Part), qsort_partstr);
	if( R_ori > 0.0 )
	{
		printf  ("Nombre de particules considéré lors des traitements : %d (au lieu de %d)\nParticule la plus lointaine : %g\n",
				NbPart,
				NbPartOri,
				posvits[NbPart-1].r
			);
	}
	else
		printf("Nouveau Maximum :: %g\n", posvits[NbPart-1].r);

	rmax        = posvits[NbPart-1].r;
	taille      = 2.0 * posvits[NbPart-1].r;
	root        = Tree_Init(NbPart, 0.0, 0.0, 0.0, taille/*2.0 * posvits[NbPart-1].r*/);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);

	root->first = posvits;

	Tree_Build2(root, NbPart, NbMin);
	Save_Part("Particule-it.dat", posvits, NbPart);

	double g_ratio = 0.0,
	       p_ratio = 0.0;

	Axial_ratio(root, taille/2.0/*R_ori*/, &g_ratio, &p_ratio);

	/********************************************************************************************************************************\
	 *					Allocation de nos histogrammes :							*
	\********************************************************************************************************************************/
	// Traitement des données sous forme d'histogramme :
	// 	(x) bin de taille constante sur r.

	// Variable pour les histogrammes :
	double  //rmax        = taille/2.0, //R_ori,
		  dr        = rmax / (double)(nb_bin),
		  Ec	    = 0.0,
		  Ep	    = 0.0,
		  Tmoy      = 0.0,
		  SAniso    = 0.0,
	         *rayon     = NULL,
		 *masse     = NULL,
		 *densite   = NULL,
		 *disp      = NULL,
		 *Jac       = NULL,
		 *distrib   = NULL,
		 *energie_c = NULL,
		 *energie_t = NULL,
		 *Deltatemp = NULL,
		 *Aniso     = NULL,
		**potentiel = NULL,
		**LogDens   = NULL;

	rayon               = double1d(NbPart);
	energie_c           = double1d(NbPart);
	energie_t           = double1d(NbPart);

	disp                = double1d(nb_bin);
	distrib             = double1d(nb_bin);

	//fprintf(stderr, "\033[35mRMAX :: %g (%g); dr :: %g ; nb_bin :: %d\033[00m\n", rmax, taille, dr, nb_bin);

	/********************************************************************************************************************************\
	 *						Calcul du potentiel								*
	\********************************************************************************************************************************/
#ifdef USE_TIMER
	start = clock();
	time(&t1);
#endif
	potentiel           = CalcPotentiel(root, theta, rsoft, header.BoxSize); // * PosFact); //3.086e16 );
#ifdef USE_TIMER
	finish = clock();
	time(&t2);
	fprintf(stderr, "\033[32mTemps d'exécution de la fonction potentiel :: \033[33m%f (%.3f) secondes\033[00m\n", difftime(t2,t1), (double)(finish - start) / (double)CLOCKS_PER_SEC);
#endif

	/********************************************************************************************************************************\
	 *				Calcul des fonctions de masses, densités et énergies						*
	\********************************************************************************************************************************/
	for(int i = 0; i<NbPart; i++)
		rayon[i]    = posvits[i].r;

	qsort(rayon, (size_t)NbPart, sizeof(double), trie_rayon);

	// Calcule de la masse :
	masse               = CalcMasse(root);
	// Calcul de la densité :
	densite             = CalcDensite(root, nb_bin, dr, rmax);
	//densite             = Dens(root, nb_bin, dr, rmax);
	LogDens             = CalcLogDensite(root, nb_bin, rayon[0], rmax, rayon[NbPart/2]);
	// Calcul de la température :
	Deltatemp           = CalcTemperature(root, nb_bin, dr, rmax, &Tmoy);
	// Calcul des énergies :
	CalcEnergie(root, energie_c, energie_t, (const double**)potentiel, &Ec, &Ep);

	double Emax         = energie_t[maxlocdouble1d(energie_t, NbPart - nb_hors_Ro)],
	       Emin         = energie_t[minlocdouble1d(energie_t, NbPart - nb_hors_Ro)],
	       dE           = (Emax - Emin) / (double)(nb_bin);

	//fprintf(stderr, "\033[31mEmax :: %g ; Emin :: %g\033[00m\n", Emax, Emin);
	Jac                 = CalcJacobien(root, nb_bin, energie_t, (const double**)potentiel, Emin, Emax, dE, G, distrib);

	Aniso               = CalcAnisotropie(root, nb_bin, dr, rmax, &SAniso);

	printf("\033[36mAxial Ratio (p, g)  : %g, %g\nRapport du Viriel   : %g\nTempérature moyenne : %g\nAnisotropie de l'objet : %g\n\033[00m",
			p_ratio,
			g_ratio,
			2.0*Ec/Ep,
			Tmoy,
			SAniso
		);

	/********************************************************************************************************************************\
	 *						Enregistrement des données							*
	\********************************************************************************************************************************/
	printf("\033[32mÉcriture du fichier : %s\n", "Masse.dat");
	if( (fich      = fopen("Masse.dat", "w")) == NULL )
	{
		perror("fopen a échouer ");
		fprintf(stderr, "\033[31m%s::%s::%d ==> Impossible d'ouvrir le fichier : %s\n", __FILE__, __func__, __LINE__, "Masse.dat");
		Part1d_libere(posvits);
		Tree_Free(root);
		double1d_libere(masse);
		double1d_libere(rayon);
		double1d_libere(Aniso);
		double1d_libere(densite);
		double2d_libere(potentiel);
		double2d_libere(LogDens);
		double1d_libere(energie_t);
		double1d_libere(energie_c);
		double1d_libere(disp);
		double1d_libere(Jac);
		double1d_libere(distrib);
		double1d_libere(Deltatemp);
		exit(EXIT_FAILURE);
	}
	fprintf(fich, "#Rayon:	Masse:\n");
	for(int i=0; i<NbPart; i++) fprintf(fich, "%.16g %.16g\n", rayon[i], masse[i]); //, potentiel[i][1]); //, densite[i]);
	fclose(fich);

	printf("Écriture du fichier : %s\n", "Densite.dat");
	if( (fich      = fopen("Densite.dat", "w")) == NULL )
	{
		perror("fopen a échouer ");
		fprintf(stderr, "\033[31m%s::%s::%d ==> Impossible d'ouvrir le fichier : %s\n", __FILE__, __func__, __LINE__, "Densite.dat");
		Part1d_libere(posvits);
		Tree_Free(root);
		double1d_libere(masse);
		double2d_libere(LogDens);
		double1d_libere(rayon);
		double1d_libere(Aniso);
		double1d_libere(densite);
		double2d_libere(potentiel);
		double1d_libere(energie_t);
		double1d_libere(energie_c);
		double1d_libere(disp);
		double1d_libere(Jac);
		double1d_libere(distrib);
		double1d_libere(Deltatemp);
		exit(EXIT_FAILURE);
	}
	fprintf(fich, "#Rayon (bin, bords gauche):	densité:	T(r):	Anisotropie(r):\n");
	for(int i=0; i<nb_bin; i++) fprintf(fich, "%.16g %.16g %.16g %.16g %.16g %.16g\n", (i+1.0)*dr, densite[i], Deltatemp[i], Aniso[i], LogDens[i][0], LogDens[i][1]);
	fclose(fich);

	printf("Écriture du fichier : %s\n", "Potentiel-tc.dat");
	if( (fich      = fopen("Potentiel-tc.dat", "w")) == NULL )
	{
		perror("fopen a échouer ");
		fprintf(stderr, "\033[31m%s::%s::%d ==> Impossible d'ouvrir le fichier : %s\n", __FILE__, __func__, __LINE__, "Potentiel-tc.dat");
		Part1d_libere(posvits);
		Tree_Free(root);
		double1d_libere(masse);
		double1d_libere(rayon);
		double2d_libere(LogDens);
		double1d_libere(Aniso);
		double1d_libere(densite);
		double2d_libere(potentiel);
		double1d_libere(energie_t);
		double1d_libere(energie_c);
		double1d_libere(disp);
		double1d_libere(Jac);
		double1d_libere(distrib);
		double1d_libere(Deltatemp);
		exit(EXIT_FAILURE);
	}
	fprintf(fich, "#Rayon:	Potentiel:\n");
	for(int i=0; i<NbPart; i++) fprintf(fich, "%.16g %.16g\n", potentiel[i][0], potentiel[i][1]);
	fclose(fich);

	printf("Écriture du fichier : %s\n", "Distribution.dat");
	if( (fich      = fopen("Distribution.dat", "w")) == NULL )
	{
		perror("fopen a échouer ");
		fprintf(stderr, "\033[31m%s::%s::%d ==> Impossible d'ouvrir le fichier : %s\n", __FILE__, __func__, __LINE__, "Distribution.dat");
		Part1d_libere(posvits);
		Tree_Free(root);
		double1d_libere(masse);
		double1d_libere(rayon);
		double1d_libere(Aniso);
		double2d_libere(LogDens);
		double1d_libere(densite);
		double2d_libere(potentiel);
		double1d_libere(energie_t);
		double1d_libere(energie_c);
		double1d_libere(disp);
		double1d_libere(Jac);
		double1d_libere(distrib);
		double1d_libere(Deltatemp);
		exit(EXIT_FAILURE);
	}
	fprintf(fich, "#Énergie (bin, bord gauche):	Jacobien:	distribution normalisé:\n");
	for(int i=0; i<nb_bin; i++) fprintf(fich, "%.16g %.16g %.16g\n", Emin + (i+1.0)*dE, Jac[i], /*posvits[0].m / dE * */ distrib[i]);
	fclose(fich);

	printf("Écriture du fichier : %s\n", "Energie.dat");
	if( (fich      = fopen("Energie.dat", "w")) == NULL )
	{
		perror("fopen a échouer ");
		fprintf(stderr, "\033[31m%s::%s::%d ==> Impossible d'ouvrir le fichier : %s\n", __FILE__, __func__, __LINE__, "Energie.dat");
		Part1d_libere(posvits);
		Tree_Free(root);
		double1d_libere(masse);
		double1d_libere(rayon);
		double2d_libere(LogDens);
		double1d_libere(Aniso);
		double1d_libere(densite);
		double2d_libere(potentiel);
		double1d_libere(energie_t);
		double1d_libere(energie_c);
		double1d_libere(disp);
		double1d_libere(Jac);
		double1d_libere(distrib);
		double1d_libere(Deltatemp);
		exit(EXIT_FAILURE);
	}
	for(int i=0; i<NbPart; i++) fprintf(fich, "%.16g %.16g %.16g %.16g\n", posvits[i].r, energie_c[i], energie_t[i], potentiel[i][1]);
	fclose(fich);

	printf("Écriture du fichier : %s\033[00m\n", timeparam);
	if( (fich      = fopen(timeparam, "a")) == NULL )
	{
		perror("fopen a échouer ");
		fprintf(stderr, "\033[31m%s::%s::%d ==> Impossible d'ouvrir le fichier : %s\n", __FILE__, __func__, __LINE__, "TimeParam.dat");
		Part1d_libere(posvits);
		Tree_Free(root);
		double1d_libere(masse);
		double1d_libere(rayon);
		double1d_libere(Aniso);
		double2d_libere(LogDens);
		double1d_libere(densite);
		double2d_libere(potentiel);
		double1d_libere(energie_t);
		double1d_libere(energie_c);
		double1d_libere(disp);
		double1d_libere(Jac);
		double1d_libere(distrib);
		double1d_libere(Deltatemp);
		exit(EXIT_FAILURE);
	}
	fprintf(fich, "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n", simu_time, p_ratio, g_ratio, 2.0*Ec/Ep, Tmoy, SAniso, rayon[(int)(NbPart * 0.1)], rayon[(int)(NbPart/2.0)], rayon[(int)(NbPart * 0.9)], TotMove.x, TotMove.y, TotMove.z, TotMove.vx, TotMove.vy, TotMove.vz);
	fclose(fich);

	Tree_Free(root);
	Part1d_libere(posvits);
	double1d_libere(Jac);
	double1d_libere(disp);
	double1d_libere(masse);
	double1d_libere(rayon);
	double1d_libere(Aniso);
	double1d_libere(distrib);
	double1d_libere(densite);
	double2d_libere(LogDens);
	double2d_libere(potentiel);
	double1d_libere(energie_t);
	double1d_libere(energie_c);
	double1d_libere(Deltatemp);

	return EXIT_SUCCESS;
}

////#define NEAREST(x, boxhalf, boxsize) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))
//Part PeriodicGravityCenter(const Part *posvits, const int NbPart, const double BoxSize)
//{
//	Part center = {	  .x  = BoxSize/2.0, //0.0,
//			  .y  = BoxSize/2.0, //0.0,
//			  .z  = BoxSize/2.0, //0.0,
//			  .r  = 0.0,
//			  .vx = 0.0,
//			  .vy = 0.0,
//			  .vz = 0.0,
//			  .v  = 0.0};
//
//	for (int i = 0; i < NbPart; i++) {
//		double dx = posvits[i].x - BoxSize/2.0,
//		       dy = posvits[i].y - BoxSize/2.0,
//		       dz = posvits[i].z - BoxSize/2.0;
//
//		center.x += NEAREST(dx, BoxSize/2.0, BoxSize);
//		center.y += NEAREST(dy, BoxSize/2.0, BoxSize);
//		center.z += NEAREST(dz, BoxSize/2.0, BoxSize);
//	}
//
//	center.x /= NbPart;
//	center.y /= NbPart;
//	center.z /= NbPart;
//
//	printf("%g %g %g\n", center.x, center.y, center.z);
//
//	return center;
//}
////#undef NEAREST



//int trie_Pot(const void *a, const void *b)
//{
////	const double **p1 = (const double**)a,
////	             **p2 = (const double**)b;
////	double * const *p1 = a;
////	double * const *p2 = b;
//	const double *p1 = *(const double **)a,
//	             *p2 = *(const double **)b;
//
//	if( p1[0] > p2[0] )//sqrt(p1->x*p1->x + p1->y*p1->y) > sqrt(p2->x*p2->x + p2->y*p2->y))
//		return 1;
//	else if( p1[0] < p2[0] )
//		return -1;
//	else
//		return 0;
//}
//
//int qsort_comp(const void *a, const void *b)
//{
//	const double *p1 = *(const double **)a,
//	             *p2 = *(const double **)b;
//
//	if(p1[0] > p2[0])
//		return 1;
//	else
//		return -1;
//}



//
//	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);
//
//	NbPartOri = NbPart;
//	printf("Chargement de %d particule de type %d effectué.\n", NbPart, type);
//
//	Save_Part("Particule-BeCorr.dat", posvits, NbPart);
//
//	/****************************************************************\
//	 *		Correction des conditions périodiques		*
//	\****************************************************************/
//	taille = 2.0  * posvits[NbPart-1].r;
//	root   = Tree_Init(NbPart, header.BoxSize/2.0, header.BoxSize/2.0, header.BoxSize/2.0, taille);
//	if( root == NULL )
//		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);
//
//	root->first = posvits;
//
//	//Center = GravityCenter(root, NbVois);
//	//Center = PeriodicGravityCenter(posvits, NbPart, header.BoxSize*3.086e16);
//	Center = MoreDenseParticule(root, NbVois);
//
//	double Bord = 1.0 * header.BoxSize;
//
//	for (int i = 0; i < NbPart; i++) {
//		posvits[i].x -= Center.x;
//		posvits[i].y -= Center.y;
//		posvits[i].z -= Center.z;
//
//		if( posvits[i].x > Bord*3.086e16 )
//			posvits[i].x -= header.BoxSize*3.086e16;
//		else if( posvits[i].x < -Bord*3.086e16 )
//			posvits[i].x += header.BoxSize*3.086e16;
//
//		if( posvits[i].y > Bord*3.086e16 )
//			posvits[i].y -= header.BoxSize*3.086e16;
//		else if( posvits[i].y < -Bord*3.086e16 )
//			posvits[i].y += header.BoxSize*3.086e16;
//
//		if( posvits[i].z > Bord*3.086e16 )
//			posvits[i].z -= header.BoxSize*3.086e16;
//		else if( posvits[i].z < -Bord*3.086e16 )
//			posvits[i].z += header.BoxSize*3.086e16;
//
//		posvits[i].r  = sqrt( pow(posvits[i].x, 2.0) + pow(posvits[i].y, 2.0) + pow(posvits[i].z, 2.0) );
//	}
//
//	TotMove = Part_add(Center, TotMove);
//
//	Save_Part("Particule-Classe.dat", posvits, NbPart);
//	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);
//
//	taille      = 2.0 * posvits[NbPart-1].r;
//	rmax        = posvits[NbPart-1].r;
//	root        = Tree_Init(NbPart, 0.0, 0.0, 0.0, taille);
//	if( root == NULL )
//		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);
//
//	root->first = posvits;
//	Tree_Build2(root, NbPart, NbMin);



/*{
	Part center = {	  .x  = 0.0,
			  .y  = 0.0,
			  .z  = 0.0,
			  .r  = 0.0,
			  .vx = 0.0,
			  .vy = 0.0,
			  .vz = 0.0,
			  .v  = 0.0},
	     *Vois  = NULL;
	double rholoc, dens = 0.0;

	Vois = Part1d(NbVois);

	for (int i = 0; i < root->N; i++)
	{
		//On initialise le tableau des voisins :
		for (int j = 0, k = 0; j < NbVois && k < root->N; j++,k++)
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
			if( k == root->N -1)
				k = 0;
		}

		qsort(Vois, (size_t)NbVois, sizeof(Part), qsort_partstr);
		//Calcul des voisins :
		Tree_Voisin(root, Vois, NbVois, &root->first[i]);

		//Calcul de la densité locale :
		rholoc     = (NbVois - 1) * root->first[0].m / ( (4.0/3.0) * acos(-1.0) * Vois[NbVois -1].r * Vois[NbVois -1].r * Vois[NbVois -1].r );

		if( rholoc > dens )
		{
			dens = rholoc;
			center = root->first[i];
		}
#ifndef __DENSITYCENTER_NOPROGRESS_P
		fprintf(stderr, "\r\033[31m%s:: %03.3f%%\033[00m", __func__, ( (double)i + 1.0)/( (double)root->N ) * 100.0);
#endif
	}
#ifndef __DENSITYCENTER_NOPROGRESS_P
	fputs("\n", stderr);
#endif

	Part1d_libere(Vois);

	return center;
}*/

