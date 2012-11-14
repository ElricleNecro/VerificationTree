/****************************************************************************************\
 *	But    : Obtenir le potentiel par brute force (calcul en N^2) afin de           *
 *	comparaison avec le tree code.							*
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
 * TODO2 : * Facteur d'anisotropie (1 +/- 2*sig^2_r / sig^2_t)				*
 *	   * Température contre rayon							*
 *	   * Fonction de distribution f(E)						*
 *	   * Rapport du Viriel								*
\****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "tree.h"
#include "utils.h"
#include "types.h"
#include "cte_phys.h"
#include "iogadget2.h"
#include "Verif_tools.h"

// Flags de debug : 1 activé, 0 sinon
#define __DBG_SORT_MYSORT_P 0
#define __DBG_SORT__QSORT_P 0

int trie_rayon(const void *a, const void *b)
{
	double *p1 = (double*)a,
	       *p2 = (double*)b;

	if( *p1 > *p2 )//sqrt(p1->x*p1->x + p1->y*p1->y) > sqrt(p2->x*p2->x + p2->y*p2->y))
		return 1;
	else if( *p1 < *p2 )
		return -1;
	else
		return 0;
}

int trie_Pot(const void *a, const void *b)
{
//	const double **p1 = (const double**)a,
//	             **p2 = (const double**)b;
//	double * const *p1 = a;
//	double * const *p2 = b;
	double *p1 = *(double **)a,
	       *p2 = *(double **)b;

	if( p1[0] > p2[0] )//sqrt(p1->x*p1->x + p1->y*p1->y) > sqrt(p2->x*p2->x + p2->y*p2->y))
		return 1;
	else if( p1[0] < p2[0] )
		return -1;
	else
		return 0;
}
int qsort_comp(const void *a, const void *b)
{
	double *p1 = *(double **)a,
	       *p2 = *(double **)b;

	if(p1[0] > p2[0])
		return 1;
	else
		return -1;
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
	{
		fprintf(fich, "%g %g %g %g %g %g %g %g %d\n", posvits[i].x, posvits[i].y, posvits[i].z, posvits[i].r, posvits[i].vx, posvits[i].vy, posvits[i].vz, posvits[i].v, posvits[i].id);
	}
	fclose(fich);
}

int main(int argc, char **argv)
{
	if( argc < 5 )
	{
		fprintf(stderr, "\033[31mUsage : (%s) (fichier de données) (Valeur de G) (Rayon d'origine) (softening en pc) [(Nombre de bin)]\n%d arguments au lieu de %d\033[00m\n", argv[0], argc - 1, 4);
		return EXIT_FAILURE;
	}

	/****************************************************************************\
	 *		Chargement des données et allocation :			    *
	\****************************************************************************/
	int     NbPart = 0,
		NbMin  = 20,
		NbVois = 20;
	double  G      = atof(argv[2]),
		R_ori  = atof(argv[3]),
		taille = 0.0,
		theta  = 0.0, //0.05,
		rsoft  =  atof(argv[4]) * 3.086e16; //0.0; // * 3.086e16;
	FILE *fich     = NULL;

	Part  *posvits = NULL,
	       Center  = {.x  = 0.0,
		          .y  = 0.0,
			  .z  = 0.0,
			  .r  = 0.0,
			  .vx = 0.0,
			  .vy = 0.0,
			  .vz = 0.0,
			  .v  = 0.0};
	TNoeud  root   = NULL;

	posvits        = Deconversion(argv[1], &NbPart);

	if( posvits == NULL )
		fprintf(stderr, "Erreur avec le tableau posvits !!!\n"),exit(EXIT_FAILURE);
	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);

	Tree_var(1, 50);
	taille      = 2.0 * posvits[NbPart-1].r;
	root        = Tree_Init(NbPart, 0.0, 0.0, 0.0, taille);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);

	root->first = posvits;
	Tree_Build2(root, NbPart, NbMin);

	Tree_Build2(root, NbPart, NbMin);

	//Center = GravityCenter(root, NbVois);
	Center = DensityCenter(root, NbVois);
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
	}
	Tree_Free(root), root = NULL;

	root        = Tree_Init(NbPart, 0.0, 0.0, 0.0, taille/*2.0 * posvits[NbPart-1].r*/);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);

	root->first = posvits;

	Tree_Build2(root, NbPart, NbMin);

	double 	**potentiel = NULL;
	potentiel  = double2d(root->N, 2);

	for (int i = 0; i < NbPart; i++)
	{
		potentiel[i][0] = posvits[i].r;
		potentiel[i][1] = 0.0;
		for (int j = 0; j < NbPart; j++)
		{
			if( i != j)
				potentiel[i][1] += -G * posvits[j].m / (
						sqrt( pow(posvits[j].x - posvits[i].x, 2.0) + pow(posvits[j].y - posvits[i].y, 2.0) + pow(posvits[j].z - posvits[i].z, 2.0) )
						+ 0.0);

		}
		fprintf(stderr, "\r\033[31mCalcul du potentiel :: %03.3f%%\033[00m", ( (double)i + 1.0)/( (double)root->N ) * 100.0);
	}
	fputs("\n", stderr);

	printf("Écriture du fichier : %s\n", "Potentiel-bf.dat");
	if( (fich      = fopen("Potentiel-bf.dat", "w")) == NULL )
	{
		perror("fopen a échouer ");
		fprintf(stderr, "\033[31m%s::%s::%d ==> Impossible d'ouvrir le fichier : %s\n", __FILE__, __func__, __LINE__, "Potentiel-bf.dat");
		Part1d_libere(posvits);
		Tree_Free(root);
		double2d_libere(potentiel);
		exit(EXIT_FAILURE);
	}
	for(int i=0; i<NbPart; i++) fprintf(fich, "%.16g %.16g\n", potentiel[i][0], potentiel[i][1]);
	fclose(fich);

	Part1d_libere(posvits);
	Tree_Free(root);

	double2d_libere(potentiel);

	return EXIT_SUCCESS;
}

