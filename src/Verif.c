/****************************************************************************************\
 *	But    : Obtenir les densité et fonction de distribution à partir de position   *
 *	         et vitesse de particule.						*
 *	Entrée : Fichier snapshot Gadget2						*
 *	Sortie : 3 fichiers : densité(r), potntiel(r) et distribution(E)		*
 * Auteur      : Guillaume Plum								*
 * Date        : Vendredi 13 Mai 2011							*
 * Cadre       : Stage de M2 sur la stabilité des systèmes auto-gravitants		*
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
#include <stdbool.h>

#ifdef USE_TIMER
#include <time.h>
#endif

#ifdef USE_MCHECK
#include <mcheck.h>
#endif

#ifdef USE_SQLITE3
#include <sqlite3.h>
#endif

#ifdef USE_HDF5
#include "HDF5.h"
#endif

#include "tree.h"
#include "utils.h"
#include "types.h"
#include "cte_phys.h"

#include "io_snapshot.h"

#include "Verif_tools.h"

// Flags de debug : 1 activé, 0 sinon
#define __DBG_SORT_MYSORT_P 0
#define __DBG_SORT__QSORT_P 0

TNoeud Create_Tree(Part *posvits, const int NbPart, const int NbMin, const Part Center, const double taille)
{
	TNoeud root = NULL;
	root        = Tree_Init(NbPart, Center.x, Center.y, Center.z, taille);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);
	root->first = posvits;

	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);

	Tree_Build2(root, NbPart, NbMin);

	return root;
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

void Save_Part(const char * const fn, const Part * const posvits, const int NbPart)
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
		"\t\033[35m-G \033[04mRéel\033[00m   : Valeur de G dans les unités du fichier (\033[36;03mdéfaut : 6.67384e-11\033[00m).\n"
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
		"\t\033[35m--timeparam \033[00m: Nom du fichier de paramètre.\n"
		"\t\033[35m--posfact   \033[00m: Unités dans lesquelles convertir les positions (\033[36;03mdéfaut : 3.086e16 (pc->m)\033[00m]]).\n"
		"\t\033[35m--vitfact   \033[00m: Unités dans lesquelles convertir les vitesses (\033[36;03mdéfaut : 1.0\033[00m]]).\n"
#ifdef USE_SQLITE3
		"\t\033[35m--database  \033[00m: Nom de la base de donnée à utiliser.\n"
#endif
#ifdef USE_HDF5
		"\t\033[35m--hdf5  \033[00m: Nom du fichier hdf5 à utiliser.\n"
#endif
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
	bool    gc         = false/*,
		periodic   = false*/;
	char   *filename   = NULL,
	       *timeparam  = NULL
#ifdef USE_HDF5
             , *hdf5file   = NULL
	     , *w_ext      = NULL
#endif
#ifdef USE_SQLITE3
	     , *database   = NULL
#endif
	       ;
	double  G          = 6.67384e-11,
		R_ori      = -1.0,
		taille     = 0.0,
		theta      = -1.0,
		simu_time  = 0.0,
		rmax       = 0.0,
		rsoft      = 0.0,
		PosFact    = 3.086e16,
		VitFact    = 1.0,
		r_norm     = 1.0;
	time_t  t1, t2;
	clock_t start, finish;

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
//				case 'p':
//					periodic = true;
//					break;
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
						if( (timeparam = malloc((strlen(argv[i+1]) + 1)*sizeof(char))) == NULL )
						{
							perror("Memory Problem");
							exit(EXIT_FAILURE);
						}
						strncpy(timeparam, argv[i+1], strlen(argv[i+1]));
						i++;
					}
#ifdef USE_HDF5
					else if( !strcmp("--hdf5", argv[i]) )
					{
						hdf5file = argv[i+1];
						i++;
					}
#endif
#ifdef USE_SQLITE3
					else if( !strcmp("--database", argv[i]) )
					{
						if( (database = malloc((strlen(argv[i+1]) + 1)*sizeof(char))) == NULL )
						{
							perror("Memory Problem");
							exit(EXIT_FAILURE);
						}
						strncpy(database, argv[i+1], strlen(argv[i+1]));
						i++;
					}
#endif
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
	if( timeparam == NULL )
	{
		if( (timeparam = malloc(20*sizeof(char))) == NULL )
		{
			perror("Memory Problem");
			exit(EXIT_FAILURE);
		}
		strncpy(timeparam, "TimeParam.dat", 20);
	}
#ifdef USE_HDF5
	if( hdf5file == NULL )
		hdf5file = "simu.h5";
#endif
#ifdef USE_SQLITE3
	if( database == NULL )
	{
		if( (database = malloc(20*sizeof(char))) == NULL )
		{
			perror("Memory Problem");
			exit(EXIT_FAILURE);
		}
		strncpy(database, "simu.db", 20);
	}
#endif
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

	Tree_SetG(G);

	if( posvits == NULL )
		fprintf(stderr, "Erreur avec le tableau de particules !!!\n"),exit(EXIT_FAILURE);
	NbPartOri = NbPart;
	printf("Chargement de %d particule de type %d effectué.\n", NbPart, type);

	Save_Part("Particule-BeCorr.dat", posvits, NbPart);

	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);
#ifdef PERIODIC
#	ifdef USE_TIMER
	start = clock();
	time(&t1);
#	endif
	TotMove = ReCentre(root, posvits, NbPart, NbVois, NbMin, header.BoxSize); // * PosFact);
#	ifdef USE_TIMER
	finish = clock();
	time(&t2);
	fprintf(stderr, "\033[32mTemps d'exécution de la fonction %s :: \033[33m%f (%.3f) secondes\033[00m\n", "ReCentre", difftime(t2,t1), (double)(finish - start) / (double)CLOCKS_PER_SEC);
#	endif
	printf("\033[31m%g\t%g\t%g\n\033[00m", TotMove.x, TotMove.y, TotMove.z);
	qsort(posvits, (size_t)NbPart, sizeof(Part), qsort_partstr);
#endif

	rmax        = posvits[NbPart-1].r;
	taille      = /*header.BoxSize; / */ 2.0 * posvits[NbPart-1].r;
	root        = Tree_Init(NbPart, 0.0, 0.0, 0.0, taille);
	if( root == NULL )
		fprintf(stderr, "Erreur avec Tree_Init !!!\n"),exit(EXIT_FAILURE);

#ifdef DBG_NEWWAY
	printf("Particule la plus lointaine : %g (%g, %g, %g, %d).\n", taille, posvits[NbPart-1].x, posvits[NbPart-1].y, posvits[NbPart-1].z, posvits[NbPart-1].id);
#endif

	Save_Part("Particule-Classe.dat", posvits, NbPart);

	root->first = posvits;
	Tree_Build2(root, NbPart, NbMin);

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
#ifdef PERIODIC
		Center = DensityCenter(root, NbVois, header.BoxSize); // * PosFact); //3.086e16 );
#else
		Center = DensityCenter(root, NbVois); // * PosFact); //3.086e16 );
#endif
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

	printf("\033[36mNombre de particule : %d de masse : %g, La plus lointaine : %g (rayon maximum = %g, taille de la boîte = %g)\nsoftening : %g, critère : %g\nNombre de particule par feuille de Tree Code : %d, Nombre de voisin à chercher : %d\033[00m\n",
					NbPart,		posvits[0].m,		  rmax,		R_ori,			taille,		rsoft,		theta,						NbMin,				NbVois);
	if( R_ori > 0.0 )
		printf("\033[36mTaille du Cut-Off : %g\033[00m\n", R_ori);

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
#ifdef PERIODIC
	potentiel           = CalcPotentiel(root, theta, rsoft, header.BoxSize); // * PosFact); //3.086e16 );
#else
	potentiel           = CalcPotentiel(root, theta, rsoft); // * PosFact); //3.086e16 );
#endif
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
	r_norm              = 1.0; //rayon[NbPart/2];
	LogDens             = CalcLogDensite(root, nb_bin, rayon[0] - 0.01*rayon[0], rmax, r_norm/*1.0*/);
	//printf("TOTO : %g\n", rayon[NbPart/2]);
	// Calcul de la température :
	Deltatemp           = CalcTemperature(root, densite, nb_bin, dr, &Tmoy);
	// Calcul des énergies :
	CalcEnergie(root, energie_c, energie_t, potentiel, &Ec, &Ep);

	double Emax         = energie_t[maxlocdouble1d(energie_t, NbPart - nb_hors_Ro)],
	       Emin         = energie_t[minlocdouble1d(energie_t, NbPart - nb_hors_Ro)],
	       dE           = (Emax - Emin) / (double)(nb_bin);

	//fprintf(stderr, "\033[31mEmax :: %g ; Emin :: %g\033[00m\n", Emax, Emin);
	Jac                 = CalcJacobien(root, nb_bin, energie_t, potentiel, Emin, Emax, dE, distrib);

	Aniso               = CalcAnisotropie(root, nb_bin, dr, &SAniso);

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
#ifdef USE_SQLITE3
	char create[][1024] = { "CREATE TABLE IF NOT EXISTS id_simu (nom TEXT PRIMARY KEY, id INT, time REAL)",
				"CREATE TABLE IF NOT EXISTS masse (id INT, type INT, r REAL, m REAL)",
				"CREATE TABLE IF NOT EXISTS densite (id INT, type INT, bin_rg REAL, rho REAL, t REAL, aniso REAL)",
				"CREATE TABLE IF NOT EXISTS densite_log (id INT, type INT, bin_rg REAL, l_densite REAL)",
				"CREATE TABLE IF NOT EXISTS distribution (id INT, type INT, e REAL, distribution REAL)",
				"CREATE TABLE IF NOT EXISTS energie (id INT, type INT, r REAL, ec REAL, epot REAL, etot REAL)",
				"CREATE TABLE IF NOT EXISTS potentiel (id INT, type INT, r REAL, pot REAL)",
				"CREATE TABLE IF NOT EXISTS Movement (id INT, type INT, x REAL, y REAL, z REAL, vx REAL, vy REAL, vz REAL)",
				"CREATE TABLE IF NOT EXISTS timeparam (id INT, type INT, time REAL, p_ratio REAL, g_ratio REAL, viriel REAL, Ec REAL, Ep REAL, tmoy REAL, aniso REAL, r10 REAL, r50 REAL, r90 REAL, x REAL, y REAL, z REAL, vx REAL, vy REAL, vz REAL)",
				//"CREATE TABLE IF NOT EXISTS particule (id INT, type INT, x REAL, y REAL, z REAL, vx REAL, vy REAL, vz REAL)",
				"BEGIN TRANSACTION",
				},
	     delete[][1024] = { "DELETE FROM id_simu WHERE id=%d",
				"DELETE FROM masse WHERE id=%d AND type=%d",
				"DELETE FROM densite WHERE id=%d AND type=%d",
				"DELETE FROM densite_log WHERE id=%d AND type=%d",
				"DELETE FROM distribution WHERE id=%d AND type=%d",
				"DELETE FROM energie WHERE id=%d AND type=%d",
				"DELETE FROM potentiel WHERE id=%d AND type=%d",
				"DELETE FROM Movement WHERE id=%d AND type=%d",
				"DELETE FROM timeparam WHERE id=%d AND type=%d",
				},
	     tampon[1024] = {0};
	int nb_table      = sizeof(create)/sizeof(create[0]),
	    id            = get_id(filename);
	sqlite3    *conn  = NULL;

	if( sqlite3_open(database, &conn) )
	{
		fprintf(stderr, "%s::%s::%d ==> Database %s connection failed.\n", __FILE__, __func__, __LINE__, database);
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
	sqlite3_busy_timeout(conn, 60000);

	/*
	 * On crée les tables n'existant pas :
	 */
	for(int i = 0; i < nb_table; i++)
		sqlite3_exec(conn, create[i], NULL, NULL, NULL);

	/*
	 * On efface les données correspondantes à l'ID du snapshot :
	 */
	for(unsigned int i = 0; i < sizeof(delete)/sizeof(delete[0]); i++)
	{
		snprintf(tampon, 1024*sizeof(char), delete[i], id, type);
		sqlite3_exec(conn, tampon, NULL, NULL, NULL);
	}

	/*
	 * Insertion des données dans la base SQLite :
	 */
	snprintf(tampon, 1024*sizeof(char), "INSERT INTO %s VALUES(\"%s\", %d, %g)", "id_simu", filename, id, simu_time);
	//printf("::%s::\n", tampon);
	sqlite3_exec(conn, tampon, NULL, NULL, NULL);

	snprintf(tampon, 1024*sizeof(char), "INSERT INTO %s VALUES(%d, %d, %.14g, %.14g, %.14g, %.14g, %.14g, %.14g)", "Movement", id, type, TotMove.x, TotMove.y, TotMove.z, TotMove.vx, TotMove.vy, TotMove.vz);
	/*printf("::%s::\n", tampon);*/
	sqlite3_exec(conn, tampon, NULL, NULL, NULL);

	snprintf(tampon, 1024*sizeof(char), "INSERT INTO %s VALUES(%d, %d, %.14g, %.14g, %.14g, %.14g, %.14g, %.14g)", "Movement", id, type, TotMove.x, TotMove.y, TotMove.z, TotMove.vx, TotMove.vy, TotMove.vz);
	sqlite3_exec(conn, tampon, NULL, NULL, NULL);

	snprintf(tampon, 1024*sizeof(char), "INSERT INTO %s VALUES (%d, %d, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g)",
			 "timeparam",
			 id,
			 type,
			 simu_time, p_ratio, g_ratio, 2.0*Ec/Ep, Ec, Ep, Tmoy, SAniso,
			 rayon[(int)(NbPart * 0.1)], rayon[(int)(NbPart/2.0)], rayon[(int)(NbPart * 0.9)],
			 TotMove.x, TotMove.y, TotMove.z, TotMove.vx, TotMove.vy, TotMove.vz);
	sqlite3_exec(conn, tampon, NULL, NULL, NULL);

	for(int i = 0; i < NbPart; i++)
	{
		snprintf(tampon, 1024*sizeof(char), "INSERT INTO %s VALUES(%d, %d, %.14g, %.14g)",
			"masse", id, type, rayon[i], masse[i]);
		sqlite3_exec(conn, tampon, NULL, NULL, NULL);

		snprintf(tampon, 1024*sizeof(char), "INSERT INTO %s VALUES(%d, %d, %.14g, %.14g, %.14g, %.14g)",
			"energie", id, type, posvits[i].r, energie_c[i], potentiel[i][1], energie_t[i]);
		sqlite3_exec(conn, tampon, NULL, NULL, NULL);

		snprintf(tampon, 1024*sizeof(char), "INSERT INTO %s VALUES(%d, %d, %.14g, %.14g)",
			"potentiel", id, type, potentiel[i][0], potentiel[i][1]);
		sqlite3_exec(conn, tampon, NULL, NULL, NULL);

		//snprintf(tampon, 1024*sizeof(char), "INSERT INTO %s VALUES(%d, %d, %.14g, %.14g)", "particule", id, type, posvits[i].x, posvits[i].y, posvits[i].z, posvits[i].vx, posvits[i].vy, posvits[i].vz);
		//sqlite3_exec(conn, tampon, NULL, NULL, NULL);
	}

	for(int i = 0; i < nb_bin; i++)
	{
		snprintf(tampon, 1024*sizeof(char), "INSERT INTO %s VALUES(%d, %d, %.14g, %.14g, %.14g, %.14g)",
			"densite", id, type, (i+1.0)*dr, densite[i], Deltatemp[i], Aniso[i]);
		sqlite3_exec(conn, tampon, NULL, NULL, NULL);

		snprintf(tampon, 1024*sizeof(char), "INSERT INTO %s VALUES(%d, %d, %.14g, %.14g)",
			"densite_log", id, type, LogDens[i][0], LogDens[i][1]);
		sqlite3_exec(conn, tampon, NULL, NULL, NULL);

		snprintf(tampon, 1024*sizeof(char), "INSERT INTO %s VALUES(%d, %d, %.14g, %.14g)",
			"distribution", id, type, Emin + (i+1.0)*dE, distrib[i]);
		sqlite3_exec(conn, tampon, NULL, NULL, NULL);
	}

	sqlite3_exec(conn, "END TRANSACTION", NULL, NULL, NULL);

	sqlite3_close(conn);
#endif
#ifdef USE_HDF5
	int str_start = -1;
	double *tmp = NULL;
	double virtmp = 2.0*Ec/Ep;
	hid_t grp, file = H5Fcreate(hdf5file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hsize_t size[2] = {1, 1};

	for(int i = strlen(filename)-1; i >= 0 && filename[i] != '/'; i--)
	{
		str_start = i;
	}

	if( str_start < 0 || str_start >= strlen(filename) )
		str_start = 0;

	w_ext = remove_ext(&filename[str_start]);

	char       *tab = NULL;
	int         tN  = strlen(w_ext) + 6 + strlen("densite_log");
	tab             = malloc(tN*sizeof(char));

	snprintf(tab, tN, "/%s", w_ext);
	printf("Hdf5 :: %s\n", tab);
	grp             = H5Gcreate(file, tab, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	/****************************************************************************************\
	 *			Saving time dependant parameter 				*
	\****************************************************************************************/
	snprintf(tab, tN, "/%s/%s", w_ext, "timeparam");

	ExtensibleDataSet eds = CreateExtensibleDS(file, tab, size);
	ExtensibleDataSet_Extend(eds, &simu_time, size);
	ExtensibleDataSet_Extend(eds, &p_ratio, size);
	ExtensibleDataSet_Extend(eds, &g_ratio, size);
	ExtensibleDataSet_Extend(eds, &virtmp, size);
	ExtensibleDataSet_Extend(eds, &Ec, size);
	ExtensibleDataSet_Extend(eds, &Ep, size);
	ExtensibleDataSet_Extend(eds, &Tmoy, size);
	ExtensibleDataSet_Extend(eds, &SAniso, size);
	ExtensibleDataSet_Extend(eds, &rayon[(int)(NbPart * 0.1)], size);
	ExtensibleDataSet_Extend(eds, &rayon[(int)(NbPart * 0.5)], size);
	ExtensibleDataSet_Extend(eds, &rayon[(int)(NbPart * 0.9)], size);
	ExtensibleDataSet_Extend(eds, &TotMove.x, size);
	ExtensibleDataSet_Extend(eds, &TotMove.y, size);
	ExtensibleDataSet_Extend(eds, &TotMove.z, size);
	ExtensibleDataSet_Extend(eds, &TotMove.vx, size);
	ExtensibleDataSet_Extend(eds, &TotMove.vy, size);
	ExtensibleDataSet_Extend(eds, &TotMove.vz, size);
	ExtensibleDataSet_Close(eds);

	/****************************************************************************************\
	 *			Saving density profile + other quantity				*
	\****************************************************************************************/
	size[0] = nb_bin;
	snprintf(tab, tN, "/%s/%s", w_ext, "densite");

	eds = CreateExtensibleDS(file, tab, size);
	if( (tmp = double1d(nb_bin)) == NULL )
	{
		perror("Allocation problem:");
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < nb_bin; i++)
	{
		tmp[i] = (i+1.0)*dr;
	}
	ExtensibleDataSet_Extend(eds, tmp, size);
	free(tmp);
	ExtensibleDataSet_Extend(eds, densite, size);
	ExtensibleDataSet_Extend(eds, Deltatemp, size);
	ExtensibleDataSet_Extend(eds, Aniso, size);
	ExtensibleDataSet_Close(eds);

	/****************************************************************************************\
	 *			Saving log density profile + other quantity			*
	\****************************************************************************************/
	size[1] = 2;
	snprintf(tab, tN, "/%s/%s", w_ext, "densite_log");

	eds = CreateExtensibleDS(file, tab, size);
	ExtensibleDataSet_Extend(eds, LogDens[0], size);
	ExtensibleDataSet_Close(eds);

	/****************************************************************************************\
	 *			Saving mass profile + other quantity				*
	\****************************************************************************************/
	size[0] = NbPart;
	size[1] = 1;
	snprintf(tab, tN, "/%s/%s", w_ext, "masse");

	eds = CreateExtensibleDS(file, tab, size);
	ExtensibleDataSet_Extend(eds, rayon, size);
	ExtensibleDataSet_Extend(eds, masse, size);
	ExtensibleDataSet_Close(eds);

	/****************************************************************************************\
	 *			Saving energy profile + other quantity				*
	\****************************************************************************************/
	snprintf(tab, tN, "/%s/%s", w_ext, "energie");

	eds = CreateExtensibleDS(file, tab, size);
	if( (tmp = double1d(NbPart)) == NULL )
	{
		perror("Allocation problem:");
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < NbPart; i++)
	{
		tmp[i] = posvits[i].r;
	}
	ExtensibleDataSet_Extend(eds, tmp, size);
	free(tmp);
	ExtensibleDataSet_Extend(eds, energie_c, size);
	ExtensibleDataSet_Extend(eds, energie_t, size);
	ExtensibleDataSet_Close(eds);

	/****************************************************************************************\
	 *			Saving potential profile + other quantity			*
	\****************************************************************************************/
	size[1] = 2;
	snprintf(tab, tN, "/%s/%s", w_ext, "potentiel");

	eds = CreateExtensibleDS(file, tab, size);
	ExtensibleDataSet_Extend(eds, potentiel[0], size);
	ExtensibleDataSet_Close(eds);

	H5Gclose(grp);
	H5Fclose(file);
	free(tab);
	free(w_ext);
#endif
#ifdef USE_FILE
	FILE   *fich   = NULL;

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
	fprintf(fich, "#Rayon (bin, bords gauche):	densité:	T(r):	Anisotropie(r):		Rayon (bin log cte) :	Densité (bin log cte) :\n");
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
#endif

	free(timeparam);
#ifdef USE_SQLITE3
	free(database);
#endif
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

	(void)G;
	return EXIT_SUCCESS;
}

