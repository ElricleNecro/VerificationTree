#ifndef TYPES_H

#define TYPES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double** PosVit;

//Flags utiles :
#define R      0
#define V      1
#define X      2
#define Y      3
#define Z      4
#define VX     5
#define VY     6
#define VZ     7
#define NB_col 8


/**
 * \struct Part
 * Contient toutes les informations relatives à une particule.
 * Sert pour le Tree-Code, mais aussi pour d'autre traitement.
 */
struct _part {
	unsigned int id;		/*< Identifiant de la particule.*/
	double x;			/*< Abscisse de la particule.*/
	double y;			/*< Ordonnée de la particule.*/
	double z;			/*< Côte de la particule.*/
	double r;			/*< Distance de la particule à un objet (ie : centre de l'amas, particule, ...).*/
	double vx;			/*< Vitesse selon x de la particule.*/
	double vy;			/*< Vitesse selon y de la particule.*/
	double vz;			/*< Vitesse selon z de la particule.*/
	double v;			/*< Vitesse de la particule.*/
	double m;			/*< Masse de la particule.*/
};
typedef struct _part Part;

/**
 * Fonction s'occupant d'allouer un tableau de type Part.
 * @param n Taille du tableau.
 * @return Adresse du premier élément du tableau.
 */
Part* Part1d(const int n) __attribute__ ((__const__));

/**
 * Fonction additionnant 2 particules.
 * @param a Particule a.
 * @param b Particule b.
 * @return La somme de a et b.
 */
Part  Part_add(Part a, Part b) __attribute ((__const__));

void Part_SortById(Part *part, const int NbPart);

/**
 * Fonction libérant le tableau alloué par \ref Part1d.
 * @param *ptf Adresse du premier élément du tableau à désallouer;
 */
void Part1d_libere(Part *ptf);

/**
 * Fonction adaptée à qsort pour le trie des tableaux de particules selon la distance à (0, 0, 0).
 * @param *a Premier objet à comparer.
 * @param *b Second objet à comparer.
 * @return -1 si a < b, 0 si a == b, 1 si a > b.
 */
int    qsort_partstr(const void *a, const void *b) __attribute__ ((__const__));

/**
 * Fonction adaptée à qsort pour le trie des tableaux de particules selon l'axe x.
 * @param *a Premier objet à comparer.
 * @param *b Second objet à comparer.
 * @return -1 si a < b, 0 si a == b, 1 si a > b.
 */
int    qsort_partaxe(const void *a, const void *b) __attribute__ ((__const__));

/**
 * Fonction se chargeant d'échanger 2 tableaux entre eux. Utilisé dans la fonction de construction de l'arbre.
 * @param *a Pointeur sur le premier élément du tableau ligne à échanger
 * @param *b Pointeur sur le premier élément du tableau ligne avec qui échanger
 */
void   Echange(Part *a, Part *b);

/**
 * \struct VolVois
 * Structure contenant tout les infos nécessaires sur les particules contenues dans la sphère de rayon `farest`.
 */
typedef struct _vois_vol {
	Part **part;	/*< Tableau de particule contenu dans la sphère.*/
	size_t size;	/*< Nombre d'éléments contenu dans la sphère.*/
	size_t cap;	/*< Nombre d'éléments pouvant être stocké avant de devoir faire une ré-allocation.*/
	double farest;	/*< Distance de la particule la plus lointaine.*/
} VolVois;

/**
 * Allocation de la structure.
 * @param pointeur à allouer.
 */
void VolVois_New(VolVois *new);

/**
 * Libération de la structure.
 * @param pointeur à libérer.
 */
void VolVois_Free(VolVois *this);

#endif /* end of include guard: TYPES_H */
