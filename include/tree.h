#ifndef __TREE_H_200112_155224
#define __TREE_H_200112_155224

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <assert.h>
#include <stdbool.h>

#include "types.h"

/**
 * \file tree.h
 *
 * Le but de ce code est de mettre en place un algorithme de Tree code
 * dans le but d'accélerer (et de faciliter) le calcul du centre de densité.
 * Par la suite, ce tree-code sera utilisé pour calculer le potentiel et le
 * rapport du Viriel avec plus de précision que par la méthode utilisant des
 * histogrammes utilisé jusqu'ici.
 *
 * Voir l'article :
 * 	A hierarchical \f$O(N \log N)\f$ force-calculation algorithm.
 * 	Josh Barnes and Piet Hut
 * 	Nature 324, 446-449, 04/12/1986
 */

/**
 * \todo Repenser \ref Tree_Save
 * \todo Ré-écrire \ref Tree_Read en fonction de \ref Tree_Save
 * \todo Compléter l'écriture du code en utilisant des structure plutôt que des
 * tableaux de double (\ref TNoeud).
 */

/**
 * \def NEAREST
 * Calcul la distance entre 2 points en prenant en compte les conditions périodiques
 */
#define NEAREST(x, boxhalf, boxsize) (((x)>(boxhalf))?((x)-(boxsize)):(((x)<-(boxhalf))?((x)+(boxsize)):(x)))

/**
 * \struct TNoeud
 * Contient toutes les informations nécessaire à chaque noeud de l'arbre (quad/oct-tree)
 * Sert à l'algorithme du tree-code
 */
struct _tnoeud {
	int    N;			/*< Nombre de particule contenue dans le noeud.*/
	int    level;			/*< Niveau du noeud dans l'arbre.*/
	double x, y, z;			/*< Coordonnées du centre géometrique du carré.*/
	double cote;			/*< Longueur du côté du carré.*/
	Part *first;			/*< Pointeur vers la particule du tableau de particule.*/
	Part cm;			/*< Position du centre de masse de la cellule.*/
	double CM;			/*< Masse totale dans la cellule.*/
	struct _tnoeud* parent;		/*< Pointeur vers le noeud parent.*/
	struct _tnoeud* frere;		/*< Pointeur vers le frére de gauche.*/
	struct _tnoeud* fils;		/*< Pointeur vers le 1er fils, le plus à droite.*/
};
typedef struct _tnoeud* TNoeud;

/**
 * Fonction initialisant les valeurs des colonnes pour le tableau.
 * Appel obligatoire pour pouvoir utiliser le Tree Code.
 * @param zp zp <= -1 : QuadTree (Espace 2D), sinon : OctTree (Espace 3D)
 * @param LMax Niveau de raffinement maximum de l'arbre.
 */
void   Tree_var(const int LMax);

/**
 * Fonction vérifiant que la particule "cherch" ne se trouve pas dans le tableau de particule "Tab".
 * @param cherch Particule recherché
 * @param *Tab Tableau de particule
 * @param NbVois Taille du tableau de particule
 * @return 0 si la particule est dans le tableau, 1 sinon.
 */
int    NotIn(Part cherch, Part *Tab, const int NbVois);

/**
 * Fonction chargé d'initialiser la racine ou un nœud qeulconque pour le Tree Code.
 * @param NbPart nombre de particules contenues dans la racine
 * @param xc abscisse du centre du cube racine
 * @param yc ordonnée du centre du cube racine
 * @param zc côte du centre du cube racine
 * @param cote taille d'un côté du cube
 * @return pointeur vers le noeud initialisé.
 */
TNoeud Tree_Init (int  NbPart, double xc, double yc, double zc, double cote);

/**
 * Fonction parcourant le tableaux de particules pour mettre celles appartenant au cube t1
 * au début du tableau t1->first (\ref TNoeud), tout en incrémentant le nombre de particule
 * t1->N pour la cellule t1.
 * Version 3d.
 * @param t1 Cellule devant recevoir les particules.
 * @param NbPart Nombre de particule à parcourir entre t1->first[0] et la fin du tableau (t1->first[NbPart-1])
 */
void Tree_Calc(TNoeud t1, const int NbPart);

/**
 * Fonction chargé de créer le Tree Code.
 * Seconde version plus efficace et plus souble que \ref Tree_Build
 * @param root noeud à partir duquel construire l'arbre
 * @param NbPart nombre de particules contenues dans le noeud
 * @param NbMin nombre de particule minimale dans un cube
 * @return 0 si succes
 */
int    Tree_Build2(TNoeud root, int NbPart, int NbMin);

/**
 * Fonction de calcul des voisins par défaut, version 3d.
 * @param *insert Tableau contenant les particules à insérer.
 * @param N Taille du tableau contenant les particules à insérer.
 * @param *Tab Tableau dans lequel on souhaite insérer les particules de *insert si possible.
 * @param NbVois Taille du second tableau
 * @param *part Particule dont on cherche les voisins
 */
void   CalcVois(Part *insert, const int N, Part *Tab, const int NbVois, const Part *part, const double BS);

/**
 * Fonction parcourant le tree-code pour en calculer les voisins.
 * \todo Généraliser cette fonction pour une utilisation plus générale,
 * comme pour le potentiel.
 * \todo Rendre périodique le calcul.
 * @param root Nœud racine à partir duquel parcourir l'arbre.
 * @param *Tab Tableau contenant les voisins (il doit être initialisé
 * avec des particules, sans critére de choix particulier).
 * @param NbVois Nombre de voisin.
 * @param *part Particule dont on veut les voisins.
 */
void   Tree_Voisin(TNoeud root, Part *Tab, int NbVois, const Part *part, const double BS);

/**
 * Fonction sauvegardant l'arbre dans un fichier.
 * @param root noeud à partir duquel commencer l'enregistrement
 * @param fich pointeur vers le fichier à écrire
 */
void   Tree_Save (TNoeud root, FILE *fich);

/**
 * Fonction lisant l'arbre dans un fichier.
 * ***À MODIFIER***
 * @param root noeud à partir duquel commencer la lecture
 * @param fich pointeur vers le fichier à écrire
 * @return 0 si tout va bien, -5 sinon
 */
int    Tree_Read (TNoeud root, FILE *fich);

/**
 * Fonction libérant la mémoire associé à l'arbre.
 * @param root le noeud à partir duquel commencer la libération
 */
void   Tree_Free (TNoeud root);

/**
 * Fonction se servant du Tree-Code pour calculer le potentiel.
 * @param root noeud à partir duquel commencer le calcul (racine de l'arbre).
 * @param *part Tableau de particule (sais plus à quoi il sert).
 * @param accept Critère d'acceptation des particules (angle d'ouverture du Tree Code).
 * @param soft Softening à appliquer lors du calcul du potentiel.
 */
double Tree_CalcPot(TNoeud root, const Part *part, const double accept, const double soft, const double BS);

void Tree_SetG(double nG);
double Tree_GetG(void);

#endif //__TREE_H_200112_155224
