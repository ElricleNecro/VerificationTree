#ifndef __VERIF_TOOLS_H_PERSO__
#define __VERIF_TOOLS_H_PERSO__

#include <math.h>

#include "tree.h"
#include "utils.h"
#include "types.h"
#include "cte_phys.h"

// Flags de debug :
//	-> 1 pour activé, 0 sinon
#define __DEBUG_QSORT_P 1

#ifdef __DEBUG_VOIS_LOG
Part MoreDenseParticule(const TNoeud root, const int NbVois, const double BS, const char * fname);
#else
#ifdef PERIODIC
Part MoreDenseParticule(const TNoeud root, const int NbVois, const double BS);
#else
Part MoreDenseParticule(const TNoeud root, const int NbVois);
#endif
#endif
Part ReCentre(TNoeud root, Part *posvits, const int NbPart, const int NbVois, const int NbMin, const double BoxSize);

//Part GravityCenter(Part *list, const int NbPart);
/**
 * Calcule le centre de gravité du système en utilisant le tree code.
 * @param root Racine du tree code.
 * @param N Dummy variable : ne sert pas réellement dans le calcul.
 * @return Coordonnée et vitesse, de type Part, du centre de gravité.
 */
Part GravityCenter(const TNoeud root, const int N) __attribute__ ((__const__));

#ifdef PERIODIC
/**
 * Calcule le centre de densité du système en utilisant le tree code.
 * @param root Racine du tree code.
 * @param N Nombre de voisin à utiliser lors du calcul du centre de densité.
 * @param BS Taille de la boîte.
 * @return Coordonnée et vitesse, de type Part, du centre de densité.
 */
Part DensityCenter(const TNoeud root, const int NbVois, const double BS) __attribute__ ((__const__));
#else
/**
 * Calcule le centre de densité du système en utilisant le tree code.
 * @param root Racine du tree code.
 * @param N Nombre de voisin à utiliser lors du calcul du centre de densité.
 * @return Coordonnée et vitesse, de type Part, du centre de densité.
 */
Part DensityCenter(const TNoeud root, const int NbVois) __attribute__ ((__const__));
#endif

/**
 * Calcul les ratios des axes principaux de la matrice d'inertie associé au système.
 * @param root Racine de l'arbre.
 * @param R_ori Rayon d'origine de l'objet, sert pour rejeter certaine particule (pas complétement implémenté).
 * @param *grand Rapport du grand axe sur l'intémerdiaire.
 * @param *petit Rapport du petit axe sur l'intermediaire.
 */
void Axial_ratio(const TNoeud root, double R_ori, double *grand, double *petit);

#ifdef PERIODIC
/**
 * Calcul le potentiel du systéme en utilisant un calcul direct à l'aide du tree code.
 * @param root Racine de l'arbre.
 * @param theta Critére de sélection des cubes.
 * @param rsoft Paramètre de lissage du potentiel.
 * @return Rayon sur la colonne 0 et potentiel sur la colonne 1, classé de la même maniére que leq particules de l'arbre.
 */
double** CalcPotentiel(const TNoeud root, const double theta, const double rsoft, const double BS) __attribute__ ((__const__));
#else
/**
 * Calcul le potentiel du systéme en utilisant un calcul direct à l'aide du tree code.
 * @param root Racine de l'arbre.
 * @param theta Critére de sélection des cubes.
 * @param rsoft Paramètre de lissage du potentiel.
 * @param BS Taille de la boîte.
 * @return Rayon sur la colonne 0 et potentiel sur la colonne 1, classé de la même maniére que leq particules de l'arbre.
 */
double** CalcPotentiel(const TNoeud root, const double theta, const double rsoft) __attribute__ ((__const__));
#endif

/**
 * Calcul la fonction de masse de l'objet.
 * @param root Racine de l'arbre.
 * @return La fonction de masse pour chaque rayon.
 */
double* CalcMasse(const TNoeud root) __attribute__ ((__const__));

/**
 * Fonction se chargeant de calculer les énergies cinétiques et totales à partir
 * du tableau de potentiel et des données de l'arbre.
 *
 * \todo Il y a des résultats bizarre à la sorti pour l'énergie cinétique totale, il faut les corriger.
 * @param[in] root Racine de l'arbre.
 * @param[out] *energie_c Énergie cinétique (\f$E_c = \frac{1}{2} m v^2\f$) classé de la même maniére que les particules dans l'arbre.
 * @param[out] *energie_t Énergie totale du sytème (\f$E_t = E_c + E_p\f$).
 * @param[in] **potentiel Potentiel de l'objet.
 * @param[out] *Ectot Énergie cinétique totale.
 * @param[out] *Eptot Énergie potentiel totale.
 */
void    CalcEnergie(const TNoeud root, double *energie_c, double *energie_t, double **potentiel, double *Ectot, double *Eptot);

/**
 * Fonction calculant le profil de densité de l'objet.
 *
 * @param[in] root Racine de l'arbre.
 * @param[in] NbBin Nombre de bin à utiliser pour calculer le profil.
 * @param[in] dr Taille d'un bin.
 * @param[in] rmax Taille maximum.
 * @return Tableau contenant le profil de Densité.
 */
double* CalcDensite(const TNoeud root, const int NbBin, const double dr, const double rmax) __attribute__ ((__const__));

/**
 * Fonction calculant le profil de densité de l'objet, en utilisant un pas constant en log.
 *
 * @param[in] root Racine de l'arbre.
 * @param[in] NbBin Nombre de bin à utiliser pour calculer le profil.
 * @param[in] rmin Borne inférieur.
 * @param[in] rmax Borne supérieur.
 * @return Tableau contenant le profil de Densité.
 */
//double*  CalcLogDensite(const TNoeud root, const int NbBin, const double cte, const double rmax);
double** CalcLogDensite(const TNoeud root, const int NbBin, const double rmin, const double rmax, const double rnorm) __attribute__ ((__const__));

/**
 * Fonction calculant le profil de température.
 *
 * @param[in] root Racine de l'arbre.
 * @param[in] NbBin Nombre de bin à utiliser pour calculer le profil.
 * @param[in] dr Taille d'un bin.
 * @param[in] rmax Taille maximum.
 * @param[out] *Tmoy Température moyenne.
 * @return Tableau contenant le profil de Température.
 */
double* CalcTemperature(const TNoeud root, const int nb_bin, const double dr, double *Tmoy);

/**
 * Fonction calculant la distribution en énergie et le Jacobien permettant la transformation de \f$f(\vec{x}, \vec{p})\f$ vers \f$f(E)\f$.
 *
 * @param[in] root Racine de l'arbre.
 * @param[in] NbBin Nombre de bin à utiliser pour calculer le profil.
 * @param[in] *energie_t Tableau contenant l'énergie totale, tel que calculé par la fonction \ref CalcEnergie.
 * @param[in] **potentiel Tableau contenant le potentiel de l'objet, calculé par la fonction \ref CalcPotentiel.
 * @param[in] Emin Énergie minimum de l'objet.
 * @param[in] Emax Énergie maximum de l'objet.
 * @param[in] dE Intervalle d'énergie entre chaque bins d'énergie.
 * @param[in] G Constante gravitationnelle.
 * @param[out] *distrib Distribution en énergie
 */
double* CalcJacobien(const TNoeud root, const int NbBin, const double *energie_t, double **potentiel, const double Emin, const double Emax, const double dE, double *distrib);

/**
 * Fonction calculant le profil de densité de l'objet.
 *
 * @param[in] root Racine de l'arbre.
 * @param[in] NbBin Nombre de bin à utiliser pour calculer le profil.
 * @param[in] dr Taille d'un bin.
 * @param[in] rmax Taille maximum.
 * @param[out] *Coeff Anisotropie totale de l'objet.
 * @return Tableau contenant le profil de Densité.
 */
double* CalcAnisotropie(const TNoeud root, const int NbBin, const double dr, double *Coeff);

/**
 * Retourne le maximum entre 3 paramètres
 * @param a 1éres valeurs à comparer
 * @param b 2nde valeurs à comparer
 * @param c 3émes valeurs à comparer
 * @return Le maximum
 */
double min(double a, double b, double c);

/**
 * Retourne le minimum entre 3 paramètres
 * @param a 1éres valeurs à comparer
 * @param b 2nde valeurs à comparer
 * @param c 3émes valeurs à comparer
 * @return Le minimum
 */
double max(double a, double b, double c);

#endif
