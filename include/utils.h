#ifndef __UTILS_H_GUI__
#define __UTILS_H_GUI__

#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/****************************************************************************************\
 * 	Code récupéré de la librairie mnitab de Mr lefrére, enseignant du MNI du M1	*
 * 	Physique Générale et du M1 SDUEE de l'université Pierre et Marie Curie Paris 6	*
\****************************************************************************************/
unsigned int* unsignedint1d(const unsigned int n);
void unsignedint1d_libere(unsigned int *ptf);
int* int1d(const int n);
void int1d_libere(int *ptf);
float* float1d(const int n);
void float1d_libere(float *ptf);
double* double1d(const int n);
void double1d_libere(double *ptf);
int maxlocdouble1d(double * x, int n);
int minlocdouble1d(double * x, int n);

double ** double2d(const int nlignes, const int mcolonnes);
void double2d_libere(double ** mat);

double maxdouble2d(double **tab, const int NbPart, const int col);

void lissage(double *tab, const int N);
char *remove_ext(const char* mystr);
int get_id(const char *filename);

#endif
