#ifndef _IOSNAPSHOT_HEADER_P
#define _IOSNAPSHOT_HEADER_P

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "types.h"

struct io_header_1
{
	int npart[6];
	double mass[6];
	double time;
	double redshift;
	int flag_sfr;
	int flag_feedback;
	int npartTotal[6];
	int flag_cooling;
	int num_files;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
};
typedef struct io_header_1 IO_Header;

struct particle_data
{
	float Pos[3];
	float Vel[3];
	float Mass;
	int Type;
	int Id;

	float Rho, U, Temp, Ne;
};
typedef struct particle_data* Particle;

/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
//int load_snapshot(char *fname, int files);
//Particule load_snapshot(char *fname, int files, int *NbPart, int *Ngas, int *uudouble *tps);
Particle load_snapshot(const char *fname, const int files, int *NbPart, int *Ngas, IO_Header *header1);

//Part* read_snapshot(const char *fname, const int files, const int type, int *NbPart);
//Part* read_snapshot(const char *fname, const int files, const int type, int *NbPart, double *time);
//Part* read_snapshot(const char *fname, const int files, const int type, const int periodic, int *NbPart, double *time, IO_Header *hea);
Part* read_snapshot(const char *fname, const int files, const int type, const double PosFact, const double VitFact, int *NbPart, double *time, IO_Header *hea, bool CorrectId);

/* this routine allocates the memory for the
 * particle data.
 */
//int allocate_memory(void);
Particle allocate_particle_memory(const int NumPart);
int* allocate_id_memory(const int NumPart);

/**
 * Convertit le fichier du format Gadget-2 vers quelque chose d'utilisable dans.
 * @param *file_name Nom du snapshot à lire.
 * @param Nb Nombre de particule lu dans le fichier
 * @return Tableau de particule vitesse de type \ref Part.
 */
Part* Deconversion(char *file_name, int *Nb);

#endif
