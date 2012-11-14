#include "io_snapshot.h"

Particle load_snapshot(const char *fname, const int files, int *NbPart, int *Ngas, double *tps, IO_Header *header2)
{
	FILE *fd;
	char buf[200];
	double Time;
	int i, k, dummy, ntot_withmasses;
	int n, pc, pc_new, pc_sph;
	int NumPart, *Id;
	Particle P = NULL;
	IO_Header header1;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

	for(i = 0, pc = 1; i < files; i++, pc = pc_new)
	{
		if(files > 1)
			sprintf(buf, "%s.%d", fname, i);
		else
			sprintf(buf, "%s", fname);

		if(!(fd = fopen(buf, "r")))
		{
			printf("can't open file `%s`\n", buf);
			exit(EXIT_FAILURE);
		}

		printf("reading `%s' ...\n", buf);
		fflush(stdout);

		fread(&dummy, sizeof(dummy), 1, fd);
		fread(&header1, sizeof(header1), 1, fd);
		fread(&dummy, sizeof(dummy), 1, fd);

		if(files == 1)
		{
			for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
				NumPart += header1.npart[k];
			*Ngas = header1.npart[0];
		}
		else
		{
			for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
				NumPart += header1.npartTotal[k];
			*Ngas = header1.npartTotal[0];
		}

		for(k = 0, ntot_withmasses = 0; k < 6; k++)
		{
			if(header1.mass[k] == 0)
				ntot_withmasses += header1.npart[k];
		}

		if(i == 0)
		{
			printf("allocating memory...\n");
			//allocate_memory();
			P = allocate_particle_memory(NumPart);
			Id = allocate_id_memory(NumPart);
			printf("allocating memory...done\n");
		}

		SKIP;
		for(k = 0, pc_new = pc; k < 6; k++)
		{
			for(n = 0; n < header1.npart[k]; n++)
			{
				fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
				pc_new++;
			}
		}
		SKIP;

		SKIP;
		for(k = 0, pc_new = pc; k < 6; k++)
		{
			for(n = 0; n < header1.npart[k]; n++)
			{
				fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
				pc_new++;
			}
		}
		SKIP;


		SKIP;
		for(k = 0, pc_new = pc; k < 6; k++)
		{
			for(n = 0; n < header1.npart[k]; n++)
			{
				fread(&Id[pc_new], sizeof(int), 1, fd);
				pc_new++;
			}
		}
		SKIP;


		if(ntot_withmasses > 0)
			SKIP;
		for(k = 0, pc_new = pc; k < 6; k++)
		{
			for(n = 0; n < header1.npart[k]; n++)
			{
				P[pc_new].Type = k;

				if(header1.mass[k] == 0)
					fread(&P[pc_new].Mass, sizeof(float), 1, fd);
				else
					P[pc_new].Mass = header1.mass[k];
				pc_new++;
			}
		}
		if(ntot_withmasses > 0)
			SKIP;


		if(header1.npart[0] > 0)
		{
			SKIP;
			for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
			{
				fread(&P[pc_sph].U, sizeof(float), 1, fd);
				pc_sph++;
			}
			SKIP;

			SKIP;
			for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
			{
				fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
				pc_sph++;
			}
			SKIP;

			if(header1.flag_cooling)
			{
				SKIP;
				for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
				{
					fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
					pc_sph++;
				}
				SKIP;
			}
			else
				for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
				{
					P[pc_sph].Ne = 1.0;
					pc_sph++;
				}
		}

		fclose(fd);
	}


	*tps = Time = header1.time;
	//Redshift = header1.time;

	for(i = 1; i < NumPart; i++)
		P[i].Id = Id[i];

	*header2 = header1;
	*NbPart = NumPart;

	Id++;
	free(Id);

	return P;
}

Part* read_snapshot(const char *fname, const int files, const int type, const double PosFact, const double VitFact, int *NbPart, double *time, IO_Header *hea)
{
	Part* part = NULL;
	Particle P = NULL;
	IO_Header header;
	int npart, NbGas;

	P = load_snapshot(fname, files, &npart, &NbGas, time, &header);

	printf("Type : %d\n", type);
	printf("\033[33m%d particules de masse %g à lire dans : %s\033[00m\n", header.npart[4], header.mass[4], fname);

	printf("\033[36mHeader du fichier Gadget (format 1) :\033[00m\n");
	printf("\033[34m\tNombre de fichier par snapshot : %d\n", header.num_files);
	printf("\tMasse et nombre d'éléments des catégories d'objet :\n");
	for(int i = 0; i < 6; i++)
		printf("\t\t%s : Masse %g, et %d élément%c (total : %d)\n", (i == 0)?"Gaz":( (i == 1)?"Halo":( (i == 2)?"Disk":( (i==3)?"Bulge":( (i==4)?"Stars":"Bndry" )))), header.mass[i], header.npart[i], (header.npart[i] > 1)?'s':' ', header.npartTotal[i]);
	puts("\033[00m");


	part    = Part1d(header.npart[type]);
	*NbPart = header.npart[type];
	*hea    = header;

	printf("%d\t%d\t%d\t%d\t%g\n", *NbPart, header.npart[type], type, npart, header.BoxSize);

	for (int i = 1, j = 0; i <= npart && j < header.npart[type]; i++)
	{
		if(P[i].Type == type)
		{
			part[j].x  = P[i].Pos[0] * PosFact;
			part[j].y  = P[i].Pos[1] * PosFact;
			part[j].z  = P[i].Pos[2] * PosFact;
			part[j].r  = sqrt(part[j].x*part[j].x + part[j].y*part[j].y + part[j].z*part[j].z);

			part[j].vx = P[i].Vel[0] * VitFact;
			part[j].vy = P[i].Vel[1] * VitFact;
			part[j].vz = P[i].Vel[2] * VitFact;
			part[j].v  = sqrt(part[j].vx*part[j].vx + part[j].vy*part[j].vy + part[j].vz*part[j].vz);

			part[j].id = P[i].Id;
			part[j].m  = P[i].Mass;
			if(part[j].m == 0.0)
			{
				fprintf(stderr, "%d:: Masse nulle, initialisation au Header : (%g, %g) -> %g\n", j, part[j].m, P[i].Mass, header.mass[type]);
				part[j].m = header.mass[type];
			}
			j++;
		}
	}

	P++;
	free(P);

	return part;
}

Part* Deconversion(char *file_name, int *Nb)
{
	return read_snapshot(file_name, 1, 4, 1.0, 1.0, Nb, NULL, NULL);
}

Particle allocate_particle_memory(const int NumPart)
{
	Particle P = NULL;

	if(!(P = malloc(NumPart * sizeof(struct particle_data))))
	{
		fprintf(stderr, "failed to allocate memory.\n");
		exit(EXIT_FAILURE);
	}

	P--;				/* start with offset 1 */

	return P;
}

int* allocate_id_memory(const int NumPart)
{
	int *Id = NULL;
	if(!(Id = malloc(NumPart * sizeof(int))))
	{
		fprintf(stderr, "failed to allocate memory.\n");
		exit(EXIT_FAILURE);
	}

	Id--;				/* start with offset 1 */

	return Id;
}


