#include <stdlib.h>

#include "tree.h"
#include "utils.h"
#include "Verif_tools.h"
#include "../Tree_Code/rand.h"

int main(int argc, char **argv)
{
	int NbPart  = atoi(argv[1]),
	    NbBin   = atoi(argv[2]),
	    seed    = 20369782;
	double rmax = 2.5e17,
	       Tmp  = 0.0,
	       *Temp;
	TNoeud root;

	printf("NbPart = %d\tNbBin = %d\n", NbPart, NbBin);

	root        = malloc(sizeof(*root));
	root->N     = NbPart;
	root->first = Part1d(root->N);

	for (int i = 0; i < NbPart; i++) {
		root->first[i].v = 300.0 * ran2(&seed) - 150.0;
		root->first[i].r = rmax  * ran2(&seed);
	}

	Temp = CalcTemperature(root, NbBin, rmax/NbBin, rmax, &Tmp);

	free(Temp);
	free(root->first);
	free(root);
}
