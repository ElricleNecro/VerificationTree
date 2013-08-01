#include "HDF5.h"

#include <stdlib.h>
#include <string.h>

ExtensibleDataSet CreateExtensibleDS(hid_t id, const char *sub_grp, hsize_t dims[2]) //const int lig, const int col)
{
	ExtensibleDataSet  new = malloc(sizeof( struct _ExtensibleDataSet ));
	hid_t              cparms;
	hid_t              ds,
			   dsp;
	hsize_t            maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
	herr_t             status;

	// Création du Chunk permettant d'insérer des colonnes ou lignes au fur et à mesure :
	cparms         = H5Pcreate(H5P_DATASET_CREATE);
	status         = H5Pset_chunk(cparms, 2, dims);

	// Création de l'espace de travail à allouer :
	dsp            = H5Screate_simple(2, dims, maxdims);

	// Création du dataset en lui-même :
	ds             = H5Dcreate(id, sub_grp, H5T_IEEE_F64LE, dsp, H5P_DEFAULT, cparms, H5P_DEFAULT);

	// Enregistrement des propriétés du Dataset :
	new->size[0]   = dims[0];
	new->size[1]   = dims[1];
	new->dataset   = ds;
	new->file      = id;
	new->offset[0] = 0;
	new->offset[1] = 0;

	H5Sclose(dsp);

	return new;
}

void ExtensibleDataSet_Extend(ExtensibleDataSet id, const double *data, hsize_t dims[2])
{
	// On étend l'espace alloué au dataset :
	id->size[1]     += dims[1];
	herr_t status    = H5Dextend(id->dataset, id->size);
	// On récupère les infos sur l'espace nouvellement alloué :
	hid_t  filespace = H5Dget_space(id->dataset);
	// On sélectionne la zone où l'on va écrire les nouvelles données :
	status           = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, id->offset, NULL, dims, NULL);
	// On crée le dataspace pour les nouvelles données :
	hid_t  dataspace = H5Screate_simple(2, dims, NULL);

	// On écrit les données :
	status           = H5Dwrite(id->dataset, H5T_NATIVE_DOUBLE, dataspace, filespace, H5P_DEFAULT, data);

	// On incrémente l'offset :
	//id->offset[0]    += dims[0];
	id->offset[1]   += dims[1];

	// On libère la mémoire :
	H5Sclose(filespace);
	H5Sclose(dataspace);
}

void ExtensibleDataSet_Close(ExtensibleDataSet id)
{
	H5Dclose(id->dataset);
}

