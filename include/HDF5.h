#ifndef VERIF_HDF5
#define VERIF_HDF5

#include <hdf5.h>
#include <stdbool.h>

typedef struct _ExtensibleDataSet {
	hsize_t size[2], offset[2];
	hid_t dataset, file, grp;
	bool incremente;
} *ExtensibleDataSet;

//ExtensibleDataSet CreateExtensibleDS(hid_t id, const char *name, const char *sub_grp, const int lig, const int col);
ExtensibleDataSet CreateExtensibleDS(hid_t id, const char *sub_grp, hsize_t dims[2]);
void ExtensibleDataSet_Extend(ExtensibleDataSet id, const double *data, hsize_t dims[2]);
void ExtensibleDataSet_Close(ExtensibleDataSet id);

#endif /* end of include guard */
