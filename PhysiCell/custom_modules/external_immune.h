#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

#ifndef __external_immune__
#define __external_immune__
	
extern Submodel_Information external_immune_info; 
	
void external_immune_model_setup( void );

void external_immune_model( double dt );

#endif 