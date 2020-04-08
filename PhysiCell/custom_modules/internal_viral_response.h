#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

#ifndef __internal_viral_response__
#define __internal_viral_response__
	
// inputs: 

// outputs: 

extern Submodel_Information internal_virus_response_model_info; 

void internal_virus_response_model_setup( void );
void internal_virus_response_model( Cell* pCell, Phenotype& phenotype, double dt );

#endif 