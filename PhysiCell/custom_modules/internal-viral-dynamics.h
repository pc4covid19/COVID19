#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#ifndef __internal_viral_dynamics__
#define __internal_viral_dynamics__
	
// inputs: 

// outputs: 

void internal_virus_model( Cell* pCell, Phenotype& phenotype, double dt );

#endif 