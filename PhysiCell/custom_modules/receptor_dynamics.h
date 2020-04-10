#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

#ifndef __receptor_dynamics__
#define __receptor_dynamics__
	
void receptor_dynamics_model_setup( void );
void receptor_dynamics_model( Cell* pCell, Phenotype& phenotype, double dt );

// this needs to be done on faster time scale; 
void receptor_dynamics_model( double dt ); 


#endif 