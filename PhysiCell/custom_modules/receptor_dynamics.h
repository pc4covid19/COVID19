#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

#ifndef __receptor_dynamics__
#define __receptor_dynamics__

extern Submodel_Information receptor_dynamics_info; 
	
void receptor_dynamics_model_setup( void );

// don't put into individual cell models -- needs to be a fast process 
void receptor_dynamics_model( Cell* pCell, Phenotype& phenotype, double dt );

// this needs to be done on faster time scale; 
void receptor_dynamics_main_model( double dt ); 


#endif 