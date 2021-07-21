#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

#ifndef __DC_history__
#define __DC_history__

extern Submodel_Information DC_history_info; 
	
void DC_history_model_setup( void );

// don't put into individual cell models
void DC_history_model( Cell* pCell, Phenotype& phenotype, double dt );

// this needs to be done on faster time scale; 
void DC_history_main_model( double dt ); 


#endif 