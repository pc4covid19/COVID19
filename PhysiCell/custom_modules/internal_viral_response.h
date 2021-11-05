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

void pyroptosis_cascade( Cell* pCell, Phenotype& phenotype, double dt );

void create_secreting_agentcallvir(double positionpass0, double positionpass1);

void create_secreting_agentvir( Cell_Definition* pCD, double positionpass0, double positionpass1);

#endif 