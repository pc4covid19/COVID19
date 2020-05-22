#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

#include "./immune_submodels.h" // done 
#include "./receptor_dynamics.h" // done 
#include "./internal_viral_dynamics.h" // done 
#include "./internal_viral_response.h"  // done 

#ifndef __epithelium_submodel__
#define __epithelium_submodel__

/*
void immune_submodel_setup( void );
void immune_submodel_model( Cell* pCell, Phenotype& phenotype, double dt );
*/

void epithelium_phenotype( Cell* pCell, Phenotype& phenotype, double dt ); 
void epithelium_mechanics( Cell* pCell, Phenotype& phenotype, double dt ); 

// this damage response will need to be added to the "infected cell response" model 
void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt ); 

void epithelium_submodel_setup( void ); 

#endif 