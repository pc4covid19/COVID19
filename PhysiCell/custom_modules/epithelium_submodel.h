#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

#include "./immune_submodels.h"  
#include "./receptor_dynamics.h"  
#include "./internal_viral_dynamics.h"  
#include "./internal_viral_response.h"   

#ifndef __epithelium_submodel__
#define __epithelium_submodel__

extern Submodel_Information epithelium_submodel_info; 

void epithelium_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt ); 
void epithelium_phenotype( Cell* pCell, Phenotype& phenotype, double dt ); 
void epithelium_mechanics( Cell* pCell, Phenotype& phenotype, double dt ); 

// this damage response will need to be added to the "infected cell response" model 
void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt ); 

void epithelium_submodel_setup( void ); 

// (Adrianne V5) Model for ROS induction of cell apoptosis
void ROS_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt ); 

void create_secreting_agentcall(double positionpass0, double positionpass1);
void create_secreting_agent( Cell_Definition* pCD, double positionpass0, double positionpass1);
// void create_secreting_agent( std::string cell_name );

#endif 