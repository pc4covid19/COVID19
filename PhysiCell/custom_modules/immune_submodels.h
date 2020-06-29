#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

#ifndef __immune_submodels__
#define __immune_submodels__

/*
void immune_submodel_setup( void );
void immune_submodel_model( Cell* pCell, Phenotype& phenotype, double dt );
*/

/* functions for checking nearby cells */
Cell* check_for_live_neighbor_for_interaction( Cell* pAttacker , double dt );
Cell* check_for_dead_neighbor_for_interaction( Cell* pAttacker , double dt ); 

/* functions for cell-cell adhesion */ 
bool attempt_immune_cell_attachment( Cell* pAttacker, Cell* pTarget , double dt );
Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt );

void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant );
void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt );

/*
void attach_cells( Cell* pCell_1, Cell* pCell_2 );
void detach_cells( Cell* pCell_1 , Cell* pCell_2 );
void remove_all_adhesions( Cell* pCell ); 
*/

std::vector<double> choose_vascularized_position( void );

//void choose_initialized_voxels( void );
//extern std::vector<int> vascularized_voxel_indices;
extern std::vector<int> vascularized_voxel_indices;
void choose_initialized_voxels( void );

extern std::vector<Cell*> cells_to_move_from_edge; 
void process_tagged_cells_on_edge( void ); 

void CD8_Tcell_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2 , double dt );
void CD8_Tcell_phenotype( Cell* pCell, Phenotype& phenotype, double dt ); 
void CD8_Tcell_mechanics( Cell* pCell, Phenotype& phenotype, double dt ); 

void immune_cell_motility_direction( Cell* pCell, Phenotype& phenotype , double dt );

void macrophage_phenotype( Cell* pCell, Phenotype& phenotype, double dt ); 
void macrophage_mechanics( Cell* pCell, Phenotype& phenotype, double dt ); 

void neutrophil_phenotype( Cell* pCell, Phenotype& phenotype, double dt ); 
void neutrophil_mechanics( Cell* pCell, Phenotype& phenotype, double dt ); 


// this damage response will need to be added to the "infected cell response" model 
// void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt ); 

// immune cell recruitment 

void create_infiltrating_immune_cell( Cell_Definition* pCD );
void create_infiltrating_immune_cell( std::string cell_name );
void create_infiltrating_immune_cell_initial( Cell_Definition* pCD );

void create_infiltrating_Tcell( void ); 
void create_infiltrating_neutrophil( void );
void create_infiltrating_macrophage( void );

void immune_cell_recruitment( double dt ); 

void initial_immune_cell_placement( void ); 

void immune_submodels_setup( void ); 

// this needs to be done on faster time scale; 
// void receptor_dynamics_model( double dt ); 

bool check_for_out_of_bounds( Cell* pC , double tolerance );
std::vector<double> set_nudge_from_edge( Cell* pC , double tolerance );

void keep_immune_cells_off_edge(void); 
void keep_immune_cells_in_bounds( double dt );

void detach_all_dead_cells( void );
void detach_all_dead_cells( double dt );

#endif 