#include "./epithelium_submodel.h" 

using namespace PhysiCell; 

std::string epithelium_submodel_version = "0.0.1"; 

Submodel_Information epithelium_submodel_info; 

void epithelium_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// receptor dynamics 
	// requires faster time scale - done in main function 
	
	// viral dynamics model 
	internal_viral_dynamics_info.phenotype_function(pCell,phenotype,dt); 
	// internal_virus_model(pCell,phenotype,dt);
	
	// viral response model 
	internal_virus_response_model_info.phenotype_function(pCell,phenotype,dt); 
	// internal_virus_response_model(pCell,phenotype,dt);	
	
	// T-cell based death
	TCell_induced_apoptosis(pCell, phenotype, dt ); 
	
	return; 
}

void epithelium_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if I'm dead, don't bother 
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		
		// Let's just fully disable now. 
		pCell->functions.custom_cell_rule = NULL; 
		return; 
	}	
	
	// if I'm adhered to something ... 
	if( pCell->state.neighbors.size() > 0 )
	{
		// add the elastic forces 
		extra_elastic_attachment_mechanics( pCell, phenotype, dt );
	}
	return; 
}

void epithelium_submodel_setup( void )
{
	Cell_Definition* pCD;
	
	// set up any submodels you need 
	// viral replication 
	
	// receptor trafficking 
	receptor_dynamics_model_setup(); // done 
	// viral replication 
	internal_virus_model_setup();	
	// single-cell response 
	internal_virus_response_model_setup(); 
 	
	// set up epithelial cells
		// set version info 
	epithelium_submodel_info.name = "epithelium model"; 
	epithelium_submodel_info.version = epithelium_submodel_version; 
		// set functions 
	epithelium_submodel_info.main_function = NULL; 
	epithelium_submodel_info.phenotype_function = epithelium_phenotype; 
	epithelium_submodel_info.mechanics_function = epithelium_mechanics; 
	
		// what microenvironment variables do you expect? 
	epithelium_submodel_info.microenvironment_variables.push_back( "virion" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "interferon 1" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "pro-inflammatory cytokine" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "chemokine" ); 
		// what custom data do I need? 
	//epithelium_submodel_info.cell_variables.push_back( "something" ); 
		// register the submodel  
	epithelium_submodel_info.register_model();	
		// set functions for the corresponding cell definition 
	pCD = find_cell_definition( "lung epithelium" ); 
	pCD->functions.update_phenotype = epithelium_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = epithelium_submodel_info.mechanics_function;
	
	return;
}

void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int apoptosis_index = phenotype.death.find_death_model_index( "apoptosis" ); 
	if( pCell->custom_data["TCell_contact_time"] > pCell->custom_data["TCell_contact_death_threshold"] )
	{
		// make sure to get rid of all adhesions! 
		// detach all attached cells 
		for( int n = 0; n < pCell->state.neighbors.size() ; n++ )
		{
			detach_cells( pCell, pCell->state.neighbors[n] ); 
		}
		
		// induce death 
		pCell->start_death( apoptosis_index ); 
		pCell->functions.update_phenotype = NULL; 
	}
	
	return; 
}

