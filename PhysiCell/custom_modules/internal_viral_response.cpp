#include "./internal_viral_response.h" 

using namespace PhysiCell; 

std::string internal_virus_response_version = "0.2.0"; 

Submodel_Information internal_virus_response_model_info; 

void internal_virus_response_model_setup( void )
{
	// set up the model 
		// set version info 
	internal_virus_response_model_info.name = "internal viral response"; 
	internal_virus_response_model_info.version = internal_virus_response_version; 
		// set functions 
	internal_virus_response_model_info.main_function = NULL; 
	internal_virus_response_model_info.phenotype_function = internal_virus_response_model; 
	internal_virus_response_model_info.mechanics_function = NULL; 
	
		// what microenvironment variables do you expect? 
		// what custom data do I need? 
	internal_virus_response_model_info.cell_variables.push_back( "max_infected_apoptosis_rate" ); 
	internal_virus_response_model_info.cell_variables.push_back( "max_apoptosis_half_max" ); 
	internal_virus_response_model_info.cell_variables.push_back( "apoptosis_hill_power" ); 	
	
		// register the submodel  
	internal_virus_response_model_info.register_model();	
		// set functions for the corresponding cell definition 
		
//	pCD = find_cell_definition( "lung epithelium" ); 
//	pCD->functions.update_phenotype = epithelium_submodel_info.phenotype_function;
//	pCD->functions.custom_cell_rule = epithelium_submodel_info.mechanics_function;
	
	return; 
}

void internal_virus_response_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	// bookkeeping -- find microenvironment variables we need

	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
	
	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
	static int nA_internal = pCell->custom_data.find_variable_index( "assembled_virion" ); 
	
	// actual model goes here 

	// now, set apoptosis rate 
	
	static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "apoptosis" );
	phenotype.death.rates[apoptosis_model_index] = 
		pCell->custom_data["max_infected_apoptosis_rate"] ; 
	
	double v = pCell->custom_data[nA_internal] / 
		pCell->custom_data["max_apoptosis_half_max"] ; 
	v = pow( v, pCell->custom_data["apoptosis_hill_power"] ); 
	
	double effect = v / (1.0+v); 
	phenotype.death.rates[apoptosis_model_index] *= effect; 
	
	return; 
}
