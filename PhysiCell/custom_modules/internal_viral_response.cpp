#include "./internal_viral_response.h" 

using namespace PhysiCell; 

Submodel_Information internal_virus_response_model_info; 

void internal_virus_response_model_setup( void )
{
	internal_virus_response_model_info.name = "internal viral response"; 
	internal_virus_response_model_info.version = "0.2.0";
	internal_virus_response_model_info.main_function= internal_virus_response_model; 
	
	// what custom data do I need? 
	
	internal_virus_response_model_info.cell_variables.push_back( "assembled virion" ); 
	
	internal_virus_response_model_info.cell_variables.push_back( "max infected apoptosis rate" ); 
	internal_virus_response_model_info.cell_variables.push_back( "max apoptosis half max" ); 
	internal_virus_response_model_info.cell_variables.push_back( "apoptosis hill power" ); 
	
	// submodel_registry.register_model( internal_virus_response_model_info ); 	
	internal_virus_response_model_info.register_model();	
	
	return; 
}

void internal_virus_response_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	// bookkeeping -- find microenvironment variables we need

	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
	
	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
	static int nA_internal = pCell->custom_data.find_variable_index( "assembled virion" ); 

	static int nUV = pCell->custom_data.find_variable_index( "uncoated virion" ); 
	static int nR  = pCell->custom_data.find_variable_index( "viral RNA" ); 
	static int nP  = pCell->custom_data.find_variable_index( "viral protein" ); 
	
	// actual model goes here 

	// now, set apoptosis rate 
	
	static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	phenotype.death.rates[apoptosis_model_index] = pCell->custom_data["max infected apoptosis rate"] ; 
	
	double v = pCell->custom_data[nA_internal] / 
		pCell->custom_data["max apoptosis half max"] ; 
	v = pow( v, pCell->custom_data["apoptosis hill power"] ); 
	
	double effect = v / (1.0+v); 
	phenotype.death.rates[apoptosis_model_index] *= effect; 
	
	
	return; 
}
