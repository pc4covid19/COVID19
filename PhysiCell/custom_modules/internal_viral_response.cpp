#include "./internal_viral_response.h" 

using namespace PhysiCell; 

std::string internal_virus_response_version = "0.3.0"; 

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
	static Cell_Definition* pCD = find_cell_definition( "lung epithelium" ); 
	
	// bookkeeping -- find microenvironment variables we need

	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
	
	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
	static int nA_internal = pCell->custom_data.find_variable_index( "assembled_virion" ); 
	
	static int nINF1 = microenvironment.find_density_index( "interferon 1" ); 
	
	// actual model goes here 

	// now, set apoptosis rate 
	
	static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "apoptosis" );
	// phenotype.death.rates[apoptosis_model_index] = 
	
	// base death rate (from cell line)
	double base_death_rate = 
		pCD->phenotype.death.rates[apoptosis_model_index]; 
	
	// additional death rate from infectoin  
	double additional_death_rate = pCell->custom_data["max_infected_apoptosis_rate"] ; 
	
	
	double v = pCell->custom_data[nA_internal] / 
		pCell->custom_data["max_apoptosis_half_max"] ; 
	v = pow( v, pCell->custom_data["apoptosis_hill_power"] ); 
	
	double effect = v / (1.0+v); 
	additional_death_rate *= effect; 
	phenotype.death.rates[apoptosis_model_index] = base_death_rate + additional_death_rate; 
	
	// if we're infected, secrete a chemokine for the immune model
	static int nAV = pCell->custom_data.find_variable_index( "assembled_virion" ); 	
	double AV = pCell->custom_data[nAV]; 
	
	static int nP = pCell->custom_data.find_variable_index( "viral_protein"); 
	double P = pCell->custom_data[nP];

	static int chemokine_index = microenvironment.find_density_index( "chemokine" ); 
	
/* old 
	if( P > 0.001 )
	{
		phenotype.secretion.secretion_rates[chemokine_index] = 
			pCell->custom_data[ "infected_cell_chemokine_secretion_rate" ];
		phenotype.secretion.saturation_densities[chemokine_index] = 1.0; 		
	}
	else
	{
		phenotype.secretion.secretion_rates[chemokine_index] = 0.0;
	}
*/	
	
	// if I am dead, make sure to still secrete the chemokine 
	
	// static int chemokine_index = microenvironment.find_density_index( "chemokine" ); 
	// static int nP = pCell->custom_data.find_variable_index( "viral_protein"); 
	// double P = pCell->custom_data[nP];
	
	// static int nAV = pCell->custom_data.find_variable_index( "assembled_virion" ); 
	// double AV = pCell->custom_data[nAV]; 
	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
		
	static int nR = pCell->custom_data.find_variable_index( "viral_RNA");
	double R = pCell->custom_data[nR];
	
	
	
	
	
	if( R >= 1.00 - 1e-16 ) 
	{
		pCell->custom_data["infected_cell_chemokine_secretion_activated"] = 1.0; 
	}

	if( pCell->custom_data["infected_cell_chemokine_secretion_activated"] > 0.1 && phenotype.death.dead == false )
	{
		double rate = AV; 
		rate /= pCell->custom_data["max_apoptosis_half_max"];
		if( rate > 1.0 )
		{ rate = 1.0; }
		rate *= pCell->custom_data[ "infected_cell_chemokine_secretion_rate" ];

		phenotype.secretion.secretion_rates[chemokine_index] = rate; 
		phenotype.secretion.saturation_densities[chemokine_index] = 1.0;

		// (Adrianne) adding pro-inflammatory cytokine secretion by infected cells
		pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = pCell->custom_data["activated_cytokine_secretion_rate"]/10;
	}
	
	// interferon signaling 

	// // approximate activation by extracellular interferons 
	// // // activation = min( 1 , extracellular_interferon / interferon_max_response_threshold ) 
	pCell->custom_data["interferon_activation"] = pCell->nearest_density_vector()[nINF1] / 
		( pCell->custom_data["interferon_max_response_threshold"] + 1e-32 ); 
	if( pCell->custom_data["interferon_activation"] > 1.0 )
	{ pCell->custom_data["interferon_activation"] = 1.0; } 
		
	// // Type-I interferon secretion 
	// // // secretion_rate = r_viral * Heaviside( RNA - 1 ) + r_paracrine * activation 
	phenotype.secretion.secretion_rates[nINF1] = pCell->custom_data["interferon_activation"]; 
	phenotype.secretion.secretion_rates[nINF1] *= pCell->custom_data["max_interferon_secretion_rate_via_paracrine"]; 
	if( R >= 1.0 - 1e-16 ) // if there is at least 1 complete set of uncoated viral RNA
	{
		double scaled_RNA = R / ( pCell->custom_data["interferon_viral_RNA_threshold"] + 1e-32 );
		if( scaled_RNA > 1 )
		{ scaled_RNA = 1.0; }
		phenotype.secretion.secretion_rates[nINF1] += pCell->custom_data["interferon_secretion_rate_via_infection"] * scaled_RNA; 
	} 
	
	// // now the interferon response 
	// // // protein_synthesis_rate = protein_synthesis_rate_0 * ( 1 - interferon_activation * interferon_max_virus_inhibition ) 
	pCell->custom_data["protein_synthesis_rate"] = pCell->custom_data["interferon_max_virus_inhibition"]; // inhibition
	pCell->custom_data["protein_synthesis_rate"] *= pCell->custom_data["interferon_activation"]; // activation*inhibition
	pCell->custom_data["protein_synthesis_rate"] *= -1; // -activation*inhibition 
	pCell->custom_data["protein_synthesis_rate"] += 1.0; // 1 - activation*inhibition 
	pCell->custom_data["protein_synthesis_rate"] *= pCD->custom_data["protein_synthesis_rate"]; 
		// protein_synthesis_rate0 * (1 - activation*inhibition)
	
	return; 
}

