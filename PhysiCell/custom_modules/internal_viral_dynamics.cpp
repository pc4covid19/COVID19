#include "./internal_viral_dynamics.h" 

using namespace PhysiCell; 

Submodel_Information internal_viral_dynamics_info; 

void internal_virus_model_setup( void )
{
	internal_viral_dynamics_info.name = "internal viral dynamics"; 
	internal_viral_dynamics_info.version = "0.2.0";
	internal_viral_dynamics_info.main_function= internal_virus_model; 

	// what custom data do I need? 
	
	internal_viral_dynamics_info.cell_variables.push_back( "virion" ); // adhered, in process of endocytosis 
	internal_viral_dynamics_info.cell_variables.push_back( "uncoated virion" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "viral RNA" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "viral protein" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "assembled virion" ); 

	internal_viral_dynamics_info.cell_variables.push_back( "virion uncoating rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "uncoated to RNA rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "protein synthesis rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "virion assembly rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "virion export rate" ); 

	// submodel_registry.register_model( internal_viral_dynamics_info ); 
	internal_viral_dynamics_info.register_model();
	
	return; 
}

void internal_virus_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	// bookkeeping -- find microenvironment variables we need

	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
	
	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
	static int nA_internal = pCell->custom_data.find_variable_index( "assembled virion" ); 

	static int nUV = pCell->custom_data.find_variable_index( "uncoated virion" ); 
	static int nR  = pCell->custom_data.find_variable_index( "viral RNA" ); 
	static int nP  = pCell->custom_data.find_variable_index( "viral protein" ); 
	
	// copy virions from "internalized variables" to "custom variables"
/*	
	pCell->custom_data[nV_internal] = 
		phenotype.molecular.internalized_total_substrates[nV_external]; 
	// this transfer is now handled in receptor dynamics 
*/		
	pCell->custom_data[nA_internal] = 
		phenotype.molecular.internalized_total_substrates[nA_external]; 
		
	// actual model goes here 

	// uncoat endocytosed virus
	double dV = dt * pCell->custom_data["virion uncoating rate"] * pCell->custom_data[nV_internal] ;
	if( dV > pCell->custom_data[nV_internal] )
	{ dV = pCell->custom_data[nV_internal]; } 
	pCell->custom_data[nV_internal] -= dV; 
	pCell->custom_data[nUV] += dV; 

	// convert uncoated virus to usable mRNA 
	double dR = dt * pCell->custom_data["uncoated to RNA rate"] * pCell->custom_data[nUV]; 
	if( dR > pCell->custom_data[nUV] )
	{ dR = pCell->custom_data[nUV]; }
	pCell->custom_data[nUV] -= dR; 
	pCell->custom_data[nR] += dR; 
	
	// use mRNA to create viral protein 
	double dP = dt * pCell->custom_data["protein synthesis rate"] * pCell->custom_data[nR];
	pCell->custom_data[nP] += dP; 

	// degrade mRNA 
	

	// degrade protein 
	
	
	// assemble virus 
	double dA = dt * pCell->custom_data["virion assembly rate"] * pCell->custom_data[nP]; 
	if( dA > pCell->custom_data[nP] )
	{ dA = pCell->custom_data[nP]; } 
	pCell->custom_data[nP] -= dA; 
	pCell->custom_data[nA_internal] += dA; 
	
	// set export rate 
	
	phenotype.secretion.net_export_rates[nA_external] = 
		pCell->custom_data["virion export rate" ] * pCell->custom_data[nA_internal]; 
 
	// copy data from custom variables to "internalized variables" 
	
/*	
	phenotype.molecular.internalized_total_substrates[nV_external] = 
		pCell->custom_data[nV_internal];
*/		
	phenotype.molecular.internalized_total_substrates[nA_external] = 
		pCell->custom_data[nA_internal];	
	
	return; 
}
	
