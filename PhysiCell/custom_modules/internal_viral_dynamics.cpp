#include "./internal_viral_dynamics.h" 

using namespace PhysiCell; 

std::string internal_virus_replication_version = "0.3.0"; 

Submodel_Information internal_viral_dynamics_info; 


void internal_virus_model_setup( void )
{
		// set version
	internal_viral_dynamics_info.name = "internal viral replication dynamics"; 
	internal_viral_dynamics_info.version = internal_virus_replication_version; 
		// set functions 
	internal_viral_dynamics_info.main_function = NULL; 
	internal_viral_dynamics_info.phenotype_function = internal_virus_model; 
	internal_viral_dynamics_info.mechanics_function = NULL; 
		// what microenvironment variables do I need? 

		// what custom data do I need? 
	internal_viral_dynamics_info.cell_variables.push_back( "virion" ); // adhered, in process of endocytosis 
	internal_viral_dynamics_info.cell_variables.push_back( "uncoated_virion" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "viral_RNA" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "viral_protein" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "assembled_virion" ); 

	internal_viral_dynamics_info.cell_variables.push_back( "virion_uncoating_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "uncoated_to_RNA_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "protein_synthesis_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "virion_assembly_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "virion_export_rate" ); 

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
	static int nA_internal = pCell->custom_data.find_variable_index( "assembled_virion" ); 

	static int nUV = pCell->custom_data.find_variable_index( "uncoated_virion" ); 
	static int nR  = pCell->custom_data.find_variable_index( "viral_RNA" ); 
	static int nP  = pCell->custom_data.find_variable_index( "viral_protein" ); 

/*	
	static bool done = false; 
	extern Cell* pInfected; 
	if( pCell == pInfected && 1 == 0 )
	{
		std::cout << std::endl << "viral dynamics : " << __LINE__ << " " 
			<< phenotype.molecular.internalized_total_substrates[ nV_external ] << " " 
			<< phenotype.molecular.internalized_total_substrates[ nA_external ] << " " 
			<< pCell->custom_data[nV_internal] << " " 
			<< pCell->custom_data[nUV] << " " 
			<< pCell->custom_data[nR] << " " 
			<< pCell->custom_data[nP] << " " 	
			<< pCell->custom_data[nA_internal] << " " 
			<< std::endl; 		
	}
*/	
	
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
	double dV = dt * pCell->custom_data["virion_uncoating_rate"] * pCell->custom_data[nV_internal] ;
	if( dV > pCell->custom_data[nV_internal] )
	{ dV = pCell->custom_data[nV_internal]; } 
	pCell->custom_data[nV_internal] -= dV; 
	pCell->custom_data[nUV] += dV; 

	// convert uncoated virus to usable mRNA 
	double dR = dt * pCell->custom_data["uncoated_to_RNA_rate"] * pCell->custom_data[nUV]; 
	if( dR > pCell->custom_data[nUV] )
	{ dR = pCell->custom_data[nUV]; }
	pCell->custom_data[nUV] -= dR; 
	pCell->custom_data[nR] += dR; 
	
	// use mRNA to create viral protein 
	double dP = dt * pCell->custom_data["protein_synthesis_rate"] * pCell->custom_data[nR];
	pCell->custom_data[nP] += dP; 

	// degrade mRNA 
	

	// degrade protein 
	
	
	// assemble virus 
	double dA = dt * pCell->custom_data["virion_assembly_rate"] * pCell->custom_data[nP]; 
	if( dA > pCell->custom_data[nP] )
	{ dA = pCell->custom_data[nP]; } 
	pCell->custom_data[nP] -= dA; 
	pCell->custom_data[nA_internal] += dA; 
	
	// set export rate 
	
	phenotype.secretion.net_export_rates[nA_external] = 
		pCell->custom_data["virion_export_rate" ] * pCell->custom_data[nA_internal]; 
 
	// copy data from custom variables to "internalized variables" 
	
/*	
	phenotype.molecular.internalized_total_substrates[nV_external] = 
		pCell->custom_data[nV_internal];
*/		
	phenotype.molecular.internalized_total_substrates[nA_external] = 
		pCell->custom_data[nA_internal];	
	
	return; 
}
	
