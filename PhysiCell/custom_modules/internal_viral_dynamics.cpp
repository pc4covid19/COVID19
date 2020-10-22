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
	internal_viral_dynamics_info.cell_variables.push_back( "max_RNA_replication_rate" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "RNA_replication_half" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "basal_RNA_degradation_rate" ); 

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

  /*
  pCell->custom_data[nA_internal] = newly assembled virion 
  phenotype.molecular.internalized_total_substrates[nA_external] = all assembled virion that is not yet exported 
  */
  // This resets the internal assembled virion count 
  // so we are commenting it out 
	// pCell->custom_data[nA_internal] = 
	// 	phenotype.molecular.internalized_total_substrates[nA_external]; 
		
	// actual model goes here 

	// uncoat endocytosed virus
	double dV = dt * pCell->custom_data["virion_uncoating_rate"] * pCell->custom_data[nV_internal] ;
	if( dV > pCell->custom_data[nV_internal] )
	{ dV = pCell->custom_data[nV_internal]; } 
	pCell->custom_data[nV_internal] -= dV; 
	pCell->custom_data[nUV] += dV; 

	// convert uncoated virus to usable mRNA 
	double dR = dt * pCell->custom_data["uncoated_to_RNA_rate"] * pCell->custom_data[nUV]; 
  // RNA replication post uncoated to RNA calc
	dR += dt * pCell->custom_data["max_RNA_replication_rate"] * pCell->custom_data[nR] /
              (pCell->custom_data[nR] + pCell->custom_data["RNA_replication_half"]);
  // RNA degradation
  dR -= dt * pCell->custom_data["basal_RNA_degradation_rate"] * pCell->custom_data[nR];
	// if( dR > pCell->custom_data[nUV] )
	// { dR = pCell->custom_data[nUV]; }
	pCell->custom_data[nUV] -= dR; 
	pCell->custom_data[nR] += dR; 
	
	// use mRNA to create viral protein 
	double dP = dt * pCell->custom_data["protein_synthesis_rate"] * pCell->custom_data[nR];
	pCell->custom_data[nP] += dP; 

	// degrade protein 
	
	// assemble virus 
	double dA = dt * pCell->custom_data["virion_assembly_rate"] * pCell->custom_data[nP]; 
	if( dA > pCell->custom_data[nP] )
	{ dA = pCell->custom_data[nP]; } 
	pCell->custom_data[nP] -= dA; 
	pCell->custom_data[nA_internal] += dA; 
	
	// set export rate 
  double nA_export_rate = dt * pCell->custom_data["virion_export_rate"] * pCell->custom_data[nA_internal];
  // next step is to get a discrete number here
  // floor + a probability of exporting an additional one
  double nA_export_count = floor( nA_export_rate );
  // TODO: Add the stochastic additional one here to nA_export_count
  // when doing so make sure the internal nA count doesn't drop below
  // zero by comparing the two
  // this is the export count to add/remove
  pCell->custom_data[nA_internal] -= nA_export_count;
  phenotype.molecular.internalized_total_substrates[nA_external] += nA_export_count;
//6. adjust phenotype.secretion.net_export_rates[nA_external] so taht dt * export_rate * phenotype.molecular.internalized_total_substrates[nA_external] = number of particles you wanted to export 
// dt * exp_rate * internal = part_to_exp
  // is this correct at all? 
  phenotype.secretion.net_export_rates[nA_external] = nA_export_count/dt;
	
	// phenotype.secretion.net_export_rates[nA_external] = 
	// 	pCell->custom_data["virion_export_rate" ] * pCell->custom_data[nA_internal]; 
 
	// copy data from custom variables to "internalized variables" 
	
/*	
	phenotype.molecular.internalized_total_substrates[nV_external] = 
		pCell->custom_data[nV_internal];
*/		

  /* Handled above now
	phenotype.molecular.internalized_total_substrates[nA_external] = 
		pCell->custom_data[nA_internal];	
  */
	
	return; 
}
