#include "./receptor_dynamics.h" 

using namespace PhysiCell; 

Submodel_Information receptor_dynamics_info; 

void receptor_dynamics_model_setup( void )
{
	receptor_dynamics_info.name = "receptor dynamics"; 
	receptor_dynamics_info.version = "0.2.0";
	receptor_dynamics_info.main_function= receptor_dynamics_model; 
	
	// what variables and parameters do you need? 
	
	receptor_dynamics_info.cell_variables.push_back( "unbound external ACE2" ); 
	receptor_dynamics_info.cell_variables.push_back( "bound external ACE2" ); 
	receptor_dynamics_info.cell_variables.push_back( "unbound internal ACE2" ); 
	receptor_dynamics_info.cell_variables.push_back( "bound internal ACE2" ); 
	
	receptor_dynamics_info.cell_variables.push_back( "ACE2 binding rate" ); 
	receptor_dynamics_info.cell_variables.push_back( "ACE2 endocytosis rate" ); 
	receptor_dynamics_info.cell_variables.push_back( "ACE2 cargo release rate" ); 	
	receptor_dynamics_info.cell_variables.push_back( "ACE2 recycling rate" ); 
	
	// submodel_registry.register_model( receptor_dynamics_info ); 
	receptor_dynamics_info.register_model(); 
	
	return; 
}

extern Cell_Definition lung_epithelium; 

void receptor_dynamics_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	// bookkeeping -- find microenvironment variables we need

	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
	
	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
	static int nA_internal = pCell->custom_data.find_variable_index( "assembled virion" ); 

	// bookkeeping -- find custom data we need 
	
	static int nR_EU = pCell->custom_data.find_variable_index( "unbound external ACE2" ); 
	static int nR_EB = pCell->custom_data.find_variable_index( "bound external ACE2" ); 
	static int nR_IU = pCell->custom_data.find_variable_index( "unbound internal ACE2" ); 
	static int nR_IB = pCell->custom_data.find_variable_index( "bound internal ACE2" ); 
	
	static int nR_bind = pCell->custom_data.find_variable_index( "ACE2 binding rate" ); 
	static int nR_endo = pCell->custom_data.find_variable_index( "ACE2 endocytosis rate" ); 
	static int nR_release = pCell->custom_data.find_variable_index( "ACE2 cargo release rate" ); 	
	static int nR_recycle = pCell->custom_data.find_variable_index( "ACE2 recycling rate" ); 
	
	// do nothing if dead 
	if( phenotype.death.dead == true )
	{ return; } 

	// if not lung epithelium, do nothing 
	if( pCell->type != lung_epithelium.type )
	{ return; } 
	
	// actual model goes here 
	
	// internalized virus tells us how many have recently bound to receptors 
	double newly_bound = phenotype.molecular.internalized_total_substrates[nV_external]; 
	// if it tried to bind to more virus than there are receptors, compensate 
	double excess_binding = newly_bound - pCell->custom_data[nR_EU]; 
	if( excess_binding > 0.0 )
	{
		// don't bring in more virus than there are receptors 
		newly_bound = pCell->custom_data[nR_EU]; 
		// dump any excess back into the microenvironment
		static double one_virion_to_density = 1.0 / microenvironment.mesh.dV; 
		// this needs omp critical because 2 cells writing to 1 voxel is not thread safe 
		#pragma omp critical 
		{
			pCell->nearest_density_vector()[nV_external] += excess_binding * one_virion_to_density; 
		}
	}
	phenotype.molecular.internalized_total_substrates[nV_external] = 0.0; 
	
	// add newly bound receptor to R_EB
	
	pCell->custom_data[nR_EB] += newly_bound; 
	
	// remove newly bound receptor from R_EU 

	pCell->custom_data[nR_EU] -= newly_bound; 
	
	// endocytosis 
	
	double dR_IB = dt*pCell->custom_data[nR_endo]*pCell->custom_data[nR_EB];
	if( dR_IB > pCell->custom_data[nR_EB] )
	{ dR_IB = pCell->custom_data[nR_EB]; }
	pCell->custom_data[nR_EB] -= dR_IB; // move from external bound
	pCell->custom_data[nR_IB] += dR_IB; // move to internal bound
	
	// viral release from endosomes 
	
	double dR_IU = dt*pCell->custom_data[nR_release]*pCell->custom_data[nR_IB];
	if( dR_IU > pCell->custom_data[nR_IB] )
	{ dR_IU = pCell->custom_data[nR_IB]; }
	pCell->custom_data[nR_IB] -= dR_IU; // move from internal bound 
	pCell->custom_data[nR_IU] += dR_IU; // move to internal unbound 
	pCell->custom_data[nV_internal] += dR_IU; // release virus into cytoplasm 
	
	// receptor recycling 
	
	double dR_EU = dt*pCell->custom_data[nR_recycle]*pCell->custom_data[nR_IU];
	if( dR_EU > pCell->custom_data[nR_IU] )
	{ dR_EU = pCell->custom_data[nR_IU]; }
	pCell->custom_data[nR_IU] -= dR_EU; // move from internal unbound 
	pCell->custom_data[nR_EU] += dR_EU; // move to external unbound 
	
	// update the virion uptake rate 
	
	phenotype.secretion.uptake_rates[nV_external] = 
		pCell->custom_data[nR_bind] * pCell->custom_data[nR_EU]; 
	
	return; 
}

void receptor_dynamics_model( double dt )
{
	#pragma omp parallel for 
	for( int n=0; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		receptor_dynamics_model( pC, pC->phenotype , dt ); 	
	}
	
	return; 
}
	
