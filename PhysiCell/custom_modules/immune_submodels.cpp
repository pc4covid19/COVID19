#include "./immune_submodels.h" 

using namespace PhysiCell; 

std::string immune_submodels_version = "0.0.1"; 

// Submodel_Information Immune_submodels_info; // not needed for now 

Submodel_Information CD8_submodel_info; 
Submodel_Information Macrophage_submodel_info; 
Submodel_Information Neutrophil_submodel_info; 

void CD8_Tcell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// for 
/*	
	if( pCell->state.neighbors.size() > 0 )
	{ std::cout << "adhered Tcell " << pCell << std::endl; } 
*/
	
	return; 
}

void CD8_Tcell_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
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
	
	// if I am not adhered to a cell, turn motility on 
	if( pCell->state.neighbors.size() == 0 )
	{ phenotype.motility.is_motile = true; }
	else
	{ phenotype.motility.is_motile = false; }	
	
	// check for contact with infected cell 
	
	// if I'm adhered to something ... 
	if( pCell->state.neighbors.size() > 0 )
	{
		// add the elastic forces 
		extra_elastic_attachment_mechanics( pCell, phenotype, dt );
		
		// induce damage to whatever we're adhered to 
		#pragma omp critical
		{
			for( int n = 0; n < pCell->state.neighbors.size() ; n++ )
			{
				pCell->state.neighbors[n]->custom_data["TCell_contact_time"] += dt; 
			}
		}

		// decide whether to detach 
		bool detach_me = false; 
		
		if( UniformRandom() < dt / ( pCell->custom_data["cell_attachment_lifetime"] + 1e-15 ) )
		{ detach_me = true; }
		
		// if I detach, go through the process 
		if( detach_me )
		{
			// detach all attached cells 
			for( int n = 0; n < pCell->state.neighbors.size() ; n++ )
			{
				detach_cells( pCell, pCell->state.neighbors[n] ); 
			}
			// resume motile behavior 
			phenotype.motility.is_motile = true; 
		}
		return; 
	}
	
	// I'm not attached, look for cells nearby and try to attach
	
	// if this returns non-NULL, we're now attached to a cell 
	if( immune_cell_check_neighbors_for_attachment( pCell , dt) )
	{
		// set motility off 
		phenotype.motility.is_motile = false; 
		return; 
	}
	phenotype.motility.is_motile = true; // I suggest eliminating this. 
	
	return; 
}

void macrophage_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
//	std::cout << __FUNCTION__ << " " << __LINE__ << std::endl; 
	static int apoptosis_index = phenotype.death.find_death_model_index( "apoptosis" ); 
	static Cell_Definition* pCD = find_cell_definition( "macrophage" ); 
	
	// make changes to volume change rate??

	// if too much debris, comit to apoptosis 	
	
//	std::cout << __FUNCTION__ << " " << __LINE__ << std::endl; 

	double ingested_debris = ( phenotype.volume.total - pCD->phenotype.volume.total ); 
	if( ingested_debris > pCell->custom_data[ "maximum_tolerated_ingested_debris" ] )
	{
		pCell->start_death( apoptosis_index ); 
//		std::cout << " I ate to much and must therefore die " << std::endl; 
//		system("pause"); 
	}

//	std::cout << __FUNCTION__ << " " << __LINE__ << std::endl; 

	// check for cells to eat 
	std::vector<Cell*> neighbors = pCell->cells_in_my_container(); 

//	std::cout << __FUNCTION__ << " " << __LINE__ << std::endl; 
	
	// at least one of the cells is pCell 
	if( neighbors.size() < 2 )
	{ return; } 
	
//	std::cout << "\t\t" << __FUNCTION__ << " " << __LINE__ << std::endl; 

	int n = 0; 
	Cell* pTestCell = neighbors[n]; 
//	std::cout << pCell << " vs " ; 
	while( n < neighbors.size() )
	{
		pTestCell = neighbors[n]; 
//		std::cout << pTestCell << " "; 
		// if it is not me and not a macrophage 
		if( pTestCell != pCell && pTestCell->phenotype.death.dead == true && 
			pTestCell->phenotype.flagged_for_removal == false )
		{
//			std::cout << std::endl; 
//			std::cout << "\t\tnom nom nom" << std::endl; 
//			std::cout << "\t\t\t" << pCell->type << " eats " << pTestCell->type << std::endl; 
//			std::cout << "\t\t\t" << pCell  << " eats " << pTestCell << std::endl; 
			pCell->ingest_cell( pTestCell ); 
			
//			system("pause");
			return; 
		}
//		else
//		{
//			std::cout << " (" << (int) pTestCell->phenotype.death.dead << " " << 
//			(int) pTestCell->phenotype.flagged_for_removal << ") " ; 
//		}
		
		n++; 
	}
//	std::cout << " " << std::endl; 
	
	return; 
}

void macrophage_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	return; 
}

void neutrophil_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	return; 
}

void neutrophil_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	return; 
}

/*
void immune_submodel_setup( void )
{
	immune_submodel_info.name = "immune submodel"; 
	immune_submodel_info.version = "0.0.1";
	immune_submodel_info.main_function= NULL; // receptor_dynamics_model; 
	
	// what variables and parameters do you need? 
	
	immune_submodel_info.cell_variables.push_back( "unbound external ACE2" ); 
	immune_submodel_info.cell_variables.push_back( "bound external ACE2" ); 
	immune_submodel_info.cell_variables.push_back( "unbound internal ACE2" ); 
	immune_submodel_info.cell_variables.push_back( "bound internal ACE2" ); 
	
	immune_submodel_info.cell_variables.push_back( "ACE2 binding rate" ); 
	immune_submodel_info.cell_variables.push_back( "ACE2 endocytosis rate" ); 
	immune_submodel_info.cell_variables.push_back( "ACE2 cargo release rate" ); 	
	immune_submodel_info.cell_variables.push_back( "ACE2 recycling rate" ); 
	
	// submodel_registry.register_model( immune_submodel_info ); 
	immune_submodel_info.register_model(); 
	
	return; 
}
*/
/*
extern Cell_Definition lung_epithelium; 
*/

/*
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
*/

/*
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
*/


void immune_submodels_setup( void )
{
	Cell_Definition* pCD;
	
	// 
	// set up CD8 Tcells
		// set version info 
	CD8_submodel_info.name = "CD8 Tcell model"; 
	CD8_submodel_info.version = immune_submodels_version; 
		// set functions 
	CD8_submodel_info.main_function = NULL; 
	CD8_submodel_info.phenotype_function = CD8_Tcell_phenotype; 
	CD8_submodel_info.mechanics_function = CD8_Tcell_mechanics; 
		// what microenvironment variables do you expect? 
	CD8_submodel_info.microenvironment_variables.push_back( "virion" ); 
	CD8_submodel_info.microenvironment_variables.push_back( "interferon 1" ); 
	CD8_submodel_info.microenvironment_variables.push_back( "pro-inflammatory cytokine" ); 
	CD8_submodel_info.microenvironment_variables.push_back( "chemokine" ); 
		// what custom data do I need? 
	//CD8_submodel_info.cell_variables.push_back( "something" ); 
		// register the submodel  
	CD8_submodel_info.register_model();	
		// set functions for the corresponding cell definition 
	pCD = find_cell_definition( "CD8 Tcell" ); 
	pCD->functions.update_phenotype = CD8_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = CD8_submodel_info.mechanics_function;
	
	// set up macrophages
	Macrophage_submodel_info = CD8_submodel_info; // much shared information 
		// set version info 
	Macrophage_submodel_info.name = "macrophage model"; 
	Macrophage_submodel_info.version = immune_submodels_version; 
		// set functions 
	Macrophage_submodel_info.main_function = NULL; 
	Macrophage_submodel_info.phenotype_function = macrophage_phenotype; 
	Macrophage_submodel_info.mechanics_function = macrophage_mechanics; 
		// what microenvironment variables do you expect? 
	// nothing unique 
		// what custom data do I need? 
	//CD8_submodel_info.cell_variables.push_back( "something" ); 
		// register the submodel  
	Macrophage_submodel_info.register_model();	
		// set functions for the corresponding cell definition 
	pCD = find_cell_definition( "macrophage" ); 
	pCD->functions.update_phenotype = Macrophage_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = Macrophage_submodel_info.mechanics_function;
	
	// set up neutrophils 
	// set up macrophages
	Neutrophil_submodel_info = CD8_submodel_info; // much shared information 
		// set version info 
	Neutrophil_submodel_info.name = "neutrophil model"; 
	Neutrophil_submodel_info.version = immune_submodels_version; 
		// set functions 
	Neutrophil_submodel_info.main_function = NULL; 
	Neutrophil_submodel_info.phenotype_function = neutrophil_phenotype; 
	Neutrophil_submodel_info.mechanics_function = neutrophil_mechanics; 
		// what microenvironment variables do you expect? 
	// nothing unique 
		// what custom data do I need? 
	//CD8_submodel_info.cell_variables.push_back( "something" ); 
		// register the submodel  
	Neutrophil_submodel_info.register_model();	
		// set functions for the corresponding cell definition 
	pCD = find_cell_definition( "neutrophil" ); 
	pCD->functions.update_phenotype = Neutrophil_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = Neutrophil_submodel_info.mechanics_function;	
	
}

Cell* check_for_live_neighbor_for_interaction( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container(); 
	int i = 0; 
	while( i < nearby.size() )
	{
		// don't try to kill yourself 
		if( nearby[i] != pAttacker && nearby[i]->phenotype.death.dead == false )
		{ return nearby[i]; }
		i++; 
	}
	return NULL; 
}

Cell* check_for_dead_neighbor_for_interaction( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container(); 
	int i = 0; 
	while( i < nearby.size() )
	{
		// don't try to kill yourself 
		if( nearby[i] != pAttacker && nearby[i]->phenotype.death.dead == true )
		{ return nearby[i]; }
		i++; 
	}
	return NULL; 
}

void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	axpy( &(pActingOn->velocity) , elastic_constant , displacement ); 
	
	return; 
}

void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	for( int i=0; i < pCell->state.neighbors.size() ; i++ )
	{
		add_elastic_velocity( pCell, pCell->state.neighbors[i], pCell->custom_data["elastic_attachment_coefficient"] ); 
	}
	return; 
}	

void attach_cells( Cell* pCell_1, Cell* pCell_2 )
{
	#pragma omp critical
	{
		bool already_attached = false; 
		for( int i=0 ; i < pCell_1->state.neighbors.size() ; i++ )
		{
			if( pCell_1->state.neighbors[i] == pCell_2 )
			{ already_attached = true; }
		}
		if( already_attached == false )
		{ pCell_1->state.neighbors.push_back( pCell_2 ); }
		
		already_attached = false; 
		for( int i=0 ; i < pCell_2->state.neighbors.size() ; i++ )
		{
			if( pCell_2->state.neighbors[i] == pCell_1 )
			{ already_attached = true; }
		}
		if( already_attached == false )
		{ pCell_2->state.neighbors.push_back( pCell_1 ); }
	}

	return; 
}

void detach_cells( Cell* pCell_1 , Cell* pCell_2 )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < pCell_1->state.neighbors.size() )
		{
			// if cell 2 is in cell 1's list, remove it
			if( pCell_1->state.neighbors[i] == pCell_2 )
			{
				int n = pCell_1->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_1->state.neighbors[i] = pCell_1->state.neighbors[n-1]; 
				// shrink by one 
				pCell_1->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}
	
		found = false; 
		i = 0; 
		while( !found && i < pCell_2->state.neighbors.size() )
		{
			// if cell 1 is in cell 2's list, remove it
			if( pCell_2->state.neighbors[i] == pCell_1 )
			{
				int n = pCell_2->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_2->state.neighbors[i] = pCell_2->state.neighbors[n-1]; 
				// shrink by one 
				pCell_2->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}
	}
	
	return; 
}



bool attempt_immune_cell_attachment( Cell* pAttacker, Cell* pTarget , double dt )
{
	// if the target is not infected, give up 
	if( pTarget->custom_data[ "assembled_virion" ] < 1 )
	{ return false; }
		
	// if the target cell is dead, give up 
	if( pTarget->phenotype.death.dead == true )
	{ return false; } 

	// if the target cell is too far away, give up 
	std::vector<double> displacement = pTarget->position - pAttacker->position;
	double distance_scale = norm( displacement ); 
	if( distance_scale > pAttacker->custom_data["max_attachment_distance"] )
	{ return false; } 

	// now, get the attachment probability 
	
	double attachment_probability = pAttacker->custom_data["cell_attachment_rate"] * dt; 

	// don't need to cap it at 1.00: if prob > 100%, 
	// then this statement always evaluates as true, 
	// just the same as capping probability at 100% 
	if( UniformRandom() <= attachment_probability )
	{
		attach_cells( pAttacker, pTarget ); 
		return true; 
	}
	
	return false; 	
}

Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container(); 
	int i = 0; 
	while( i < nearby.size() )
	{
		// don't try to kill yourself 
		if( nearby[i] != pAttacker )
		{
			if( attempt_immune_cell_attachment( pAttacker, nearby[i] , dt ) )
			{
				return nearby[i]; 
			}
		}
		i++; 
	}
	
	return NULL; 
}

void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int apoptosis_index = phenotype.death.find_death_model_index( "apoptosis" ); 
	if( pCell->custom_data["TCell_contact_time"] > pCell->custom_data["TCell_contact_death_threshold"] )
	{
		// std::cout << "I die now" << std::endl; 
		
		// make sure to get rid of all adhesions! 
		// detach all attached cells 
		for( int n = 0; n < pCell->state.neighbors.size() ; n++ )
		{
			detach_cells( pCell, pCell->state.neighbors[n] ); 
		}
		
		pCell->start_death( apoptosis_index ); 
		pCell->functions.update_phenotype = NULL; 
	}
	
	return; 
}

