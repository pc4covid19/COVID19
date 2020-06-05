#include "./immune_submodels.h" 

using namespace PhysiCell; 

std::string immune_submodels_version = "0.0.1"; 

// Submodel_Information Immune_submodels_info; // not needed for now 

Submodel_Information CD8_submodel_info; 
Submodel_Information Macrophage_submodel_info; 
Submodel_Information Neutrophil_submodel_info; 

void remove_all_adhesions( Cell* pCell )
{
	// detach all attached cells 
	for( int n = 0; n < pCell->state.neighbors.size() ; n++ )
	{ detach_cells( pCell, pCell->state.neighbors[n] ); }		
	
	return; 
}

void create_infiltrating_immune_cell( Cell_Definition* pCD )
{
	static double Xmin = microenvironment.mesh.bounding_box[0]; 
	static double Ymin = microenvironment.mesh.bounding_box[1]; 
	static double Zmin = microenvironment.mesh.bounding_box[2]; 

	static double Xmax = microenvironment.mesh.bounding_box[3]; 
	static double Ymax = microenvironment.mesh.bounding_box[4]; 
	static double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	static bool setup_done = false; 
	
	// warning hardcoded
	static double relative_edge_margin = 0.0; // 0.1; 
	static double relative_interior = 1 - 2 * relative_edge_margin; 
	
	if( default_microenvironment_options.simulate_2D == true && setup_done == false )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	static double Xrange = (Xmax - Xmin); 
	static double Yrange = (Ymax - Ymin); 
	static double Zrange = (Zmax - Zmin); 
	
	// keep cells away from the outer edge 
	
	if( setup_done == false )
	{
		Xmin += relative_edge_margin*Xrange; 
		Ymin += relative_edge_margin*Yrange; 
		Zmin += relative_edge_margin*Zrange;
		
		Xrange *= relative_interior;
		Yrange *= relative_interior;
		Zrange *= relative_interior; 
		setup_done = true; 
	}
	
	std::vector<double> position = {0,0,0}; 
	position[0] = Xmin + UniformRandom()*Xrange; 
	position[1] = Ymin + UniformRandom()*Yrange; 
	position[2] = Zmin + UniformRandom()*Zrange + parameters.doubles("immune_z_offset"); 
		
	Cell* pC = create_cell( *pCD ); 
	pC->assign_position( position );
	
	return; 
}

void create_infiltrating_immune_cell( std::string cell_name )
{
	create_infiltrating_immune_cell( find_cell_definition( cell_name ) ); 
	return;
}

void create_infiltrating_neutrophil(void)
{
	static Cell_Definition* pCD = find_cell_definition( "neutrophil" );
	create_infiltrating_immune_cell( pCD ); 
	return;
}

void create_infiltrating_Tcell(void)
{
	static Cell_Definition* pCD = find_cell_definition( "CD8 Tcell" );
	create_infiltrating_immune_cell( pCD ); 

	return; 
}

void CD8_Tcell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	return; 
}

void CD8_Tcell_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if I'm dead, don't bother 
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		
		// detach all attached cells 
		remove_all_adhesions( pCell ); 
		
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
		#pragma omp critical(track_contact_time)
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

void keep_immune_cells_off_edge( void )
{
	static double Xmin = microenvironment.mesh.bounding_box[0]; 
	static double Ymin = microenvironment.mesh.bounding_box[1]; 
	static double Zmin = microenvironment.mesh.bounding_box[2]; 

	static double Xmax = microenvironment.mesh.bounding_box[3]; 
	static double Ymax = microenvironment.mesh.bounding_box[4]; 
	static double Zmax = microenvironment.mesh.bounding_box[5]; 

	static bool setup_done = false; 
	if( default_microenvironment_options.simulate_2D == true && setup_done == false )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	static double Xrange = (Xmax - Xmin); 
	static double Yrange = (Ymax - Ymin); 
	static double Zrange = (Zmax - Zmin); 
	
	// warning hardcoded
	static double relative_edge_margin = 0; // 0.1; 
	static double relative_interior = 1 - 2 * relative_edge_margin; 
	
	if( setup_done == false )
	{
		Xmin += relative_edge_margin*Xrange; 
		Ymin += relative_edge_margin*Yrange; 
		Zmin += relative_edge_margin*Zrange;
		
		Xrange *= relative_interior;
		Yrange *= relative_interior;
		Zrange *= relative_interior;  
		setup_done = true; 
	}
	
	static int epithelial_type = get_cell_definition( "lung epithelium" ).type; 

	for( int n=0 ; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if( pC->phenotype.death.dead == false && pC->is_out_of_domain && pC->type != epithelial_type )
		{
			
			pC->is_out_of_domain = false; 
			pC->is_active = true; 
			pC->is_movable = true; 			
			
			
			std::vector<double> position = pC->position; 
			position[0] = Xmin + Xrange * UniformRandom(); 
			position[1] = Ymin + Yrange * UniformRandom(); 
			position[2] = Zmin + Zrange * UniformRandom() + parameters.doubles("immune_z_offset"); 

			#pragma omp critical(move_from_edge)
			{
				std::cout << " moving cell " << pC << " of type " << pC->type_name << std::endl; 
				std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
				pC->assign_position( position ); 	
				std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
				// pC->update_voxel_in_container(); // cut this? 
				std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
			}
		}
	}
	return; 
/*	
	// keep cells away from the outer edge 
	
	// check for out of bounds 
	std::vector<double> position = pCell->position; 
	static std::vector<double>* pBB = &(microenvironment.mesh.bounding_box); 
	if( position[0] < (*pBB)[0] || position[0] > (*pBB)[3] || 
		position[1] < (*pBB)[1] || position[1] > (*pBB)[4] || 
		position[2] < (*pBB)[2] || position[2] > (*pBB)[5] )
	{
		position[0] = Xmin + Xrange * UniformRandom(); 
		position[1] = Ymin + Yrange * UniformRandom(); 
		position[2] = Zmin + Zrange * UniformRandom(); 

		pCell->assign_position( position ); 
		return; 
	}			
*/
	return;
}

void keep_immune_cells_in_bounds( double dt )
{
	static double dt_bounds = 5; 
	static double next_time = 0.0; 

	static double t_bounds = 0.0; 
	static double t_last_bounds = 0.0; 
	static double t_next_bounds = 0.0; 
	
	static double tolerance = 0.1 * diffusion_dt; 
	
	// is it time for the next immune recruitment? 
	if( t_bounds > t_next_bounds- tolerance )
	{
		double elapsed_time = (t_bounds - t_last_bounds );
		
		keep_immune_cells_off_edge(); 
		
		t_last_bounds = t_bounds; 
		t_next_bounds = t_bounds + dt_bounds; 
	}
	t_bounds += dt; 

	return; 
}

void macrophage_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	
/*	
	// check for out of bounds 
	if( pCell->is_out_of_domain == true ) 
	{ 
		// std::cout << pCell->type_name << " is out of bounds" << std::endl; 
		
		std::vector<double> position = {0,0,0}; 
		position[0] = Xmin + Xrange * UniformRandom(); 
		position[1] = Ymin + Yrange * UniformRandom(); 
		position[2] = Zmin + Zrange * UniformRandom(); 

		pCell->assign_position( position ); 
		return; 
	}
*/	
	
//	std::cout << __FUNCTION__ << " " << __LINE__ << std::endl; 
	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" ); 
	static Cell_Definition* pCD = find_cell_definition( "macrophage" ); 
	
	// make changes to volume change rate??

	// if too much debris, comit to apoptosis 	
	
//	std::cout << __FUNCTION__ << " " << __LINE__ << std::endl; 

	double ingested_debris = ( phenotype.volume.total - pCD->phenotype.volume.total ); 
	if( ingested_debris > pCell->custom_data[ "maximum_tolerated_ingested_debris" ] )
	{
		pCell->start_death( apoptosis_index ); 
		return; 
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

	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");

	int n = 0; 
	Cell* pTestCell = neighbors[n]; 
//	std::cout << pCell << " vs " ; 
	while( n < neighbors.size() )
	{
		pTestCell = neighbors[n]; 
//		std::cout << pTestCell << " "; 
		// if it is not me and not a macrophage 
		if( pTestCell != pCell && pTestCell->phenotype.death.dead == true && 
			pTestCell->phenotype.volume.total > 1e-15 )
//			pTestCell->phenotype.flagged_for_removal == false )
		{
//			std::cout << std::endl; 
			#pragma omp critical(macrophage_eat) 
			{
				std::cout << "\t\tnom nom nom" << std::endl; 
				std::cout << "\t\t\t" << pCell->type_name << " eats " << pTestCell->type_name << std::endl; 
				std::cout << "\t\t\t" << pCell  << " eats " << pTestCell << std::endl; 
				pCell->ingest_cell( pTestCell ); 
				remove_all_adhesions( pTestCell ); // debug 
				
				std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
				
				phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 
					pCell->custom_data["activated_macrophage_secretion_rate"]; // 10;

				std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;

				phenotype.motility.migration_speed = 
					pCell->custom_data[ "activated_macrophage_speed" ]; 
					
				std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;

				// warning : hardcoded 
				phenotype.motility.migration_bias = 0.5; 
				
				std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
			}
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
	#pragma omp critical(attach)
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
	#pragma omp critical(detach)
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

void immune_cell_recruitment( double dt )
{
	static int proinflammatory_cytokine_index = 
		microenvironment.find_density_index("pro-inflammatory cytokine");
	
	static double dt_immune = parameters.doubles( "immune_dt" ); 
	static double t_immune = 0.0; 
	static double t_last_immune = 0.0; 
	static double t_next_immune = 0.0; 
	
	static double tolerance = 0.1 * diffusion_dt; 
	
	// is it time for the next immune recruitment? 
	if( t_immune > t_next_immune- tolerance )
	{
		double elapsed_time = (t_immune - t_last_immune );
//		std::cout << "Immune time! " << t_immune << " (elapsed: " << elapsed_time << ") " << std::endl; 
		
		// neutrophil recruitment 
		
		static double neutrophil_recruitment_rate = parameters.doubles( "neutrophil_max_recruitment_rate" ); 
		static double NR_min_signal = parameters.doubles( "neutrophil_recruitment_min_signal" ); 
		static double NR_sat_signal = parameters.doubles( "neutrophil_recruitment_saturation_signal" ); 
		static double NR_max_minus_min = NR_sat_signal - NR_min_signal; 
		
		double total_rate = 0;
		// integrate \int_domain r_max * (signal-signal_min)/(signal_max-signal_min) * dV 
		double total_scaled_signal= 0.0;
		for( int n=0; n<microenvironment.mesh.voxels.size(); n++ )
		{
			// (signal(x)-signal_min)/(signal_max/signal_min)
			double dRate = ( microenvironment(n)[proinflammatory_cytokine_index] - NR_min_signal ); 
			dRate /= NR_max_minus_min; 
			// crop to [0,1] 
			if( dRate > 1 ) 
			{ dRate = 1; } 
			if( dRate < 0 )
			{ dRate = 0; }
			total_rate += dRate; 
		}	
		// multiply by dV and rate_max 
		total_scaled_signal = total_rate; 
		
		total_rate *= microenvironment.mesh.dV; 
		total_rate *= neutrophil_recruitment_rate; 
		
		// expected number of new neutrophils 
		int number_of_new_cells = (int) round( total_rate * elapsed_time ); 
		if( number_of_new_cells )
		{
			std::cout << "\tRecruiting " << number_of_new_cells << " neutrophils ... " << std::endl; 
			
//			std::cout << "\tTotal signal/dV : " << total_scaled_signal << std::endl;
//			std::cout << "\tTotal signa : " << total_scaled_signal * microenvironment.mesh.dV << std::endl; 
//			double total_volume = microenvironment.mesh.dV * microenvironment.mesh.voxels.size() ; 
//			std::cout << "\tmean signal : " << total_scaled_signal * microenvironment.mesh.dV / total_volume << std::endl; 
			for( int n = 0; n < number_of_new_cells ; n++ )
			{ create_infiltrating_neutrophil(); }
		}
		
		// CD8 T cell recruitment 
		
		static double CD8_Tcell_recruitment_rate = parameters.doubles( "CD8_Tcell_max_recruitment_rate" ); 
		static double TC_min_signal = parameters.doubles( "CD8_Tcell_recruitment_min_signal" ); 
		static double TC_sat_signal = parameters.doubles( "CD8_Tcell_recruitment_saturation_signal" ); 
		static double TC_max_minus_min = TC_sat_signal - TC_min_signal; 
		
		total_rate = 0;
		// integrate \int_domain r_max * (signal-signal_min)/(signal_max-signal_min) * dV 
		total_scaled_signal= 0.0;
		for( int n=0; n<microenvironment.mesh.voxels.size(); n++ )
		{
			// (signal(x)-signal_min)/(signal_max/signal_min)
			double dRate = ( microenvironment(n)[proinflammatory_cytokine_index] - TC_min_signal ); 
			dRate /= TC_max_minus_min; 
			// crop to [0,1] 
			if( dRate > 1 ) 
			{ dRate = 1; } 
			if( dRate < 0 )
			{ dRate = 0; }
			total_rate += dRate; 
		}	
		// multiply by dV and rate_max 
		total_scaled_signal = total_rate; 
		
		total_rate *= microenvironment.mesh.dV; 
		total_rate *= CD8_Tcell_recruitment_rate; 
		
		// expected number of new neutrophils 
		number_of_new_cells = (int) round( total_rate * elapsed_time ); 
		if( number_of_new_cells )
		{
			std::cout << "\tRecruiting " << number_of_new_cells << " CD8 T cells ... " << std::endl; 
			
//			std::cout << "\tTotal signal/dV : " << total_scaled_signal << std::endl;
//			std::cout << "\tTotal signa : " << total_scaled_signal * microenvironment.mesh.dV << std::endl; 
//			double total_volume = microenvironment.mesh.dV * microenvironment.mesh.voxels.size() ; 
//			std::cout << "\tmean signal : " << total_scaled_signal * microenvironment.mesh.dV / total_volume << std::endl; 
			for( int n = 0; n < number_of_new_cells ; n++ )
			{ create_infiltrating_Tcell(); }
		}
		
		t_last_immune = t_immune; 
		t_next_immune = t_immune + dt_immune; 
		
		std::cout << "\t\tnext immune time: " << t_next_immune << std::endl;  
	}
	t_immune += dt; 
	return; 
}

void initial_immune_cell_placement( void )
{
	Cell_Definition* pCD8 = find_cell_definition( "CD8 Tcell" ); 
	Cell_Definition* pMF = find_cell_definition( "macrophage" ); 
	Cell_Definition* pN = find_cell_definition( "neutrophil" ); 

	// CD8+ T cells; 
	for( int n = 0 ; n < parameters.ints("number_of_CD8_Tcells") ; n++ )
	{ create_infiltrating_immune_cell( pCD8 ); }		

	// macrophages 
	for( int n = 0 ; n < parameters.ints("number_of_macrophages") ; n++ )
	{ create_infiltrating_immune_cell( pMF ); }		

	// neutrophils 
	for( int n = 0 ; n < parameters.ints("number_of_neutrophils") ; n++ )
	{ create_infiltrating_immune_cell( pN ); }		

	return;
}