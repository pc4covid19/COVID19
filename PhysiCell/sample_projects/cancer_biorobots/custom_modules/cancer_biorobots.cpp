/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./cancer_biorobots.h"

Cell_Definition* cargo_cell; 
Cell_Definition* worker_cell; 

void create_cargo_cell_type( void ) 
{
	cargo_cell = find_cell_definition( "cargo cell" ); 

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	
	int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0 
	int attract_ID = microenvironment.find_density_index( "chemoattractant" ); // 1
	int therapy_ID = microenvironment.find_density_index( "therapeutic" ); // 2
	
	// reduce o2 uptake 
	
	cargo_cell->phenotype.secretion.uptake_rates[oxygen_ID] *= 
		parameters.doubles("cargo_o2_relative_uptake");   
	
	cargo_cell->phenotype.mechanics.cell_cell_adhesion_strength *= 
		parameters.doubles("cargo_relative_adhesion"); // 0.0;
	cargo_cell->phenotype.mechanics.cell_cell_repulsion_strength *= 
		parameters.doubles("cargo_relative_repulsion"); // 5.0;
		
	// figure out mechanics parameters 
	
	cargo_cell->phenotype.mechanics.relative_maximum_attachment_distance 
		= parameters.doubles( "max_attachment_distance" ) / cargo_cell->phenotype.geometry.radius ; 

	cargo_cell->phenotype.mechanics.relative_detachment_distance 
		= parameters.doubles( "max_elastic_displacement" ) / cargo_cell->phenotype.geometry.radius ; 
		
	cargo_cell->phenotype.mechanics.attachment_elastic_constant 
		= parameters.doubles("elastic_coefficient"); 
	
	// set functions 
	
	cargo_cell->functions.update_phenotype = cargo_cell_phenotype_rule; 
	cargo_cell->functions.custom_cell_rule = cargo_cell_rule; 
	cargo_cell->functions.contact_function = biorobots_contact_function; 
	cargo_cell->functions.update_migration_bias = NULL;	
	
	// set custom data values 
	
	return;
}	

void create_worker_cell_type( void )
{
	worker_cell = find_cell_definition( "worker cell" ); 

	int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0 
	int attract_ID = microenvironment.find_density_index( "chemoattractant" ); // 1
	int therapy_ID = microenvironment.find_density_index( "therapeutic" ); // 2	
	
	// reduce o2 uptake 
	
	worker_cell->phenotype.secretion.uptake_rates[oxygen_ID] *= 
		parameters.doubles("worker_o2_relative_uptake"); // 0.1; 
	
	worker_cell->phenotype.mechanics.cell_cell_adhesion_strength *= 
		parameters.doubles("worker_relative_adhesion"); // 0.0;
	worker_cell->phenotype.mechanics.cell_cell_repulsion_strength *= 
		parameters.doubles("worker_relative_repulsion"); // 5.0;
		
	// figure out mechanics parameters 
	
	worker_cell->phenotype.mechanics.relative_maximum_attachment_distance 
		= parameters.doubles( "max_attachment_distance" ) / worker_cell->phenotype.geometry.radius ; 
		
	worker_cell->phenotype.mechanics.attachment_elastic_constant 
		= parameters.doubles("elastic_coefficient"); 		
	
	worker_cell->phenotype.mechanics.relative_detachment_distance 
		= parameters.doubles( "max_elastic_displacement" ) / worker_cell->phenotype.geometry.radius ; 

	// set functions 
	
	worker_cell->functions.update_phenotype = NULL; // worker_cell_rule; 
	worker_cell->functions.custom_cell_rule = worker_cell_rule;  
	worker_cell->functions.contact_function = biorobots_contact_function; 
	worker_cell->functions.update_migration_bias = worker_cell_motility;	
	
	// set custom data values 
	
	return; 
}

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	SeedRandom( parameters.ints("random_seed") ); 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	
	cell_defaults.parameters.o2_proliferation_saturation = 38.0;  
	cell_defaults.parameters.o2_reference = 38.0; 
	
	int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0 
	int attract_ID = microenvironment.find_density_index( "chemoattractant" ); // 1
	int therapy_ID = microenvironment.find_density_index( "therapeutic" ); // 2	
	
	// set the default cell type to o2-based proliferation with the effect of the 
	// on oncoprotein, and secretion of the immunostimulatory factor 
	
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_therapy; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	// change the max cell-cell adhesion distance 
	cell_defaults.phenotype.mechanics.set_relative_maximum_adhesion_distance(parameters.doubles("max_relative_cell_adhesion_distance") );
	
	cell_defaults.phenotype.mechanics.attachment_elastic_constant = parameters.doubles( "elastic_coefficient" );

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	// create the biorobot types 
	create_cargo_cell_type(); 
	create_worker_cell_type(); 
	
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters
	
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "WARNING: overriding from 3-D to 2-D" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
			
	initialize_microenvironment(); 	

	return; 
}	

//keep 
void introduce_biorobots( void )
{
	// idea: we'll "inject" them in a little column
		
	static double worker_fraction = 
		parameters.doubles("worker_fraction"); // 0.10; /* param */ 
	static int number_of_injected_cells = 
		parameters.ints("number_of_injected_cells"); // 500; /* param */ 
	
	// make these vary with domain size 
	double left_coordinate = default_microenvironment_options.X_range[1] - 150.0; // 600.0; 
	double right_cooridnate = default_microenvironment_options.X_range[1] - 50.0; // 700.0;

	double bottom_coordinate = default_microenvironment_options.Y_range[0] + 50.0; // -700; 
	double top_coordinate = default_microenvironment_options.Y_range[1] - 50.0; // 700; 
		
	for( int i=0 ;i < number_of_injected_cells ; i++ )
	{
		std::vector<double> position = {0,0,0}; 
		position[0] = left_coordinate + (right_cooridnate-left_coordinate)*UniformRandom(); 
		position[1] = bottom_coordinate + (top_coordinate-bottom_coordinate)*UniformRandom(); 
		
		Cell* pCell;  
		if( UniformRandom() <= worker_fraction )
		{ pCell = create_cell( *worker_cell ); }
		else
		{ pCell = create_cell( *cargo_cell ); }
		pCell->assign_position( position ); 
	}
	
	return; 
}

void setup_tissue( void )
{
	// place a cluster of tumor cells at the center 
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = parameters.doubles("tumor_radius"); // 200.0; 
	
	Cell* pCell = NULL; 
	
	double x = 0.0; 
	double x_outer = tumor_radius; 
	double y = 0.0; 
	
	int n = 0; 
	while( y < tumor_radius )
	{
		x = 0.0; 
		if( n % 2 == 1 )
		{ x = 0.5*cell_spacing; }
		x_outer = sqrt( tumor_radius*tumor_radius - y*y ); 
		
		while( x < x_outer )
		{
			pCell = create_cell(); // tumor cell 
			pCell->assign_position( x , y , 0.0 );
			
			if( fabs( y ) > 0.01 )
			{
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( x , -y , 0.0 );				
			}
			
			if( fabs( x ) > 0.01 )
			{ 
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( -x , y , 0.0 );
				
				if( fabs( y ) > 0.01 )
				{
					pCell = create_cell(); // tumor cell 
					pCell->assign_position( -x , -y , 0.0 );
				}
			}
			x += cell_spacing; 
			
		}
		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 		
	
	return; 
}

// done 
std::vector<std::string> cancer_biorobots_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4, "black" ); 
	
	static int damage_i = pCell->custom_data.find_variable_index( "damage" ); 
	static double max_damage = 1.0 * cell_defaults.custom_data["damage_rate"] / (1e-16 + cell_defaults.custom_data[ "repair_rate" ] );
	
	// cargo cell 
	if( pCell->type == cargo_cell->type )
	{
		output[0] = "blue";
		output[1] = "blue";
		output[2] = "blue"; 
		output[3] = "none"; // no nuclear outline color 
		return output;
	}
	
	// worker cell 
	if( pCell->type == worker_cell->type )
	{
		output[0] = "red";
		output[1] = "red";
		output[2] = "red"; 
		output[3] = "none"; // no nuclear outline color 
		return output;
	}
	
	// apoptotic tumor - cyan 
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - cyan
	{
		output[0] = "cyan";
		output[2] = "darkcyan"; 
		return output; 
	}	
	
	// Necrotic tumor - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
		return output; 
	}		
	
	// live tumor -- shade by level of damage 
	
	
	// if live: color by damage 
	if( pCell->phenotype.death.dead == false )
	{
		int damage = (int) round( pCell->custom_data[damage_i] * 255.0 / max_damage ); 
		
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)" , damage , 255-damage , damage );
		output[0].assign( szTempString );
		output[1].assign( szTempString );
		sprintf( szTempString , "rgb(%u,%u,%u)" , damage/4 , (255-damage)/4 , damage/4 );
		output[2].assign( szTempString );
	}
	return output; 
}

// keep! 
Cell* worker_cell_check_neighbors_for_attachment( Cell* pWorker , double dt )
{
	std::vector<Cell*> nearby = pWorker->cells_in_my_container(); 
	int i = 0; 
	while( i < nearby.size() )
	{
		// don't try to attach to yourself 
		if( nearby[i] != pWorker )
		{
			if( worker_cell_attempt_attachment( pWorker, nearby[i] , dt ) )
			{ return nearby[i]; }
		}
		i++; 
	}
	
	return NULL; 
}

// keep! 
bool worker_cell_attempt_attachment( Cell* pWorker, Cell* pCargo , double dt )
{
	static int receptor_i = pCargo->custom_data.find_variable_index( "receptor" ); 

	static double receptor_threshold = 
		parameters.doubles("attachment_receptor_threshold"); // 0.1; 
	
	static double max_attachment_distance = 
		parameters.doubles("max_attachment_distance"); // 18.0; 
	static double min_attachment_distance = 
		parameters.doubles("min_attachment_distance"); // 14.0; 
	static double attachment_difference = max_attachment_distance - min_attachment_distance; 
	
	if( pCargo->custom_data[receptor_i] > receptor_threshold )
	{
		std::vector<double> displacement = pCargo->position - pWorker->position;
		double distance = norm( displacement ); 
		if( distance > max_attachment_distance )
		{ return false; } 
	
		if( distance < min_attachment_distance )
		{ 
			attach_cells( pWorker, pCargo );
			return true; 
		}
		
		return true; 
	}
	
	return false; 
}

void worker_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if I am dead, don't bother

	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		
		// Let's just fully disable now. 
		pCell->functions.custom_cell_rule = NULL; 
		return; 
	}
	
	// am I searching for cargo? if so, see if I've found it
	if( pCell->state.number_of_attached_cells() == 0 )
	{
		std::vector<Cell*> nearby = pCell->cells_in_my_container(); 
		bool attached = false; // want to limit to one attachment 
		int i =0;
		while( i < nearby.size() && attached == false )
		{
			// if it is expressing the receptor, dock with it 
			if( nearby[i]->custom_data["receptor"] > 0.5 && attached == false )
			{
				attach_cells( pCell, nearby[i] ); 
				// nearby[i]->custom_data["receptor"] = 0.0; // put into cargo cell rule instead? 
				// nearby[i]->phenotype.secretion.set_all_secretion_to_zero(); // put into cargo rule instead? 
				attached = true; 
			}
			i++; 
		}
		
	}
	
	return; 
}

void worker_cell_motility( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int o2_index = microenvironment.find_density_index( "oxygen" ); 
	static int signal_index = microenvironment.find_density_index( "chemoattractant" ); 

	static double detection_threshold = 
		parameters.doubles("motility_shutdown_detection_threshold"); // 0.001; 
	
	// if attached, biased motility towards director chemoattractant 
	// otherwise, biased motility towards cargo chemoattractant 
	
	static double attached_worker_migration_bias = 
		parameters.doubles("attached_worker_migration_bias"); 
	static double unattached_worker_migration_bias = 
		parameters.doubles("unattached_worker_migration_bias"); 
	
	if( pCell->state.number_of_attached_cells() > 0 )
	{
		phenotype.motility.migration_bias = attached_worker_migration_bias; 

		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(o2_index);	
		phenotype.motility.migration_bias_direction *= -1.0; 
		normalize( &( phenotype.motility.migration_bias_direction ) );			
	}
	else
	{
		// if there is no detectable signal, shut down motility (permanently)
		if( pCell->nearest_density_vector()[signal_index] < detection_threshold )
		{
			phenotype.motility.is_motile = false; 
			pCell->functions.update_migration_bias = NULL; 
		}
		
		phenotype.motility.migration_bias = unattached_worker_migration_bias; 
		
		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(signal_index);	
		normalize( &( phenotype.motility.migration_bias_direction ) );			
	}
	
	return; 
}

/*  CARGO CELL RULES */ 

void cargo_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		
		// Let's just fully disable now. 
		pCell->functions.custom_cell_rule = NULL; 
		return; 
	}
	
	// if I'm docked
	if( pCell->state.number_of_attached_cells() > 0 )
	{
		phenotype.motility.is_motile = false; 
		return; 
	}
	
	return; 
}

// phenotype rule 
void cargo_cell_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int o2_index = microenvironment.find_density_index( "oxygen" ); 
	static int signal_index = microenvironment.find_density_index( "chemoattractant" ); 
	static int drug_index = microenvironment.find_density_index( "therapeutic" ); 
	
	static int drop_index = pCell->custom_data.find_variable_index( "cargo_release_o2_threshold" ); 
	static int receptor_index = pCell->custom_data.find_variable_index( "receptor" ); 
	
	static int apoptosis_model_index = phenotype.death.find_death_model_index( "apoptosis" );
	
	// if dettached and receptor on, secrete signal 
	
	// if dettached and receptor off, secrete chemo
	
	if( pCell->state.number_of_attached_cells() == 0 )
	{
		if( pCell->custom_data[receptor_index] > 0.1 )
		{
			phenotype.secretion.secretion_rates[signal_index] = 10.0; 
			phenotype.secretion.secretion_rates[drug_index] = 0.0; 
		}
		else
		{
			phenotype.secretion.secretion_rates[signal_index] = 0.0; 
			phenotype.secretion.secretion_rates[drug_index] = 10.0; 
		}
		return; 
	}
	
	// if you reach this point of the code, the cell is attached 
	

	// if attached and oxygen high, secrete nothing, receptor off 
	
	// if attached and oxygen low, dettach, start secreting chemo, receptor off   
	
	if( pCell->nearest_density_vector()[o2_index] > pCell->custom_data[drop_index] )
	{
		phenotype.secretion.secretion_rates[signal_index] = 0.0; 
		phenotype.secretion.secretion_rates[drug_index] = 0.0; 
		pCell->custom_data[receptor_index] = 0.0; 
	}
	else
	{
		phenotype.secretion.secretion_rates[signal_index] = 0.0; 
		phenotype.secretion.secretion_rates[drug_index] = 10.0; 
		pCell->custom_data[receptor_index] = 0.0; 		
		
		pCell->remove_all_attached_cells(); 
	}
	
	return; 
}

/* TUMOR CELL RULES */ 

void tumor_cell_phenotype_with_therapy( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int damage_i = pCell->custom_data.find_variable_index( "damage" ); 

	static int damage_rate_i = pCell->custom_data.find_variable_index( "damage_rate" ); 
	static int repair_rate_i = pCell->custom_data.find_variable_index( "repair_rate" ); 
	static int death_rate_i = pCell->custom_data.find_variable_index( "drug_death_rate" ); 
	
	static int chemo_i = microenvironment.find_density_index( "therapeutic" ); 
	
	static int apoptosis_model_index = phenotype.death.find_death_model_index( "apoptosis" );	
	
	static double max_damage = 1.0 * cell_defaults.custom_data["damage_rate"] / (1e-16 + cell_defaults.custom_data[ "repair_rate" ] );
	
	// if I'm dead, don't bother. disable my phenotype rule
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 
		return; 
	}
	
	// first, vary the cell birth and death rates with oxygenation
	
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// the update the cell damage 
	
	// dD/dt = alpha*c - beta-D by implicit scheme 
	
	double temp = pCell->nearest_density_vector()[chemo_i];
	
	// reuse temp as much as possible to reduce memory allocations etc. 
	temp *= dt; 
	temp *= pCell->custom_data[damage_rate_i]; 
	
	pCell->custom_data[damage_i] += temp; // d_prev + dt*chemo*damage_rate 
	
	temp = pCell->custom_data[repair_rate_i];
	temp *= dt; 
	temp += 1.0; 
	pCell->custom_data[damage_i] /= temp;  // (d_prev + dt*chemo*damage_rate)/(1 + dt*repair_rate)
	
	// then, see if the cell undergoes death from the therapy 
	
	temp = dt; 
	temp *= pCell->custom_data[damage_i]; 
	temp *= pCell->custom_data[death_rate_i]; 
	temp /= max_damage; // dt*(damage/max_damage)*death_rate 

	if( UniformRandom() <= temp )
	{
		pCell->start_death( apoptosis_model_index );
		pCell->functions.update_phenotype = NULL; 		
		pCell->functions.custom_cell_rule = NULL; 
	}

	return; 
}

void biorobots_contact_function( Cell* pActingOn, Phenotype& pao, Cell* pAttachedTo, Phenotype& pat , double dt )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	
	static double max_elastic_displacement = pao.geometry.radius * pao.mechanics.relative_detachment_distance; 
	static double max_displacement_squared = max_elastic_displacement*max_elastic_displacement; 
	
	// detach cells if too far apart 
	
	if( norm_squared( displacement ) > max_displacement_squared )
	{
		detach_cells( pActingOn , pAttachedTo );
		return; 
	}
	
	axpy( &(pActingOn->velocity) , pao.mechanics.attachment_elastic_constant , displacement ); 
	
	return; 
}

