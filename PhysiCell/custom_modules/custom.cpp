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
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
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

#include "./custom.h"

// declare cell definitions here 

Cell_Definition lung_epithelium; 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "default"; 
	
	// set default cell cycle model 

	cell_defaults.functions.cycle_model = live; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = viral_dynamics; // NULL; 
	cell_defaults.functions.custom_cell_rule = NULL; // receptor_dynamics_model; 
	
	// needed for a 2-D simulation: 
	
	/* grab code from heterogeneity */ 
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );

	int live_phase_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::live );

	// initially no necrosis or apoptosis 
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0; 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 
	
	// set cycle rate to zero 
	cell_defaults.phenotype.cycle.data.transition_rate(live_phase_index,live_phase_index) = 0.0; 

	// set all secretion, uptake, and export rates to zero 
	// set all to be fully released by apoptotic cells 
	for( int n = 0; n < microenvironment.number_of_densities(); n++ )
	{
		cell_defaults.phenotype.secretion.uptake_rates[n] = 0; 
		cell_defaults.phenotype.secretion.secretion_rates[n] = 0; 
		cell_defaults.phenotype.secretion.saturation_densities[n] = 0; 
		cell_defaults.phenotype.secretion.net_export_rates[n] = 0; 
		
		cell_defaults.phenotype.molecular.fraction_released_at_death[n] = 
			parameters.doubles("virus_fraction_released_at_death"); 
	}
	
	// register the submodels 
	// (which ensures that the cells have all the internal variables they need) 
	
	receptor_dynamics_model_setup(); 
	internal_virus_model_setup();
	internal_virus_response_model_setup();
	submodel_registry.display( std::cout ); 
	
	// set uptake rate of virions 
	static int nV = microenvironment.find_density_index( "virion" ); 
	// cell_defaults.phenotype.secretion.uptake_rates[nV] = parameters.doubles( "virion_uptake_rate" ); 
	
	// disable motility 
	cell_defaults.phenotype.motility.is_motile = false; 	
	
	// add custom data here
	// set variable values here 

	// viral dynamics parameters 
	Parameter<double> paramD; 	
	paramD = parameters.doubles["virion_uncoating_rate"]; 
//	cell_defaults.custom_data.add_variable( "virion_uncoating_rate" , paramD.units, paramD.value ); 
	cell_defaults.custom_data[ "virion uncoating rate" ] = paramD.value; 
	
	paramD = parameters.doubles["uncoated_to_RNA_rate"]; 
	cell_defaults.custom_data[ "uncoated to RNA rate" ] = paramD.value; 
	
	paramD = parameters.doubles["protein_synthesis_rate"]; 
	cell_defaults.custom_data[ "protein synthesis rate" ] = paramD.value; 
	
	paramD = parameters.doubles["virion_assembly_rate"]; 
	cell_defaults.custom_data[ "virion assembly rate" ] = paramD.value; 
	
	paramD = parameters.doubles["virion_export_rate"]; 
	cell_defaults.custom_data[ "virion export rate" ] = paramD.value; 

	// viral response parameters. 
	paramD = parameters.doubles["max_infected_apoptosis_rate"]; 
	cell_defaults.custom_data[ "max infected apoptosis rate" ] = paramD.value; 

	paramD = parameters.doubles["max_apoptosis_half_max"]; 
	cell_defaults.custom_data[ "max apoptosis half max" ] = paramD.value; 

	paramD = parameters.doubles["apoptosis_hill_power"]; 
	cell_defaults.custom_data[ "apoptosis hill power" ] = paramD.value; 
	
	// receptor dynamics parameters. 

	paramD = parameters.doubles["ACE2_binding_rate"]; 
	cell_defaults.custom_data[ "ACE2 binding rate" ] = paramD.value; 

	paramD = parameters.doubles["ACE2_endocytosis_rate"]; 
	cell_defaults.custom_data[ "ACE2 endocytosis rate" ] = paramD.value; 

	paramD = parameters.doubles["ACE2_cargo_release_rate"]; 
	cell_defaults.custom_data[ "ACE2 cargo release rate" ] = paramD.value; 

	paramD = parameters.doubles["ACE2_recycling_rate"]; 
	cell_defaults.custom_data[ "ACE2 recycling rate" ] = paramD.value; 

	paramD = parameters.doubles["ACE2_receptors_per_cell"]; 
	cell_defaults.custom_data[ "unbound external ACE2" ] = paramD.value; 
	
	
	// Now, let's define another cell type. 
	// It's best to just copy the default and modify it. 
	
	// make this cell type randomly motile, less adhesive, greater survival, 
	// and less proliferative 
	
	lung_epithelium = cell_defaults; 
	lung_epithelium.type = 1; 
	lung_epithelium.name = "lung epithelium"; 
	
	// make sure the new cell type has its own reference phenotype
	
	lung_epithelium.parameters.pReference_live_phenotype = &( lung_epithelium.phenotype ); 
	
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
/* now this is in XML 
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
*/
	
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
/* now this is in XML 	
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = false; 
	
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0 }; 
*/
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	static int nV = microenvironment.find_density_index( "virion" ); 
	
	
	// create some cells near the origin
	
	Cell* pC;
	
	// hexagonal cell packing 
	
	double cell_radius = lung_epithelium.phenotype.geometry.radius; 
	double spacing = 0.95 * cell_radius * 2.0; 
	
	double x_min = microenvironment.mesh.bounding_box[0] + cell_radius; 
	double x_max = microenvironment.mesh.bounding_box[3] - cell_radius; 

	double y_min = microenvironment.mesh.bounding_box[1] + cell_radius; 
	double y_max = microenvironment.mesh.bounding_box[4] - cell_radius; 
	
	double x = x_min; 
	double y = y_min; 
	
	double center_x = 0.5*( x_min + x_max ); 
	double center_y = 0.5*( y_min + y_max ); 
	
	double triangle_stagger = sqrt(3.0) * spacing * 0.5; 
	
	// find hte cell nearest to the center 
	double nearest_distance_squared = 9e99; 
	Cell* pNearestCell = NULL; 
	
	int n = 0; 
	while( y < y_max )
	{
		while( x < x_max )
		{
			pC = create_cell( lung_epithelium ); 
			pC->assign_position( x,y, 0.0 );
			
			double dx = x - center_x;
			double dy = y - center_y; 
			
			double temp = dx*dx + dy*dy; 
			if( temp < nearest_distance_squared )
			{
				nearest_distance_squared = temp;
				pNearestCell = pC; 
			}
			
/*			
			// if this cell is at (0,0,0), insert one virion
			
			if( fabs( x-5 ) < 5 && fabs( y-5 ) < 5 )
			{
				pC->phenotype.molecular.internalized_total_substrates[ nV ] = 1.0; 
			}
*/			
			
			x += spacing; 
		}
		x = x_min; 
		
		n++; 
		y += triangle_stagger; 
		// in odd rows, shift 
		if( n % 2 == 1 )
		{
			x += 0.5 * spacing; 
		}
	}
	
	int number_of_virions = (int) ( parameters.doubles("multiplicity_of_infection") * 
		(*all_cells).size() ); 
	double single_virion_density_change = 1.0 / microenvironment.mesh.dV; 
	
	// infect the cell closest to the center  

	if( parameters.bools( "use_single_infected_cell" ) == true )
	{
		std::cout << "Infecting center cell with one virion ... " << std::endl; 
		pNearestCell->phenotype.molecular.internalized_total_substrates[ nV ] = 1.0; 
	}
	else
	{
		std::cout << "Placing " << number_of_virions << " virions ... " << std::endl; 
		for( int n=0 ; n < number_of_virions ; n++ )
		{
			// pick a random voxel 
			std::vector<double> position = {0,0,0}; 
			position[0] = x_min + (x_max-x_min)*UniformRandom(); 
			position[1] = y_min + (y_max-y_min)*UniformRandom(); 
			
			int m = microenvironment.nearest_voxel_index( position ); 
			
			// int n = (int) ( ( microenvironment.number_of_voxels()-1.0 ) * UniformRandom() ); 
			// microenvironment(i,j)[nV] += single_virion_density_change; 
			microenvironment(m)[nV] += single_virion_density_change; 
		}
	}
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector< std::string> output( 4, "black" ); 

	// static int color_index = cell_defaults.custom_data.find_variable_index( "assembled virion" ); 
	static int color_index = cell_defaults.custom_data.find_variable_index( parameters.strings["color_variable"].value ); 
	static int nV = cell_defaults.custom_data.find_variable_index( "virion" ); 
	
	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int nR_EB = cell_defaults.custom_data.find_variable_index( "bound external ACE2" ); 
	static int nR_IB = cell_defaults.custom_data.find_variable_index( "bound internal ACE2" ); 
	
	
	// color by assembled virion 
	
//	static double my_max = -9e9; 

	if( pCell->phenotype.death.dead == false )
	{
		// find fraction of max viral load 
		double v = pCell->custom_data[ color_index ] ; 
		
		double interpolation = 0; 
		if( v < 1 )
		{ interpolation = 0; } 
		if( v >= 1.0 && v < 10 )
		{ interpolation = 0.25; } 
		if( v >= 10.0 && v < 100 )
		{ interpolation = 0.5; } 
		if( v >= 100.0 && v < 1000 )
		{ interpolation = 0.75; } 
		if( v >= 1000.0 )
		{ interpolation = 1.0; } 

		int red = (int) floor( 255.0 * interpolation ) ; 
		int green = red; 
		int blue = 255 - red; 

		char color [1024]; 
		sprintf( color, "rgb(%u,%u,%u)" , red,green,blue ); 

		output[0] = color; 
		output[2] = color; 
		output[3] = color; 
	}
/*
		// color boundary by bound ACE2 receptor on 
		// surface or inside cell 
		
		v = pCell->custom_data[ nR_EB ] + pCell->custom_data[ nR_IB ] + 
			pCell->phenotype.molecular.internalized_total_substrates[ nV_external ]; 

		interpolation = 0; 
		if( v < 1 )
		{ interpolation = 0; } 
		if( v >= 1.0 && v < 10 )
		{ interpolation = 0.25; } 
		if( v >= 10.0 && v < 100 )
		{ interpolation = 0.5; } 
		if( v >= 100.0 && v < 1000 )
		{ interpolation = 0.75; } 
		if( v >= 1000.0 )
		{ interpolation = 1.0; } 

		red = (int) floor( 255.0 * interpolation ) ; 
		green = red; 
		blue = 255 - red; 

		sprintf( color, "rgb(%u,%u,%u)" , red,green,blue ); 
		
		output[1] = color;			
	}
*/
	
	return output; 
}


void viral_dynamics( Cell* pCell, Phenotype& phenotype, double dt )
{
	// viral dynamics model 
	internal_virus_model(pCell,phenotype,dt);
	
	// viral response model 
	
	internal_virus_response_model(pCell,phenotype,dt);

	return; 
}


void move_exported_to_viral_field( void )
{
	static int nV = microenvironment.find_density_index( "virion" ); 
	static int nA = microenvironment.find_density_index( "assembled virion" ); 
	
	#pragma omp parallel for 
	for( int n = 0 ; n < microenvironment.number_of_voxels() ; n++ )
	{
		microenvironment(n)[nV] += microenvironment(n)[nA]; 
		microenvironment(n)[nA] = 0; 
	}
	
	return;
}

