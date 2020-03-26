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
	
	// set uptake rate of virions 
	static int nV = microenvironment.find_density_index( "virion" ); 
	cell_defaults.phenotype.secretion.uptake_rates[ nV] = parameters.doubles( "virion_uptake_rate" ); 
	
	// disable motility 
	cell_defaults.phenotype.motility.is_motile = false; 	
	
	// add custom data here

	Parameter<double> paramD; 	
	paramD = parameters.doubles["virion_uncoating_rate"]; 
	cell_defaults.custom_data.add_variable( "virion_uncoating_rate" , paramD.units, paramD.value ); 
	
	paramD = parameters.doubles["uncoated_to_RNA_rate"]; 
	cell_defaults.custom_data.add_variable( "uncoated_to_RNA_rate" , paramD.units, paramD.value ); 
	
	paramD = parameters.doubles["protein_synthesis_rate"]; 
	cell_defaults.custom_data.add_variable( "protein_synthesis_rate" , paramD.units, paramD.value ); 
	
	paramD = parameters.doubles["virion_assembly_rate"]; 
	cell_defaults.custom_data.add_variable( "virion_assembly_rate" , paramD.units, paramD.value ); 
	
	paramD = parameters.doubles["virion_export_rate"]; 
	cell_defaults.custom_data.add_variable( "virion_export_rate" , paramD.units, paramD.value ); 


	paramD = parameters.doubles["max_infected_apoptosis_rate"]; 
	cell_defaults.custom_data.add_variable( "max_infected_apoptosis_rate" , paramD.units, paramD.value ); 

	paramD = parameters.doubles["max_apoptosis_half_max"]; 
	cell_defaults.custom_data.add_variable( "max_apoptosis_half_max" , paramD.units, paramD.value ); 

	paramD = parameters.doubles["apoptosis_hill_power"]; 
	cell_defaults.custom_data.add_variable( "apoptosis_hill_power" , paramD.units, paramD.value ); 
	

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
	
	double triangle_stagger = sqrt(3.0) * spacing * 0.5; 

	int n = 0; 
	while( y < y_max )
	{
		while( x < x_max )
		{
			pC = create_cell( lung_epithelium ); 
			pC->assign_position( x,y, 0.0 );
			
			
			// if this cell is at (0,0,0), insert one virion
			
			if( fabs( x-5 ) < 5 && fabs( y-5 ) < 5 )
			{
				pC->phenotype.molecular.internalized_total_substrates[ nV ] = 1.0; 
			}
			
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
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector< std::string> output( 4, "black" ); 
	// std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
	
	static int color_index = microenvironment.find_density_index( "assembled virion" ); 
	
	// color by assembled virion 
	
	static double max_virion = 1.0 * pCell->custom_data[ "max_apoptosis_half_max" ] + 1e-16; 

	if( pCell->phenotype.death.dead == false )
	{
		// find fraction of max viral load 
		double v = pCell->phenotype.molecular.internalized_total_substrates[ color_index ] ; 
		
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
	
	return output; 
}


void viral_dynamics( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int nE = microenvironment.find_density_index( "virion" ); 
	static int nUV = microenvironment.find_density_index( "uncoated virion" ); 
	static int nR = microenvironment.find_density_index( "viral RNA" ); 
	static int nP = microenvironment.find_density_index( "viral protein" ); 
	static int nA = microenvironment.find_density_index( "assembled virion" ); 
	
	// if dead, stop this business, 
	if( phenotype.death.dead == true )
	{
		phenotype.secretion.set_all_uptake_to_zero();
		pCell->functions.update_phenotype = NULL; 
		return;
	}
	
	// uncoat endocytosed virus
	double dE = dt * pCell->custom_data["virion_uncoating_rate"] *  ( phenotype.molecular.internalized_total_substrates[nE] ); 
	if( dE > phenotype.molecular.internalized_total_substrates[nE] )
	{ dE = phenotype.molecular.internalized_total_substrates[nE]; } 
	phenotype.molecular.internalized_total_substrates[nE] -= dE; 
	phenotype.molecular.internalized_total_substrates[nUV] += dE; 
	
	// convert uncoated virus to usable mRNA 
	double dR = dt * pCell->custom_data["uncoated_to_RNA_rate"] *  ( phenotype.molecular.internalized_total_substrates[nUV] ); 
	if( dR > phenotype.molecular.internalized_total_substrates[nUV] )
	{ dR = phenotype.molecular.internalized_total_substrates[nUV]; } 
	phenotype.molecular.internalized_total_substrates[nUV] -= dR; 
	phenotype.molecular.internalized_total_substrates[nR] += dR; 
	
	// use mRNA to create viral protein 
	double dP = dt * pCell->custom_data["protein_synthesis_rate"] *  ( phenotype.molecular.internalized_total_substrates[nR] ); 
	phenotype.molecular.internalized_total_substrates[nP] += dP; 
	
	// degrade mRNA 


	// degrade protein 
	
	
	// assemble virus 
	double dA = dt * pCell->custom_data["virion_assembly_rate"] *  ( phenotype.molecular.internalized_total_substrates[nP] ); 
	if( dA > phenotype.molecular.internalized_total_substrates[nP] )
	{ dA = phenotype.molecular.internalized_total_substrates[nP]; } 
	phenotype.molecular.internalized_total_substrates[nP] -= dA; 
	phenotype.molecular.internalized_total_substrates[nA] += dA; 
	
	
	// set export rate 
	
	phenotype.secretion.net_export_rates[nA] = pCell->custom_data["virion_export_rate" ] *  ( phenotype.molecular.internalized_total_substrates[nA] ); 
	

	// now, set apoptosis rate 
	
	static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	phenotype.death.rates[apoptosis_model_index] = pCell->custom_data["max_infected_apoptosis_rate"] ; 
	
	double v = phenotype.molecular.internalized_total_substrates[nA] /
		pCell->custom_data["max_apoptosis_half_max"] ; 
	v = pow( v, pCell->custom_data["apoptosis_hill_power"] ); 
	
	double effect = v / (1.0+v); 
	phenotype.death.rates[apoptosis_model_index] *= effect; 
	
/*
	if( fabs( pCell->position[0]-5 ) < 5 && fabs( pCell->position[1]-5 ) < 5 )
	{
		std::cout << phenotype.molecular.internalized_total_substrates << std::endl; 
		std::cout << phenotype.secretion.net_export_rates << std::endl; 
		std::cout << v << " : " << effect << " :: " << phenotype.death.rates[apoptosis_model_index] << std::endl; 
	}
*/	
	
	/*
	
		cell_defaults.custom_data.add_variable( "assembled_virion" , "virion" , 0.0 ); 

	
	
	// add custom data here

 

	paramD = parameters.doubles["max_infected_apoptosis_rate"]; 
	cell_defaults.custom_data.add_variable( "max_infected_apoptosis_rate" , paramD.units, paramD.value ); 

	paramD = parameters.doubles["max_apoptosis_half_max"]; 
	cell_defaults.custom_data.add_variable( "max_apoptosis_half_max" , paramD.units, paramD.value ); 

	paramD = parameters.doubles["apoptosis_hill_power"]; 
	cell_defaults.custom_data.add_variable( "apoptosis_hill_power" , paramD.units, paramD.value ); 
	
*/

	
	
	
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

