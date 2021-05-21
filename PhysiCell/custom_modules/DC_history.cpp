#include "./DC_history.h" 
#include <algorithm> 

using namespace PhysiCell; 

std::string DC_history_version = "0.5.0"; 

Submodel_Information DC_history_info; 

void DC_history_model_setup( void )
{
		// set version 
	DC_history_info.name = "DC history"; 
	DC_history_info.version = DC_history_version; 
		// set functions 
	DC_history_info.main_function = DC_history_main_model; 
	DC_history_info.phenotype_function = NULL; // pushed into the "main" model  
	DC_history_info.mechanics_function = NULL; 	
		// what microenvironment variables do you need 
		
		// what cell variables and parameters do you need? 
	
	// submodel_registry.register_model( receptor_dynamics_info ); 
	DC_history_info.register_model(); 
	
	return; 
}

void DC_history_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int DC_type = get_cell_definition( "DC" ).type; 
	
	// bookkeeping -- find microenvironment variables we need

	// bookkeeping -- find custom data we need 
	static double DCprob = parameters.doubles( "DC_leave_prob" ); 
	extern double DCAMOUNT; //declare existance of counter
	// do nothing if dead 
	if( phenotype.death.dead == true )
	{ return; } 

	// if not DC, do nothing 
	if( pCell->type != DC_type )
	{ return; } 
	
	// (Adrianne) if DC is already activated, then check whether it leaves the tissue
	if( pCell->custom_data["activated_immune_cell"] >  0.5 && UniformRandom() < DCprob)
	{
		// (Adrianne) DC leaves the tissue and so we lyse that DC
		std::cout<<"DC leaves tissue"<<std::endl;
		pCell->lyse_cell(); 
		#pragma omp critical 
		{ DCAMOUNT++; } // add one	
		return;
		
	}
	
	return; 
}

void DC_history_main_model( double dt )
{
	extern double DCAMOUNT;
	extern std::vector<int>history;
	DCAMOUNT=0;
	
	#pragma omp parallel for 
	for( int n=0; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if( pC->phenotype.death.dead == false )
		{ DC_history_model( pC, pC->phenotype , dt ); }
	}
	std::rotate(history.rbegin(),history.rbegin()+1,history.rend());
	history.front() = DCAMOUNT;
	
	/* std::copy(history.begin(), history.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl; 
	 */
	return; 
}
