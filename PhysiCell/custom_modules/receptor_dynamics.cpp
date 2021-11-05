#include "./receptor_dynamics.h" 

using namespace PhysiCell; 

std::string receptor_model_version = "0.5.0"; 

Submodel_Information receptor_dynamics_info; 

void receptor_dynamics_model_setup( void )
{
		// set version 
	receptor_dynamics_info.name = "receptor dynamics"; 
	receptor_dynamics_info.version = receptor_model_version; 
		// set functions 
	receptor_dynamics_info.main_function = receptor_dynamics_main_model; 
	receptor_dynamics_info.phenotype_function = NULL; // pushed into the "main" model  
	receptor_dynamics_info.mechanics_function = NULL; 	
		// what microenvironment variables do you need 
	receptor_dynamics_info.microenvironment_variables.push_back( "virion" ); 		
	receptor_dynamics_info.microenvironment_variables.push_back( "assembled virion" ); 		
		// what cell variables and parameters do you need? 
	receptor_dynamics_info.cell_variables.push_back( "unbound_external_ACE2" ); 
	receptor_dynamics_info.cell_variables.push_back( "bound_external_ACE2" ); 
	receptor_dynamics_info.cell_variables.push_back( "unbound_internal_ACE2" ); 
	receptor_dynamics_info.cell_variables.push_back( "bound_internal_ACE2" ); 
	
	receptor_dynamics_info.cell_variables.push_back( "ACE2_binding_rate" ); 
	receptor_dynamics_info.cell_variables.push_back( "ACE2_endocytosis_rate" ); 
	receptor_dynamics_info.cell_variables.push_back( "ACE2_cargo_release_rate" ); 	
	receptor_dynamics_info.cell_variables.push_back( "ACE2_recycling_rate" ); 
	
	// submodel_registry.register_model( receptor_dynamics_info ); 
	receptor_dynamics_info.register_model(); 
	
	return; 
}

void receptor_dynamics_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int lung_epithelial_type = get_cell_definition( "lung epithelium" ).type; 
	
	// bookkeeping -- find microenvironment variables we need

	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
	
	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
	static int nA_internal = pCell->custom_data.find_variable_index( "assembled_virion" ); 

	// bookkeeping -- find custom data we need 
	
	static int nR_EU = pCell->custom_data.find_variable_index( "unbound_external_ACE2" ); 
	static int nR_EB = pCell->custom_data.find_variable_index( "bound_external_ACE2" ); 
	static int nR_IU = pCell->custom_data.find_variable_index( "unbound_internal_ACE2" ); 
	static int nR_IB = pCell->custom_data.find_variable_index( "bound_internal_ACE2" ); 
	
	static int nR_bind = pCell->custom_data.find_variable_index( "ACE2_binding_rate" ); 
	static int nR_endo = pCell->custom_data.find_variable_index( "ACE2_endocytosis_rate" ); 
	static int nR_release = pCell->custom_data.find_variable_index( "ACE2_cargo_release_rate" ); 	
	static int nR_recycle = pCell->custom_data.find_variable_index( "ACE2_recycling_rate" ); 
	

	// do nothing if dead 
	if( phenotype.death.dead == true )
	{ return; } 

	// if not lung epithelium, do nothing 
	if( pCell->type != lung_epithelial_type )
	{ return; } 
	
	// actual model goes here 
	// reaction set
	
	double x[4][6]={{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};//initialize x
	double f[4][6]={{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};//initialize f
	int j;//initialize counter
	
	//initial values for RK4
	x[0][0] = pCell->custom_data[nR_EU]; 
	x[0][1] = pCell->custom_data[nR_EB];
	x[0][2] = pCell->custom_data[nR_IB]; 
	x[0][3] = pCell->custom_data[nR_IU]; 
	x[0][4] = pCell->custom_data[nV_internal]; 
	x[0][5] = 0; 
	
	
	//SOLVE ODE BEFORE STOCHASTIC PORTION
	for(j = 0; j < 4; j++){
			f[j][0] = {pCell->custom_data[nR_recycle]*x[j][3]}; //define SPECIAL function
			f[j][1] = {-pCell->custom_data[nR_endo]*x[j][1]-0*x[j][1]}; //define SPECIAL function
			f[j][2] = {pCell->custom_data[nR_endo]*x[j][1]-pCell->custom_data[nR_release]*x[j][2]}; //define function
			f[j][3] = {pCell->custom_data[nR_release]*x[j][2]-pCell->custom_data[nR_recycle]*x[j][3]}; //define function
			f[j][4] = {pCell->custom_data[nR_release]*x[j][2]}; //define function
			f[j][5] = {0*x[j][1]}; //counter for export
			if (j== 0 || j==1){
				x[j+1][0]=x[0][0]+dt/2*f[j][0]; //first and second x approximations
				x[j+1][1]=x[0][1]+dt/2*f[j][1]; //first and second x approximations
				x[j+1][2]=x[0][2]+dt/2*f[j][2]; //first and second x approximations
				x[j+1][3]=x[0][3]+dt/2*f[j][3]; //first and second x approximations
				x[j+1][4]=x[0][4]+dt/2*f[j][4]; //first and second x approximations
				x[j+1][5]=x[0][5]+dt/2*f[j][5]; //first and second x approximations
			}
			if (j== 2){
				x[j+1][0]=x[0][0]+dt*f[j][0]; //third approximation
				x[j+1][1]=x[0][1]+dt*f[j][1]; //third approximation
				x[j+1][2]=x[0][2]+dt*f[j][2]; //third approximation
				x[j+1][3]=x[0][3]+dt*f[j][3]; //third approximation
				x[j+1][4]=x[0][4]+dt*f[j][4]; //third approximation
				x[j+1][5]=x[0][5]+dt*f[j][5]; //third approximation
			}
		}

		pCell->custom_data[nR_EU]=(x[0][0]+dt*(f[0][0]/6+f[1][0]/3+f[2][0]/3+f[3][0]/6));
		pCell->custom_data[nR_EB]=x[0][1]+dt*(f[0][1]/6+f[1][1]/3+f[2][1]/3+f[3][1]/6);
		pCell->custom_data[nR_IB]=x[0][2]+dt*(f[0][2]/6+f[1][2]/3+f[2][2]/3+f[3][2]/6); //detirmine n+1
		pCell->custom_data[nR_IU]=x[0][3]+dt*(f[0][3]/6+f[1][3]/3+f[2][3]/3+f[3][3]/6); //detirmine n+1
		pCell->custom_data[nV_internal]=x[0][4]+dt*(f[0][4]/6+f[1][4]/3+f[2][4]/3+f[3][4]/6); //detirmine n+1
	if( immune_cell_check_neighbors_for_attachmentv( pCell , dt) )
	{
		pCell->custom_data[nR_EU] -= 1;
		pCell->custom_data[nR_EB] += 1;
		//cell is removed in called function if passed
	}
	return; 
}

Cell* immune_cell_check_neighbors_for_attachmentv( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container(); //check on this function, can we search nearby containers also?
	int i = 0; 
	while( i < nearby.size() )
	{
		// don't try to kill yourself 
		if( nearby[i] != pAttacker )
		{
			if( attempt_immune_cell_attachmentv( pAttacker, nearby[i] , dt ) )
			{
				return nearby[i]; 
			}
		}
		i++; 
	}
	
	return NULL; 
}

bool attempt_immune_cell_attachmentv( Cell* pAttacker, Cell* pTarget , double dt )
{
	static int nR_EU = pAttacker->custom_data.find_variable_index( "unbound_external_ACE2" ); 
	static int v_type = get_cell_definition( "virion" ).type;
	static int nR_bind = pAttacker->custom_data.find_variable_index( "ACE2_binding_rate" ); 	
	// if the target is not a virion, give up
	if( pTarget->type != v_type )
	{ return false; }
	// if the target cell is dead, give up 
	if( pTarget->phenotype.death.dead == true )
	{ return false; } 

	// if the target cell is too far away, give up 
	std::vector<double> displacement = pTarget->position - pAttacker->position;
	double distance_scale = norm( displacement ); 

	double cell_radius = pAttacker->phenotype.geometry.radius;
	// is the virion close to the cell?
	if( distance_scale > cell_radius - 0.1) //add a small amount for virion 100nm size (not radius, which is increased for visualization)
	{ return false; } 

	// now, get the uptake probability 
	// check binding
	double dt_bind = dt* pAttacker->custom_data[nR_bind]* pAttacker->custom_data[nR_EU]*0.01;

	// don't need to cap it at 1.00: if prob > 100%, 
	// then this statement always evaluates as true, 
	// just the same as capping probability at 100% 
	if( UniformRandom()<= dt_bind )
		{
			// don't attach more virus than the available unbound receptor
			if(pAttacker->custom_data[nR_EU] >= 1)
			{
				pTarget->lyse_cell(); 
				return true; 
			}
	}
	return false; 	
}

void receptor_dynamics_main_model( double dt )
{
	#pragma omp parallel for 
	for( int n=0; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if( pC->phenotype.death.dead == false )
		{ receptor_dynamics_model( pC, pC->phenotype , dt ); }
	}
	
	return; 
}
	
