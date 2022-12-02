#include "./receptor_dynamics.h" 

using namespace PhysiCell; 

std::string receptor_model_version = "0.6.0"; 

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
	
	double x[4][5]={{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};//initialize x
	double f[4][5]={{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};//initialize f
	int j;//initialize counter
	
	//initial values for RK4
	x[0][0] = pCell->custom_data[nR_EU]; 
	x[0][1] = pCell->custom_data[nR_EB];
	x[0][2] = pCell->custom_data[nR_IB]; 
	x[0][3] = pCell->custom_data[nR_IU]; 
	x[0][4] = pCell->custom_data[nV_internal];
	
/* 	// internalize
	double dR_IB = pCell->custom_data[nR_endo]*pCell->custom_data[nR_EB];	
	// viral release from endosomes 	
	double dR_IU = pCell->custom_data[nR_release]*pCell->custom_data[nR_IB];	
	// receptor recycling 	
	double dR_EU = pCell->custom_data[nR_recycle]*pCell->custom_data[nR_IU];
 */
	
    double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]* pCell->custom_data[nR_EU]; //use FE to find what loop to enter
	if( dt_bind<1 )
	{
		//SOLVE ODE BEFORE STOCHASTIC PORTION
	
		for(j = 0; j < 4; j++){
			f[j][0] = {pCell->custom_data[nR_recycle]*x[j][3]}; //define SPECIAL function
			f[j][1] = {-pCell->custom_data[nR_endo]*x[j][1]}; //define SPECIAL function
			f[j][2] = {pCell->custom_data[nR_endo]*x[j][1]-pCell->custom_data[nR_release]*x[j][2]}; //define function
			f[j][3] = {pCell->custom_data[nR_release]*x[j][2]-pCell->custom_data[nR_recycle]*x[j][3]}; //define function
			f[j][4] = {pCell->custom_data[nR_release]*x[j][2]}; //define function
			if (j== 0 || j==1){
				x[j+1][0]=x[0][0]+dt/2*f[j][0]; //first and second x approximations
				x[j+1][1]=x[0][1]+dt/2*f[j][1]; //first and second x approximations
				x[j+1][2]=x[0][2]+dt/2*f[j][2]; //first and second x approximations
				x[j+1][3]=x[0][3]+dt/2*f[j][3]; //first and second x approximations
				x[j+1][4]=x[0][4]+dt/2*f[j][4]; //first and second x approximations
			}
			else if (j== 2){
				x[j+1][0]=x[0][0]+dt*f[j][0]; //third approximation
				x[j+1][1]=x[0][1]+dt*f[j][1]; //third approximation
				x[j+1][2]=x[0][2]+dt*f[j][2]; //third approximation
				x[j+1][3]=x[0][3]+dt*f[j][3]; //third approximation
				x[j+1][4]=x[0][4]+dt*f[j][4]; //third approximation
			}
		}

		pCell->custom_data[nR_EU]=(x[0][0]+dt*(f[0][0]/6+f[1][0]/3+f[2][0]/3+f[3][0]/6));
		pCell->custom_data[nR_EB]=x[0][1]+dt*(f[0][1]/6+f[1][1]/3+f[2][1]/3+f[3][1]/6);
		pCell->custom_data[nR_IB]=x[0][2]+dt*(f[0][2]/6+f[1][2]/3+f[2][2]/3+f[3][2]/6); //detirmine n+1
		pCell->custom_data[nR_IU]=x[0][3]+dt*(f[0][3]/6+f[1][3]/3+f[2][3]/3+f[3][3]/6); //detirmine n+1
		pCell->custom_data[nV_internal]=x[0][4]+dt*(f[0][4]/6+f[1][4]/3+f[2][4]/3+f[3][4]/6); //detirmine n+1
		
		//START STOCHASTIC PORTION
		if( dt_bind>0 && UniformRandom()<= dt_bind )
		{
			// don't attach more virus than the available unbound receptor
			if(pCell->custom_data[nR_EU] >= 1)
			{
				pCell->custom_data[nR_EU] -= 1;
				pCell->custom_data[nR_EB] += 1;

				#pragma omp critical
				{pCell->nearest_density_vector()[nV_external] -= 1.0 / microenvironment.mesh.dV;}
			}
		}
		
	}

	if( dt_bind>1 )
	{
		//SOLVE ODE BEFORE STOCHASTIC PORTION
		//THIS RK4 METHOD SOLVES THE BINDING IN A DETERMINISTIC FASHON
	
		for(j = 0; j < 4; j++){
			f[j][0] = {pCell->custom_data[nR_recycle]*x[j][3]}; //define SPECIAL function
			f[j][1] = {-pCell->custom_data[nR_endo]*x[j][1]}; //define SPECIAL function
			f[j][2] = {pCell->custom_data[nR_endo]*x[j][1]-pCell->custom_data[nR_release]*x[j][2]}; //define function
			f[j][3] = {pCell->custom_data[nR_release]*x[j][2]-pCell->custom_data[nR_recycle]*x[j][3]}; //define function
			f[j][4] = {pCell->custom_data[nR_release]*x[j][2]}; //define function
			if (j== 0 || j==1){
				x[j+1][0]=x[0][0]+dt/2*f[j][0]; //first and second x approximations
				x[j+1][1]=x[0][1]+dt/2*f[j][1]; //first and second x approximations
				x[j+1][2]=x[0][2]+dt/2*f[j][2]; //first and second x approximations
				x[j+1][3]=x[0][3]+dt/2*f[j][3]; //first and second x approximations
				x[j+1][4]=x[0][4]+dt/2*f[j][4]; //first and second x approximations
			}
			else if (j== 2){
				x[j+1][0]=x[0][0]+dt*f[j][0]; //third approximation
				x[j+1][1]=x[0][1]+dt*f[j][1]; //third approximation
				x[j+1][2]=x[0][2]+dt*f[j][2]; //third approximation
				x[j+1][3]=x[0][3]+dt*f[j][3]; //third approximation
				x[j+1][4]=x[0][4]+dt*f[j][4]; //third approximation
			}
		}

		pCell->custom_data[nR_EU]=(x[0][0]+dt*(f[0][0]/6+f[1][0]/3+f[2][0]/3+f[3][0]/6));
		pCell->custom_data[nR_EB]=x[0][1]+dt*(f[0][1]/6+f[1][1]/3+f[2][1]/3+f[3][1]/6);
		pCell->custom_data[nR_IB]=x[0][2]+dt*(f[0][2]/6+f[1][2]/3+f[2][2]/3+f[3][2]/6); //detirmine n+1
		pCell->custom_data[nR_IU]=x[0][3]+dt*(f[0][3]/6+f[1][3]/3+f[2][3]/3+f[3][3]/6); //detirmine n+1
		pCell->custom_data[nV_internal]=x[0][4]+dt*(f[0][4]/6+f[1][4]/3+f[2][4]/3+f[3][4]/6); //detirmine n+1
		
		//attempting proper integration to stochastic portion, still needs some thought
		double alpha = dt_bind;
		double n_virion = pCell->nearest_density_vector()[nV_external]* microenvironment.mesh.dV;
		
		//limit to number of virons in a voxel
		if(alpha > n_virion)
        {
        	alpha = n_virion;
        }
	    double alpha1 = floor( alpha );
	    double alpha2 = alpha - alpha1;
		
		//limit to number of unbound receptor on cell surface
        if(alpha1 > pCell->custom_data[nR_EU])
        {
        	alpha1 = pCell->custom_data[nR_EU];
		}	
		
		//STOCHASTIC PORTION

        pCell->custom_data[nR_EU] -= alpha1;
		pCell->custom_data[nR_EB] += alpha1;
			
		#pragma omp critical
		{pCell->nearest_density_vector()[nV_external] -= alpha1 / microenvironment.mesh.dV;}
	    
	    if( UniformRandom()<= alpha2 )
	    {
			// don't attach more virus than the available unbound receptor
			if(pCell->custom_data[nR_EU] >= 1)
			{
				pCell->custom_data[nR_EU] -= 1;
				pCell->custom_data[nR_EB] += 1;
				
				#pragma omp critical
				{pCell->nearest_density_vector()[nV_external] -= 1.0 / microenvironment.mesh.dV;}
			}
				
	    }
	} 
	
	return; 
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
	
