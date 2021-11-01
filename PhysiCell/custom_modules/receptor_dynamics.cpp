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
	
/* 	// internalize
	double dR_IB = pCell->custom_data[nR_endo]*pCell->custom_data[nR_EB];	
	// viral release from endosomes 	
	double dR_IU = pCell->custom_data[nR_release]*pCell->custom_data[nR_IB];	
	// receptor recycling 	
	double dR_EU = pCell->custom_data[nR_recycle]*pCell->custom_data[nR_IU];
 */
	//int ignore_smoothing_flag=1;
	static int ignore_smoothing_flag = parameters.ints( "ignore_smoothing_flag" ); 
	double x_min = microenvironment.mesh.bounding_box[0]; 
	double x_max = microenvironment.mesh.bounding_box[3]; 
	double y_min = microenvironment.mesh.bounding_box[1]; 
	double y_max = microenvironment.mesh.bounding_box[4]; 
	
    double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*
     					phenotype.volume.total* pCell->custom_data[nR_EU]; //use FE to find what loop to enter
	const double PI  =3.141592653589793238463;
	std::vector<double> cellpos(3, 0);
	std::vector<double> dummypos_y(3, 0);
	std::vector<double> dummypos_x(3, 0);
	std::vector<double> dummypos_xy(3, 0);
	double xc;
	double yc;
	double zc;
	double a;
	double b;
	double c;
	double ze;
	double s;
	double xmesh;
	double ymesh;
	int dummy_voxel_index;
	int dummy_voxel_index_x;
	int dummy_voxel_index_y;
	int dummy_voxel_index_xy;
	double cell_radius;
	double cell_volume;
	double h;
	double VSCapy;
	double fractout;
	double fractout_y;
	double fractout_x;
	double fractout_xy;
	double VSCapx;
	double Vxy;
	double circx;
	double circy;
	double xang;
	double yang;
	
	if (ignore_smoothing_flag > 0.5){
		double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*
		phenotype.volume.total* pCell->custom_data[nR_EU]; //use FE to find what loop to enter
	}
	else {
		cellpos[0]=pCell->position[0]; //x
		cellpos[1]=pCell->position[1]; //y
		cellpos[2]=pCell->position[2]; //z
		
		dummy_voxel_index= microenvironment.nearest_voxel_index( cellpos );
		
		xc = microenvironment.mesh.voxels[dummy_voxel_index].center[0];
		yc = microenvironment.mesh.voxels[dummy_voxel_index].center[1];
		zc = microenvironment.mesh.voxels[dummy_voxel_index].center[2];
		xmesh=microenvironment.mesh.dx;
		ymesh=microenvironment.mesh.dy;
		cell_radius = pCell->phenotype.geometry.radius;
		cell_volume = pCell->phenotype.volume.total;
		
		// current looping assumes that voxel>cell size
		if (cellpos[0] > xc){
			if (cellpos[0] + cell_radius > xc + xmesh/2) { 
				//enter x greater than loop
				dummypos_x[0]=cellpos[0]+xmesh;
				dummypos_x[1]=cellpos[1];
				dummypos_x[2]=cellpos[2];
				dummy_voxel_index_x= microenvironment.nearest_voxel_index( dummypos_x );
				if (cellpos[1] > yc){
					if (cellpos[1] + cell_radius > yc + ymesh/2) { 
						//enter x and y greater than loop
						dummypos_y[0]=cellpos[0];
						dummypos_y[1]=cellpos[1]+ymesh;
						dummypos_y[2]=cellpos[2];
						dummy_voxel_index_y= microenvironment.nearest_voxel_index( dummypos_y );
						//solve for percent inside voxels in y
						h = abs((cellpos[1] + cell_radius)-(yc + ymesh/2));
						VSCapy=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout_y=VSCapy/cell_volume;
						//solve for percent inside voxels in x
						h = abs((cellpos[0] + cell_radius)-(xc + xmesh/2));
						VSCapx=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout_x=VSCapx/cell_volume;
						xang = acos((xc + xmesh/2-cellpos[0])/cell_radius);
						yang = acos((yc + ymesh/2-cellpos[1])/cell_radius);
						if(PI/2<xang+yang){
							//find corner
							dummypos_xy[0]=cellpos[0]+xmesh;
							dummypos_xy[1]=cellpos[1]+ymesh;
							dummypos_xy[2]=cellpos[2];
							dummy_voxel_index_xy= microenvironment.nearest_voxel_index( dummypos_xy );
							//solve for corner
							circx=sqrt(pow(cell_radius,2)-pow(yc + ymesh/2 - cellpos[1],2))+cellpos[0];
							circy=sqrt(pow(cell_radius,2)-pow(xc + xmesh/2 - cellpos[0],2))+cellpos[1];
/* 							a=sqrt(pow(circx-(xc + xmesh/2),2)+pow(yc + ymesh/2-circy,2));
							b=sqrt(pow(yc + ymesh/2-circy,2));
							c=sqrt(pow(circx-(xc + xmesh/2),2));
							s=(pow(a,2)+pow(b,2)+pow(c,2))/2; */
							a=circy - (yc + ymesh/2);
							b=circx - (xc + xmesh/2);
							ze=sqrt(pow(cell_radius,2)-pow(xc + xmesh/2-cellpos[0],2)-pow(yc + ymesh/2-cellpos[1],2));
							//Vxy=sqrt(s*(s-a)*(s-b)*(s-c))*2/3*ze;
							Vxy=a*b*1/3*ze;
							fractout_xy=Vxy/cell_volume;
							if (fractout_xy>0.25) {
								std::cout<<"ISSUE"<<std::endl;
							}
							double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout_x-fractout_y+fractout_xy) *
							microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external]*fractout_x *
							microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external]*fractout_y *
							microenvironment.nearest_density_vector( dummy_voxel_index_xy )[nV_external]*fractout_xy; //use FE to find what loop to enter
						}
						else {
							double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout_x-fractout_y) *
							microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external]*fractout_x *
							microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external]*fractout_y; //use FE to find what loop to enter
						}
					}
					else {
						//solve for percent inside voxels in x
						h = abs((cellpos[0] + cell_radius)-(xc + 10));
						VSCapx=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout=VSCapx/cell_volume;
						dt_bind = dt * pCell->custom_data[nR_bind] * pCell->nearest_density_vector()[nV_external]*(1.0-fractout) *
						microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external]*fractout; //use FE to find what loop to enter
					}
				}
				else {
					if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
						//enter x greater and y less than loop
						dummypos_y[0]=cellpos[0];
						dummypos_y[1]=cellpos[1]-20;
						dummypos_y[2]=cellpos[2];
						dummy_voxel_index_y= microenvironment.nearest_voxel_index( dummypos_y );
						//solve for percent inside voxels in y
						h = abs((yc - 10)-(cellpos[1] + cell_radius));
						VSCapy=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout_y=VSCapy/cell_volume;
						//solve for percent inside voxels in x
						h = abs((cellpos[0] + cell_radius)-(xc + 10));
						VSCapx=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout_x=VSCapx/cell_volume;
						xang = acos((-cellpos[0] + (xc + xmesh/2))/cell_radius);
						yang = acos((-cellpos[1] + (yc - ymesh/2))/cell_radius)-PI/2;
						if(PI/2<xang+yang){
/* 							//find corner
							dummypos_xy[0]=cellpos[0]+xmesh;
							dummypos_xy[1]=cellpos[1]-ymesh;
							dummypos_xy[2]=cellpos[2];
							dummy_voxel_index_xy= microenvironment.nearest_voxel_index( dummypos_xy );
							//solve for corner
							circx=sqrt(pow(cell_radius,2)-pow(yc - ymesh/2 - cellpos[0],2))+cellpos[1];
							circy=sqrt(pow(cell_radius,2)-pow(xc + xmesh/2 - cellpos[1],2))+cellpos[0];
							a=sqrt(pow(circx-(xc + xmesh/2),2)+pow(yc - ymesh/2-circy,2));
							b=sqrt(pow(yc - ymesh/2-circy,2));
							c=sqrt(pow(circx-(xc + xmesh/2),2));
							s=(pow(a,2)+pow(b,2)+pow(c,2))/2;
							ze=sqrt(pow(cell_radius,2)-pow(circx-cellpos[0],2)-pow(circy-cellpos[1],2));
							Vxy=sqrt(s*(s-a)*(s-b)*(s-c))*2/3*ze;
							fractout_xy=Vxy/cell_volume; */
							fractout_xy=0.01;
							double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout_x-fractout_y+fractout_xy) *
							microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external]*fractout_x *
							microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external]*fractout_y *
							microenvironment.nearest_density_vector( dummy_voxel_index_xy )[nV_external]*fractout_xy; //use FE to find what loop to enter
						}
						else {
							double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout_x-fractout_y) *
							microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external]*fractout_x *
							microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external]*fractout_y; //use FE to find what loop to enter
						}
					}
					else {
						//solve for percent inside voxels in x
						h = abs((cellpos[0] + cell_radius)-(xc + 10));
						VSCapx=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout=VSCapx/cell_volume;
						double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout) *
						microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external]*fractout; //use FE to find what loop to enter
					}
				}
			}
			else {
				if (cellpos[1] > yc){
					if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
						//enter y greater than loop
						dummypos_y[0]=cellpos[0];
						dummypos_y[1]=cellpos[1]+20;
						dummypos_y[2]=cellpos[2];
						dummy_voxel_index_y= microenvironment.nearest_voxel_index( dummypos_y );
						
						//solve for percent inside voxels in y
						h = abs((cellpos[1] + cell_radius)-(yc + 10));
						VSCapy=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout=VSCapy/cell_volume;
						double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout) *
						microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external]*fractout; //use FE to find what loop to enter
					}
					else {
						double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*
						phenotype.volume.total* pCell->custom_data[nR_EU]; //use FE to find what loop to enter
					}
				}
				else {
					if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
						//enter y less than loop
						dummypos_y[0]=cellpos[0];
						dummypos_y[1]=cellpos[1]-20;
						dummypos_y[2]=cellpos[2];
						dummy_voxel_index_y= microenvironment.nearest_voxel_index( dummypos_y );
						
						//solve for percent inside voxels in y
						h = abs((yc - 10)-(cellpos[1] + cell_radius));
						VSCapy=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout=VSCapy/cell_volume;
						double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout) *
						microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external]*fractout; //use FE to find what loop to enter
					}
					else {
						double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*
						phenotype.volume.total* pCell->custom_data[nR_EU]; //use FE to find what loop to enter
					}
				}
			}
		}
		else {
			if (cellpos[0] - cell_radius < xc - xmesh/2) { //FIND RELIABLE MESH SIZING REF
				//enter x less than loop
				dummypos_x[0]=cellpos[0]-xmesh;
				dummypos_x[1]=cellpos[1];
				dummypos_x[2]=cellpos[2];
				dummy_voxel_index_x= microenvironment.nearest_voxel_index( dummypos_x );
				if (cellpos[1] > yc){
					if (cellpos[1] + cell_radius > yc + ymesh/2) { //FIND RELIABLE MESH SIZING REF
						//enter y greater and x less than loop
						dummypos_y[0]=cellpos[0];
						dummypos_y[1]=cellpos[1]+ymesh;
						dummypos_y[2]=cellpos[2];
						dummy_voxel_index_y= microenvironment.nearest_voxel_index( dummypos_y );
						
						//solve for percent inside voxels in y
						h = abs((cellpos[1] + cell_radius)-(yc + 10));
						VSCapy=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout_y=VSCapy/cell_volume;
						//solve percent in x
						h = abs(xc - 10 - (cellpos[0] + cell_radius));
						VSCapx=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout_x=VSCapx/cell_volume;
						xang = acos((abs(xc - xmesh/2-cellpos[0]))/cell_radius);
						yang = acos((yc + ymesh/2-cellpos[1])/cell_radius);
						if(PI/2<xang+yang){
 							//find corner
							dummypos_xy[0]=cellpos[0]-xmesh;
							dummypos_xy[1]=cellpos[1]+ymesh;
							dummypos_xy[2]=cellpos[2];
							dummy_voxel_index_xy= microenvironment.nearest_voxel_index( dummypos_xy );
							//solve for corner
							circx=-sqrt(pow(cell_radius,2)-pow(yc + ymesh/2 - cellpos[1],2))+cellpos[0];
							circy=sqrt(pow(cell_radius,2)-pow(xc - xmesh/2 - cellpos[0],2))+cellpos[1];

							a=circy - (yc + ymesh/2);
							b=abs(circx - (xc - xmesh/2));
							ze=sqrt(pow(cell_radius,2)- pow(xc - xmesh/2 - cellpos[0],2)- pow(yc + ymesh/2 - cellpos[1],2));
							Vxy=a*b*1/3*ze;
							//fractout_xy=Vxy/cell_volume; 
							fractout_xy=0.01;
							//std::cout<<ze<<std::endl;
							double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout_x-fractout_y+fractout_xy) *
							microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external]*fractout_x *
							microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external]*fractout_y *
							microenvironment.nearest_density_vector( dummy_voxel_index_xy )[nV_external]*fractout_xy; //use FE to find what loop to enter
						}
						else {
							double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout_x-fractout_y) *
							microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external]*fractout_x *
							microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external]*fractout_y; //use FE to find what loop to enter
						}
					}
					else {
						//solve percent in x
						h = abs(xc - 10 - (cellpos[0] + cell_radius));
						VSCapx=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout=VSCapx/cell_volume;
						double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout) *
						microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external]*fractout; //use FE to find what loop to enter
					}
				}
				else {
					if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
						//enter x and y less than loop
						dummypos_y[0]=cellpos[0];
						dummypos_y[1]=cellpos[1]-20;
						dummypos_y[2]=cellpos[2];
						dummy_voxel_index_y= microenvironment.nearest_voxel_index( dummypos_y );
						
						//solve for percent inside voxels in y
						h = abs((yc - 10)-(cellpos[1] + cell_radius));
						VSCapy=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout_y=VSCapy/cell_volume;
						//solve percent in x
						h = abs(xc - 10 - (cellpos[0] + cell_radius));
						VSCapx=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout_x=VSCapx/cell_volume;
						xang = acos((-cellpos[0] + (xc - xmesh/2))/cell_radius)-PI/2;
						yang = acos((-cellpos[1] + (yc - ymesh/2))/cell_radius)-PI/2;
						if(PI/2<xang+yang){
/* 							//find corner
							dummypos_xy[0]=cellpos[0]-xmesh;
							dummypos_xy[1]=cellpos[1]-ymesh;
							dummypos_xy[2]=cellpos[2];
							dummy_voxel_index_xy= microenvironment.nearest_voxel_index( dummypos_xy );
							//solve for corner
							circx=sqrt(pow(cell_radius,2)-pow(yc - ymesh/2 - cellpos[0],2))+cellpos[1];
							circy=sqrt(pow(cell_radius,2)-pow(xc - xmesh/2 - cellpos[1],2))+cellpos[0];
							a=sqrt(pow(circx-(xc - xmesh/2),2)+pow(yc - ymesh/2-circy,2));
							b=sqrt(pow(yc - ymesh/2-circy,2));
							c=sqrt(pow(circx-(xc - xmesh/2),2));
							s=(pow(a,2)+pow(b,2)+pow(c,2))/2;
							ze=sqrt(pow(cell_radius,2)-pow(circx-cellpos[0],2)-pow(circy-cellpos[1],2));
							Vxy=sqrt(s*(s-a)*(s-b)*(s-c))*2/3*ze;
							fractout_xy=Vxy/cell_volume; */
							fractout_xy=0.01;
							double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout_x-fractout_y+fractout_xy) *
							microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external]*fractout_x *
							microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external]*fractout_y *
							microenvironment.nearest_density_vector( dummy_voxel_index_xy )[nV_external]*fractout_xy; //use FE to find what loop to enter
						}
						else {
							double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout_x-fractout_y) *
							microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external]*fractout_x *
							microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external]*fractout_y; //use FE to find what loop to enter
						}
					}
					else {
						//solve percent in x
						h = abs(xc - 10 - (cellpos[0] + cell_radius));
						VSCapx=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout=VSCapx/cell_volume;
						double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout) *
						microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external]*fractout; //use FE to find what loop to enter
					}
				}
			}
			else {
				if (cellpos[1] > yc){
					if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
						//enter y greater than loop
						dummypos_y[0]=cellpos[0];
						dummypos_y[1]=cellpos[1]+20;
						dummypos_y[2]=cellpos[2];
						dummy_voxel_index_y= microenvironment.nearest_voxel_index( dummypos_y );
						
						//solve for percent inside voxels in y
						h = abs((cellpos[1] + cell_radius)-(yc + 10));
						VSCapy=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout=VSCapy/cell_volume;
						double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout) *
						microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external]*fractout; //use FE to find what loop to enter
					}
					else {
						double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*
						phenotype.volume.total* pCell->custom_data[nR_EU]; //use FE to find what loop to enter
					}
				}
				else {
					if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
						//enter y less than loop
						dummypos_y[0]=cellpos[0];
						dummypos_y[1]=cellpos[1]-20;
						dummypos_y[2]=cellpos[2];
						dummy_voxel_index_y= microenvironment.nearest_voxel_index( dummypos_y );
						
						//solve for percent inside voxels in y
						h = abs((yc - 10)-(cellpos[1] + cell_radius));
						VSCapy=1/3*PI*pow(h,2)*(3*cell_radius-h);
						fractout=VSCapy/cell_volume;
						double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*(1.0-fractout) *
						microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external]*fractout; //use FE to find what loop to enter
					}
					else {
						double dt_bind = dt* pCell->custom_data[nR_bind]* pCell->nearest_density_vector()[nV_external]*
						phenotype.volume.total* pCell->custom_data[nR_EU]; //use FE to find what loop to enter
					}
				}
			}
		}	
	}
	if( dt_bind<1 )
	{
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
		
		#pragma omp critical
		{ pCell->nearest_density_vector()[nV_external] += dt*(f[0][5]/6+f[1][5]/3+f[2][5]/3+f[3][5]/6) / microenvironment.mesh.dV; }
		
		//START STOCHASTIC PORTION
		if( dt_bind>0 && UniformRandom()<= dt_bind )
		{
			// don't attach more virus than the available unbound receptor
			if(pCell->custom_data[nR_EU] >= 1)
			{
				pCell->custom_data[nR_EU] -= 1;
				pCell->custom_data[nR_EB] += 1;

				if (ignore_smoothing_flag > 0.5){
					#pragma omp critical
					{pCell->nearest_density_vector()[nV_external] -= 1.0 / microenvironment.mesh.dV;}
				}
				else {
					// current looping assumes that voxel>cell size
					if (cellpos[0] > xc){
						if (cellpos[0] + cell_radius > xc + 10) { //FIND RELIABLE MESH SIZING REF
							//enter x greater than loop
							if (cellpos[1] > yc){
								if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
									//enter x and y greater than loop
									if(PI/2<xang+yang){
										#pragma omp critical
										{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout_x-fractout_y+fractout_xy) / microenvironment.mesh.dV;
										microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout_x / microenvironment.mesh.dV;
										microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout_y / microenvironment.mesh.dV;
										microenvironment.nearest_density_vector( dummy_voxel_index_xy )[nV_external] -= fractout_xy / microenvironment.mesh.dV;}
									}
									else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout_x-fractout_y) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout_x / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout_y / microenvironment.mesh.dV;}
									}
								}
								else {
									//solve for percent inside voxels in x
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
							}
							else {
								if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
									//enter x greater and y less than loop
									if(PI/2<xang+yang){
										#pragma omp critical
										{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout_x-fractout_y+fractout_xy) / microenvironment.mesh.dV;
										microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout_x / microenvironment.mesh.dV;
										microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout_y / microenvironment.mesh.dV;
										microenvironment.nearest_density_vector( dummy_voxel_index_xy )[nV_external] -= fractout_xy / microenvironment.mesh.dV;}
									}
									else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout_x-fractout_y) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout_x / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout_y / microenvironment.mesh.dV;}
									}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
							}
						}
						else {
							if (cellpos[1] > yc){
								if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
									//enter y greater than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= 1.0 / microenvironment.mesh.dV;}
								}
							}
							else {
								if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
									//enter y less than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= 1.0 / microenvironment.mesh.dV;}
								}
							}
						}
					}
					else {
						if (cellpos[0] - cell_radius < xc - 10) { //FIND RELIABLE MESH SIZING REF
							//enter x less than loop
							if (cellpos[1] > yc){
								if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
									//enter y greater and x less than loop
									if(PI/2<xang+yang){
										#pragma omp critical
										{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout_x-fractout_y+fractout_xy) / microenvironment.mesh.dV;
										microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout_x / microenvironment.mesh.dV;
										microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout_y / microenvironment.mesh.dV;
										microenvironment.nearest_density_vector( dummy_voxel_index_xy )[nV_external] -= fractout_xy / microenvironment.mesh.dV;}
									}
									else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout_x-fractout_y) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout_x / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout_y / microenvironment.mesh.dV;}
									}								}
								else {
									//solve percent in x
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
							}
							else {
								if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
									//enter x and y less than loop
									if(PI/2<xang+yang){
										#pragma omp critical
										{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout_x-fractout_y+fractout_xy) / microenvironment.mesh.dV;
										microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout_x / microenvironment.mesh.dV;
										microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout_y / microenvironment.mesh.dV;
										microenvironment.nearest_density_vector( dummy_voxel_index_xy )[nV_external] -= fractout_xy / microenvironment.mesh.dV;}
									}
									else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout_x-fractout_y) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout_x / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout_y / microenvironment.mesh.dV;}
									}								}
								else {
									//solve percent in x
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
							}
						}
						else {
							if (cellpos[1] > yc){
								if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
									//enter y greater than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= 1.0 / microenvironment.mesh.dV;}
								}
							}
							else {
								if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
									//enter y less than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= 1.0 / microenvironment.mesh.dV;}
								}
							}
						}
					}	
				}
			}
		}
		
	}

	if( dt_bind>1 )
	{
		//SOLVE ODE BEFORE STOCHASTIC PORTION
		//THIS RK4 METHOD SOLVES THE BINDING IN A DETERMINISTIC FASHON
	
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
		
		//attempting proper integration to stochastic portion, still needs some thought
		double alpha = dt_bind;
		double n_virion = pCell->nearest_density_vector()[nV_external]* microenvironment.mesh.dV;
		
		//not used
		#pragma omp critical
		{ pCell->nearest_density_vector()[nV_external] += dt*(f[0][5]/6+f[1][5]/3+f[2][5]/3+f[3][5]/6) / microenvironment.mesh.dV; }
		
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
		if (ignore_smoothing_flag > 0.5){
			#pragma omp critical
			{pCell->nearest_density_vector()[nV_external] -= alpha1 / microenvironment.mesh.dV;}
		}
		else {
			// current looping assumes that voxel>cell size
			if (cellpos[0] > xc){
						if (cellpos[0] + cell_radius > xc + 10) { //FIND RELIABLE MESH SIZING REF
							//enter x greater than loop
							if (cellpos[1] > yc){
								if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
									//enter x and y greater than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1*(1.0-fractout_x-fractout_y) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= alpha1*fractout_x / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= alpha1*fractout_y / microenvironment.mesh.dV;}
								}
								else {
									//solve for percent inside voxels in x
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1*(1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= alpha1*fractout / microenvironment.mesh.dV;}
								}
							}
							else {
								if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
									//enter x greater and y less than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1*(1.0-fractout_x-fractout_y) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= alpha1*fractout_x / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= alpha1*fractout_y / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1*(1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= alpha1*fractout / microenvironment.mesh.dV;}
								}
							}
						}
						else {
							if (cellpos[1] > yc){
								if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
									//enter y greater than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1*(1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= alpha1*fractout / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1 / microenvironment.mesh.dV;}
								}
							}
							else {
								if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
									//enter y less than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1*(1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= alpha1*fractout / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1 / microenvironment.mesh.dV;}
								}
							}
						}
					}
					else {
						if (cellpos[0] - cell_radius < xc - 10) { //FIND RELIABLE MESH SIZING REF
							//enter x less than loop
							if (cellpos[1] > yc){
								if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
									//enter y greater and x less than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1*(1.0-fractout_x-fractout_y) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= alpha1*fractout_x / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= alpha1*fractout_y / microenvironment.mesh.dV;}
								}
								else {
									//solve percent in x
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1*(1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= alpha1*fractout / microenvironment.mesh.dV;}
								}
							}
							else {
								if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
									//enter x and y less than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1*(1.0-fractout_x-fractout_y) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= alpha1*fractout_x / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= alpha1*fractout_y / microenvironment.mesh.dV;}
								}
								else {
									//solve percent in x
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1*(1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= alpha1*fractout / microenvironment.mesh.dV;}
								}
							}
						}
						else {
							if (cellpos[1] > yc){
								if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
									//enter y greater than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1*(1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= alpha1*fractout / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1 / microenvironment.mesh.dV;}
								}
							}
							else {
								if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
									//enter y less than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1*(1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= alpha1*fractout / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= alpha1 / microenvironment.mesh.dV;}
								}
							}
						}
					}	
				}
	    
	    if( UniformRandom()<= alpha2 )
	    {
			// don't attach more virus than the available unbound receptor
			if(pCell->custom_data[nR_EU] >= 1)
			{
				pCell->custom_data[nR_EU] -= 1;
				pCell->custom_data[nR_EB] += 1;
				
				if (ignore_smoothing_flag > 0.5){
					#pragma omp critical
					{pCell->nearest_density_vector()[nV_external] -= 1.0 / microenvironment.mesh.dV;}
				}
				else {
					// current looping assumes that voxel>cell size
					if (cellpos[0] > xc){
						if (cellpos[0] + cell_radius > xc + 10) { //FIND RELIABLE MESH SIZING REF
							//enter x greater than loop
							if (cellpos[1] > yc){
								if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
									//enter x and y greater than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout_x-fractout_y) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout_x / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout_y / microenvironment.mesh.dV;}
								}
								else {
									//solve for percent inside voxels in x
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
							}
							else {
								if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
									//enter x greater and y less than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout_x-fractout_y) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout_x / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout_y / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
							}
						}
						else {
							if (cellpos[1] > yc){
								if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
									//enter y greater than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= 1.0 / microenvironment.mesh.dV;}
								}
							}
							else {
								if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
									//enter y less than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= 1.0 / microenvironment.mesh.dV;}
								}
							}
						}
					}
					else {
						if (cellpos[0] - cell_radius < xc - 10) { //FIND RELIABLE MESH SIZING REF
							//enter x less than loop
							if (cellpos[1] > yc){
								if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
									//enter y greater and x less than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout_x-fractout_y) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout_x / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout_y / microenvironment.mesh.dV;}
								}
								else {
									//solve percent in x
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
							}
							else {
								if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
									//enter x and y less than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout_x-fractout_y) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout_x / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout_y / microenvironment.mesh.dV;}
								}
								else {
									//solve percent in x
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_x )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
							}
						}
						else {
							if (cellpos[1] > yc){
								if (cellpos[1] + cell_radius > yc + 10) { //FIND RELIABLE MESH SIZING REF
									//enter y greater than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= 1.0 / microenvironment.mesh.dV;}
								}
							}
							else {
								if (cellpos[1] - cell_radius < yc - 10) { //FIND RELIABLE MESH SIZING REF
									//enter y less than loop
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= (1.0-fractout) / microenvironment.mesh.dV;
									microenvironment.nearest_density_vector( dummy_voxel_index_y )[nV_external] -= fractout / microenvironment.mesh.dV;}
								}
								else {
									#pragma omp critical
									{pCell->nearest_density_vector()[nV_external] -= 1.0 / microenvironment.mesh.dV;}
								}
							}
						}
					}	
				}
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
	
