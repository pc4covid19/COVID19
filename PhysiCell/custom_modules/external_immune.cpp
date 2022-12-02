#include "./external_immune.h" 
#include <algorithm> 

using namespace PhysiCell; 

std::string external_immune_version = "0.6.0"; 

Submodel_Information external_immune_info; 


void external_immune_model_setup( void )
{
		// set version
	external_immune_info.name = "external immune"; 
	external_immune_info.version = external_immune_version; 
		// set functions 
	external_immune_info.main_function = external_immune_model; 
	external_immune_info.phenotype_function = NULL; 
	external_immune_info.mechanics_function = NULL; 
		// what microenvironment variables do I need? 

		// what custom data do I need? 
	//external_immune_info.parameters.doubles.push_back( "DM" );
	//external_immune_info.parameters.doubles.push_back( "TC" );

	// submodel_registry.register_model( internal_viral_dynamics_info ); 
	external_immune_info.register_model();
	
	return; 
}

void external_immune_model( double dt )
{
	// bookkeeping -- find microenvironment variables we need

	extern double DM;
	extern double DL;
	extern double TC;
	extern double TH1;
	extern double TH2;
	extern double TCt;
	extern double Tht;
	extern double Bc;
	extern double Ps;
	extern double Ig;
	extern double TCN;
	extern double THN;
	extern double BN;
	extern double EPICOUNT;
	extern double tissueCD8; 
	extern double tissueCD4; 
	
	static double dC = parameters.doubles( "TC_death_rate" ); 
	static double rT1 = parameters.doubles( "max_activation_TC" ); 
	static double rT2 = parameters.doubles( "half_max_activation_TC" ); 
	static double dT1 = parameters.doubles( "max_clearance_TC" ); 
	static double dT2 = parameters.doubles( "half_max_clearance_TC" ); 
	static double Tc0 = parameters.doubles( "TC_population_threshold" ); 
	static double dDm = parameters.doubles( "DM_decay" );
	static double sTh1 = parameters.doubles( "Th1_max_activation" );
	static double pTh1 = parameters.doubles( "Th1_damping" );
	static double dTh1 = parameters.doubles( "Th1_decay" );
	static double mTh = parameters.doubles( "Th_base_decay" );
	static double sTh2 = parameters.doubles( "Th2_self_feeback" );
	static double pTh2 = parameters.doubles( "Th2_max_conversion" );
	static double CD8_Tcell_recruitment_rate = parameters.doubles( "T_Cell_Recruitment" ); 
	static double dB = parameters.doubles( "BCell_removal" );
	static double Bh = parameters.doubles( "BCell_activation" );
	static double Bh2 = parameters.doubles( "BCell_activation_half" );
	static double rB1 = parameters.doubles( "BCell_DC_proliferation" );
	static double h = parameters.doubles( "BCell_Th2_wieght_function" );
	static double rB2 = parameters.doubles( "BCell_damping" );
	static double pSc = parameters.doubles( "PCell_recuitment" );
	static double dS = parameters.doubles( "PCell_degradation" );
	static double pAS = parameters.doubles( "Ig_recuitment" );
	static double dMc = parameters.doubles( "Ig_degradation" );
	static double pT = parameters.doubles( "TC_activation" );
	static double pT2 = parameters.doubles( "TC_activation_half" );
	static double kh = parameters.doubles( "Th_half" );
	static double r = parameters.doubles( "transport_ratio" )*(1-0.95)/0.95;
		
	double lypmh_scale = EPICOUNT / 500000;
	// actual model goes here 
	
	double x[4][13]={{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};//initialize x
	double f[4][13]={{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};//initialize f
	int j;
	
	extern std::vector<int>history;
	
	static double timedelay = parameters.doubles( "Lymph_node_td" ); 
	double td_l = round(timedelay*1440/dt);
	//std::cout<<history[td_l]<<std::endl;
	
	x[0][0] = (DM+history[td_l])/lypmh_scale; 
	x[0][1] = TC; //initial values
	x[0][2] = TH1; //initial values
	x[0][3] = TH2; //initial values
	x[0][4] = TCt/lypmh_scale;
	x[0][5] = Tht/lypmh_scale;
	x[0][6] = Bc;
	x[0][7] = Ps;
	x[0][8] = Ig/lypmh_scale;
	x[0][9] = TCN;
	x[0][10] = THN;
	x[0][11] = BN;
	x[0][12] = (DL+(THN+TH1+TH2)/(THN+TH1+TH2+5000)*history[td_l])/lypmh_scale;
	//(pT/(y(5)+pT2) added
	//rT1->pT1
	//rT2->pT2
	static double ro=pTh1/(10*pTh2);
	
	
	
    for(j = 0; j < 4; j++){
		f[j][0] = {-dDm*x[j][0]}; //DM
        	f[j][1] = {pT/(x[j][12]+pT2)*x[j][9]*x[j][12]+rT1*x[j][1]*(1E6-x[j][1])/1E6*x[j][12]/(x[j][12]+rT2) - dT1*x[j][1]*x[j][12]/(x[j][12]+dT2)}; //Tc
		f[j][2] = {sTh1*x[j][0]*x[j][10]*(0.1/kh*(1+kh*x[j][3])+x[j][2])/((1+kh*x[j][3])*(1+kh*x[j][3]))+(pTh1*kh*x[j][0]*x[j][2]*x[j][2])/((1+kh*x[j][3])*(1+kh*x[j][3]))-(dTh1*(kh*kh)*x[j][0]*x[j][2]*x[j][2]*x[j][2])/(1+kh*x[j][3])-mTh*x[j][2]}; //Th1
		f[j][3] = {sTh2*x[j][0]*x[j][10]*(0.1/kh+x[j][3])/(1+kh*x[j][3])+(pTh2*(kh*kh)*(ro+x[j][2]/(1+kh*x[j][3]))*x[j][0]*x[j][3]*x[j][3])/((1+kh*x[j][2]+kh*x[j][3]))-mTh*x[j][3]}; //Th2
		f[j][4] = {CD8_Tcell_recruitment_rate*(1-(tissueCD8/lypmh_scale)/(r*x[j][1]+(tissueCD8/lypmh_scale)+1E-4))*x[j][1]}; //CD8 export
		f[j][5] = {CD8_Tcell_recruitment_rate*(1-(tissueCD4/lypmh_scale)/(r*(x[j][2]+x[j][3])+(tissueCD4/lypmh_scale)+1E-4))*(x[j][2]+x[j][3])}; //CD4 export
		f[j][6] = {Bh*x[j][0]/(x[j][0]+Bh2)*x[j][11]+rB1*x[j][6]*(x[j][0]+h*x[j][3])/(x[j][0]+h*x[j][3]+rB2)-dB*x[j][6]-2*pSc*x[j][6]}; //B-Cell
		f[j][7] = {pSc*x[j][6]-dS*x[j][7]}; //P-Cell
		f[j][8] = {pAS*x[j][7]-dMc*x[j][8]}; //Ig
		f[j][9] = {dC*(Tc0-x[j][9]) - pT/(x[j][0]+pT2)*x[j][9]*x[j][0]}; //TcN
		f[j][10] = {dC*(Tc0-x[j][10]) - sTh1*x[j][0]*x[j][10]*(0.1/kh*(1+kh*x[j][3])+x[j][2])/((1+kh*x[j][3])*(1+kh*x[j][3])) - sTh2*x[j][0]*x[j][10]*(0.1/kh+x[j][3])/(1+kh*x[j][3])}; //ThN
		f[j][11] = {dC*(Tc0-x[j][11])-Bh*x[j][0]/(x[j][0]+Bh2)*x[j][11]}; //bN
		f[j][12] = {-dDm*x[j][12]}; //DM
        if (j== 0 || j==1){
            x[j+1][0]=x[0][0]+dt/2*f[j][0]; //first and second x approximations
			x[j+1][1]=x[0][1]+dt/2*f[j][1]; //first and second x approximations
			x[j+1][2]=x[0][2]+dt/2*f[j][2]; //first and second x approximations
			x[j+1][3]=x[0][3]+dt/2*f[j][3]; //first and second x approximations
			x[j+1][4]=x[0][4]+dt/2*f[j][4]; //first and second x approximations
			x[j+1][5]=x[0][5]+dt/2*f[j][5]; //first and second x approximations
			x[j+1][6]=x[0][6]+dt/2*f[j][6]; //first and second x approximations
			x[j+1][7]=x[0][7]+dt/2*f[j][7]; //first and second x approximations
			x[j+1][8]=x[0][8]+dt/2*f[j][8]; //first and second x approximations
			x[j+1][9]=x[0][9]+dt/2*f[j][9]; //first and second x approximations
			x[j+1][10]=x[0][10]+dt/2*f[j][10]; //first and second x approximations
			x[j+1][11]=x[0][11]+dt/2*f[j][11]; //first and second x approximations
			x[j+1][12]=x[0][12]+dt/2*f[j][12]; //first and second x approximations
		}
        if (j== 2){
            x[j+1][0]=x[0][0]+dt*f[j][0]; //third approximation
			x[j+1][1]=x[0][1]+dt*f[j][1]; //third approximation
			x[j+1][2]=x[0][2]+dt*f[j][2]; //third approximation
			x[j+1][3]=x[0][3]+dt*f[j][3]; //third approximation
			x[j+1][4]=x[0][4]+dt*f[j][4]; //third approximation
			x[j+1][5]=x[0][5]+dt*f[j][5]; //third approximation
			x[j+1][6]=x[0][6]+dt*f[j][6]; //third approximation
			x[j+1][7]=x[0][7]+dt*f[j][7]; //third approximation
			x[j+1][8]=x[0][8]+dt*f[j][8]; //third approximation
			x[j+1][9]=x[0][9]+dt*f[j][9]; //third approximation
			x[j+1][10]=x[0][10]+dt*f[j][10]; //third approximation
			x[j+1][11]=x[0][11]+dt*f[j][11]; //third approximation
			x[j+1][12]=x[0][12]+dt*f[j][12]; //third approximation
		}
    } 

	DM=(x[0][0]+dt*(f[0][0]/6+f[1][0]/3+f[2][0]/3+f[3][0]/6))*lypmh_scale;
	TC=x[0][1]+dt*(f[0][1]/6+f[1][1]/3+f[2][1]/3+f[3][1]/6);
	TH1=x[0][2]+dt*(f[0][2]/6+f[1][2]/3+f[2][2]/3+f[3][2]/6); //detirmine n+1
	TH2=x[0][3]+dt*(f[0][3]/6+f[1][3]/3+f[2][3]/3+f[3][3]/6); //detirmine n+1
	TCt=(x[0][4]+dt*(f[0][4]/6+f[1][4]/3+f[2][4]/3+f[3][4]/6))*lypmh_scale;
	Tht=(x[0][5]+dt*(f[0][5]/6+f[1][5]/3+f[2][5]/3+f[3][5]/6))*lypmh_scale;
	Bc=x[0][6]+dt*(f[0][6]/6+f[1][6]/3+f[2][6]/3+f[3][6]/6);
	Ps=x[0][7]+dt*(f[0][7]/6+f[1][7]/3+f[2][7]/3+f[3][7]/6);
	Ig=(x[0][8]+dt*(f[0][8]/6+f[1][8]/3+f[2][8]/3+f[3][8]/6))*lypmh_scale;
	TCN=(x[0][9]+dt*(f[0][9]/6+f[1][9]/3+f[2][9]/3+f[3][9]/6));
	THN=(x[0][10]+dt*(f[0][10]/6+f[1][10]/3+f[2][10]/3+f[3][10]/6));
	BN=(x[0][11]+dt*(f[0][11]/6+f[1][11]/3+f[2][11]/3+f[3][11]/6));
	DL=(x[0][12]+dt*(f[0][12]/6+f[1][12]/3+f[2][12]/3+f[3][12]/6))*lypmh_scale;

	double x_min = microenvironment.mesh.bounding_box[0] + 1e-6; 
	double x_max = microenvironment.mesh.bounding_box[3] - 1e-6; 
	double y_min = microenvironment.mesh.bounding_box[1] + 1e-6; 
	double y_max = microenvironment.mesh.bounding_box[4] - 1e-6; 
	
	double number_of_Ig=floor(Ig);
	Ig -= number_of_Ig;
	
	
	static int nAb = microenvironment.find_density_index( "Ig" ); 
	static int nV = microenvironment.find_density_index( "virion" ); 

	static double kv = parameters.doubles("Ig_neutralization_rate");
	
	for( int n=0 ; n < number_of_Ig ; n++ )
		{
			// pick a random voxel 
			std::vector<double> position = {0,0,0}; 
			position[0] = x_min + (x_max-x_min)*UniformRandom(); 
			position[1] = y_min + (y_max-y_min)*UniformRandom(); 
			
			int m = microenvironment.nearest_voxel_index( position ); 

			microenvironment(m)[nAb] += 1.0 / microenvironment.mesh.dV; 
		}
	#pragma omp parallel for
	for( int n=0 ; n < microenvironment.number_of_voxels() ; n++ )
		{
			if (microenvironment(n)[nV]>0 && microenvironment(n)[nAb]>0) {
			double rate = kv * microenvironment(n)[nAb] * microenvironment(n)[nV] * dt; 
			if (rate < 0) {
				rate=0;
			}
			if (rate > microenvironment(n)[nAb]) {
				rate = microenvironment(n)[nAb];
			}
			if (rate > microenvironment(n)[nV]) {
				rate = microenvironment(n)[nV];
			}
			microenvironment(n)[nAb] -= rate;
			microenvironment(n)[nV] -= rate;
			}
		}
	
	return; 
}
	
