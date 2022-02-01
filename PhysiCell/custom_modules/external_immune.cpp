#include "./external_immune.h" 
#include <algorithm> 

using namespace PhysiCell; 

std::string external_immune_version = "0.5.0"; 

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
	static double ro = parameters.doubles( "Th1_Th2_conversion_wieght" );
	static double CD8_Tcell_recruitment_rate = parameters.doubles( "T_Cell_Recruitment" ); 
	static double dB = parameters.doubles( "BCell_base_rate" );
	static double B0 = parameters.doubles( "BCell_base_value" );
	static double rB1 = parameters.doubles( "BCell_DC_proliferation" );
	static double h = parameters.doubles( "BCell_Th2_wieght_function" );
	static double rB2 = parameters.doubles( "BCell_damping" );
	static double pSc = parameters.doubles( "PCell_recuitment" );
	static double dS = parameters.doubles( "PCell_degradation" );
	static double pAS = parameters.doubles( "Ig_recuitment" );
	static double dMc = parameters.doubles( "Ig_degradation" );
		
	double lypmh_scale = EPICOUNT / 500000;
	// actual model goes here 
	
	double x[4][13]={{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};//initialize x
	double f[4][13]={{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};//initialize f
	int j;
	
	extern std::vector<int>history;
	
	x[0][0] = (DM+history.back())/lypmh_scale; 
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
	x[0][12] = (DL+(THN+TH1+TH2)/(THN+TH1+TH2+5000)*history.back())/lypmh_scale;
    
	double pT=0.002;
	double pT2=100;
	//(pT/(y(5)+pT2) added
	//rT1->pT1
	//rT2->pT2

    
/*     g*y(6) - (y(8)+y(9)+y(15))/(y(8)+y(9)+y(15)+5000)*dT1*y(7)*y(4)/(y(4)+dT2) % y(7) = Tct
    
    sTh1*y(5)*y(15)*(sig1*(1+y(9))+y(8))/((1+y(9))^2) + pTh1*y(5)*(y(8)^2)/((1+y(9))^2) - dTh1*y(5)*(y(8)^3)/(KTh1+y(9)) - mTh*y(8)  % y(8) = Th1
    
    sTh2*y(5)*y(15)*(sig2+y(9))/(1+y(9)) + pTh2*(ro+y(8))*y(5)*(y(9)^2)/((1+y(9))*(1+y(8)+y(9))) - mTh*y(9)  % y(9) = Th2
    
    g*(y(8)+y(9)) - mTh*y(10) % y(10) = Tct
    
    3*y(5)/(y(5)+10000)*y(17) + rB1*(y(11))*(y(5) + h*y(9))/(rB2 + y(5) + h*y(9)) - pS*y(11) - pL*y(11)*y(9) - db*y(11)   % y(11) = B
    
    pS*y(11) - dS*y(12)  % y(12) = pS
    
    pAS*y(12)+pAL*y(18) - dM*y(13)  - eV*y(3)*y(13) - eV*y(2)*y(13)      %    y(13) = Ig

    dD0*(D0-y(14)) - bD*y(14)*y(3)  %    y(14) = Dn
    
    dHn*(1000-y(15)) - sTh1*y(5)*y(15)*(sig1*(1+y(9))+y(8))/((1+y(9))^2) - sTh2*y(5)*y(15)*(sig2+y(9))/(1+y(9)) %    y(15) = Thn
    
    dC*(Tc0-y(16)) - pT/(y(5)+pT2)*y(16)*y(5)  %    y(16) = Tcn
    
    1E-3*(1000-y(17))-3*y(5)/(y(5)+10000)*y(17) %    y(17) = Bn
    
    pL*y(11)*y(9) - 3E-2*y(18)];%    y(18) = pL */
	
	
	
    for(j = 0; j < 4; j++){
		f[j][0] = {-dDm*x[j][0]}; //DM
        f[j][1] = {pT/(x[j][12]+pT2)*x[j][9]*x[j][12]+rT1*x[j][1]*x[j][12]/(x[j][12]+rT2) - dT1*x[j][1]*x[j][12]/(x[j][12]+dT2)}; //Tc
		f[j][2] = {sTh1*x[j][0]*x[j][10]*(0.1*(1+x[j][3])+x[j][2])/((1+x[j][3])*(1+x[j][3]))+(pTh1*x[j][0]*x[j][2]*x[j][2])/((1+x[j][3])*(1+x[j][3]))-(dTh1*x[j][0]*x[j][2]*x[j][2]*x[j][2])/(500+x[j][3])-mTh*x[j][2]}; //Th1
		f[j][3] = {sTh2*x[j][0]*x[j][10]*(0.1+x[j][3])/(1+x[j][3])+(pTh2*(ro+x[j][2])*x[j][0]*x[j][3]*x[j][3])/((1+x[j][3])*(1+x[j][2]+x[j][3]))-mTh*x[j][3]}; //Th2
		f[j][4] = {CD8_Tcell_recruitment_rate*x[j][1]}; //CD8 export
		f[j][5] = {CD8_Tcell_recruitment_rate*(x[j][2]+x[j][3])}; //CD4 export
		f[j][6] = {0.0021*x[j][0]/(x[j][0]+10000)*x[j][11]+rB1*x[j][6]*(x[j][0]+h*x[j][3])/(x[j][0]+h*x[j][3]+rB2)-dB*x[j][6]-2*pSc*x[j][6]}; //B-Cell
		f[j][7] = {pSc*x[j][6]-dS*x[j][7]}; //P-Cell
		f[j][8] = {pAS*x[j][7]-dMc*x[j][8]}; //Ig
		f[j][9] = {dC*(Tc0-x[j][9]) - pT/(x[j][0]+pT2)*x[j][9]*x[j][0]}; //TcN
		f[j][10] = {dC*(Tc0-x[j][10]) - sTh1*x[j][0]*x[j][10]*(0.1*(1+x[j][3])+x[j][2])/((1+x[j][3])*(1+x[j][3])) - sTh2*x[j][0]*x[j][10]*(0.1+x[j][3])/(1+x[j][3])}; //ThN
		f[j][11] = {1E-3*(Tc0-x[j][11])-0.0021*x[j][0]/(x[j][0]+10000)*x[j][11]}; //bN
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
	
	//std::cout << "Placing " << number_of_Ig << " Ig ... " << std::endl; 
	if( number_of_Ig > 1000)
	{
		number_of_Ig=1000;
	}
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
			double rate = 1.5 * microenvironment(n)[nAb] * microenvironment(n)[nV] * dt; //rate is 1.5 after conversions for now - set to zero in no Ig cases
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
	
