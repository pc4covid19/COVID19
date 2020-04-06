#include "./submodel_data_structures.h" 

using namespace PhysiCell; 

Submodel_Registry submodel_registry; 

Submodel_Information::Submodel_Information( void )
{
	name = "none";
	version = "-1"; 
	main_function = NULL; 

	return; 
}

void Submodel_Information::display( std::ostream& os )
{
	os << "Submodel: " << name << " (Version " << version << ")" << std::endl 
	<< "\tfunction: " ; 
	if( main_function )
	{ os << (long long int) main_function; }
	else
	{ os << "NULL"; }
	os << std::endl; 
	
	return; 
}

void Submodel_Registry::register_model( Submodel_Information& model )
{
	#pragma omp critical 
	{
		// already registered? 
		bool found = false; 
		for( int n = 0; n < submodels.size() ; n++ )
		{
			if( submodels[n] == &model )
			{ found = true; } 
		} 
		
		if( found == false )
		{
			submodels.push_back( &model ); 
//			add_software_citation(); 
// void add_software_citation( std::string name , std::string version, std::string DOI , std::string URL )
		}
	}
	
	return;
}

void Submodel_Registry::display( std::ostream& os )
{
	os << "The following submodels are registered: " << std::endl; 
	os << "=======================================" << std::endl; 
	for( int n = 0 ; n < submodels.size(); n++ )
	{
		submodels[n]->display( os ); 
	}
	os << std::endl; 
	
	return;
}

	

