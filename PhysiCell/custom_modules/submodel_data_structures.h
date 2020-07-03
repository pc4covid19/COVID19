#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM;
using namespace PhysiCell;

#ifndef __submodel_data__
#define __submodel_data__

class Submodel_Information
{
 private:
 public:
	std::string name; 
	std::string version;
	
	// for use in the main program call 
	void(*main_function)(double); 
	
	// for use in cell phenotype functions 
	void(*phenotype_function)(Cell*,Phenotype&,double); 
	// for use in cell custom / mechanics functions 
	void(*mechanics_function)(Cell*,Phenotype&,double); 
	
	
	std::vector< std::string > microenvironment_variables; 
	std::vector< std::string > cell_variables; // custom data and parameters  
	
	Submodel_Information( void ); 
	
	void register_model( void ); 
	
	void display( std::ostream& os ); 
}; 

class Submodel_Registry
{
 private:
 public:  
	std::vector<Submodel_Information*> submodels; 
	void register_model( Submodel_Information& model ); 
	void display( std::ostream& os ); 
	
};

void execute_all_submodel_main_functions( double dt );

extern Submodel_Registry submodel_registry; 

#endif