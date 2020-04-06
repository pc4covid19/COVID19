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
	void(*main_function)(Cell*,Phenotype&,double); 
	
	Submodel_Information( void ); 
	
	void display( std::ostream& os ); 
}; 

class Submodel_Registry
{
 private:
	std::vector<Submodel_Information*> submodels; 
 public:  
	void register_model( Submodel_Information& model ); 
	void display( std::ostream& os ); 
	
};

extern Submodel_Registry submodel_registry; 

#endif