# COVID19 tissue simulator
**Version:** 0.6.0

**Release date:** 30 January 2023

## Overview
This model simulates replication dynamics of SARS-CoV-2 (coronavirus / COVID19) in a layer of epithelium with an initial immune reaction. It is being rapidly prototyped and refined with community support (see below).

In this model, SARS-CoV-2 (coronavirus / COVID19) infects a single cell, or a solution of virions is administered to the extracellular space. The virus is uncoated to explose viral RNA, which synthesizes viral proteins that are assembled into a virion. Assembled virions are exported to the environment, where they can diffuse and infect other cells. In the extracellular space, virions adhere to ACE2 receptors and get internalized through endocytosis. Internalized ACE2 receptors release their virus cargo and are recycled back to the surface. 

Resident macrophages ingest apototic cells and release a pro-inflammatory cytokine that recruits additional macrophages, neutrophils, and CD8+ T cells. CD8+ T cells chemotax towards cytokines released by infected cells and adhere. Cumulative CD8+ T cell contact time can induce apoptosis in infected cells. Activated macrophages and neutrophils chemotaxis chemotax along chemokine and debris gradients and continue to phagocytose dead cells. Neutrophils also absorb free (extracellular) virus. Immunoglobulin can bind and remove Virion in the external field and bind to infected Cells.

The model includes a basic pharmacodynamic response (to assembled virions) to cause cell apoptosis. Apoptosed cells release some or all of their internal contents, notably including virions.

### Caveats and disclaimers: 
**This model is under active development**
* It has not been peer reviewed. 
* It is intended to drive basic scientific research and public education at this stage. 
* **It cannot be used for public policy decisions.**
* **It cannot be used for individual medical decisions.**

### Key makefile rules:

**make**               : compiles the project.
 
**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

**make reset**         : reset to default settings (restores config file)

### More references 

**Preprint:**      https://doi.org/10.1101/2020.04.02.019075 

**Model details:** https://github.com/MathCancer/COVID19/wiki/About 

**Homepage:**      http://covid19.PhysiCell.org

**Support:**       https://sourceforge.net/p/physicell/tickets/

**Latest info:**   follow [@PhysiCell](https://twitter.com/PhysiCell) and [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)

* * * 

## Release summary: 
### 0.6.0:
This release incorporates major v4 model feedback and adds: 
+ 

**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.

### New features and changes:
#### 0.6.0:

+ Removed multiple hard coded values to xml including lymph node values, the standard deviation of initial viral placement, and mechanics voxel size

+ Changed macrophage activation so it needs to find viral RNA in the dead cell to activate.

+ Fixed viral uptake to observed rates, removed cell volume from the calculation

+ Added function to pre dose the system with Ig in the xml including time delays

+ Moved immune_dt to the phenotype timescale

+ Updated to PhysiCell 1.10.2

### Bugfixes 
#### 0.6.0: 

+ fixed a division by 0 that can rarely occur when multithreading

+ multiple small fixes to xml

+ fix to rare crash when a out of bounds cell interacts

### Notices for intended changes that may affect backwards compatibility:
 
+ added params that are needed for the proper functioning of the system 

+ removed unused params that are not in the current version

+ changed how ACE2 binding is calculated

### Planned future improvements: 
 
+ Continue to vet model biology with collaborators.

* * * 
**Version:** 0.5.0

**Release date:** 21 July 2021

## Overview
This model simulates replication dynamics of SARS-CoV-2 (coronavirus / COVID19) in a layer of epithelium with an initial immune reaction. It is being rapidly prototyped and refined with community support (see below).

In this model, SARS-CoV-2 (coronavirus / COVID19) infects a single cell, or a solution of virions is administered to the extracellular space. The virus is uncoated to explose viral RNA, which synthesizes viral proteins that are assembled into a virion. Assembled virions are exported to the environment, where they can diffuse and infect other cells. In the extracellular space, virions adhere to ACE2 receptors and get internalized through endocytosis. Internalized ACE2 receptors release their virus cargo and are recycled back to the surface. 

Resident macrophages ingest apototic cells and release a pro-inflammatory cytokine that recruits additional macrophages, neutrophils, and CD8+ T cells. CD8+ T cells chemotax towards cytokines released by infected cells and adhere. Cumulative CD8+ T cell contact time can induce apoptosis in infected cells. Activated macrophages and neutrophils chemotaxis chemotax along chemokine and debris gradients and continue to phagocytose dead cells. Neutrophils also absorb free (extracellular) virus. Immunoglobulin can bind and remove Virion in the external field and bind to infected Cells.

The model includes a basic pharmacodynamic response (to assembled virions) to cause cell apoptosis. Apoptosed cells release some or all of their internal contents, notably including virions.

### Caveats and disclaimers: 
**This model is under active development using rapid prototyping:**
* It has not been peer reviewed. 
* It is intended to drive basic scientific research and public education at this stage. 
* **It cannot be used for public policy decisions.**
* **It cannot be used for individual medical decisions.**

**This model will be continually refined with input from the community, particularly experts in infectious diseases. The validation state will be updated as this progresses.**

### Key makefile rules:

**make**               : compiles the project.
 
**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

**make reset**         : reset to default settings (restores config file)

### More references 

**Preprint:**      https://doi.org/10.1101/2020.04.02.019075 

**Model details:** https://github.com/MathCancer/COVID19/wiki/About 

**Homepage:**      http://covid19.PhysiCell.org

**Support:**       https://sourceforge.net/p/physicell/tickets/

**Latest info:**   follow [@PhysiCell](https://twitter.com/PhysiCell) and [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)

See changes.md for the full change log. 

* * * 

## Release summary: 
### 0.5.0:
This release incorporates major v4 model feedback and adds: 
+ An improved systemtic immune model (v0.2.0) that adds B-Cells, Plasma Cells, and Ig. Time delays can be set for transport of Cells.

+ An updated immune submodel (immune_submodels v0.3.0) that includes more macrophage states, ROS secretion, and T Cell generation changes.

+ An improved fibroblast/ collagen model with improved dynamics and finite anti-inflammatory secretion. 

**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.

### New features and changes:
#### 0.5.0:

+ Lymph node delay and history to simulate transport to lymph node.

+ ROS induced apoptosis; neutrophils secrete ROS on phagocytosis

+ fibroblast collagen changes, anti-inflammatory secretion after CD8+ T Cell induced death is now temporary 

+ BCell,PCell,Ig addition; Ig can travel to tissue to both bind to cells and remove virion in external environment

+ CD8+ and CD4+ T Cell tissue phenotype changes; generation counter kills CD8s faster once "generation" passes

+ CD8 contact to active macrophage turns off proinflammatory and turns on antiinflammatory response.

+ General parameter changes.

### Bugfixes 
#### 0.5.0: 
+ None. 

### Notices for intended changes that may affect backwards compatibility:
 
+ None.  

### Planned future improvements: 
 
+ Continue to vet model biology with collaborators. 

+ Improve lymph node module time scales to better reflect bilogical data. (v6) 

+ Integrate SBML support for submodels (v6)

+ Refine viral binding model (v6)

+ Treat virion as an agent in the tissue space (v6)

+ Refine immune model (including more cell types and improved parameter estimates) (in process v5-v6)

* * * 
**Version:** 0.4.0

**Release date:** 30 September 2020

## Overview
This model simulates replication dynamics of SARS-CoV-2 (coronavirus / COVID19) in a layer of epithelium with an initial immune reaction. It is being rapidly prototyped and refined with community support (see below).

In this model, SARS-CoV-2 (coronavirus / COVID19) infects a single cell, or a solution of virions is administered to the extracellular space. The virus is uncoated to explose viral RNA, which synthesizes viral proteins that are assembled into a virion. Assembled virions are exported to the environment, where they can diffuse and infect other cells. In the extracellular space, virions adhere to ACE2 receptors and get internalized through endocytosis. Internalized ACE2 receptors release their virus cargo and are recycled back to the surface. 

Resident macrophages ingest apototic cells and release a pro-inflammatory cytokine that recruits additional macrophages, neutrophils, and CD8+ T cells. CD8+ T cells chemotax towards cytokines released by infected cells and adhere. Cumulative CD8+ T cell contact time can induce apoptosis in infected cells. Activated macrophages and neutrophils chemotaxis chemotax along chemokine and debris gradients and continue to phagocytose dead cells. Neutrophils also absorb free (extracellular) virus. 

The model includes a basic pharmacodynamic response (to assembled virions) to cause cell apoptosis. Apoptosed cells release some or all of their internal contents, notably including virions.

### Caveats and disclaimers: 
**This model is under active development using rapid prototyping:**
* It has not been peer reviewed. 
* It is intended to drive basic scientific research and public education at this stage. 
* **It cannot be used for public policy decisions.**
* **It cannot be used for individual medical decisions.**

**This model will be continually refined with input from the community, particularly experts in infectious diseases. The validation state will be updated as this progresses.**

### Key makefile rules:

**make**               : compiles the project.
 
**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

**make reset**         : reset to default settings (restores config file)

### More references 

**Preprint:**      https://doi.org/10.1101/2020.04.02.019075 

**Model details:** https://github.com/MathCancer/COVID19/wiki/About 

**Homepage:**      http://covid19.PhysiCell.org

**Support:**       https://sourceforge.net/p/physicell/tickets/

**Latest info:**   follow [@PhysiCell](https://twitter.com/PhysiCell) and [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)

See changes.md for the full change log. 

* * * 

## Release summary: 
### 0.4.0:
This release incorporates major v3 model feedback and adds: 
+ The first systemtic immune model (v0.1.0) that models immune cell expansion and trafficking to/from the local tissue model. 

+ An updated immune submodel (immune_submodels v0.2.0) that includes more macrophage states, dendritic cells, and immune cell trafficking to/from a systemic immune model 

+ An improved virus-ACE2 receptor binding model (receptor_dynamics v0.3.0) that disallows partial cell infection (only integer virus-receptor binding pairs are allowed). 

+ The first models of interferon signaling (internal_viral_response v0.3.0)

+ The first detailed model of infected cell pyroptosis (internal_viral_response v0.3.0) 

**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.

### New features and changes:
#### 0.4.0:
+ Updated immune submodel to v0.2.0

+ Added systemic immune model v0.1.0 to include Dendritic Cell (DC) dependent immune reponse in the lymph node. DC can leave the spatial system to enter an ODE system (lymph node), to inform the recuitment of CD4/CD8.

+ Updated infected cell responses to v0.3.0 to include interferon signaling: infected cells secrete interferon. Any epithelial cell secretes interferon if it senses interferon. Reading interferon causes cells to inhibit viral protein synthesis 

+ SVG plots now have user-specified opacity for epithelial and non-epithelial cells. 



### Bugfixes 
#### 0.4.0: 
+ None. 

### Notices for intended changes that may affect backwards compatibility:
 
+ None.  

### Planned future improvements: 
 
+ Continue to vet model biology with collaborators. 

+ Add lymph node module. (in process v4) 

+ Add tissue damage models. (in process v4-v5)

+ Integrate SBML support for submodels (v5?)

+ Refine viral replication model (v5?) 

+ Refine immune model (including more cell types and improved parameter estimates) (in process v4)

+ Add interferon response model (in process v4-v5)

* * * 
**Version:** 0.3.2

**Release date:** 15 July 2020

## Overview
This model simulates replication dynamics of SARS-CoV-2 (coronavirus / COVID19) in a layer of epithelium with an initial immune reaction. It is being rapidly prototyped and refined with community support (see below).

In this model, SARS-CoV-2 (coronavirus / COVID19) infects a single cell, or a solution of virions is administered to the extracellular space. The virus is uncoated to explose viral RNA, which synthesizes viral proteins that are assembled into a virion. Assembled virions are exported to the environment, where they can diffuse and infect other cells. In the extracellular space, virions adhere to ACE2 receptors and get internalized through endocytosis. Internalized ACE2 receptors release their virus cargo and are recycled back to the surface. 

Resident macrophages ingest apototic cells and release a pro-inflammatory cytokine that recruits additional macrophages, neutrophils, and CD8+ T cells. CD8+ T cells chemotax towards cytokines released by infected cells and adhere. Cumulative CD8+ T cell contact time can induce apoptosis in infectd cells. Activated macrophages and neutrophils chemotaxis chemotax along chemokine and debris gradients and continue to phagocytose dead cells. Neutrophils also absorb free (extracellular) virus. 

The model includes a basic pharmacodynamic response (to assembled virions) to cause cell apoptosis. Apoptosed cells release some or all of their internal contents, notably including virions.

### Caveats and disclaimers: 
**This model is under active development using rapid prototyping:**
* It has not been peer reviewed. 
* It is intended to drive basic scientific research and public education at this stage. 
* **It cannot be used for public policy decisions.**
* **It cannot be used for individual medical decisions.**

**This model will be continually refined with input from the community, particularly experts in infectious diseases. The validation state will be updated as this progresses.**

### Key makefile rules:

**make**               : compiles the project.
 
**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

**make reset**         : reset to default settings (restores config file)

### More references 

**Preprint:**      https://doi.org/10.1101/2020.04.02.019075 

**Model details:** https://github.com/MathCancer/COVID19/wiki/About 

**Homepage:**      http://covid19.PhysiCell.org

**Support:**       https://sourceforge.net/p/physicell/tickets/

**Latest info:**   follow [@PhysiCell](https://twitter.com/PhysiCell) and [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)

See changes.md for the full change log. 

* * * 

## Release summary: 
### 0.3.2:
This release removes rules from the macrophages where (1) they cannot phagocytose over a certain size and (2) they apoptose if exceeding a certain size. These rules will require additional development in v4. 

### 0.3.1:
This release improves parameter estimates for digestion of phagocytosed material and has an immune model refinement to prevent runaway macrophage death. 

### 0.3.0:
This release incorporates major v2 model feedback and adds the first immune submodel. 

**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.

### New features and changes:
#### 0.3.2: 
+ Removed macrophage rules to (1) check size of target cells prior to phagocytosis and (2) trigger apoptosis if exceeding a tolerance. (Upgrade from immune submodel 0.1.1 to 0.1.2.)

* Added "make reset" rule to restore the default configuration file (from PhysiCell_settings.xml-default)

#### 0.3.1: 
+ Refined macrophage and neutrophil models and parameters for phagocytosis. (Upgrade from immune submodel 0.1.0 to 0.1.1.) 

#### 0.3.0:
+ Refactored modular design to include refinements from immune model. 

+ First integration of new immune submodel. 

+ Upgrade to PhysiCell Version 1.7.1, allowing use of XML-based cell definitions to define the behavior of immune cell types. 

+ Upgrade to PhysiCell Version 1.7.2beta to improve multithreaded performance, add new cell-cell interaction features, and fix concurrency issues on some platforms. 

### Bugfixes 
#### 0.3.0: 
+ None. 

### Notices for intended changes that may affect backwards compatibility:
 
+ None.  

### Planned future improvements: 
 
+ Continue to vet model biology with collaborators. 

+ Add lymph node module. 

+ Add tissue damage models. 

+ Integrate SBML support for submodels.  

+ Refine viral replication model. 

+ Refine immune model (including more cell types and improved parameter estimates).

+ Add interferon response model. 

* * * 

# COVID19 tissue simulator 
**Version:** 0.3.1

**Release date:** 3 July 2020

## Overview
This model simulates replication dynamics of SARS-CoV-2 (coronavirus / COVID19) in a layer of epithelium with an initial immune reaction. It is being rapidly prototyped and refined with community support (see below).

In this model, SARS-CoV-2 (coronavirus / COVID19) infects a single cell, or a solution of virions is administered to the extracellular space. The virus is uncoated to explose viral RNA, which synthesizes viral proteins that are assembled into a virion. Assembled virions are exported to the environment, where they can diffuse and infect other cells. In the extracellular space, virions adhere to ACE2 receptors and get internalized through endocytosis. Internalized ACE2 receptors release their virus cargo and are recycled back to the surface. 

Resident macrophages ingest apototic cells and release a pro-inflammatory cytokine that recruits additional macrophages, neutrophils, and CD8+ T cells. CD8+ T cells chemotax towards cytokines released by infected cells and adhere. Cumulative CD8+ T cell contact time can induce apoptosis in infectd cells. Activated macrophages and neutrophils chemotaxis chemotax along chemokine and debris gradients and continue to phagocytose dead cells. Neutrophils also absorb free (extracellular) virus. 

The model includes a basic pharmacodynamic response (to assembled virions) to cause cell apoptosis. Apoptosed cells release some or all of their internal contents, notably including virions.

### Caveats and disclaimers: 
**This model is under active development using rapid prototyping:**
* It has not been peer reviewed. 
* It is intended to drive basic scientific research and public education at this stage. 
* **It cannot be used for public policy decisions.**
* **It cannot be used for individual medical decisions.**

**This model will be continually refined with input from the community, particularly experts in infectious diseases. The validation state will be updated as this progresses.**

### Key makefile rules:

**make**               : compiles the project.
 
**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

### More references 

**Preprint:**      https://doi.org/10.1101/2020.04.02.019075 

**Model details:** https://github.com/MathCancer/COVID19/wiki/About 

**Homepage:**      http://covid19.PhysiCell.org

**Support:**       https://sourceforge.net/p/physicell/tickets/

**Latest info:**   follow [@PhysiCell](https://twitter.com/PhysiCell) and [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)

See changes.md for the full change log. 

* * * 

## Release summary: 
### 0.3.1:
This release improves parameter estimates for digestion of phagocytosed material and has an immune model refinement to prevent runaway macrophage death. 

### 0.3.0:
This release incorporates major v2 model feedback and adds the first immune submodel. 

**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.

### New features and changes:
#### 0.3.1: 
+ Refined macrophage and neutrophil models and parameters for phagocytosis. (Upgrade from immune submodel 0.1.0 to 0.1.1.) 

#### 0.3.0:
+ Refactored modular design to include refinements from immune model. 

+ First integration of new immune submodel. 

+ Upgrade to PhysiCell Version 1.7.1, allowing use of XML-based cell definitions to define the behavior of immune cell types. 

+ Upgrade to PhysiCell Version 1.7.2beta to improve multithreaded performance, add new cell-cell interaction features, and fix concurrency issues on some platforms. 

### Bugfixes 
#### 0.3.0: 
+ None. 

### Notices for intended changes that may affect backwards compatibility:
 
+ None.  

### Planned future improvements: 
 
+ Continue to vet model biology with collaborators. 

+ Add lymph node module. 

+ Add tissue damage models. 

+ Integrate SBML support for submodels.  

+ Refine viral replication model. 

+ Refine immune model (including more cell types and improved parameter estimates).

+ Add interferon response model. 

* * * 

# COVID19 tissue simulator 
**Version:** 0.3.0

**Release date:** 3 July 2020

## Overview
This model simulates replication dynamics of SARS-CoV-2 (coronavirus / COVID19) in a layer of epithelium with an initial immune reaction. It is being rapidly prototyped and refined with community support (see below).

In this model, SARS-CoV-2 (coronavirus / COVID19) infects a single cell, or a solution of virions is administered to the extracellular space. The virus is uncoated to explose viral RNA, which synthesizes viral proteins that are assembled into a virion. Assembled virions are exported to the environment, where they can diffuse and infect other cells. In the extracellular space, virions adhere to ACE2 receptors and get internalized through endocytosis. Internalized ACE2 receptors release their virus cargo and are recycled back to the surface. 

Resident macrophages ingest apototic cells and release a pro-inflammatory cytokine that recruits additional macrophages, neutrophils, and CD8+ T cells. CD8+ T cells chemotax towards cytokines released by infected cells and adhere. Cumulative CD8+ T cell contact time can induce apoptosis in infectd cells. Activated macrophages and neutrophils chemotaxis chemotax along chemokine and debris gradients and continue to phagocytose dead cells. Neutrophils also absorb free (extracellular) virus. 

The model includes a basic pharmacodynamic response (to assembled virions) to cause cell apoptosis. Apoptosed cells release some or all of their internal contents, notably including virions.

### Caveats and disclaimers: 
**This model is under active development using rapid prototyping:**
* It has not been peer reviewed. 
* It is intended to drive basic scientific research and public education at this stage. 
* **It cannot be used for public policy decisions.**
* **It cannot be used for individual medical decisions.**

**This model will be continually refined with input from the community, particularly experts in infectious diseases. The validation state will be updated as this progresses.**

### Key makefile rules:

**make**               : compiles the project.
 
**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

### More references 

**Preprint:**      https://doi.org/10.1101/2020.04.02.019075 

**Model details:** https://github.com/MathCancer/COVID19/wiki/About 

**Homepage:**      http://covid19.PhysiCell.org

**Support:**       https://sourceforge.net/p/physicell/tickets/

**Latest info:**   follow [@PhysiCell](https://twitter.com/PhysiCell) and [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)

See changes.md for the full change log. 

* * * 

## Release summary: 
### 0.3.0:
This release incorporates major v2 model feedback and adds the first immune submoel. 




**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.

### New features and changes:
#### 0.3.0:
+ Refactored modular design to include refinements from immune model. 

+ First integration of new immune submodel. 

+ Upgrade to PhysiCell Version 1.7.1, allowing use of XML-based cell definitions to define the behavior of immune cell types. 

+ Upgrade to PhysiCell Version 1.7.2beta to improve multithreaded performance, add new cell-cell interaction features, and fix concurrency issues on some platforms. 

### Bugfixes 
#### 0.2.0: 
+ None. 

### Notices for intended changes that may affect backwards compatibility:
 
+ None.  

### Planned future improvements: 
 
+ Continue to vet model biology with collaborators. 

+ Add lymph node module. 

+ Add tissue damage models. 

+ Integrate SBML support for submodels.  

+ Refine viral replication model. 

+ Refine immune model (including more cell types and improved parameter estimates).

+ Add interferon response model. 

* * * 

# COVID19 tissue simulator 
**Version:** 0.2.1

**Release date:** 10 April 2020

## Overview
This model simulates replication dynamics of SARS-CoV-2 (coronavirus / COVID19) in a layer of epithelium. It is being rapidly prototyped and refined with community support (see below).

In this model, SARS-CoV-2 (coronavirus / COVID19) infects a single cell, or a solution of virions is administered to the extracellular space. The virus is uncoated to explose viral RNA, which synthesizes viral proteins that are assembled into a virion. Assembled virions are exported to the environment, where they can diffuse and infect other cells. In the extracellular space, virions adhere to ACE2 receptors and get internalized through endocytosis. Internalized ACE2 receptors release their virus cargo and are recycled back to the surface. 

The model includes a basic pharmacodynamic response (to assembled virions) to cause cell apoptosis. Apoptosed cells release some or all of their internal contents, notably including virions.

### Caveats and disclaimers: 
**This model is under active development using rapid prototyping:**
* It has not been peer reviewed. 
* It is intended to drive basic scientific research and public education at this stage. 
* **It cannot be used for public policy decisions.**
* **It cannot be used for individual medical decisions.**

**This model will be continually refined with input from the community, particularly experts in infectious diseases. The validation state will be updated as this progresses.**

### Key makefile rules:

**make**               : compiles the project.
 
**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

### More references 

**Preprint:**      https://doi.org/10.1101/2020.04.02.019075 

**Model details:** https://github.com/MathCancer/COVID19/wiki/About 

**Homepage:**      http://covid19.PhysiCell.org

**Support:**       https://sourceforge.net/p/physicell/tickets/

**Latest info:**   follow [@PhysiCell](https://twitter.com/PhysiCell) and [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)

See changes.md for the full change log. 

* * * 

## Release summary: 
### 0.2.1:
This is a bug release fix, which fixes a rendering bug where all cells were visualized as black. 

### 0.2.0:
This release incorporates major v1 model feedback, particularly a refactoring into a more modular architecture with submodels, a placeholder ACE2 receptor traffickign model, and receptor-modulated endocytosis. 

**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.

### New features and changes:
#### 0.2.1:
+ None. 

#### 0.2.0:
+ Refactored into modular design based on v1 preprint feedback. 

+ Set default max time to 7200 minutes

+ Set default diffusion coefficient to 90 micron^2/min based on v1 feedback 

+ Added a basic functionality to "register" all submodels into a list with basic information and automatically creation of custom cell variables. 

+ Cells now automatically record their internal virus variables in output data, in response to v1 feedback. 

+ Added ACE2 receptor trafficking model based on v1 feedback. This submodel is now responsible for delivering virions to the cell cytoplasm. Notably, ACE2 receptor availability modulates the virus uptake rate. 

+ Added ability to specify the MOI (multiplicity of infection) at the simulation start, in response to v1 feedback. 

+ Simplified the diffusing fields. 

+ Worked to improve parameter estimates based on ACE2 papers. 

### Bugfixes 
#### 0.2.1: 

#### 0.2.0: 
+ A single cell is infected at the center of the domain, rather than a hard-coded test that fails for some domain sizes. 

### Notices for intended changes that may affect backwards compatibility:
 
+ None.  

### Planned future improvements: 
 
+ Improved parameter estimates. 

+ Continue to vet model biology with collaborators. 

+ Add inflammatory response, and potentially link to ARDS. 

+ Add tissue damage models.  

* * * 

# COVID19 tissue simulator 
**Version:** 0.2.0

**Release date:** 9 April 2020

## Overview
This model simulates replication dynamics of SARS-CoV-2 (coronavirus / COVID19) in a layer of epithelium. It is being rapidly prototyped and refined with community support (see below).

In this model, SARS-CoV-2 (coronavirus / COVID19) infects a single cell, or a solution of virions is administered to the extracellular space. The virus is uncoated to explose viral RNA, which synthesizes viral proteins that are assembled into a virion. Assembled virions are exported to the environment, where they can diffuse and infect other cells. In the extracellular space, virions adhere to ACE2 receptors and get internalized through endocytosis. Internalized ACE2 receptors release their virus cargo and are recycled back to the surface. 

The model includes a basic pharmacodynamic response (to assembled virions) to cause cell apoptosis. Apoptosed cells release some or all of their internal contents, notably including virions.

### Caveats and disclaimers: 
**This model is under active development using rapid prototyping:**
* It has not been peer reviewed. 
* It is intended to drive basic scientific research and public education at this stage. 
* **It cannot be used for public policy decisions.**
* **It cannot be used for individual medical decisions.**

**This model will be continually refined with input from the community, particularly experts in infectious diseases. The validation state will be updated as this progresses.**

### Key makefile rules:

**make**               : compiles the project.
 
**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

### More references 

**Preprint:**      https://doi.org/10.1101/2020.04.02.019075 

**Model details:** https://github.com/MathCancer/COVID19/wiki/About 

**Homepage:**      http://covid19.PhysiCell.org

**Support:**       https://sourceforge.net/p/physicell/tickets/

**Latest info:**   follow [@PhysiCell](https://twitter.com/PhysiCell) and [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)

See changes.md for the full change log. 

* * * 

## Release summary: 
**This release is part of the v2 prototyping iteration.**

This release incorporates major v1 model feedback, particularly a refactoring into a more modular architecture with submodels, a placeholder ACE2 receptor trafficking model, and receptor-modulated endocytosis. 

**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.

### New features and changes:

+ Refactored into modular design based on v1 preprint feedback. 

+ Set default max time to 7200 minutes

+ Set default diffusion coefficient to 90 micron^2/min based on v1 feedback 

+ Added a basic functionality to "register" all submodels into a list with basic information and automatically creation of custom cell variables. 

+ Cells now automatically record their internal virus variables in output data. 

+ Added ACE2 receptor trafficking model based on v1 feedback. This submodel is now responsible for delivering virions to the cell cytoplasm.

+ Added ability to specify the MOI (multiplicity of infection) at the simulation start based on v1 feedback. 

+ Simplified the diffusing fields. 

+ Worked to improve parameter estimates based on ACE2 papers. 

+ Cell surface is colored baced on bound ACE2 receptor elevel. 

### Bugfixes 

+ A single cell is infected at the center of the domain, rather than a hard-coded test that fails for some domain sizes. 

### Notices for intended changes that may affect backwards compatibility:
 
+ None.  

### Planned future improvements: 
 
+ Improved parameter estimates. 

+ Continue to vet model biology with collaborators. 

+ Add inflammatory response, and potentially link to ARDS. 

+ Add tissue damage models.  

* * * 

# COVID19 tissue simulator 
**Version:** 0.1.1, 0.1.2, 0.1.3 

**Release date:** 26 March 2020 

## Overview
This model simulates replication dynamics of SARS-CoV-2 (coronavirus / COVID19) in a layer of epithelium. It is being rapidly prototyped and refined with community support (see below).

In this model, SARS-CoV-2 (coronavirus / COVID19) infects a single cell. The virus is uncoated to explose viral RNA, which synthesizes viral proteins that are assembled into a virion. Assembled virions are exported to the environment. Virions diffuse and can infect other cells (including the original cell).

The model includes a basic pharmacodynamic response (to assembled virions) to cause cell apoptosis. Apoptosed cells release some or all of their internal contents, notably including virions.

Released 0.1.2 and 0.1.3 have no changes to 0.1.1, other than enabling Zenodo snapshots. 

### Caveats and disclaimers: 
**This model is under active development using rapid prototyping:**
* It has not been peer reviewed. 
* It is intended to drive basic scientific research and public education at this stage. 
* **It cannot be used for public policy decisions.**
* **It cannot be used for individual medical decisions.**

**This model will be continually refined with input from the community, particularly experts in infectious diseases. The validation state will be updated as this progresses.**

### Key makefile rules:

**make**               : compiles the project.
 
**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

### More references 

**Model details:** https://github.com/MathCancer/COVID19/wiki/About 

**Homepage:**     http://PhysiCell.org

**Downloads:**    http://PhysiCell.sf.net

**Support:**      https://sourceforge.net/p/physicell/tickets/

**Quick Start:**  Look at QuickStart.pdf in the documentation folder. 

**User Guide:**   Look at UserGuide.pdf in the documentation folder. 
 
**Tutorials:**    http://www.mathcancer.org/blog/physicell-tutorials/

**Latest info:**  follow [@PhysiCell](https://twitter.com/PhysiCell) and [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)

See changes.md for the full change log. 

* * * 

## Release summary: 

This releases first documentation on the math. See the /math directory. 

**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.
 
### Major new features and changes:

+ Math and equations first documented. 

+ More information on the model assumptions will be updated and refined at https://github.com/MathCancer/COVID19/wiki/About 

### Minor new features and changes: 
 
+ Refined parameter file for better xml2jupyter compatibility.  

### Bugfixes: 

+ None. 
 
### Notices for intended changes that may affect backwards compatibility:
 
+ None.  

### Planned future improvements: 
 
+ Improved parameter estimates. 

+ Vet model biology with collaborators. 

+ Add molecular-scale model of virus-cell interactions via ACE2, and potentially receptor trafficking that modulates virus uptake rate. 

+ Add inflammatory response, and potentially link to ARDS. 

* * * 

# COVID19 tissue simulator 

**Version:** 0.1.0

**Release date:** 26 March 2020 

## Overview
This model simulates replication dynamics of SARS-CoV-2 (coronavirus / COVID19) in a layer of epithelium. It is being rapidly prototyped and refined with community support (see below).

In this model, SARS-CoV-2 (coronavirus / COVID19) infects a single cell. The virus is uncoated to explose viral RNA, which synthesizes viral proteins that are assembled into a virion. Assembled virions are exported to the environment. Virions diffuse and can infect other cells (including the original cell).

The model includes a basic pharmacodynamic response (to assembled virions) to cause cell apoptosis. Apoptosed cells release some or all of their internal contents, notably including virions.


### Caveats and disclaimers: 
**This model is under active development using rapid prototyping:**
* It has not been peer reviewed. 
* It is intended to drive basic scientific research and public education at this stage. 
* **It cannot be used for public policy decisions.**
* **It cannot be used for individual medical decisions.**

**This model will be continually refined with input from the community, particularly experts in infectious diseases. The validation state will be updated as this progresses.**

### Key makefile rules:

**make**               : compiles the project.
 
**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

### More references 

**Model details:** https://github.com/MathCancer/COVID19/wiki/About 

**Homepage:**     http://PhysiCell.org

**Downloads:**    http://PhysiCell.sf.net

**Support:**      https://sourceforge.net/p/physicell/tickets/

**Quick Start:**  Look at QuickStart.pdf in the documentation folder. 

**User Guide:**   Look at UserGuide.pdf in the documentation folder. 
 
**Tutorials:**    http://www.mathcancer.org/blog/physicell-tutorials/

**Latest info:**  follow [@PhysiCell](https://twitter.com/PhysiCell) and [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)

See changes.md for the full change log. 

* * * 

## Release summary: 

This is the initial release, based on rapid prototyping on March 25, 2020.

**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.
 
### Major new features and changes:

+ This initial release includes RNA virus replication in a single layer of (lung) epithelium. 

+ More information on the model assumptions will be updated and refined at https://github.com/MathCancer/COVID19/wiki/About 

### Minor new features and changes: 
 
+ First release. 

### Bugfixes: 

+ None. 
 
### Notices for intended changes that may affect backwards compatibility:
 
+ None.  

### Planned future improvements: 
 
+ Improved parameter estimates. 

+ Vet model biology with collaborators. 

* * * 
