# COVID19 tissue simulator 
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
