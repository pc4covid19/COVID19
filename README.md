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

This release incorporates major v1 model feedback, particularly a refactoring into a more modular architecture with submodels, a placeholder ACE2 receptor traffickign model, and receptor-modulated endocytosis. 

**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.

### New features and changes:

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

+ A single cell is infected at the center of the domain, rather than a hard-coded test that fails for some domain sizes. 

### Notices for intended changes that may affect backwards compatibility:
 
+ None.  

### Planned future improvements: 
 
+ Improved parameter estimates. 

+ Continue to vet model biology with collaborators. 

+ Add inflammatory response, and potentially link to ARDS. 

+ Add tissue damage models.  

* * * 