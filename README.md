# COVID19 tissue simulator 
**Version:** 0.3.0

**Release date:** ?? May 2020

## Overview
This model simulates replication dynamics of SARS-CoV-2 (coronavirus / COVID19) in a layer of epithelium. It is being rapidly prototyped and refined with community support (see below).

In this model, SARS-CoV-2 (coronavirus / COVID19) infects a single cell, or a solution of virions is administered to the extracellular space. The virus is uncoated to explose viral RNA, which synthesizes viral proteins that are assembled into a virion. Assembled virions are exported to the environment, where they can diffuse and infect other cells. In the extracellular space, virions adhere to ACE2 receptors and get internalized through endocytosis. Internalized ACE2 receptors release their virus cargo and are recycled back to the surface. 

Macrophages ingest apototic cells and release a pro-inflammatory cytokine that recruits CD8+ T cells and neutrophils. CD8+ T cells chemotax towards cytokines released by infected cells, adhere, and cause damage that can induce apoptosis. Neutrophils also chemotax and can absorb free (extracellular) virus. 

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
This release incorporates major v2 model feedback, ... 




**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.

### New features and changes:
#### 0.3.0:
+ Refactored modular design to include refinements from immune model. 

+ First integration of new immune submodel. 

### Bugfixes 
#### 0.2.0: 
+ None. 

### Notices for intended changes that may affect backwards compatibility:
 
+ None.  

### Planned future improvements: 
 
+ Continue to vet model biology with collaborators. 

+ Add lymph node module for 

+ Add tissue damage models. 

+ Integrate SBML support for submodels.  

* * * 