# COVID19 tissue simulator 
**Version:** 0.1.0

**Release date:** 26 March 2020 

## Overview

### Key makefile rules:

**make**               : compiles the current project. If no 
                     project has been defined, it first 
                     populates the cancer heterogeneity 2D 
                     sample project and compiles it 
   
**make \[project-name\]**: populates the indicated sample project. 
                     Use "make" to compile it. 

  \[project-name\] choices:
    template2D 
    template3D
    biorobots-sample
    cancer-biorobots-sample
    heterogeneity-sample
    cancer-immune-sample 
    virus-macrophage-sample
 
**make list-projects** : list all available sample projects 

**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

**make reset**         : de-populates the sample project and returns to the original PhysiCell state. Use this when switching to a new PhysiCell sample project. 


**Homepage:**     http://PhysiCell.MathCancer.org

**Downloads:**    http://PhysiCell.sf.net

**Support:**      https://sourceforge.net/p/physicell/tickets/

**Quick Start:**  Look at QuickStart.pdf in the documentation folder. 

**User Guide:**   Look at UserGuide.pdf in the documentation folder. 
 
**Tutorials:**    http://www.mathcancer.org/blog/physicell-tutorials/

**Latest info:**  follow [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)

See changes.md for the full change log. 

* * * 

## Release summary: 

This release fixes minor bugs and improves the documentation. It also adds some minor new capabilities, such as setting time step sizes in the XML configuration file. 

**NOTE:** OSX users must now define PHYSICELL_CPP system variable. See the documentation.
 
### Major new features and changes:

+ List here. 
 
### Minor new features and changes: 
 
+ List here. 
  
### Bugfixes: 

+ List here. 
 
### Notices for intended changes that may affect backwards compatibility:
 
+ List here. 

### Planned future improvements: 
 
+ List here. 

* * * 
