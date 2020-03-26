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
 
**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 


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
