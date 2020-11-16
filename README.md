# COVID19 tissue simulator 
**Version:** 0.4.0

**Release date:** 16 November 2020

## Overview
This model simulates replication dynamics of SARS-CoV-2 (coronavirus / COVID19) in a layer of epithelium and several submodels (such as single-cell response, pyroptosis death model, tissue-damage model, lymph node model and immune response). It is being rapidly prototyped and refined with community support (see below).

This multiscale simulator combines several model components:
* **Tissue:** Virus, cell debris, and chemokines diffuse within the extracellular space. They may also “decay” to reflect removal by interstitial flow into nearby blood vessels or airways.
* **ACE2 receptor dynamics:** Virions bind to ACE2 receptors on the surface, which are internalized (endocytosed) into cells. After virions are released from internalized receptors, they can return to the surface.
* **Viral replication:** Internalized virus is uncoated to expose viral RNA, which synthesizes viral proteins that are assembled into virions. Assembled virions are transported to the cell surface to be exported to the tissue (exocytosed).
* **Single-cell response:** Infected cells secrete a chemokine that may attract immune cells. In a simple pharmacodynamics response (to assembled virions), infected cells can undergo apoptosis. In version `0.4.0`, infected cells also secrete Type-I interferon which induces the inhibition of viral protein synthesis. Apoptosed cells release some or allof their internal contents, notably including virions.
* **Pyroptosis death model:** In version `0.4.0`, we introduced a model of cell death via pyroptosis. Viral RNA levels within the cell act as a DAMP/PAMP that initiates the pyroptosis cascade. The intracellular processes result in the secretion of cytokines IL-1β and IL-18, as well as cell swelling and rupture.
* **Tissue-damage model:** In version `0.4.0`, fibroblast-mediated collagen deposition was included to account for the fibrosis at the damaged sitein response to immune response-induced tissue injury, in which fibroblast cells are recruited into the tissue by following the gradient of anti-inflammatory cytokine and deposit collagen in the place where infected cells are killed by CD8+T cells.
* **Lymph node (LN) model:** In version `0.4.0`, the external immune system was incorporated in the whole model. The presence of virus and infected cells in the tissue induces dendritic cells to activate and egress out of tissue to lymph nodes, where they present antigen to induce activation and proliferation of virus-specific CD4+T cells and CD8+T cells.
* **Immune response:**
  * **Resident (and recruited) macrophages** seek apoptotic cells. They phagocytose (ingest) dead cells upon contact and activate. They also break down ("digest") ingested materials.
  * **Activated macrophages** release a pro-inflammatory cytokine to recruit other immune cells, while seeking both apoptotic and infected cells by chemotaxis. Activated macrophages can “wear out” and apoptose after phagocytosing too much material.
  * **Exhausted macrophages** are active macrophages that internalised debris is above a threshold. They would stop phagocytosing in this stage.
  * **Hyperactivated macrophages** are able to phagocytose infected cells with at least one intracellular viral protein when CD4+T cells induce macrophages to a hyperactive state.
  * **Neutrophils** are recruited by accumulated pro-inflammatory cytokine. They seek apoptotic cells, phagocytose them, and activate. Activated neutrophils seek both apoptotic and infected cells. Neutrophils also capture extracellular virions.
  * **CD8+T cells** are recruited by accumulated pro-inflammatory cytokine. They seek and adhere to infected cells. After sufficient contact time with one or more CD8+T cells, infectedcells undergo apoptosis. They would increase proliferation rate and killing efficacy after activated DCs present antigen tothem.
  * **CD4+T cells** are recruited by accumulated pro-inflammatory cytokine. They apoptose naturally and become dead cells. They are activated in the lymph node by three signals: (1) antigenic presentation by the DCs, (2) direct activation by cytokines secreted by DCs, (3) direct activation by cyto-kines secreted by CD4+T cells.
  * **Dendritic cells** residential exit in the tissue and are activated by infected cells and/or virus. Portion of activated DCs leave the tissue to travel to the lymph node. Activated DCs present antigen to CD8+T cells, which benefits them from improving proliferation rate and killing efficacy.
  * **Fibroblast** are recruited into the tissue by anti-inflammatory cytokine. Fibroblast cells apoptose naturally and become dead cells. They move locally in the tissue along up gradients of anti-inflammatory cytokine and deposit collagen continuously in the location where infected cells are killed by CD8+T cells.

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
