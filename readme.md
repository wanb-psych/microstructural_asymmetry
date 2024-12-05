# Microstructural asymmetry in the human cortex

This repo is for documentation of the work "Microstructural asymmetry in the human cortex"
---
Paper: https://www.nature.com/articles/s41467-024-54243-9
It includes 5 visualization ipython notebooks, 2 built-in python scripts, and 1 bash script. 

---
### $ BASH
`run_solar.sh`  
This script helps compute twin-based heritability using SOLAR (https://solar-eclipse-genetics.org/) 

### $ built-in python scripts
`func_plot.py` and `func_stats.py` help plot brain surfaces and bar charts, and statistics about peemutation using variogram and cohen's d, FDR correction, etc.

### $ ipython notebooks
- **Vis_BigBrain.ipynb** includes BigBrain main and supplementary figures: *post-mortem* cytoarchitectural asymmetry

- **Vis_HCP.ipynb** includes main and supplementary results and figures using Human Connectome Project (HCP) data: *in vivo* T1w/T2w asymmetry, sex and age effects, heritability, and relation to functional connectome

- **Vis_CCA.ipynb** includes behavioral relevence results and figures using canonical correlation analysis (CCA) in HCP: language and mental health association

- **Vis_MICs.ipynb** includes the supplementary figure using MICA-MICs sample (N=50): *in vivo* quantitative T1 asymmetry

- **Vis_NSPN.ipynb** includes the supplementary figure using NSPN sample (N=286): *in vivo* magnetization transfer (MT) asymmetry

---
### Key toolboxes
- Micapipe (https://micapipe.readthedocs.io/en/latest/)
- BigBrainWrap (https://bigbrainwarp.readthedocs.io/en/latest/)
- Brainspace (https://brainspace.readthedocs.io/en/latest/index.html)
- Brainstat (https://brainstat.readthedocs.io/en/master/)
- BrainSMASH (https://brainsmash.readthedocs.io/en/latest/)
- SOLAR (https://solar-eclipse-genetics.org/)
- Scikit-learn (https://scikit-learn.org/stable/)


---

### Correspondence
Bin Wan (binwan.academic@gmail.com)  
Sofie L. Valk (valk@cbs.mpg.de)  

Cognitive Neurogenetics (CNG) lab,   
Max Planck Institute for Human Cognitive and Brain Sciences &  
Institute of Neuroscience and Medicine (INM-7), Research Centre Jülich 
