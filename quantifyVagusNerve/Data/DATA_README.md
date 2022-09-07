# Data 
The directory **spatial-neuro/quantifyVagusNerve/Data/** contains the following subdirectories:

* **U-Net predictions:** This subdirectory contains images (*.png*) of the segmented unmyelinated axons in the vagus and pelvic nerve cross-sections under consideration.  

* **Inputs:** There are *.csv* and *.rds* files in this subdirectory. The *.csv* files contain several morphometric characteristics of the segmented unmyelinated axons in each nerve cross-section, along with their centroid locations, which we use in particular to construct the spatial point patterns. The morphometric characteristics are extracted using an open source image processing package Fiji (@schindelin2012fiji). The *.rds* files, drawn using R graphics tools, store the outer boundary and inner holes of the cross-sections that are used for spatial point pattern construction also. These are the main input files for the analysis.

* **Demo Inputs:** This subdirectory contains a small subset of the inputs for quick demonstration purpose.

* **Supporting Files:** Some additional information of the nerve samples, such as sample ID, nerve location and sex, are saved in this subdirectory. The pair-wise orientation of the cross-sections are computed at a certain stage of the pipeline. Those orientations are also stored inside this subdirectory.

