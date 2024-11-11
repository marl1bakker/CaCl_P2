Will need Labeo package Umit-master: https://github.com/LabeoTech/ and matnwb-master: https://github.com/NeurodataWithoutBorders/matnwb/releases

Example data is available on dandiarchive.org: https://doi.org/10.48324/dandi.001210/0.241111.1757 

Try "Example_data_pipeline" to go through the steps of a single acquisition. The first two steps of the pipeline are skipped if you get the data from DANDI. The data you get will be in .nwb format, and you will have to run a small function to transform it into .dat and .mat files, which the rest of the pipeline will work with. This is all explained in the example pipeline as well.  

For the pipeline with all mice, do Preprocessing_Pipeline_CaCl_GCaMP first, then Umit_Pipeline_CaCl_GCaMP and lastly Pipeline_CaCl_GCaMP
