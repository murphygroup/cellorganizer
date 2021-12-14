slml2info
********
This function generates a report from information extracted from a generative model file.

**Example:**

*slml2info(filenames, options)*

=============================  ===============================================================
        Inputs                                             Outputs
=============================  ===============================================================
  filenames                     List of files
  options                       Options structure
=============================  ===============================================================

General Options
================

Generic Options
^^^^^^^^^^^^^^^

*options.output_directory* (mandatory) **['pwd/report']**
  * Name of directory that the resulting report is saved as

*options.labels* (optional) **[{}]**
  * List of labels used for shape space plots
  
 T-Cell Model
 ^^^^^^^^^^^^^^^
 
*options.model.tcell.region_2_radius* (optional) **[8, 8/3]**
    * Define the cylinder for the larger region (numerator) the radius of the ring.
     
*options.model.tcell.region_2_thickness* (optional) **[4, 4]**
    * The height of the cylinder.
    
*options.model.tcell.region_2_start_ind* (optional) **[1, 1]**
    * The index of the top of the cylinder (indices start from the top of the cell, i.e., 1 means from the topmost slice, nearest the synapse).

*options.model.tcell.enrichment_region_type* (optional) **['ring']**
    * Define the enrichment region as the ring

*options.model.tcell.should_use_global_enrichment_region* (optional) **[false]**
    * Don’t use the top n% fluorescence to create the region
    
*options.model.tcell.use_user_defined_enrichment_region* (optional) **[true]**
    * Use the user defined enrichment region
    
*options.model.tcell.enrichment_over_certain_region* (optional) **[true]**
    * Calculate enrichment over a certain region 

*options.model.tcell.enrichment_bottom_region_type* (optional) **['top_fluorescence']**
    * Use intensity in this region as the denominator for calculating enrichment.	

*options.model.tcell.error_bar_type* (optional) **['sem']**
    * Set the errorbar type, either SEM or SD
    
*options.model.tcell.save_result_filename* (optional) **[N/A]**
    * Specify filename in which to save enrichment results
