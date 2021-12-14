slml2slml
********
This function combines multiple models into a single model file.

**Example:**

*slml2report(files, options)*

=============================  ===============================================================
        Inputs                                             Outputs
=============================  ===============================================================
  files                         list of paths of models that need to be combined
  options                       List of options
=============================  ===============================================================

*options.selection* (mandatory)
  * a matrix used to specify what submodels should be used from each file.

*options.output_filename* (optional) **["model.mat"]**
  * the file name of output model.
