img2SPHARMparameterization
********
This function takes an image array and turns it into SPHARM parameters that can then be used to recreate the image or a mesh counterpart.

A CellOrganizer SPHARM parameters consists of five components,

1) Face Vectors
2) Vertices
3) Faces
4) Spherical Vectors
5) Cost Matrix

**Example Call:**

*img2SPHARMparameterization(image, options)*

=============================  ===============================================================
        Inputs                                            Description
=============================  ===============================================================
  image array                   2D/3D array
  options                       (optional) List of options
=============================  ===============================================================


General Options
================

Generic Options
^^^^^^^^^^^^^^^

*options.NMfirsttry_maxiter* (optional) **default: [300]**
    * Maximum iterations of optimization of topological geometry in spherical model in the first run.

*options.NMretry_maxiter* (optional) **default: [100]**
    * Maximum iterations of optimization of topological geometry in spherical model if the first try fails.

*options.NMretry_maxiterbig* (optional) **default: [300]**
    * Maximum iterations of optimization of topological geometry in spherical model if the second try fails.

*options.NMcost_tol* (optional) **default: [1e-7]**
    * The minimum cost of lagrangian optimization if the cost function is less than this the optimization completes. Decreasing this value will reduce compute time but potentially will also reduce model quality.

*options.NMlargr_tol* (optional) **default: [1e-7]**
    * The absolute difference between two iterations of lagrangian optimization, if smaller than this value the optimization completes. Decreasing this value will reduce compute time but potentially will also reduce model quality.

*options.maxDeg* (optional) **default: [31]**
    * Degree of spherical harmonic descriptor

*options.hd_thresh* (optional) **default: [10]**
    * Threshold for error tolerance for a given cell. If above this parameter the cell is discarded.
