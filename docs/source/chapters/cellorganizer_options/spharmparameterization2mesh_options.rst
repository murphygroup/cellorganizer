SPHARMparameterization2mesh
********
This function takes SPHARM parameters and reconstructs a mesh from them.

A CellOrganizer SPHARM parameters consists of five components,

1) Face Vectors
2) Vertices
3) Faces
4) Spherical Vectors
5) Cost Matrix

**Example Call:**

*SPHARMparameterization2mesh(params, options)*

=============================  ===============================================================
        Inputs                                            Description
=============================  ===============================================================
  params                        SPHARM descriptors (.mat file)
  options                       (optional) List of options
=============================  ===============================================================


General Options
================

Generic Options
^^^^^^^^^^^^^^^

*options.debug* (optional) **default: ['false']**
    * Debug flag to troubleshoot specific issues.

*options.cropping* (optional) **default: ['tight']**
    * Method in how to rasterize a mesh.

*options.oversampling_scale* (optional) **default: [1]**
    * If cropping is not 'tight', it is the degree in which the mesh is scaled.

*options.meshtype.type* (optional) **default: ['even']**
    * Type of 3D mesh that is generated. Options are even or triangular.

*options.meshtype.nPhi* (optional) **default: [64]**
    * Number of meridians if mesh type is 'even'.

*options.meshtype.nTheta* (optional) **default: [32]**
    * Number of parallels if mesh type is 'even'.

*options.meshtype.nVerticies* (optional) **default: [4002]**
    * Number of verticies if mesh type is 'triangular'.

*options.figtitle* (optional) **default: [[]]**
    * Title to put on figure that is generated.

*options.plot* (optional) **default: ['false']**
    * Generate figure or not.

*options.filename* (optional) **default: [[]]**
    * Filename to save figure as.

*options.dpi* (optional) **default: [150]**
    * resolution (in pixels per inch) to export bitmap outputs at, keeping the dimensions of the on-screen figure.