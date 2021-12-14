Virtual Cell
------------

The software supports saving generated geometries as valid VCML (Virtual Cell Markup Language) files that can be imported into `Virtual Cell <https://vcell.org/>`_, a platform to perform biological system modeling and simulations. To do this, use the output flag

.. code:: matlab

	options.output.writeVCML = true;

For simple examples, investigate and run:

* ``demo3D58``: Generate a single framework and vesicles from a 3D HeLa ratio model and write to VCML.
* ``demo3D60``: Generate a single framework and vesicles from a 3D HeLa SPHARM model and write to VCML.
* ``demo3D63``: Generate a single framework from a 3D HeLa SPHARM model and write to VCML with a reaction network and spatial and compartmental simulations ready to be run in Virtual Cell.

For an example application that generates 100 synthetic geometries with reaction networks (as in ``demo3D63``), writes them to individual VCML files, and combines them into one VCML file for easy import into Virtual Cell, see ``generate_simulation_instances_SarmaGhosh2012``.

Options
^^^^^^^

The table below describes all VCML-related options. All options are fields within the ``options.output.VCML`` structure.

Options are listed as required or optional in the case that ``options.output.writeVCML = true``; otherwise, they have no effect.

=============================   ========    ===========
Option                          Required    Description
=============================   ========    ===========
writeVCML                       required    boolean flag specifying whether to write out VCML files for use with Virtual Cell. Default is false.
input_filename                  optional    string specifying Virtual Cell VCML file with biochemistry which will be combined with generated geometry in output file. Default is empty string.
downsampling                    optional    downsampling fraction for the creation of object files (1 means no downsampling, 1/5 means 1/5 the size). Default is 1.
translations                    optional    N x 2 cell array of strings (first column) to be replaced by other strings (second column).
=============================   ========    ===========

Tutorials
^^^^^^^^

``demo3D63``
~~~~~~~~~~~~

``demo3D63`` generates a single framework from a 3D HeLa SPHARM model and writes it to VCML with a reaction network and spatial and compartmental simulations ready to be run in Virtual Cell.

Key options:

* ``options.synthesis = 'framework';``: The reaction network occupies only the cytoplasm and the nucleus.
* ``options.output.VCML.writeVCML = true;``: Instructs CellOrganizer to create a VCML output file.
* ``options.VCML.downsampling = 1/2;``: Reduces the size of the 3D image holding the cell, which makes the simulation run faster and reduces the storage requirements for the simulation's results.
* ``options.VCML.translations = {'cell', 'CP'; 'nuc', 'NU'; 'nucleus', 'NU'; 'lamp2_mat_tfr_mat', 'EN'; 'CP_EC', 'PM'; 'CP_EN', 'EM'; 'CP_NU', 'NM'};``: Instructs CellOrganizer to translate ``'cell'`` into ``'CP'``, etc. so the geometry's compartment names match those of the input VCML file.
* ``reaction_network_file = '../../../data/SarmaGhosh2012ForCOdiff1e-2.vcml';``, ``options.VCML.input_filename = reaction_network_file;``: Instructs CellOrganizer to incorporate the input VCML file's reaction network and simulations into the output VCML file.
* ``options.VCML.end_time = 4000;``, etc.: Set Virtual Cell simulation parameters.

Expected outputs:

* ``demo3D63/cell1/cell.vcml``: The VCML file to be opened in Virtual Cell.

Using the generated VCML in Virtual Cell:

1. In Virtual Cell, select ``File > Open > Local...`` and select this file. You should see no errors or warnings listed in the ``Problems`` tab.
2. Click ``img_Spatial`` under ``Applications``, then ``Simulations``, then ``img_Spatial_Simulation_1``, then click the green triangle button to run the simulation on the Virtual Cell server. The simulation will take several hours to complete.
