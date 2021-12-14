MCell
-----

The software supports saving generated geometries as valid MCell MDL (Model Description Language) files that can be run using `MCell <https://mcell.org/>`_, a program for simulating reaction networks in realistic geometry. To do this, use the output flag

.. code:: matlab

	options.output.writeMCellMDL = true;

For an example application that generates multiple synthetic geometries and combines them with reaction networks and writes them to modified copies of input MCell files, see ``generate_simulation_instances_SarmaGhosh2012``.

Options
^^^^^^^

The table below describes all MCell-related options. All options are fields within the ``options.output.VCML`` structure.

Options are listed as required or optional in the case that ``options.output.writeMCellMDL = true``; otherwise, they have no effect.

=============================   ========    ===========
Option                          Required    Description
=============================   ========    ===========
output.writeMCellMDL            required    boolean flag specifying whether to write out MCell files for use with MCell. Default is false.
input_filename_pattern          optional    string specifying pattern matching a set of MCell MDL files to be combined with generated MDL files. This should be empty or `[path][prefix].*.mdl`. If not empty, CellOrganizer will only generate the geometry file and will copy the other files matching the pattern to the output directory, and it is the user's responsibility to ensure compatibility between the input and CellOrganizer's output. Default is `''`.
=============================   ========    ===========

Tutorials
^^^^^^^^

``generate_simulation_instances_SarmaGhosh2012``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``generate_simulation_instances_SarmaGhosh2012`` generates a single framework from a 3D HeLa SPHARM model and writes it to VCML with a reaction network and spatial and compartmental simulations ready to be run in Virtual Cell.

Key options:

* ...

Expected outputs:

* ``applications/generate_simulation_instances_SarmaGhosh2012/img_lotvol5r3surfcatal_0.50_*_4000sec/cell1/cell.*.mdl``: The MDL files to be run using MCell.

Using the generated MDL in MCell:

1. ...

