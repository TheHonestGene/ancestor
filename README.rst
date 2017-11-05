**********************************
Ancestor: Simple ancestry analysis
**********************************

Ancestor is a very simple library for doing ancestry analysis
of a genotype dataset (i.e. 23andMe)

Users can run Ancestor as either a standalone command line tool
or import it into their python library


Requirements
***********

The ancestry relies on the 1000 Genomes dataset.
The dataset has to be provided as an HDF5 file and be in a specific format

To play around pcs files, a test genotype and the weights in HDF5 (calculated based on the above 1000 genomes dataset) format can be found in the https://github.com/TheHonestGene/ancestor/tree/master/tests/data folder.


How-To
***********

There are 2 steps to finish in the ancestry pipeline.
Some of the steps only have to be done once or once for each genotyping plattform and dataset respectively
Each step maps to a subcommand of the command line script.
To get information about the various subcommands the user can run:

.. code:: bash

      ancestor -h

Step 1 (Optional) - Convert weights file from CSV format to HDF5 format
===============================================================
This step is not required but can speed up the rest of the pipeline a bit.
It also only have to be done once.

.. code:: bash

      ancestor convert weights_file weights_file.hdf5



Step 2 - Calculating PC projections and admixture decomposition information for a reference panel
===============================================================
For the admixture analysis and membership testing of an individual genotype in a specific population, the PC projections
and admixture decomposition for a reference genotype panel have to be calculated.
This has to be done once per genotyping platform/version and weights file.
It is possible to use any reference genotype panel but the most comprehensive one is the 1000 genomes reference panel.

.. code:: bash

      ancestor prepare 1001genomes.hdf5 weights.hdf5 1000_ref_pcs_file.hdf5 --ntmap 23andme_v4_nt_map.pickled

The --ntmap argument specifies the nucleotide map which is specific to the genotyping platform/version.
This file can be created with the `imputor library <https://github.com/TheHonestGene/imputor>`_


Step 3 - Membership test and admixture analysis of individual genotype & plotting
===============================================================
To calculate the PC projections and admixture decomposition of an individual genotype
and optionally test membership and plot the PCs the individual imputed (using the imputor library) genotype
as well as the PCs file for the reference genotype panel that was generated in Step 2 and the weights file have to be provided.


.. code:: bash

      ancestor pcs genome_imputed.hdf5 weights.hdf5 1000_ref_pcs_file.hdf5 --plot pc_plot.png --check GBR

The paramters --check and --plot are optional and used for testing membership in a population and plotting

Test
-------------

The test suite can be run with::

      $ python setup.py test

Installation
--------------

Of course, the recommended installation method is pip::

    $ pip install ancestor

Thank You
-----------

Thanks for checking this library out! We hope you find it useful.

Of course, there's always room for improvement. Feel free to `open an issue <https://github.com/TheHonestGene/ancestor/issues>`_ so we can make it better.
