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

The command line can run ancestry analysis for a given genotype and optionally create a plot of the PC space.::

      $ ancestor --hapmap hapmap.hdf5 --pcs hapmap_pcs.hdf5 --plot ancestry.png genotype.hdf5 weights.csv

If the --pcs parameter is specified it will store the calculated PCs for all individuals in the Hapmap dataset.
The next time ancestry is run for another genotype, the cached version can be used::

      $ ancestor --pcs hapmap_pcs.hdf5 --plot ancestry.png genotype.hdf5 weights.csv 

Alternatively the weights can also be provided as in an hdf5 format::

      $ ancestor --pcs hapmap_pcs.hdf5 --plot ancestry.png genotype.hdf5 weights.hdf5


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
