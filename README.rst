**********************************
Ancestor: Simple ancestry analysis
**********************************

Ancestor is a very simple library for doing ancestry analysis
of a genotype dataset (i.e. 23andMe)

Users can run Ancestor as either a standalone command line tool
or import it into their python library


Requirements
***********

The ancestry relies on the 1000 Genomes HapMap dataset.
The dataset has to be provided as an HDF5 file and be in a specific format
A version can be downloaded here.


How-To
***********

The command line can run ancestry analysis for a given genotype and optionally create a plot of the PC space.::

      $ ancestor --weights weights.csv --hapmap hapmap.hdf5 --pcs hapmap_pcs.hdf5 --plot ancestry.png genotype.hdf5

If the --pcs parameter is specified it will store the calculated PCs for all individuals in the Hapmap dataset.
The next time ancestry is run for another genotype, the cached version can be used:

      $ ancestor --weights weights.csv --pcs hapmap_pcs.hdf5 --plot ancestry.png genotype.hdf5


Installation
--------------

Of course, the recommended installation method is pip::

    $ pip install ancestor

Thank You
-----------

Thanks for checking this library out! We hope you find it useful.

Of course, there's always room for improvement. Feel free to `open an issue <https://github.com/TheHonestGene/ancestor/issues>`_ so we can make Imputor better.
