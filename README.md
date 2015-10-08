FAST2HDF5 Converter program
===========================
![alt text](https://img.shields.io/badge/Python-2.7-green.svg 'Python version')
![alt text](https://img.shields.io/badge/Linux-supported-green.svg 'Supported OS')

Purpose
-------
Contributors:
* Leendert Hayen (leendert.hayen@fys.kuleuven.be)

Uses the fasterac library, written by LPC Caen, to extract both ADC and Counter data from .faster files generated by a MOSAHR motherboard (http://www.faster.in2p3.fr). It then writes this data to an HDF5 file. For each FASTER channel, a Detector group is constructed, container two tables: ADC and Counter. The group has the FASTER label as an attribute.

Compilation options
-------------------
h5cc -o fast2hdf5 -lfasterac -std=c11 fast2hdf5.c

Dependencies
------------
* [HDF5] (https://www.hdfgroup.org/HDF5/)
* [fasterac] (http://faster.in2p3.fr/index.php/installation/fasterac-package/13-divers)

