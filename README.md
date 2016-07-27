SYNOPSIS
========

The program `sff-adaptersplitter` utilizes source from `sff2fastq` written by Indraniel Das (indraniel@gmail.com or idas@wustl.edu)from The Genome Institute at Washington University, St. Louis, MO.

It takes a multiplexed SFF file, produced by the 454 genome sequencer, and demultiplexes the reads according to given barcodes into separate SFF files per adapter.

USAGE
=====

Given an SFF file, `file.sff`, and a barcode file, `barcode.txt`,you can simply run:

    sff-adaptersplitter barcode.txt file.sff


Options
-------

`sff-adaptersplitter` can search for adapters on the entire read sequence or only using the quality trimmed sequence. Trimming is turned on by default, but can be turned off like this:

    sff-adaptersplitter barcode.txt file.sff <trimming 0|1 (default: 1)>
    
Barcode files
-------------

The barcode files supplied to `sff-adaptersplitter` are simple tab-separated text files with barcode-ID in the first column and barcode-sequence in the second column, e.g.

    IonXpress_001	CTAAGGTAAC
    IonXpress_002	TAAGGAGAAC
    IonXpress_003	AAGAGGATTC


INSTALLATION
============

In order to build NextGenMap only "cMake":http://www.cmake.org/ (>=2.8) and g++ are required. Typically these tools should be already available on your computer. If not, please install them.
The installation process utilizes "cMake":http://www.cmake.org/ (>=2.8) to create a simple Makefile, which can then build with standard "make".

    mkdir build
    cd build
    cmake ..
    make
    ../bin/sffsplit

NOTES
=====

This has been successfully compiled on Linux/Ubuntu 14.04 LTS and Linux/Ubuntu CentOS 7
workstations (64-bit machines), and on Mac OS X (10.11).  Compiling on other types of operating systems and architectures
has not been experimented upon.

AUTHORS
=======

Tobias Neumann (tobias.neumann.at@gmail.com)

DISCLAIMER
==========

This software is provided "as is" without warranty of any kind.


August 9, 2012
