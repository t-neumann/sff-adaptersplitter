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

The installation process currently consists of a very simple Makefile.

Just do the following:

    git clone git://github.com/indraniel/sff2fastq.git;
    cd sff2fastq;
    make; # try 'make genome' if at the Genome Center at Washington University
          # or on a Linux distribution from 2008 or earlier

The `sff2fastq` executable should be in the working directory.
Afterwards, you can move the executable to wherever you wish.

NOTES
=====

This has been successfully compiled on Linux/Ubuntu 8.04 & 9.10
workstations (both 32-bit and 64-bit machines), and on Mac OS X (version
10.5).  Compiling on other types of operating systems and architectures
has not been experimented upon.

The FASTQ output produced is of the Sanger FASTQ format as described
here (http://maq.sourceforge.net/fastq.shtml).

Without any given options the default approach is to output trimmed
sequence and quality values.  This is similar in nature to the sff tools
produced by 454 Life Sciences/Roche.

AUTHORS
=======

Indraniel Das (indraniel@gmail.com or idas@wustl.edu)
The Genome Institute at Washington University

Contributors
------------

* [Björn Winckler](https://github.com/b4winckler)
* [James Casbon](https://github.com/jamescasbon)

ACKNOWLEDGEMENTS
================

This software was developed at The Genome Institute at Washington
University, St. Louis, MO.

DISCLAIMER
==========

This software is provided "as is" without warranty of any kind.


March 23, 2010
