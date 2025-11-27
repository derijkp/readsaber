Readsaber
========= 
A program for analysis of read structures in NGS data
Copyright VIB and University of Antwerp

Readsaber
---------

Readsaber was originally developed to detect what kind of artefacts were abundant when optimizing
library preparation and sequencing protocols for enrichment 10x long read sequencing
It will annotate reads in given fastq files or unaligned bam files with
several structural elements such as adapters, TSO sequences, genomic sequence, etc., create a schema or structure
for each each read, and summarizes the results by outputing all schemes present and their number and percentage.
This can be used more general for all kinds of optimization as well as QC.

Installation
------------
Binary packages for Linux can be downloaded from github
([https://github.com/derijkp/readsaber](https://github.com/derijkp/readsaber))

readsaber is distributed as a portable application directory: A
self-contained directory with the readsaber executable (readsaber) and all
needed depencies compiled in a way that they should work on all (except very
ancient) Linux systems.

Installation of the package is as simple as downloading the
[distribution](https://github.com/derijkp/readsaber/releases/download/0.1.0/readsaber-0.1.0-linux-x86_64.tar.gz)
from github
([https://github.com/derijkp/readsaber](https://github.com/derijkp/readsaber))
and unpacking it, e.g.: 
```
cd ~/bin
wget https://github.com/derijkp/readsaber/releases/download/0.1.0/readsaber-0.1.0-linux-x86_64.tar.gz
tar xvzf readsaber-0.1.0-linux-x86_64.tar.gz
rm readsaber-0.1.0-linux-x86_64.tar.gz
```

You can call the executables (readsaber, cg) directly from the directory
using the path (e.g. `~/bin/readsaber-0.1.0-linux-x86_64/readsaber ..`) 
or by placing the directory in the PATH environment variable (e.g. using 
`export PATH=~/bin/:$PATH`)
You can also place soft-links to the executables in a directory already in
the PATH. (remark: The executable itself needs to stay in the application
directory to find it's dependencies), e.g.
```
cd ~/bin
ln -s readsaber-0.1.0-linux-x86_64/readsaber .
ln -s readsaber-0.1.0-linux-x86_64/cg .
```

readsaber is largely implemented within [genomecomb](https://github.com/derijkp/genomecomb), 
and its distribution comes with an appropriate full version of genomecomb,
which can be run using the cg executable, which also provides
multiple usefull extra tools for querying tsv files, etc.

Example/test run
----------------
As an example/test, the following code shows you how to download an example data set and run readsaber on it:
```
# download and unpack test data
wget https://github.com/derijkp/readsaber/releases/download/v0.1.0/readsaber_test.tar.gz
tar xvzf readsaber_test.tar.gz

cd readsaber_test

# make refdir; This test data is human, but limited to chromosome 17, so we can use only a partial genome reference.
readsaber_makerefdir g17 genome.fa

# Run readsaber
readsaber -refseq g17 annotations.fa readannots.tsv fastq/*

# Check summary results
cg viz readannots_summary.tsv

```

Reference data
--------------
A readsaber analysis can use a reference genome, with different types of indexes and supporting
files. These must be provided in a reference directory.

You can use the included command readsaber_makerefdir to make a reference
directory starting from a fasta.
```
readsaber_makerefdir refdir genomesequence.fasta
```
where 
* `genomesequence.fasta` is a multifasta file with the genomesequence and transcripts.gtf

You can also use genomecomb reference directories for this; These can be downloaded from
the genomecomb website for a number of species (or created new) as
described in the [genomecomb installation
documentation](https://derijkp.github.io/genomecomb/install.html). 

When downloadinga reference, be sure to also download and install the matching minimap2 indexes, e.g.
```
wget https://genomecomb.bioinf.be/download/refdb_hg38-0.110.0.tar.gz
tar xvzf refdb_hg38-0.110.0.tar.gz
wget https://genomecomb.bioinf.be/download/refdb_hg38-minimap2-0.110.0.tar.gz
tar xvzf refdb_hg38-minimap2-0.110.0.tar.gz
```
### PacBio indexes
By default readsaber_makerefdir will only create an ont minimap2 index for the genome.
If you want to analyze PacBio data, the appropriate index can be added by giving the 
option '-pacbioindex 1' to readsaber_makerefdir

Running readsaber
-----------------
You can run readsaber using the following command
```
readsaber ?options? annotationfile resultfile fastq ?fastq? ...
```
This will analyse the data in the given fastq (or unaligned bam) files
and write the results per read in **result**
and a summary (counts and percentage of each structure detected) of the 
results (which you mostly look at) in <result.base>_summary.tsv


Possible options are:
`-refseq`
    reference directory with genomesequence, etc. as described previously.
    Alternatively, you can also give an alignment (bam or cram file) of the reads to the genomic reference if you have already have this
    You can give this option more than once; Annotations from subsequent refdirs are indicated with transcript2, transcript3, etc.

`-refseqannot`
    by default annotations from subsequent refdirs are indicated with transcript2, transcript3, etc. ('transcript$postfix')
    Using -refseqannot you can give an alternative value, e.g. 
        '$chromosome' to show the chromosome the transcript/genomic sequence was found on,
        'genomic' for allways giving the string genomic (without the numeric postfix)
    You can give this option more than once.

`-addsequences 0/1`
    set to 1 to add sequence data to the per read output file

`-completeness number`
    annotations/assignments matching less then **number** percent of the annotation sequence are ignored (default 50)

`-ignoreN number`
    stretches with >= **number** unassigned bases will be indicated with N in the structure (default 3)

`-polyT number`
    A polyT will only be called by the specific polyT caller if at least **number** Ts (or As) were seen in a stretch (default 7).
    Set to 0 to not run the specific polyT caller

`-alimethod string`
    define the alignment method used for the annotation alignment. A preset can be added using a "_",
    e.g. the default is minimap2_ontshort for using the minimap2 aligner with the ontshort (minimap2 ont 
    preset with special options optimized for short matches) preset.
    For short read data, the method "bwa_short" is advised.

`-refalimethod string`
    define the alignment method used for the reference alignments. A preset can be added using a "_",
    e.g. the default is minimap2_splicesens for using the minimap2 aligner with the "splicesens" 
    preset (minimap2 splice preset with options for extra sensitivity).
    For short read data, the method "bwa_short" is advised.

`-keepintermediate 0/1`
    set to 1 to keep intermediate files for development/debugging

`-d`
    By default the command is run using a single core (=slow). Use the `-d` option to specify the
    manner of job distribution/parallelisation. Use a number to specify distribution
    over (max) the given number of cores on the local machine, while specifying "sge" or
    "slurm" here is one way (see next for other) to distribute jobs over a Grid Engine or SLURM cluster.
    On a cluster the command will finish after submitting all jobs (with dependencies).
    For distributed runs a tab separated log file is created (in the projectdir) named
    process_project_<projectname>.<starttime>.running
    This contains information on all started jobs, and when all jobs are finished, 
    this log file will be renamed to process_project_<projectname>.<starttime>.finished on success
    or to process_project_<projectname>.<starttime>.error when there was an error
    encountered. In this case, specific jobs that had errors can be found in the 
    logfile.
    More information on options for distribution options can be found in the
    [genomecomb joboptions help](https://derijkp.github.io/genomecomb/joboptions.html)

`-dsubmit`
    Some clusters limit the number of concurrent submitted jobs. This option provides an alternative
    method for distribution that works on those: You specify the (maximum) number of jobs
    that can be submitted using the -d <number> option, and add sge or slurm here to indicate
    which job management system will be submitted to. Using this method, the command will have keep
    running until all jobs are finished.

`-dmaxmem`
    The maximum memory to be used (reserved) when running distributed on a local machine.
    This will stop only jobs from starting as long as they request more memory than currently 
    available (based on requested memory by running jobs). Jobs that use more memory than requested
    will not be stopped by this. This option will have no effect when running on a cluster


Various **other run options** are

`-v`
    default 0, increase (up to 2) to increase the verbosity level, i.e. how much information 
    starting up jobs, dependencies, etc. is displayed

`-stack`
    set to 1 (default 0) to show an extended stack trace on error (mainly for debugging)

Resultfile
----------
The resultfile given to the commmand is a tab-separated value file with the detected schemas for all (tested) reads.
This file is mainly for checking in detail. We will typically check the summary files (see next) first (or only).
Resultfile has the following fields

`rname`
    read name

`readsize`
    size of the read

`schema`
    schema/structure of the read: it consists of a list alternating strand with detected element,e.g.
    `+ Read1 ~ N + polyT ~ transcript + TSO` where 
        "Read1" and "TSO" are matches (forward strand) to sequences given in the annotation file.
        "~ transcript" is a match (spliced alignment) to anywhaer on the first reference genome given (would be "transcript2" for the second, etc.).
        "~ N" is added when there are >= 'ignoreN' unassigned bases in the stretch

`shortschema`
    short version of the schema: same there are no 'N's recorded in the structure

`schema2`
    long version of the schema: all recorded elements have added _<size_of_element>, e.g.
    `+ Read1_33 ~ N_28 + polyT_9 ~ transcript_91 + TSO_30`
    
`sequences`
    optional field added when `-addsequences 1` is used: the same as schema, but after each element the sequence of the element is added.

Summary files
--------------
Readsaber also creates 2 summary files giving a breakdown of the abundance of the different structures/schemas. 
Their name is based on the the root of the resultfile, with added _summary (before the extension), so e.g. if
resultfile has a file name "readannot-test.tsv", the summary file will be named "readannot-test_summary.tsv"
The summary file is a tab-separated value file with the following fields:

`schema`
    the schema this line is about

`count`
    the number of reads that follows this schema

`percent`
    the percentage of reads that follows this schema

`q1_readsize`
    quartile 1 of readsize (of reads with this schema)

`avg_readsize`
    average readsize (of reads with this schema)

`q3_readsize`
    quartile 3 of readsize (of reads with this schema)

The file "readannot-test_shortsummary.tsv" is also made, and gives the same information bases on the "shortschema" field

License
-------
The use of this application is governed by the GPL (license.txt).

How to contact me
-----------------

Peter De Rijk
VIB - UAntwerp Center for Molecular Neurology, Neuromics Support Facility - Bioinformatics
University of Antwerp
Universiteitsplein 1
B-2610 Antwerpen, Belgium

tel.: +32-03-265.10.40
E-mail: Peter.DeRijk@uantwerpen.vib.be

