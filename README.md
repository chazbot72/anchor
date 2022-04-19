# anchor
A consensus sequence exploration tool.

## Overview
Anchor is a simple python tool for interrogating amino acid sequence databases in reference to a single sequence.  It works by preparing a database of pairwise alignments to your reference sequence keyed on position that can be queried for homology.

### What?

Anchor parses any number of pairwise alignments to a reference amino acid sequence and numbers them in reference to the original sequence.  Gaps are labelled starting with "a" and progress to "aa", etc.  This positionally keyed database can then be queried in an easy way to check for conservation of residues at certain positions.

## Installing anchor
Anchor is written in Python 3 and requires only two additional packages via pip:

- tqdm - https://github.com/tqdm/tqdm
- BioPython - https://biopython.org/

You also will need access to Clustal Omega binaries for alignment:

- Clustal Omega - http://www.clustal.org/omega/

Once you have these installed on your system, you only need a reference amino acid sequence and amino acid sequence database in FASTA format.

## Using anchor
Example data (HXB2.fasta / HIV1_SFL_2018_env_PRO.fasta) is provided for your convenience in testing the system.  You will need to create the directory `alignments/` to run the examples below.  anchor runs in 3 steps:

#### General Arguments 
`--action align/analyze/query` tells anchor what you want to do
`--db_root PATH` tells anchor where you want to do it

### Alignment
To align all of your database sequences to the reference, run anchor with the following options:

`python anchor.py --action align --anchor HXB2.fasta --ali_dbfile HIV1_SFL_2018_env_PRO.fasta`

This will run alignments using ClustalO and save the results in the default folder structure (alignments/).  Note that a temporary file is created for every alignment (see Other Arguments to change it from the defaults)

#### Other Arguments
`--ali_root STR` Alignment root (default is <db_root>/alignments)
`--ali_infile STR` Temp file for ClustalO alignment (default is <db_root>/inFile.fasta)


### Parsing
This command parses alignments and builds the database .pickle file that can later be loaded and queried:

`python anchor.py --action analyze`

once the alignments are parsed, a file (`db.pickle`) will be output.

#### Other Arguments
`--db_pickle STR` DB pickle file name (default is $PWD/db.pickle)

### Querying
Now that we have a database, we can query it for information.  A query has two components:  anchor and reporting regions.  Anchor regions are denoted in the following manner:

`--query "240:T, 232:T/K"`

Separate positions by commas, and add the residues you want to anchor on at the positions as backslash (`/`) delimited lists preceded by a colon (`:`).

If you don't designate reporting residues, the program returns results for all positions in the database.  Otherwise, you can specify reported positions like so:

`--report "240, 241, 242"`

The following command works on the example dataset as a test so you can verify the program is working:

`python anchor.py --action query --query "240:T,232:T/K" --report "240, 241, 242"`

This will output a .csv file keyed with timestamp into anchor's working directory.

#### Other Arguments
`--db_pickle STR` DB pickle file name (default is $PWD/db.pickle)

### Results
The resulting csv file is a breakdown of every position (or requested positions) relative to the reference sequence by amino acid (or gap).  We use this information for tracking prevalence of certain sequences in known sequence databases.

## Example Data
The example dataset provided here was compiled from the LANL HIV sequence database (http://www.hiv.lanl.gov/) 2018 Group M super filtered amino acid web alignment.  The reference sequence used is HxB2 Env. (https://www.ncbi.nlm.nih.gov/nuccore/K03455)

## Citing anchor
If you use anchor in your research, please cite our paper here:

- https://doi.org/10.1371/journal.ppat.1008753

## Future Plans
- Add anchor as a web functionality to ward.scripps.edu/gld
- Add options to report the IDs of sequences containing desired sequences.
