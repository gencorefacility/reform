# <i>ref</i>orm

[*ref*orm](https://gencore.bio.nyu.edu/reform/) is a python-based command line tool that allows for fast, easy and robust editing of reference genome sequence and annotation files.

Execution of *ref*orm requires a reference sequence (fasta), reference annotation (GFF or GTF), the novel sequences to be added (fasta), and corresponding novel annotations (GFF or GTF). A user provides as arguments the name of the modified chromosome and either the position at which the novel sequence is inserted, or the upstream and downstream sequences flanking the novel sequences. This results in the addition and/or deletion of sequence from the reference in the modified fasta file. In addition to the novel annotations, any changes to the reference annotations that result from deleted or interrupted sequence are incorporated into the modified gff.  Importantly, modified gff and fasta files include a record of the modifications.

In addition to the editing functionality described above, *ref*orm also supports the addition of novel chromosomes. This requires a reference sequence (FASTA) and annotation (GFF or GTF), along with the novel chromosome to be added. Users must provide the novel sequence (FASTA) and its corresponding annotation (GFF or GTF). The new chromosome is appended to the reference files, and all associated annotations are incorporated accordingly.

Learn more at https://gencore.bio.nyu.edu/reform/

## Usage

*ref*orm requires Python3 and Biopython v1.78 or higher. 

Install biopython if you don't already have it:

`pip install biopython`

Reform supports reading and writing .gz files using gzip. To accelerate compression and decompression, it optionally supports pgzip, a parallel implementation of gzip. Users must install pgzip separately to enable this feature.

Install pgzip if you don't already have it:

`pip install pgzip`   

Invoke the python script:

```bash
### Edit a sequence within position
python3 reform.py 
  --chrom=<chrom> \
  --position=<pos1>,<pos2>,<pos2> \ 
  --in_fasta=<in_fasta1>,<in_fasta2>,<in_fasta3> \
  --in_gff=<in_gff1>,<in_gff2>,<in_gff3> \
  --ref_fasta=<ref_fasta> \
  --ref_gff=<ref_gff>
```

```bash
### Edit a sequence within upstream & downstream
python3 reform.py 
  --chrom=<chrom> \
  --upstream_fasta=<upstream1>, <upstream2>, <upstream3> \
  --downstream_fasta=<downstream1>,<downstream2>,<downstream3> \
  --in_fasta=<in_fasta1>,<in_fasta2>,<in_fasta3> \
  --in_gff=<in_gff1>,<in_gff2>,<in_gff3> \
  --ref_fasta=<ref_fasta> \
  --ref_gff=<ref_gff>
```

```bash
### Append a novel chromosome sequence
python3 reform.py 
  --new_chrom=<chrom> \
  --in_fasta=<in_fasta1>,<in_fasta2>,<in_fasta3> \
  --in_gff=<in_gff1>,<in_gff2>,<in_gff3> \
  --ref_fasta=<ref_fasta> \
  --ref_gff=<ref_gff>
```

## Parameters

`chrom` ID of the chromsome to modify

`new_chrom` ID of the novel chromsome to append

`position` Position in chromosome at which to insert <in_fasta>. Can use `-1` to add to end of chromosome. Note: Either position, or upstream AND downstream sequence must be provided. To perform multiple edits in one run, provide multiple positions separated by commas (e.g., 0,5,-1). **Note: Position is 0-based**

`upstream_fasta` Paths to Fasta file with upstream sequence. Note: Either position, or upstream AND downstream sequence must be provided.

`downstream_fasta` Paths to Fasta file with downstream sequence. Note: Either position, or upstream AND downstream sequence must be provided.

`in_fasta` Paths to new sequence to be inserted into reference genome in fasta format.

`in_gff` Paths to GFF file describing new fasta sequence to be inserted.

`ref_fasta` Path to reference fasta file.

`ref_gff` Path to reference gff file.

## Example

```
python3 reform.py 
  --chrom="I" \
  --upstream_fasta="data/up1.fa,data/up2.fa,data/up3.fa" \
  --downstream_fasta="data/down1.fa,data/down2.fa,data/down3.fa" \
  --in_fasta="data/new1.fa,data/new2.fa,data/new3.fa" \
  --in_gff="data/new1.gff,data/new2.gff,data/new3.gff" \
  --ref_fasta="data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa" \
  --ref_gff="data/Saccharomyces_cerevisiae.R64-1-1.34.gff3"
```

## Output

`reformed.fa` Modified fasta file.

`reformed.gff3` Modified GFF file.

## Test
After local deployment or modification, you can use `test_reform.py` to verify the functionality of Reform.  
This script contains an automated test suite using the Python `unittest` framework. It verifies the correctness of Reform across a variety of genome editing scenarios.

To run all tests:

```bash
python3 test_reform.py
```