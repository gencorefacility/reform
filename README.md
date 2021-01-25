# <i>ref</i>orm

[*ref*orm](https://gencore.bio.nyu.edu/reform/) is a python-based command line tool that allows for fast, easy and robust editing of reference genome sequence and annotation files.

Execution of *ref*orm requires a reference sequence (fasta), reference annotation (GFF or GTF), the novel sequences to be added (fasta), and corresponding novel annotations (GFF or GTF). A user provides as arguments the name of the modified chromosome and either the position at which the novel sequence is inserted, or the upstream and downstream sequences flanking the novel sequences. This results in the addition and/or deletion of sequence from the reference in the modified fasta file. In addition to the novel annotations, any changes to the reference annotations that result from deleted or interrupted sequence are incorporated into the modified gff.  Importantly, modified gff and fasta files include a record of the modifications.

Learn more at https://gencore.bio.nyu.edu/reform/

## Usage

*ref*orm requires Python3 and Biopython v1.78 or higher. 

Install biopython if you don't already have it:

`pip install biopython`

Invoke the python script:

```
python3 reform.py 
  --chrom=<chrom> \
  --position=<pos> \ 
  --in_fasta=<in_fasta> \
  --in_gff=<in_gff> \
  --ref_fasta=<ref_fasta> \
  --ref_gff=<ref_gff>
```

## Parameters

`chrom` ID of the chromsome to modify

`position` Position in chromosome at which to insert <in_fasta>. Can use `-1` to add to end of chromosome. Note: Either position, or upstream AND downstream sequence must be provided. **Note: Position is 0-based**

`upsteam_fasta` Path to Fasta file with upstream sequence. Note: Either position, or upstream AND downstream sequence must be provided.

`downsteam_fasta` Path to Fasta file with downstream sequence. Note: Either position, or upstream AND downstream sequence must be provided.

`in_fasta` Path to new sequence to be inserted into reference genome in fasta format.

`in_gff` Path to GFF file describing new fasta sequence to be inserted.

`ref_fasta` Path to reference fasta file.

`ref_gff` Path to reference gff file.

## Example

```
python3 reform.py 
  --chrom="I" \
  --upstream_fasta="data/up.fa" \
  --downstream_fasta="data/down.fa" \
  --in_fasta="data/new.fa" \
  --in_gff="data/new.gff" \
  --ref_fasta="data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa" \
  --ref_gff="data/Saccharomyces_cerevisiae.R64-1-1.34.gff3"
```

## Output

`reformed.fa` Modified fasta file.

`reformed.gff3` Modified GFF file.

