# <i>ref</i>orm

[*ref*orm](https://gencore.bio.nyu.edu//) is a Python-based command-line tool for fast, robust, and flexible editing of reference genome sequence and annotation files.

To perform an edit, *ref*orm requires a reference genome (FASTA), its annotation file (GFF or GTF), a novel sequence to be inserted (FASTA), and the corresponding annotation (GFF or GTF). The user specifies either:

- the chromosome and the position at which to insert the novel sequence, or
- the chromosome along with the upstream and downstream flanking sequences.

The result is a modified reference genome (FASTA) and annotation file (GFF), incorporating the novel sequence and its annotations. Any reference annotations affected by the insertion or deletion are automatically updated. All modifications are documented within the output files.

In addition to modifying existing chromosomes, *ref*orm also supports appending entirely new chromosomes. In this mode, users provide the novel chromosome’s sequence and annotations, which are added to the reference genome and integrated into the annotation file.

Learn more at https://gencore.bio.nyu.edu/reform/

## Usage

*ref*orm requires Python3 and Biopython v1.78 or higher. 

Install biopython if you don't already have it:

`pip install biopython>=1.78`

*ref*orm supports reading and writing .gz files using gzip. To accelerate compression and decompression, it optionally supports pgzip, a parallel implementation of gzip. Users must install pgzip separately to enable this feature.

*Optional:* Install pgzip if you don't already have it:

`pip install pgzip`   

Invoke the python script:

```
### Minimal Example (Single Edit)
python3 reform.py \
  --chrom=<chrom> \
  --position=<position> \
  --in_fasta=<input_fasta.fa> \
  --in_gff=<input_annotations.gff> \
  --ref_fasta=<reference_genome.fa> \
  --ref_gff=<reference_annotations.gff3>
```

## Parameters

- `chrom`: ID of the chromosome to **modify**. **Required** unless `new_chrom` is specified. Cannot be used together with `new_chrom`.

- `new_chrom`: ID of the novel chromosome to **append**. **Required** if you're adding a new chromosome. Cannot be used together with `chrom`.

- `position`: 0-based insertion position(s) in the reference chromosome where `in_fasta` should be inserted. Use `-1` to insert at the end of the chromosome. For **multiple edits**, provide a comma-separated list (e.g., `0,5,-1`). **Note:** Either `position`, or both `upstream_fasta` and `downstream_fasta`, must be provided.

- `upstream_fasta`: Path(s) to FASTA file(s) containing the upstream flanking sequence(s) for insertion. For **multiple edits**, provide a comma-separated list (e.g., `up1.fa,up2.fa,up3.fa`). Must be used with `downstream_fasta`. Cannot be used together with `position`.

- `downstream_fasta`: Path(s) to FASTA file(s) containing the downstream flanking sequence(s) for insertion. For **multiple edits**, provide a comma-separated list (e.g., `down1.fa,down2.fa,down3.fa`). Must be used with `upstream_fasta`. Cannot be used together with `position`.

- `in_fasta`: Path(s) to FASTA file(s) containing the new sequence(s) to insert. For multiple edits, provide a comma-separated list. **The number of entries must match the number of `position` values or the number of upstream/downstream pairs.**

- `in_gff`: Path(s) to GFF3 file(s) describing the `in_fasta` sequence(s). For multiple edits, provide a comma-separated list. **The number of entries must match the number of `in_fasta` files.**

- `ref_fasta` Path to the reference genome FASTA file.

- `ref_gff` Path to the reference genome annotation (GFF3 or GTF) file.

## Examples

### Single Edit by Position

```
python3 reform.py \
  --chrom="I" \
  --position=1500 \
  --in_fasta="data/edit.fa" \
  --in_gff="data/edit.gff" \
  --ref_fasta="data/ref.fa" \
  --ref_gff="data/ref.gff3"
```

### Single Edit with Upstream/Downstream Flanks

```
python3 reform.py \
  --chrom="I" \
  --upstream_fasta="data/up.fa" \
  --downstream_fasta="data/down.fa" \
  --in_fasta="data/edit.fa" \
  --in_gff="data/edit.gff" \
  --ref_fasta="data/ref.fa" \
  --ref_gff="data/ref.gff3"
```

### Batch Edits (Multiple Positions)

```
python3 reform.py \
  --chrom="I" \
  --position=1000,2500,3000 \
  --in_fasta="data/edit1.fa,data/edit2.fa,data/edit3.fa" \
  --in_gff="data/edit1.gff,data/edit2.gff,data/edit3.gff" \
  --ref_fasta="data/ref.fa" \
  --ref_gff="data/ref.gff3"
```

### Append a Novel Chromosome

```
python3 reform.py \
  --new_chrom="new_chr1" \
  --in_fasta="data/new1.fa" \
  --in_gff="data/new1.gff" \
  --ref_fasta="data/ref.fa" \
  --ref_gff="data/ref.gff3"
```

## Output

`reformed.fa` Modified fasta file.

`reformed.gff3` Modified GFF file.

## Tests
After local deployment or modification, you can run `test_reform.py` to verify the functionality of *ref*orm. This script contains an automated test suite built with Python’s `unittest` framework and validates *ref*orm across a range of genome editing scenarios.

To run all tests:

```bash
python3 test_reform.py
```

## How to Cite

If you use *ref*orm in your research, please cite the GitHub repository:

> *ref*orm: https://github.com/gencorefacility/reform

You may also cite our article:

> Mohammed Khalfan\*, Eric Borenstein, Pieter Spealman, Farah Abdul-Rahman, and David Gresham\* (2021).  
> *Modifying Reference Sequence and Annotation Files Quickly and Reproducibly with reform*.

This article was published — to our knowledge — as the first scientific research article minted as a non-fungible token (NFT):

- 📄 [Read the full PDF](https://gencore.bio.nyu.edu/wp-content/uploads/2021/07/reform.pdf)  
- 🖼️ [View NFT on OpenSea](https://opensea.io/item/ethereum/0x495f947276749ce646f68ac8c248420045cb7b5e/89295771465272658208657695219245348516590738176651091797615877953749424013313)


