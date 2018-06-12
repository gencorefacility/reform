from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from bcbb import gff

chrom = "I"
inserted_seq = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
inserted_seq_name = "insertedSeqNameTest"
left_strand = "CATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTT"
right_strand = "ACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAAC"
ref_fasta = "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
ref_gff = "Saccharomyces_cerevisiae.R64-1-1.34.gff3"

def main():
	## Generate index of sequences from reference fasta
	chrom_seqs = SeqIO.index(ref_fasta,'fasta')
	
	## Obtain the sequence of the chromosome we wish to modify
	seq = chrom_seqs[chrom]
	seq_str = str(seq.seq)
	
	## Obtain the starting position of the left_strand
	index = seq_str.find(left_strand + right_strand)
	
	## Build the new chromosome sequence with the inserted_seq 
	new_seq = seq_str[:index + len(left_strand)] + inserted_seq + seq_str[index + (len(left_strand + right_strand)):]
	new_record = SeqRecord(Seq(new_seq, generic_dna), id=seq.id, description=seq.description)
	
	## Create new fasta file with modified chromosome 
	with open(ref_fasta.replace('.fa', '_reform.fa'), "w") as f:
		for s in chrom_seqs:
			if s == seq.id:
				SeqIO.write([new_record], f, "fasta")
			else:
				SeqIO.write([chrom_seqs[s]], f, "fasta")

	## Make new gff record for inserted_seq
	new_record = SeqRecord(Seq(inserted_seq, generic_dna), seq.id)
	qualifiers = {"source": "Custom", "ID": inserted_seq_name}
	top_feature = SeqFeature(FeatureLocation(index + len(left_strand), index + len(left_strand) + len(inserted_seq)), type="gene", strand=1, qualifiers=qualifiers)
	new_record.features = [top_feature]

	## Create new gff file
	gff_out = open(ref_gff.replace('.gff', '_reform.gff'), "w")
	new_record_created = False
	with open(ref_gff, "r") as f:
		for rec in gff.parse(f, target_lines=1):
			if rec.gff_id != seq_id or rec.start < index + len(left_strand):
				gff.write([rec], gff_out)
			else:
				if not new_record_created:
					gff.write([new_record], gff_out)
					new_record_created = True
				modified_record = rec
				modified_record.start = rec.start + len(inserted_seq)
				modified_record.end = rec.end + len(inserted_seq)
				gff.write([rec], gff_out)
	gff_out.close()

#call to main function to run the program
if __name__ == "__main__":
    main()

