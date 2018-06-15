#module load biopython/intel/1.70
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

#chrom = "I"
#inserted_fasta = "../reform_files/new.fa"
#upstream = "CATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTT"
#downstream = "ACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAAC"
#ref_fasta = "../reform_files/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
#ref_gff = "../reform_files/Saccharomyces_cerevisiae.R64-1-1.34.gff3"
#position = 100

def main():
	## Retrieve command line arguments
	in_arg = get_input_args()

	## Read the new fasta (to be inserted into the ref genome)
	record = list(SeqIO.parse(in_arg.in_fasta, "fasta"))[0]
	
	## Generate index of sequences from ref reference fasta
	chrom_seqs = SeqIO.index(in_arg.ref_fasta,'fasta')
	
	## Obtain the sequence of the chromosome we wish to modify
	seq = chrom_seqs[in_arg.chrom]
	seq_str = str(seq.seq)
	
	## Get the position to insert the new sequence
	position = get_position(in_arg.position, in_arg.upstream, in_arg.downstream, in_arg.chrom, seq_str)
	print("Proceeding to insert {} at position {}".format("inserted_fasta", position))
	
	## Build the new chromosome sequence with the inserted_seq 
	new_seq = seq_str[:position] + str(record.seq) + seq_str[position:]
	new_record = SeqRecord(Seq(new_seq, generic_dna), id=seq.id, description=seq.description)
	
	## Create new fasta file with modified chromosome 
	with open(in_arg.ref_fasta.replace('.fa', '_reform.fa'), "w") as f:
		for s in chrom_seqs:
			if s == seq.id:
				SeqIO.write([new_record], f, "fasta")
			else:
				SeqIO.write([chrom_seqs[s]], f, "fasta")
	
	i='''
	## Make new gff record for inserted_seq
	new_record = SeqRecord(Seq(inserted_seq, generic_dna), seq.id)
	qualifiers = {"source": "Custom", "ID": inserted_seq_name}
	top_feature = SeqFeature(FeatureLocation(index + len(upstream), index + len(upstream) + len(inserted_seq)), type="gene", strand=1, qualifiers=qualifiers)
	new_record.features = [top_feature]
	
	
	## Create new gff file
	gff_out = open(ref_gff.replace('.gff', '_reform.gff'), "w")
	new_record_created = False
	with open(ref_gff, "r") as f:
		for rec in gff.parse(f, target_lines=1):
			if rec.gff_id != seq_id or rec.start < index + len(upstream):
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
	'''
	
def get_input_args():
	parser = argparse.ArgumentParser()
	
	parser.add_argument('--chrom', type = str, required = True,
					help = "Chromosome name (String)") 
	parser.add_argument('--in_fasta', type = str, 
					default = "../reform_files/new.fa", 
					help = "Path to new sequence to be inserted into reference genome in fasta format") 
	parser.add_argument('--upstream', type = str, default = None, 
					help = "Upstream sequence. Either position, or upstream AND downstream sequence must be provided.")
	parser.add_argument('--downstream', type = str, default = None, 
					help = "Downstream sequence. Either position, or upstream AND downstream sequence must be provided.")
	parser.add_argument('--position', type = int, default = None,
					help = "Position at which to insert new sequence. Either position, or upstream AND downstream sequence must be provided.") 
	parser.add_argument('--ref_fasta', type = str, default = "../reform_files/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa",
					help = "Path to reference fasta file")
	parser.add_argument('--ref_gff', type = str, 
					default = "../reform_files/Saccharomyces_cerevisiae.R64-1-1.34.gff3",
					help = "Path to reference gff file") 
					
	in_args = parser.parse_args()
	if in_args.position is None and (in_args.upstream is None or in_args.downstream is None):
		print("** Error: You must provide either the position, or the upstream + downstream sequences.")
		exit()
		
	return in_args

	
def get_position(position, upstream, downstream, chrom, seq_str):
	''' 
	Determine the position in seq_str to insert the new sequence given 
	the position, upstream, downstream, and chrom arguments.
	Note that either position, or upstream AND downstream sequences must
	be provided.
	'''
	if position is not None and position >= 0:
		print("Checking position validity")
		if position > len(seq_str):
				print("** ERROR: Position greater than length of chromosome.")
				print("Chromosome: {}\Chromosome length: {}\nPosition: \n{}".format(chrom, len(seq_str), position))
				exit()
		else:
			print("Position valid")
	else:
		print("No valid position specified, checking for upstream + downstream sequence")
		if upstream is not None and downstream is not None:
			# Ensure the (upstream + downstream) target sequence exists
			# once in the selected chromosome, else die
			count = seq_str.count(upstream + downstream)
			if count == 1:
				## Obtain the starting position of the left_strand
				index = seq_str.find(upstream + downstream)
				position = index + len(upstream)
			else:
				print("** ERROR: The upstream + downstream target sequence must be present and unique in the specified chromosome.")
				print("Chromosome: {}\nNumber of occurrences: {}\nUpstream + Downstream Sequence: \n{}".format(chrom, count, upstream + downstream))
				exit()
		else:
			print("** ERROR: You must specify a valid position or upstream + downstream sequences.")
			exit()
	return position
	
#call to main function to run the program
if __name__ == "__main__":
	main()

