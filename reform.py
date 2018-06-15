#module load biopython/intel/1.70
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

#chrom = "I"
#inserted_fasta = "../reform_files/new.fa"
#inserted_gff = "../reform_files/new.gff"
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
	new_fasta = in_arg.ref_fasta.replace('.fa', '_reform.fa')
	with open(new_fasta, "w") as f:
		for s in chrom_seqs:
			if s == seq.id:
				SeqIO.write([new_record], f, "fasta")
			else:
				SeqIO.write([chrom_seqs[s]], f, "fasta")
				
	print("New fasta file created: ", new_fasta)
	print("Preparing to create new GFF file")
	
	## Create new gff file
	with open(in_arg.in_gff, "r") as f:
		in_gff_lines = []
		for line in f:
			# convert spaces to tabs
			line = re.sub("\s\s+" , "\t", line)
			line_elements = line.split('\t')
			if len(line_elements) != 9:
				print("** ERROR: GFF file does not have 9 columns, it has", len(line_elements))
				print(line_elements)
				exit()
			in_gff_lines.append(line_elements)
	
	gff_out = open(in_arg.ref_gff.replace('.gff', '_reform.gff'), "w")
	in_gff_lines_appended = False
	with open(in_arg.ref_gff, "r") as f:
		for line in f:
			if line.startswith("#"):
				## TODO: need to update the chromosome coordinates too (to include new sequence in chrom)
				## both in the header, and first row/feature of that chrom
				gff_out.write(line)
			else:
				line = re.sub("\s\s+" , "\t", line)
				line_elements = line.split('\t')
				gff_chrom_id = line_elements[0]
				gff_feat_start = int(line_elements[3])
				if gff_chrom_id != seq.id or gff_feat_start < position:
					gff_out.write(line)
				else:
					if not in_gff_lines_appended:
						for l in in_gff_lines:
							gff_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(l[0], l[1], l[2], int(l[3]) + position - 1, int(l[4]) + position - 1, l[5], l[6], l[7], l[8]))
						in_gff_lines_appended = True
					gff_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(line_elements[0], line_elements[1], line_elements[2], int(line_elements[3]) + len(str(record.seq)), int(line_elements[4]) + len(str(record.seq)), line_elements[5], line_elements[6], line_elements[7], line_elements[8]))
	gff_out.close()			
	print("New GFF file created: ", gff_out.name)
	
def get_input_args():
	parser = argparse.ArgumentParser()
	
	parser.add_argument('--chrom', type = str, required = True,
					help = "Chromosome name (String)") 
	parser.add_argument('--in_fasta', type = str, required = True,
					help = "Path to new sequence to be inserted into reference genome in fasta format") 
	parser.add_argument('--in_gff', type = str, required = True,
					help = "Path to GFF file describing new fasta sequence to be inserted") 
	parser.add_argument('--upstream', type = str, default = None, 
					help = "Upstream sequence. Either position, or upstream AND downstream sequence must be provided.")
	parser.add_argument('--downstream', type = str, default = None, 
					help = "Downstream sequence. Either position, or upstream AND downstream sequence must be provided.")
	parser.add_argument('--position', type = int, default = None,
					help = "Position at which to insert new sequence. Either position, or upstream AND downstream sequence must be provided.") 
	parser.add_argument('--ref_fasta', type = str, required = True,
					help = "Path to reference fasta file")
	parser.add_argument('--ref_gff', type = str, required = True,
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

