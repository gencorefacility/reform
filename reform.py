#module load biopython/intel/python3.6/1.72
#python3 reform.py --chrom="I" --upstream_fasta="data/up.fa" --downstream_fasta="data/down.fa" --in_fasta="data/new.fa" --in_gff="data/new.gff" --ref_fasta="data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa" --ref_gff="data/Saccharomyces_cerevisiae.R64-1-1.34.gff3"

import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

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
	positions = get_position(in_arg.position, in_arg.upstream_fasta, in_arg.downstream_fasta, in_arg.chrom, seq_str)
	position = positions['position']
	down_position = positions['down_position']
	if down_position < position:
		print("** ERROR: Upstream and Downstream sequences must not overlap. Exiting.")
		exit()
	if position != down_position:
		print("Removing nucleotodes from position {} - {}".format(position, down_position - 1))
	print("Proceeding to insert sequence '{}' from {} at position {}".format(record.description, in_arg.in_fasta, position))
	
	## Build the new chromosome sequence with the inserted_seq 
	## If the chromosome sequence length is in the header, replace it with new length
	new_seq = seq_str[:position] + str(record.seq) + seq_str[down_position:]
	chrom_length = str(len(seq_str))
	new_length = str(len(new_seq))
	new_record = SeqRecord(Seq(new_seq, generic_dna), id=seq.id, description=seq.description.replace(chrom_length, new_length))
	
	## Create new fasta file with modified chromosome 
	new_fasta = 'reformed.fa'
	with open(new_fasta, "w") as f:
		for s in chrom_seqs:
			if s == seq.id:
				SeqIO.write([new_record], f, "fasta")
			else:
				SeqIO.write([chrom_seqs[s]], f, "fasta")
				
	print("New fasta file created: ", new_fasta)
	print("Preparing to create new GFF file")
	
	## Read in new GFF features from in_gff
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
	
	## Create new gff file
	new_gff = 'reformed.gff'
	gff_out = open(new_gff, "w")
	in_gff_lines_appended = False
	split_features = []
	with open(in_arg.ref_gff, "r") as f:
		for line in f:
			if line.startswith("#"):
				line_elements = line.split()
				if line_elements[0] == "##sequence-region" and line_elements[1] == seq.id:
					original_length = int(line_elements[3])
					new_length = original_length - (down_position - position) + len(str(record.seq))
					line = line.replace(str(original_length), str(new_length))
				gff_out.write(line)
			else:
				line = re.sub("\s\s+" , "\t", line)
				line_elements = line.split('\t')
				gff_chrom_id = line_elements[0]
				gff_feat_start = int(line_elements[3])
				gff_feat_end = int(line_elements[4])
				
				if gff_chrom_id != seq.id or gff_feat_end < position:
					gff_out.write(line)
				elif line_elements[2] in ['chromosome', 'region']:
					# modify length of chromosome
					original_length = int(line_elements[4])
					new_length = original_length - (down_position - position) + len(str(record.seq))
					line = line.replace(str(original_length), str(new_length))
					gff_out.write(line)
				elif gff_feat_start < position and gff_feat_end > down_position:
					print("Feature split")
					# split feature into 2
					gff_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(line_elements[0], line_elements[1], line_elements[2], int(line_elements[3]), position, line_elements[5], line_elements[6], line_elements[7], line_elements[8] + ";reform_comment=original feature split by inserted sequence, this is the 5' end"))
					# the downstream flank(s) will be added immediately after the in_gff features are written
					renamed_id_attributes = rename_id(line)
					split_features.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(line_elements[0], line_elements[1], line_elements[2], position + len(str(record.seq)) + 1, int(line_elements[4]) + len(str(record.seq)), line_elements[5], line_elements[6], line_elements[7], renamed_id_attributes + ";reform_comment=original feature split by inserted sequence, this is the 3' end"))
				elif gff_feat_start < position and gff_feat_end > position and gff_feat_end < down_position:
					print("Feature cut off - 3' side (downstream side) of feature cut off")
					# change end position of feature to cut off point (position)
					gff_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(line_elements[0], line_elements[1], line_elements[2], int(line_elements[3]), position, line_elements[5], line_elements[6], line_elements[7], line_elements[8] + ";reform_comment=3' side of feature cut-off by inserted sequence"))
				elif gff_feat_start > position and gff_feat_end < down_position:
					print("Skip feature (removed from sequence)")
					# skip this feature
					continue
				else:
					if not in_gff_lines_appended:
						for l in in_gff_lines:
							gff_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(l[0], l[1], l[2], int(l[3]) + position, int(l[4]) + position, l[5], l[6], l[7], l[8]))
						in_gff_lines_appended = True
						for sf in split_features:
							gff_out.write(sf)
					if gff_feat_start > position and gff_feat_start < down_position and gff_feat_end > down_position:
						print("Feature cut off - 5' side (upstream side) of feature cut off")
						# change start position of feature to after cutoff point
						gff_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(line_elements[0], line_elements[1], line_elements[2], position + len(str(record.seq)) + 1, int(line_elements[4]) + len(str(record.seq)) - (down_position - position), line_elements[5], line_elements[6], line_elements[7], line_elements[8] + ";reform_comment=5' side of feature cut-off by inserted sequence"))
					elif gff_feat_start > down_position:
						gff_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(line_elements[0], line_elements[1], line_elements[2], int(line_elements[3]) + len(str(record.seq)) - (down_position - position), int(line_elements[4]) + len(str(record.seq)) - (down_position - position), line_elements[5], line_elements[6], line_elements[7], line_elements[8]))
					else:
						print("** Error: Unknown case for GFF modification. Exiting")
						exit()
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
	parser.add_argument('--upstream_fasta', type = str, default = None, 
					help = "Path to Fasta file with upstream sequence. Either position, or upstream AND downstream sequence must be provided.")
	parser.add_argument('--downstream_fasta', type = str, default = None, 
					help = "Path to Fasta file with downstream sequence. Either position, or upstream AND downstream sequence must be provided.")
	parser.add_argument('--position', type = int, default = None,
					help = "Position at which to insert new sequence. Either position, or upstream AND downstream sequence must be provided.") 
	parser.add_argument('--ref_fasta', type = str, required = True,
					help = "Path to reference fasta file")
	parser.add_argument('--ref_gff', type = str, required = True,
					help = "Path to reference gff file") 
					
	in_args = parser.parse_args()
	if in_args.position is None and (in_args.upstream_fasta is None or in_args.downstream_fasta is None):
		print("** Error: You must provide either the position, or the upstream and downstream sequences.")
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
			down_position = position
			print("Position valid")
	else:
		print("No valid position specified, checking for upstream and downstream sequence")
		if upstream is not None and downstream is not None:
			seq_str = seq_str.upper()
			upstream_fasta = list(SeqIO.parse(upstream, "fasta"))
			upstream_seq = str(upstream_fasta[0].seq).upper()
			downstream_fasta = list(SeqIO.parse(downstream, "fasta"))
			downstream_seq = str(downstream_fasta[0].seq).upper()
			# Ensure the upstream and downstream target sequences exists
			# once in the selected chromosome, else die
			upstream_seq_count = seq_str.count(upstream_seq)
			downstream_seq_count = seq_str.count(downstream_seq)
			if upstream_seq_count == 1 and downstream_seq_count == 1:
				## Obtain the starting position of the left_strand
				index = seq_str.find(upstream_seq)
				position = index + len(upstream_seq)
				down_position = seq_str.find(downstream_seq)
			else:
				print("** ERROR: The upstream and downstream target sequences must be present and unique in the specified chromosome.")
				print("Chromosome: {}\n".format(chrom))
				print("Upstream sequence found {} times".format(upstream_seq_count))
				print("Downstream sequence found {} times".format(downstream_seq_count))
				exit()
		else:
			print("** ERROR: You must specify a valid position or upstream and downstream sequences.")
			exit()
	return {'position': position, 'down_position': down_position}
	
def rename_id(line):
	attributes = line.split('\t')[8]
	elements = attributes.split(';')
	if elements[0].startswith("ID="):
		print("Renaming split feature {} --> {}, {}_split".format(elements[0], elements[0], elements[0]))
		return ("{}_split;{}".format(elements[0], ';'.join(elements[1:])))
	else:
		print("This feature will not be renamed because it does not has an ID attribute:\n", line)
		return attributes
	
#call to main function to run the program
if __name__ == "__main__":
	main()

