#!/bin/env python
import argparse
import re
import os
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
try:
    import pgzip as gzip_module
    print("Using pgzip for gzip operations.")
except ImportError:
    import gzip as gzip_module
    print("pgzip not found, falling back to gzip.")


def main():
	## Retrieve command line arguments and number of iterations
	in_arg, iterations = get_input_args()

	## List for previous postion and modification length
	prev_modifications = []

	## Path for the files generated in sequential processing, mainly managed by tempfile.
	prev_fasta_path = None
	prev_gff_path = None

	## Sequential processing
	for index in range(iterations):
		## Read the new fasta (to be inserted into the ref genome)
		record = list(SeqIO.parse(in_arg.in_fasta[index], "fasta"))[0]
		
		## Generate index of sequences from ref reference fasta
		if prev_fasta_path:
			chrom_seqs = index_fasta(prev_fasta_path)
			os.remove(prev_fasta_path)
		else:	
			chrom_seqs = index_fasta(in_arg.ref_fasta)
		
		## Obtain the sequence of the chromosome we want to modify
		seq = chrom_seqs[in_arg.chrom]
		seq_str = str(seq.seq)
		
		## Get the position to insert the new sequence
		positions = get_position(index, in_arg.position, in_arg.upstream_fasta, in_arg.downstream_fasta, in_arg.chrom, seq_str, prev_modifications)
		position = positions['position']
		down_position = positions['down_position']
		
		## Save current modification which include position(index) and length changed.
		if position == down_position:
			length_changed = len(str(record.seq))
		else:
			length_changed = len(str(record.seq)) - (down_position - position - 1)
		prev_modifications.append((position,length_changed))
		
		if position != down_position:
			print("Removing nucleotides from position {} - {}".format(position, down_position))
		print("Proceeding to insert sequence '{}' from {} at position {} on chromsome {}"
			.format(record.description, in_arg.in_fasta[index], position, in_arg.chrom))
		
		## Build the new chromosome sequence with the inserted_seq 
		## If the chromosome sequence length is in the header, replace it with new length
		new_seq = seq_str[:position] + str(record.seq) + seq_str[down_position:]
		chrom_length = str(len(seq_str))
		new_length = str(len(new_seq))
		new_record = SeqRecord(
			Seq(new_seq), 
			id=seq.id, 
			description=seq.description.replace(chrom_length, new_length)
		)
		
		## Create new fasta file with modified chromosome 
		if index < iterations - 1:
			new_fasta_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fa')
			new_fasta_file.close()
			new_fasta = new_fasta_file.name
			prev_fasta_path = new_fasta
		else:
			ref_basename, _ = get_ref_basename(in_arg.ref_fasta)
			ref_name = ref_basename
			new_fasta = ref_name + '_reformed.fa'
		with open(new_fasta, "w") as f:
			for s in chrom_seqs:
				if s == seq.id:
					SeqIO.write([new_record], f, "fasta")
				else:
					SeqIO.write([chrom_seqs[s]], f, "fasta")
		print("New fasta file created: ", new_fasta)
		
		print("Preparing to create new annotation file")
		
		## Read in new GFF features from in_gff
		in_gff_lines = get_in_gff_lines(in_arg.in_gff[index])
		
		## Create a temp file for gff, if index is not equal to last iteration
		annotation_name, annotation_ext = get_ref_basename(in_arg.ref_gff)
		print(in_arg.ref_gff)
		if index < iterations - 1:
			temp_gff = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=annotation_ext)
			temp_gff_name = temp_gff.name
			temp_gff.close()
			if prev_gff_path:
				new_gff_path = create_new_gff(temp_gff_name, prev_gff_path, in_gff_lines, position, down_position, seq.id, len(str(record.seq)))
				os.remove(prev_gff_path)
			else:
				new_gff_path = create_new_gff(temp_gff_name, in_arg.ref_gff, in_gff_lines, position, down_position, seq.id, len(str(record.seq)))
		else:
			new_gff_name = annotation_name + '_reformed' + annotation_ext
			if prev_gff_path: 
				new_gff_path = create_new_gff(new_gff_name, prev_gff_path, in_gff_lines, position, down_position, seq.id, len(str(record.seq)))
			else:
				new_gff_path = create_new_gff(new_gff_name, in_arg.ref_gff, in_gff_lines, position, down_position, seq.id, len(str(record.seq)))
		prev_gff_path = new_gff_path
		print("New {} file created: {} ".format(annotation_ext.upper(), prev_gff_path))


def index_fasta(fasta_path):
	'''
	Process incoming FASTA files, supporting both decompressed and uncompressed 
	formats. It takes a file path (fasta_path), decompresses the file into a 
	temporary file if it is a compressed file ending in .gz, and then indexes 
	the temporary file. Finally, the indexing result is returned.
	'''
	if fasta_path.endswith('.gz'):
		# Create a tempfile to store uncompressde content
		with tempfile.NamedTemporaryFile(delete=False, mode='w') as tmp_f:
			tmp_f_path = tmp_f.name
			# Use pgzip or gzip to decompress parallely. Set thread=None means use all cores
			with gzip_module.open(fasta_path, 'rt', thread=None) as f:
				tmp_f.write(f.read())
		chrom_seqs = SeqIO.index(tmp_f_path, 'fasta')
        	# remove temp file
		os.remove(tmp_f_path)
	else:
		chrom_seqs = SeqIO.index(fasta_path, 'fasta')
	return chrom_seqs

def get_ref_basename(filepath):
	'''
	Takes a filepath and returns the filename and extension minus the .gz.
	'''
	base = os.path.basename(filepath)
	if base.endswith('.gz'):
		base = base[:-3]  # remove .gz
	name, ext = os.path.splitext(base)
	return name, ext


def modify_gff_line(elements, start=None, end=None, comment=None):
	'''
	Modifies an existing GFF line and returns the modified line. Currently, you can 
	override the start position, end position, and comment column. 
	
	gff_out: an open file handle for writing new gff lines to
	elements: a list containing each column (9 in total) of a single feature line in a gff file
	start: if provided, overwrite the start position in elements with this start position
	end: if provided, overwrite the end position in elements with this end position
	comment: if provided, overwrite the comments column in elements with this comment
	'''
	if start == None: 
		start = int(elements[3])
	if end == None:
		end = int(elements[4])
	if comment == None:
		comment = elements[8]
	if not comment.endswith('\n'):
		comment += '\n'

	## Return the modified line
	return("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(elements[0], elements[1], elements[2], start, end, elements[5], elements[6], elements[7], comment))
	
def get_in_gff_lines(in_gff):
	'''
	Takes a gff file and returns a list of lists where 
	each parent list item is a single line of the gff file
	and the child elements are the columns of the line
	'''
	with open(in_gff, "r") as f:
		in_gff_lines = []
		for line in f:
			line_elements = line.split('\t')
			# Skip comment lines
			if line.startswith("#"):
				continue
			# Ensure lines have 9 columns
			if len(line_elements) != 9:
				print("** ERROR: in_gff file does not have 9 columns, it has", len(line_elements))
				print(line_elements)
				exit()
			in_gff_lines.append(line_elements)
	return in_gff_lines
	
def get_position(index, positions, upstream, downstream, chrom, seq_str, prev_modifications):
	''' 
	Determine the position in seq_str to insert the new sequence given 
	the position, upstream, downstream, and chrom arguments.
	Note that either position, or upstream AND downstream sequences must
	be provided.
	'''
	if positions and index < len(positions) and positions[index] >= -1:
		print("Checking position validity")
		position = positions[index]
		## Update current postion based on previous modification
		for pos, lc in prev_modifications:
			## Also ignore when postion == -1
			if position >= pos:
				position += lc
		if position > len(seq_str):
			print("** ERROR: Position greater than length of chromosome.")
			print("Chromosome: {}\Chromosome length: {}\nPosition: \n{}".format(chrom, len(seq_str), position))
			exit()
		elif position == -1:
			position = len(seq_str)
		# At the moment we don't accept down position as a param
		# so set down_position = position
		down_position = position
		print("Position valid")
	else:
		print("No valid position specified, checking for upstream and downstream sequence")
		if index < len(upstream) and index < len(downstream):
			seq_str = seq_str.upper()
			upstream_fasta = list(SeqIO.parse(upstream[index], "fasta"))
			upstream_seq = str(upstream_fasta[0].seq).upper()
			downstream_fasta = list(SeqIO.parse(downstream[index], "fasta"))
			downstream_seq = str(downstream_fasta[0].seq).upper()
			# Ensure the upstream and downstream target sequences exists once in the selected chromosome, else die
			upstream_seq_count = seq_str.count(upstream_seq)
			downstream_seq_count = seq_str.count(downstream_seq)
			### TODO: Update postion based on previous modifications
			if upstream_seq_count == 1 and downstream_seq_count == 1:
				## Obtain the starting position of the left_strand
				new_index = seq_str.find(upstream_seq)
				position = new_index + len(upstream_seq)
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
	if down_position < position:
		print("** ERROR: Upstream and Downstream sequences must not overlap. Exiting.")
		exit()
	return {'position': position, 'down_position': down_position}

def write_in_gff_lines(gff_out, in_gff_lines, position, split_features, sequence_length):
	'''
	in_gff_lines: a list of lists where each nested list is a list of 
		columns (in gff format) associated with each new feature to insert
	split_features: contains information about features in the original GFF 
		file that were split due to the insertion of the new sequence.
	sequence_length: length of the inserted sequence, used to determine 
		the new end positions in the GFF file.
	'''
  for l in in_gff_lines:
		# Replace the original chromosome ID from in_gtf with chrom (seq.id)
		l[0] = chrom
  # Handling of single-line comments
	if len(in_gff_lines) == 1:
		l = in_gff_lines[0]
		## Check length
		## l[3] is start position of fasta in in.gtf and l[4] is end position
		seq_id = l[0]
		if int(l[4]) - int(l[3]) + 1 != sequence_length:
			print(f"** WARNING: Inconsistent length for {seq_id}. Correcting start position to 1 and end position to {sequence_length}.")
		## Correct start(l[3]) to 1 and end(l[4]) to length of insert fasta
		new_gff_line = modify_gff_line(
			l, start=1 + position, end=sequence_length + position)
		gff_out.write(new_gff_line)
	# Handling of multiple-line comments
	else:
		### Step1: extract all start and end into corresponding set()
		start_positions = set()
		end_positions = set()
		for l in in_gff_lines:
			start_positions.add(int(l[3]))
			end_positions.add(int(l[4]))
		### Step2: Find min start and max end
		min_start = min(start_positions)
		max_end = max(end_positions)
		### Step3: Check sequence length, validness of min_start and max_end
		if max_end - min_start + 1 != sequence_length:
			raise ValueError(f"Error: Annotation length does not match sequence length. "
							f"Expected {sequence_length}, but got {max_end - min_start + 1}")
		if min_start < 1:
			raise ValueError(f"Error: Invalid min_start value: {min_start}. It must be a positive integer.")
		if min_start > max_end:
			raise ValueError(f"Error: min_start ({min_start}) cannot be greater than max_end ({max_end}).")
		### Step4: Adjust start and end. Offset will be 0 if no adjust need
		offset = min_start - 1
		for l in in_gff_lines:
			### Step5: Correct start(l[3]) and end(l[4]) by minus offset
			new_gff_line = modify_gff_line(
				l, start = int(l[3]) - offset + position, end = int(l[4]) - offset + position)
			gff_out.write(new_gff_line)		

	## If insertion caused any existing features to be split, add
	## the split features now immediately after adding the new features
	for sf in split_features:
		modified_line = modify_gff_line(
			sf[0], start = sf[1], end = sf[2], comment = sf[3])
		gff_out.write(modified_line)
		
	## Return True after writing the new GFF lines
	return True

def create_new_gff(new_gff_name, ref_gff, in_gff_lines, position, down_position, chrom_id, new_seq_length):
	'''
	Goes line by line through a gff file to remove existing features 
	(or parts of existing features) and/or insert new features. 
	In the process of adding new features, existing features 
	may be cut-off on either end or split into two. 
	This attempts to handle all cases.
	new_gff_name: the name of the new gff file to create
	ref_gff: the reference gff file to modify
	in_gff_lines: a list of lists where each nested list is a list of 
		columns (in gff format) associated with each new feature to insert
	position: start position of removal of existing sequence
	down_position: end position of removal of existing sequence
	chrom_id: the ID of the chromosome to modify, which is seq.id
	new_seq_length: the length of the new sequence being added to the chromosome
	'''
	with open(new_gff_name, "w") as gff_out:
		in_gff_lines_appended = False
		split_features = []
		last_seen_chrom_id = None
		gff_ext = new_gff_name.split('.')[-1]
		ref_gff_path = ref_gff
		if ref_gff.endswith('.gz'):
			with gzip_module.open(ref_gff, 'rt') as f:
						# Create a tempfile to store uncompressde content
				with tempfile.NamedTemporaryFile(delete=False, mode='w') as tmp_f:
					tmp_f.write(f.read())
					ref_gff_path = tmp_f.name
		with open(ref_gff_path, "r") as f:
			for line in f:
				line_elements = line.split('\t')
				if line.startswith("#"):
					if line_elements[0] == "##sequence-region" and line_elements[1] == chrom_id:
						## Edit the length of the chromosome 
						original_length = int(line_elements[3])
						new_length = original_length - (down_position - position) + new_seq_length
						line = line.replace(str(original_length), str(new_length))
					gff_out.write(line)
				else:
					gff_chrom_id = line_elements[0]
					gff_feat_start = int(line_elements[3])
					gff_feat_end = int(line_elements[4])
					gff_feat_type = line_elements[2]
					gff_feat_strand = line_elements[6]
					gff_comments = line_elements[8].strip()
					
					# If we've seen at least one chromosome
					# and the last chromosome seen was the chromosome of interest (i.e. chrom_id)
					# and now we're on a new chromosome in the original gff file
					# and the in_gff_lines have not yet been appended:
					# we assume they need to be appended to the end of the chromosome
					# so append them before proceeding
					if (last_seen_chrom_id is not None
						and last_seen_chrom_id == chrom_id
						and gff_chrom_id != last_seen_chrom_id 
						and not in_gff_lines_appended):
						in_gff_lines_appended = write_in_gff_lines(
							gff_out, in_gff_lines, position, split_features, new_seq_length)
					
					last_seen_chrom_id = gff_chrom_id
					
					# If this is not the chromosome of interest
					# Or if the current feature ends before any 
					# modification (which occurs at position)
					# Then simply write the feature as is (no modification)
					# (remember gff co-ordinates are 1 based and position 
					# is 0 based, therefor check less than *or equal to*)
					if gff_chrom_id != chrom_id or gff_feat_end <= position:
						gff_out.write(line)
						
					# Modify chromosome feature length (if feature is 
					# "chromosome" or "region")
					elif gff_feat_type in ['chromosome', 'region']:
						original_length = gff_feat_end
						new_length = original_length - (down_position - position) + new_seq_length
						line = line.replace(str(original_length), str(new_length))
						gff_out.write(line)
						
					# Split feature into 2 if feature starts before position
					# and ends after down_position
					elif gff_feat_start <= position and gff_feat_end > down_position:
						print("Feature split")
						# Which side of the feature depends on the strand (we add this as a comment)
						(x, y) = ("5", "3") if gff_feat_strand == "+" else ("3", "5")
						
						new_comment = format_comment(
							"original feature split by inserted sequence, this is the {} prime end".format(x),
							gff_ext
						)
						# Modified feature ends at 'position'
						modified_line = modify_gff_line( 
							line_elements, 
							end = position, 
							comment = gff_comments + new_comment
						)
						gff_out.write(modified_line)
						
						# The downstream flank(s) will be added immediately after the 
						# in_gff (new) features are written.
						# First, attempt to rename IDs (to indicate split)
						renamed_id_attributes = rename_id(line)
						new_comment = format_comment(
							"original feature split by inserted sequence, this is the {} prime end".format(y),
							gff_ext
						)
						split_features.append(
							(
								line_elements, 
								position + new_seq_length + 1, 
								gff_feat_end + new_seq_length - (down_position - position), 
								renamed_id_attributes + new_comment 
							)
						)
					
					# Change end position of feature to cut off point (position) if the
					# feature ends within the deletion (between position & down_position)
					elif gff_feat_start <= position and gff_feat_end <= down_position:
						# Which side of the feature depends on the strand (we add this as a comment)
						x = "3" if gff_feat_strand == "+" else "5"
						print("Feature cut off - {} prime side of feature cut off ({} strand)"
							.format(x, gff_feat_strand))
						new_comment = format_comment(
							"{} prime side of feature cut-off by inserted sequence".format(x),
							gff_ext
						)
						modified_line = modify_gff_line(
							line_elements, 
							end = position, 
							comment = gff_comments + new_comment 
						)
						gff_out.write(modified_line)
					
					# Skip this feature if it falls entirely within the deletion 
					elif gff_feat_start > position and gff_feat_end <= down_position:
						print("Skip feature (this feature was removed from sequence)")
						continue
						
					else:
						if not in_gff_lines_appended:
							in_gff_lines_appended = write_in_gff_lines(
								gff_out, in_gff_lines, position, split_features, new_seq_length)
							
						# Change start position of feature to after cutoff point if
						# the feature starts within the deletion
						if (gff_feat_start > position 
							and gff_feat_start <= down_position 
							and gff_feat_end > down_position):
							x = "5" if gff_feat_strand == "+" else "3"
							print("Feature cut off - {} prime side of feature cut off ({} strand)"
								.format(x, gff_feat_strand))
							new_comment = format_comment(
								"{} prime side of feature cut-off by inserted sequence".format(x),
								gff_ext
							)
							modified_line = modify_gff_line(
								line_elements, 
								start = position + new_seq_length + 1, 
								end = gff_feat_end + new_seq_length - (down_position - position), 
								comment = gff_comments + new_comment
							)
							gff_out.write(modified_line)
							
						# Offset all downstream feature positions by offset length
						elif gff_feat_start > down_position:
							offset_length = new_seq_length - (down_position - position)
							modified_line = modify_gff_line(
								line_elements, 
								start = gff_feat_start + offset_length, 
								end = gff_feat_end + offset_length
							)
							gff_out.write(modified_line)
							
						else:
							print("** Error: Unknown case for GFF modification. Exiting " + str(line_elements))
							exit()
							
			# If we've iterated over the entire original gff
			# (i.e. out of the for loop which iterates over it)
			# and still haven't written in_gff_lines
			# and the last chromosome seen was the chromosome of interest (i.e. chrom_id):
			# we assume they need to be appended to the end of the genome
			# (i.e. end of the last chromosome)
			# so append them now
			if (last_seen_chrom_id is not None 
				and last_seen_chrom_id == chrom_id
				and not in_gff_lines_appended):
				in_gff_lines_appended = write_in_gff_lines(
					gff_out, in_gff_lines, position, split_features, new_seq_length)
			
			# Checking to ensure in_gff_lines written
			if not in_gff_lines_appended:
				print("** Error: Something went wrong, in_gff not added to reference gff. Exiting")
				exit()
		# Remove temp gtf file
		if ref_gff.endswith('.gz'):
			os.remove(ref_gff_path)
	return new_gff_name

def format_comment(comment, ext):
	'''
	Format comment according to ext (GFF or GTF) and return
	'''
	if ext.lower() == 'gtf':
		new_comment = 'reform_comment "{}";'.format(comment)
	elif ext.lower().startswith('gff'):
		new_comment = "; reform_comment={}".format(comment)
	else:
		print("** Error: Unrecognized extension {} in format_comment(). Exiting".format(ext))
		exit()
	return new_comment
	
def rename_id(line):
	'''
	Given a gff or gtf  line, this function will append the string "_split" 
	to the end of the ID/gene_id attribute
	'''
	attributes = line.split('\t')[8].strip()
	elements = attributes.split(';')
	if elements[0].startswith("ID="):
		print("Renaming split feature {} --> {}_split".format(elements[0], elements[0]))
		return ("{}_split;{}".format(elements[0], ';'.join(elements[1:])))
	elif elements[0].startswith("gene_id "):
		gene_id = re.match(r'gene_id \"(.+)\"', elements[0])[1]
		print('Renaming split feature {} --> {}_split'.format(gene_id, gene_id))
		return ('gene _id "{}_split";{}'.format(gene_id, ';'.join(elements[1:])))

	else:
		print("This feature will not be renamed because it does not has an ID/gene_id attribute:\n", line)
		return attributes
		
def get_input_args():
	parser = argparse.ArgumentParser()
	
	parser.add_argument('--chrom', type = str, required = True,
					help = "Chromosome name (String)") 
	parser.add_argument('--in_fasta', type=str, required=True,
                    help="Path(s) to new sequence(s) to be inserted into reference genome in fasta format") 
	parser.add_argument('--in_gff', type=str, required=True,
                    help="Path(s) to GFF file(s) describing new fasta sequence(s) to be inserted") 
	parser.add_argument('--upstream_fasta', type=str, default=None, 
                    help="Path(s) to Fasta file(s) with upstream sequence. Either position, or upstream AND downstream sequence must be provided.")
	parser.add_argument('--downstream_fasta', type=str, default=None, 
                    help="Path(s) to Fasta file(s) with downstream sequence. Either position, or upstream AND downstream sequence must be provided.")
	parser.add_argument('--position', type=str, default=None,
                    help="Comma-separated positions at which to insert new sequence. Note: Position is 0-based, no space between each comma. Either position, or upstream AND downstream sequence must be provided.")
	parser.add_argument('--ref_fasta', type = str, required = True,
					help = "Path to reference fasta file")
	parser.add_argument('--ref_gff', type = str, required = True,
					help = "Path to reference gff file") 
					
	in_args = parser.parse_args()
	in_args.in_fasta = in_args.in_fasta.split(',')
	in_args.in_gff = in_args.in_gff.split(',')
	if in_args.upstream_fasta:
		in_args.upstream_fasta = in_args.upstream_fasta.split(',')
	if in_args.downstream_fasta:
		in_args.downstream_fasta = in_args.downstream_fasta.split(',')

	if in_args.position is None and (in_args.upstream_fasta is None or in_args.downstream_fasta is None):
		print("** Error: You must provide either the position, or the upstream and downstream sequences.")
		exit()
	
	if not in_args.position and len(in_args.upstream_fasta) != len(in_args.downstream_fasta):
		print("** Error: The number of upstream_fasta and downstream_fasta files does not match.")
		exit()

	if in_args.position is not None:
		try:
			in_args.position = list(map(int, in_args.position.split(',')))
			iterations = iterations = len(in_args.position)
		except ValueError:
			print("** Error: Position must be a comma-separated list of integers, like 1,5,-1.")
			exit()
	else:
		iterations = len(in_args.upstream_fasta)

	if (len(in_args.in_fasta) != len(in_args.in_gff)) or (len(in_args.in_fasta) != iterations):
		print("** Error: The number of inserted FASTA files does not match the number of GTF files, or their counts and positions do not align.")
		exit()
		
	return in_args, iterations
	
if __name__ == "__main__":
	main()
	
