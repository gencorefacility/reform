#!/bin/env python
import argparse
import re
import os
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

## Importing gzip or pgzip module for file compression
print("------------------------------------------")
print(f"Compression Library Use:")
print("------------------------------------------")
try:
	import pgzip as gzip_module
	print(f"Using pgzip for gzip operations.")
except ImportError:
	import gzip as gzip_module
	print(f"pgzip not found, falling back to gzip.")


def main():
	## Retrieve command line arguments and number of iterations
	in_arg, iterations = get_input_args()
	
	## Print reference file paths at start of process
	print("------------------------------------------")
	print(f"Path of Reference Files:")
	print("------------------------------------------")
	print(f"Reference FASTA: {os.path.realpath(in_arg.ref_fasta)}")
	print(f"Reference Annotation: {os.path.realpath(in_arg.ref_gff)}")

	## List for previous postion and modification length
	prev_modifications = []

	## Path for the files generated in sequential processing, mainly managed by tempfile.
	prev_fasta_path = None
	prev_gff_path = None

	## Sequential processing
	for index in range(iterations):
		# Start interation
		if hasattr(in_arg, 'chrom') and in_arg.chrom is not None:
			## Modify existing chrom seq
			print("-------------------------------------------")
			print(f"Begin modification from in{index+1}.fa")
			print("-------------------------------------------")
			new_fasta, annotation_ext, new_gff_path, prev_fasta_path, prev_gff_path = \
				modify_existing_chrom_seq(in_arg, index, prev_fasta_path, prev_modifications, \
				iterations, prev_gff_path)
		else:
			## Add new chrom seq
			print("-------------------------------------------")
			print(f"Begin adding a new chromosome from in{index+1}.fa")
			print("-------------------------------------------")
			new_fasta, annotation_ext, new_gff_path, prev_fasta_path, prev_gff_path = \
				add_new_chrom_seq(in_arg, index, prev_fasta_path, prev_gff_path, iterations)

	print("------------------------------------------")
	print(f"Reform Complete")
	print("------------------------------------------")
	print(f"New .fa file created:  {os.path.realpath(new_fasta)}")
	print(f"New {annotation_ext} file created: {os.path.realpath(new_gff_path)}")

def modify_existing_chrom_seq(in_arg, index, prev_fasta_path, prev_modifications, iterations, prev_gff_path):
	"""
	Modifies a specified and existing chromosome sequence by inserting/replacing a new sequence at a given 
 	position and updates the corresponding FASTA and GFF files.
 	"""
 	## Read FASTA
	record, chrom_seqs= read_fasta(in_arg, index, prev_fasta_path)
	## Read annotation (gff/gtf)
	check_gff(in_arg, index)
	## Obtain the sequence of the chromosome we want to modify
	existing_seq = chrom_seqs[in_arg.chrom]
	existing_seq_str = str(existing_seq.seq)

	## Get the position to insert the new sequence
	positions = get_position(index, in_arg.position, in_arg.upstream_fasta, in_arg.downstream_fasta, in_arg.chrom, existing_seq_str, prev_modifications)
	position = positions['position']
	down_position = positions['down_position']
	## Save current modification which include position(index) and length changed.
	if position == down_position:
		length_changed = len(str(record.seq))
	else:
		length_changed = len(str(record.seq)) - (down_position - position - 1)
	prev_modifications.append((position,length_changed))
	if position != down_position:
		print(f"Removing nucleotides from position {position} - {down_position}")
	print(f"Proceeding to insert sequence '{record.description}' from {in_arg.in_fasta[index]} at position {position} on chromosome {in_arg.chrom}")
	## Build the new chromosome sequence with the inserted_seq 
	## If the chromosome sequence length is in the header, replace it with new length
	new_seq = existing_seq_str[:position] + str(record.seq) + existing_seq_str[down_position:]
	chrom_length = str(len(existing_seq_str))
	new_length = str(len(new_seq))
	new_record = SeqRecord(
		Seq(new_seq), 
		id=existing_seq.id, 
		description=existing_seq.description.replace(chrom_length, new_length)
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
			if s == existing_seq.id:
				SeqIO.write([new_record], f, "fasta")
			else:
				SeqIO.write([chrom_seqs[s]], f, "fasta")
	## Read in new GFF features from in_gff, False means modify existing chrom
	in_gff_lines = get_in_gff_lines(in_gff=in_arg.in_gff[index], existing_chrom=in_arg.chrom, new_chrom=None)
	## Create a temp file for gff, if index is not equal to last iteration
	annotation_name, annotation_ext = get_ref_basename(in_arg.ref_gff)
	if index < iterations - 1:
		temp_gff = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=annotation_ext)
		temp_gff_name = temp_gff.name
		temp_gff.close()
		if prev_gff_path:
			new_gff_path = create_new_gff(temp_gff_name, prev_gff_path, in_gff_lines, position, down_position, existing_seq.id, len(str(record.seq)))
			os.remove(prev_gff_path)
		else:
			new_gff_path = create_new_gff(temp_gff_name, in_arg.ref_gff, in_gff_lines, position, down_position, existing_seq.id, len(str(record.seq)))
	else:
		new_gff_name = annotation_name + '_reformed' + annotation_ext
		if prev_gff_path: 
			new_gff_path = create_new_gff(new_gff_name, prev_gff_path, in_gff_lines, position, down_position, existing_seq.id, len(str(record.seq)))
		else:
			new_gff_path = create_new_gff(new_gff_name, in_arg.ref_gff, in_gff_lines, position, down_position, existing_seq.id, len(str(record.seq)))
	
	prev_gff_path = new_gff_path
	
	return new_fasta, annotation_ext, new_gff_path, prev_fasta_path, prev_gff_path

def add_new_chrom_seq(in_arg, index, prev_fasta_path, prev_gff_path, iterations):
    ## Read FASTA
	record, chrom_seqs= read_fasta(in_arg, index, prev_fasta_path)
	## Read annotation (gff/gtf)
	check_gff(in_arg, index)
	## Check if new chromosome existed in sequence
	if in_arg.new_chrom[index] in chrom_seqs:
		raise ValueError(f"Chromosome {in_arg.new_chrom[index]} already exists in the FASTA file.")
 	## Build the new chromosome sequence by append new sequence below
	## Using new_chrom as the id of the new chromosome
	new_seq = str(record.seq)
	new_seq_length = len(new_seq)
	new_record = SeqRecord(
		Seq(new_seq), 
		id=in_arg.new_chrom[index], # Use new_chrom
		description=f"{in_arg.new_chrom[index]}"
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
	## Write all the original chromosomes first, then new chromosomes.
	with open(new_fasta, "w") as f:
		for s in chrom_seqs:
			SeqIO.write([chrom_seqs[s]], f, "fasta")
		SeqIO.write([new_record], f, "fasta")
	## Read in new GFF features from in_gff
	## Pass the new_chrom name from command line and the length of the new sequence to correct ##sequence-region line
	in_gff_lines = get_in_gff_lines(in_gff=in_arg.in_gff[index], new_chrom=in_arg.new_chrom[index], sequence_length=new_seq_length)
	## Create a temp file for gff, if index is not equal to last iteration
	annotation_name, annotation_ext = get_ref_basename(in_arg.ref_gff)
	if index < iterations - 1:
		temp_gff = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=annotation_ext)
		temp_gff_name = temp_gff.name
		temp_gff.close()
		if prev_gff_path:
			new_gff_path = create_new_gff_for_existing_gff(temp_gff_name, prev_gff_path, in_gff_lines, in_arg.new_chrom[index], new_seq_length)
			os.remove(prev_gff_path)
		else:
			new_gff_path = create_new_gff_for_existing_gff(temp_gff_name, in_arg.ref_gff, in_gff_lines, in_arg.new_chrom[index], new_seq_length)
	else:
		new_gff_name = annotation_name + '_reformed' + annotation_ext
		if prev_gff_path: 
			new_gff_path = create_new_gff_for_existing_gff(new_gff_name, prev_gff_path, in_gff_lines, in_arg.new_chrom[index], new_seq_length)
		else:
			new_gff_path = create_new_gff_for_existing_gff(new_gff_name, in_arg.ref_gff, in_gff_lines, in_arg.new_chrom[index], new_seq_length)
	prev_gff_path = new_gff_path
	return new_fasta, annotation_ext, new_gff_path, prev_fasta_path, prev_gff_path

def read_fasta(in_arg, index, prev_fasta_path):
	"""
	Reads the FASTA file, extracts the sequence to be inserted,
	and indexes the reference genome chromosome sequences.
	"""
	## Read the new fasta (to be inserted into the ref genome)
	try:
		filename_fa = in_arg.in_fasta[index]
		if not os.path.exists(filename_fa):
			raise FileNotFoundError(f"Error: File {filename_fa} does not exist.")
		real_path_fa = os.path.realpath(filename_fa)
		record = list(SeqIO.parse(in_arg.in_fasta[index], "fasta"))[0]
		# Check for mismatch between FASTA record ID and command line chromosome name
		if hasattr(in_arg, 'new_chrom') and in_arg.new_chrom is not None:
			if record.id != in_arg.new_chrom[index]:
				print(f"** WARNING: Mismatch detected between chromosome name in input FASTA ({record.id}) "
                      f"and command line parameter ({in_arg.new_chrom[index]}).")
				print(f"Using command line chromosome name: {in_arg.new_chrom[index]}")
                # The actual override happens in add_new_chrom_seq where a new SeqRecord is created

	except IndexError:
		raise ValueError(f"Error: {filename_fa} is not a valid FASTA file.")
	except Exception as e:
		raise ValueError(f"Error parsing FASTA file: {str(e)}")
	print(f"Preparing to create new FASTA file")
	print(f"Original Input FASTA: {real_path_fa}")
	## Generate index of sequences from ref reference fasta
	if prev_fasta_path:
		chrom_seqs = index_fasta(prev_fasta_path)
		os.remove(prev_fasta_path)
	else:	
		chrom_seqs = index_fasta(in_arg.ref_fasta)
	return record, chrom_seqs
		
def check_gff(in_arg, index):
	"""
	Reads the GFF file, verifies its existence, and prints its file path information.
 	"""
	filename_gff = in_arg.in_gff[index]
	if not os.path.exists(filename_gff):
		raise FileNotFoundError(f"Error: File {filename_gff} does not exist.")
	real_path_gff = os.path.realpath(filename_gff)
	print("Preparing to create new annotation file")
	print(f"Original Input Annotation: {real_path_gff}")
	print() ### print new line

def index_fasta(fasta_path):
	'''
	Process incoming FASTA files, supporting both decompressed and uncompressed 
	formats. It takes a file path (fasta_path), decompresses the file into a 
	temporary file if it is a compressed file ending in .gz, and then indexes 
	the temporary file. Finally, the indexing result is returned.
	'''
	if fasta_path.endswith('.gz'):
		## Create a tempfile to store uncompressde content
		with tempfile.NamedTemporaryFile(delete=False, mode='w') as tmp_f:
			tmp_f_path = tmp_f.name
			## Use pgzip or gzip to decompress parallely. Set thread=None means use all cores
			with gzip_module.open(fasta_path, 'rt', thread=None) as f:
				tmp_f.write(f.read())
		chrom_seqs = SeqIO.index(tmp_f_path, 'fasta')
        ## remove temp file
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

def valid_gff_line(line_elements):
	'''
	Checks if the splited line is a valid GFF line. 
	Returns True if valid, False otherwise.
	'''
	if not line_elements[0].startswith("##sequence-region"):
		if len(line_elements) != 9:
			print(f"** ERROR: in_gff file does not have 9 columns, it has {len(line_elements)}")
			print(line_elements)
			return False
	else:
		## Check if ##sequence-region line has 4 columns, the reason why use 5 here is because last element is
		## spliting format indicator.
		if len(line_elements) != 5:
			print(f"** ERROR: ##sequence-region line does not have 4 columns, it has {len(line_elements) - 1}")
			print(line_elements)
			return False
	return True

def get_in_gff_lines(in_gff=None, existing_chrom=None, new_chrom=None, sequence_length=None):
	'''
	Takes a gff file and returns a list of lists where 
	each parent list item is a single line of the gff file
	and the child elements are the columns of the line
	'''
	with open(in_gff, "r") as f:
		in_gff_lines = []
		for line in f:
			# Skip empty lines
			if not line.strip():
				continue
			
			# Handle differently based on whether we're adding a new chromosome or modifying existing
			if line.startswith("##sequence-region") and new_chrom is not None:
				## Paste ##sequence-region line which only exists in gtf/gff for adding new chromosome.
				## Select user used delimiter based on content
				if '\t' in line:
					line_elements = line.split('\t')
					line_elements.append('\t') ## Add tab to the end as a format indicator
				else:
					line_elements = line.split()
					line_elements.append(' ') ## Add whitespace to the end as a format indicator
				if not valid_gff_line(line_elements):
					exit()
				# Validate new_chrom value and correct if needed
				original_chrom = line_elements[1]
				if original_chrom != new_chrom:
					print(f"** INFO: Updating chromosome name from {original_chrom} to {new_chrom} to fit the\
        				input from new_chrom parameter in command-line.")
					line_elements[1] = new_chrom
				# Validate sequence_length value and correct if needed
				original_start, original_end = line_elements[2], line_elements[3]
				if original_start != "1":
					print(f"** INFO: Updating start position from {original_start} to 1 to fit the\
         				format requirement of annotation file.")
					line_elements[2] = "1"
				if original_end != str(sequence_length):
					print(f"** INFO: Updating sequence length from {original_end} to {sequence_length} to fit the\
        				length of sequence in input FASTA file.")
				line_elements[3] = str(sequence_length)
			elif line.startswith("#"):
				## Ignore other comment lines
				continue
			else:
				## Split, check and add feature lines
				line_elements = line.split('\t')
				chorme_id = existing_chrom if existing_chrom else new_chrom
				if line_elements[0] != chorme_id:
					print("** Warning: The chromosome name in the GFF file does not match the new chromosome name.")
					print(f"Correct the chromosome name {line_elements[0]} to {chorme_id}")
					line_elements[0] = chorme_id
				if not valid_gff_line(line_elements):
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
			print(f"Chromosome: {chrom}\nChromosome length: {len(seq_str)}\nPosition: {position}")
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
			if upstream_seq_count == 1 and downstream_seq_count == 1:
				## Obtain the starting position of the left_strand
				new_index = seq_str.find(upstream_seq)
				position = new_index + len(upstream_seq)
				down_position = seq_str.find(downstream_seq)
			else:
				print("** ERROR: The upstream and downstream target sequences must be present and unique in the specified chromosome.")
				print("Chromosome: {}\n".format(chrom))
				print(f"Upstream sequence found {upstream_seq_count} times")
				print(f"Downstream sequence found {downstream_seq_count} times")
				exit()
		else:
			print("** ERROR: You must specify a valid position or upstream and downstream sequences.")
			exit()
	if down_position < position:
		print("** ERROR: Upstream and Downstream sequences must not overlap. Exiting.")
		exit()
	return {'position': position, 'down_position': down_position}

def calculate_new_length_for_in_gff(in_gff_lines, position, sequence_length):
	'''
	Calculate the new length of the chromosome after modification.
	This function checks the start and end positions of the features in the GFF file
	and adjusts them based on the insertion position and the length of the inserted sequence.
	'''
	## List to store new gff lines
	new_gff_lines = []
	# Handling of single-line comments
	if len(in_gff_lines) == 1:
		## Check if the line isn't a comment line
		if in_gff_lines[0][0].startswith("##sequence-region"):
			new_gff_lines.append(in_gff_lines[0])
		else:
			l = in_gff_lines[0]
			## Check length
			## l[3] is start position of fasta in in.gtf and l[4] is end position
			seq_id = l[0]
			if int(l[4]) - int(l[3]) + 1 != sequence_length:
				print(f"** WARNING: Inconsistent length for {seq_id}. Correcting start position to 1 and end position to {sequence_length}.")
			## Correct start(l[3]) to 1 and end(l[4]) to length of insert fasta
			new_gff_line = modify_gff_line(
				l, start=1 + position, end=sequence_length + position)
			new_gff_lines.append(new_gff_line)
	# Handling of multiple-line comments
	else:
		### Step1: extract all start and end into corresponding set()
		start_positions = set()
		end_positions = set()
		for l in in_gff_lines:
			## Check length and ignore comment lines
			if l[0].startswith("##sequence-region"):
				continue
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
			## Update length and ignore comment lines for lenght calculation
			if l[0].startswith("##sequence-region"):
				new_gff_lines.append(l)
				continue
			### Step5: Correct start(l[3]) and end(l[4]) by minus offset
			new_gff_line = modify_gff_line(
				l, start = int(l[3]) - offset + position, end = int(l[4]) - offset + position)
			new_gff_lines.append(new_gff_line)
	return new_gff_lines

def write_in_gff_lines(gff_out, in_gff_lines, position, split_features, sequence_length, chrom):
	'''
	in_gff_lines: a list of lists where each nested list is a list of 
		columns (in gff format) associated with each new feature to insert
	split_features: contains information about features in the original GFF 
		file that were split due to the insertion of the new sequence.
	sequence_length: length of the inserted sequence, used to determine 
		the new end positions in the GFF file.
	'''
 	## Replace the chromosome ID from in_gff with the correct chromosome ID
	for l in in_gff_lines:
		l[0] = chrom		
	## Check if the lenght in_gff_lines are valid, and correct if needed
	new_gff_lines = calculate_new_length_for_in_gff(in_gff_lines, position, sequence_length)
	for gff_line in new_gff_lines:
		gff_out.write(gff_line)
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
	chrom_id: the ID of the chromosome to modify
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
				## Create a tempfile to store uncompressde content
				with tempfile.NamedTemporaryFile(delete=False, mode='w') as tmp_f:
					tmp_f.write(f.read())
					ref_gff_path = tmp_f.name
		with open(ref_gff_path, "r") as f:
			for line in f:
				# For header lines, handle both tab and space-delimited formats
				if line.startswith("#"):
					if '\t' in line:
						line_elements = line.split('\t')
					else:
						line_elements = line.split()
					
					if line.startswith("##sequence-region") and line_elements[1] == chrom_id:
						## Edit the length of the chromosome 
						original_length = int(line_elements[3])
						new_length = calculate_new_length(original_length, position, down_position, new_seq_length)
						line = line.replace(str(original_length), str(new_length))
					gff_out.write(line)
				else:
					# Regular feature lines are always tab-delimited
					line_elements = line.split('\t')
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
							gff_out, in_gff_lines, position, split_features, new_seq_length, chrom_id)
					
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
						new_length = calculate_new_length(original_length, position, down_position, new_seq_length)
						line = line.replace(str(original_length), str(new_length))
						gff_out.write(line)
						
					# Split feature into 2 if feature starts before position
					# and ends after down_position
					elif gff_feat_start <= position and gff_feat_end > down_position:
						print("Feature split")
						print(line)
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
						print(f"Feature cut off - {x} prime side of feature cut off ({gff_feat_strand} strand)")
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
								gff_out, in_gff_lines, position, split_features, new_seq_length, chrom_id)
							
						# Change start position of feature to after cutoff point if
						# the feature starts within the deletion
						if (gff_feat_start > position 
							and gff_feat_start <= down_position 
							and gff_feat_end > down_position):
							x = "5" if gff_feat_strand == "+" else "3"
							print(f"Feature cut off - {x} prime side of feature cut off ({gff_feat_strand} strand)")
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
							print(f"** Error: Unknown case for GFF modification. Exiting {line_elements}")
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
					gff_out, in_gff_lines, position, split_features, new_seq_length, chrom_id)
			
			# Checking to ensure in_gff_lines written
			if not in_gff_lines_appended:
				print("** Error: Something went wrong, in_gff not added to reference gff. Exiting")
				exit()
		# Remove temp gtf file
		if ref_gff.endswith('.gz'):
			os.remove(ref_gff_path)
	return new_gff_name

def create_new_gff_for_existing_gff(new_gff_name, ref_gff, in_gff_lines, chrom_id, sequence_length):
	"""
	Appends new annotations to an existing GFF file without modifying existing features.
	"""
	gff_splitor = ''
	ref_gff_path = ref_gff
	## Handle compressed .gz GFF files
	if ref_gff.endswith('.gz'):
		with gzip_module.open(ref_gff, 'rt') as f:
			with tempfile.NamedTemporaryFile(delete=False, mode='w') as tmp_f:
				for line in f:
					tmp_f.write(line)
				tmp_f.flush()
				ref_gff_path = tmp_f.name
    ## Open new GFF file for writing
	with open(new_gff_name, "w") as gff_out:
        ## Copy all existing annotations to new GFF file
		with open(ref_gff_path, "r", encoding="utf-8") as f:
			for line in f:
				if line.startswith("##sequence-region") and gff_splitor == '':
					## Select user used delimiter based on content
					if '\t' in line:
						gff_splitor = ('\t')
					else:
						gff_splitor = (' ')
				gff_out.write(line)
		## Append new annotations if present
		if in_gff_lines:
			## Call calculate_new_length_for_in_gff to correct the length of the chromosome
			## Since we are not modifying existing features, we can set position to 0
			new_gff_lines = calculate_new_length_for_in_gff(in_gff_lines, 0, sequence_length)
			print(f"Appending {len(in_gff_lines)} new annotations to chromosome {chrom_id}.")
			for new_annotation in new_gff_lines:
				if new_annotation[0] == "##sequence-region":
					## Use predefined format for sequence-region line
					if gff_splitor == '':
						gff_splitor = new_annotation[-1]
					## Remove format indicater, and add new line
					gff_out.write(gff_splitor.join(new_annotation[:-1])+'\n')
				elif new_annotation:
					gff_out.write(new_annotation)
	
 	## Cleanup temp file if needed
	if ref_gff.endswith('.gz') and ref_gff_path != ref_gff:
		try:
			os.remove(ref_gff_path)
		except Exception as e:
			print(f"Warning: Failed to delete temp file {ref_gff_path}: {e}")
	
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
		print(f"** Error: Unrecognized extension {ext} in format_comment(). Exiting")
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
		print(f"Renaming split feature {elements[0]} --> {elements[0]}_split")
		return ("{}_split;{}".format(elements[0], ';'.join(elements[1:])))
	elif elements[0].startswith("gene_id "):
		gene_id = re.match(r'gene_id \"(.+)\"', elements[0])[1]
		print(f"Renaming split feature {gene_id} --> {gene_id}_split")
		return ('gene _id "{}_split";{}'.format(gene_id, ';'.join(elements[1:])))

	else:
		print(f"This feature will not be renamed because it does not have an ID/gene_id attribute:\n{line}")
		return attributes

def calculate_new_length(original_length, position, down_position, new_seq_length):
	"""
	Calculate the new chromosome/region length after a modification.
	Returns: The new calculated length
	"""
	return original_length - (down_position - position) + new_seq_length
	
def get_input_args():
	parser = argparse.ArgumentParser()
	chrom_group = parser.add_mutually_exclusive_group(required=True)
	chrom_group.add_argument('--chrom', type=str, 
                             help="Chromosome name (String)")
	chrom_group.add_argument('--new_chrom', type=str, 
                             help="Comma-separated new chromosome name (String)")
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
	if (len(in_args.in_fasta) != len(in_args.in_gff)):
		print("** Error: The number of inserted FASTA files does not match the number of GTF files, or their counts and positions do not align.")
		exit()
	else:
		iterations = len(in_args.in_fasta)
	## Modify existing chrom
	if in_args.chrom:
		if in_args.upstream_fasta:
			in_args.upstream_fasta = in_args.upstream_fasta.split(',')
		if in_args.downstream_fasta:
			in_args.downstream_fasta = in_args.downstream_fasta.split(',')
		if in_args.position is None and (in_args.upstream_fasta is None or in_args.downstream_fasta is None):
			print("** Error: You must provide either the position, or the upstream and downstream sequences.")
			exit()
		if in_args.position is not None:
			try:
				in_args.position = list(map(int, in_args.position.split(',')))
			except ValueError:
				print("** Error: Position must be a comma-separated list of integers, like 1,5,-1.")
				exit()
			if iterations != len(in_args.position):
				print("** Error: Position must be a equal to number of input FASTA")
				exit()
		else:
			if iterations != len(in_args.upstream_fasta):
				print("** Error: Upstream FASTA must be a equal to number of input FASTA")
				exit()
		if not in_args.position and len(in_args.upstream_fasta) != len(in_args.downstream_fasta):
			print("** Error: The number of upstream_fasta and downstream_fasta files does not match.")
			exit()
	## Add new chrom
	else:
		if in_args.position or in_args.upstream_fasta or in_args.downstream_fasta:
			parser.error("** Error: When using --new_chrom, you cannot provide --position, --upstream_fasta, or --downstream_fasta.")
			exit()
		## Convert new_chrom from string to list
		in_args.new_chrom = in_args.new_chrom.split(',')
	
	return in_args, iterations
	
if __name__ == "__main__":
	main()

