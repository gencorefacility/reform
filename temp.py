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
						new_length = original_length - (down_position - position) + new_seq_length
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
								gff_out, in_gff_lines, position, split_features, new_seq_length, chrom_id)
							
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
					gff_out, in_gff_lines, position, split_features, new_seq_length, chrom_id)
			
			# Checking to ensure in_gff_lines written
			if not in_gff_lines_appended:
				print("** Error: Something went wrong, in_gff not added to reference gff. Exiting")
				exit()
		# Remove temp gtf file
		if ref_gff.endswith('.gz'):
			os.remove(ref_gff_path)
	return new_gff_name


"""
python3 ../../reform.py \
		--new_chrom="Y" \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf
"""