import unittest
import filecmp
import subprocess
import os

class TestReform(unittest.TestCase):
	def setUp(self):
		self.wd = os.getcwd()
		self.cleanup()  # Clean up before each test case

	def tearDown(self):
		os.chdir(self.wd)
		self.cleanup()  # Clean up after each test case

	def cleanup(self):
		print("Clean Up")
		# Cleanup command to remove ref_reformed.fa and ref_reformed.gtf files
		cleanup_command = 'find test_data/ -type f \( -name "ref_reformed.fa" -o -name "ref_reformed.gtf" -o -name "ref_reformed.gff3" \) -exec rm -f {} +'
		os.system(cleanup_command)
		print("Done")
	
	def test_case_1(self):
		"""
		Case 1:
		Indel within a feature resulting in deleted sequence,
		inserted sequence, and splitting of existing features
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/1/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--upstream_fasta=up.fa \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
		
	def test_case_2(self):
		"""
		Case 2:
		Indel resulting in truncating 3' side of feature.
		Indel starts at position 2 of feature and ends at last position 
		of feature. All but first position of original feature remains. 
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/2/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--upstream_fasta=up.fa \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
		
	def test_case_3(self):
		"""
		Case 3:
		Deletion resulting in removal of entire feature.
		Deletion starts one base upstream of the feature and ends 
		at last position of feature. 
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/3/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--upstream_fasta=up.fa \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
		
	def test_case_4(self):
		"""
		Case 4:
		Indel causing feature split. 
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/4/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--upstream_fasta=up.fa \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
		
	def test_case_5(self):
		"""
		Case 5:
		1-base deletion, 10 base insertion. 
		Deletion is position 1 of existing feature causing a 1 bp 
		truncation upstream side (5' side) of feature (first base), and the 
		remainder of the feature to be offset 9 bases downstream. 
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/5/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--upstream_fasta=up.fa \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
		
	def test_case_6(self):
		"""
		Case 6:
		2 bp deletion, 10 bp insertion. 
		Indel is upstream of existing features.
		Result is all existing features offset (down)
		by 8 bases. 
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/6/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--upstream_fasta=up.fa \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
		
	def test_case_7(self):
		"""
		Case 7:
		Simple insertion using the position argument.
		Insertion upstream of existing feature, will get
		offset by length of insertion. 
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/7/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--position=3
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)

	def test_case_8(self):
		"""
		Case 8:
		Testing truncating features on reverse strand
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/8/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--upstream_fasta=up.fa \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)

	def test_case_9(self):
		"""
		Case 9:
		Insertion using the position argument at position 0
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/9/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--position=0
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
		
	def test_case_10(self):
		"""
		Case 10:
		Insertion using the position argument at position -1 
		(end of chromosome)
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/10/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--position=-1
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
		
	def test_case_11(self):
		"""
		Case 11:
		Testing GTF3 comment format (all above are GTF)
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/11/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--in_fasta=in.fa \
		--in_gff=in.gff3 \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gff3 \
		--position=5
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gff3', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gff3', 'r') as f:
			new_gff = f.read()
		print("Testing GFF3")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
	

	def test_case_12(self):
		"""
		Case 12:
		Testing Case 2 with .gz ref sequence input
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/12/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--upstream_fasta=up.fa \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa.gz \
		--ref_gff=ref.gtf.gz \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
	
	def test_case_13(self):
		"""
		Case 13:
		Testing Case 10 with .gz ref sequence input
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/13/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa.gz \
		--ref_gff=ref.gtf.gz \
		--position=-1
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
	
	def test_case_14(self):
		"""
		Case 14:
		Testing Sequential Processing which use multiple up.fa and down.fa files
		"""

		wd = os.getcwd()
		os.chdir('test_data/14/')

		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--upstream_fasta=up1.fa,up2.fa,up3.fa \
		--in_fasta=in1.fa,in2.fa,in3.fa \
		--in_gff=in1.gtf,in2.gtf,in3.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--downstream_fasta=down1.fa,down2.fa,down3.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)


	def test_case_15(self):
		"""
		Case 15:
		Testing Sequential Processing which combine test case 9, 10, 11
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/15/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--in_fasta=in1.fa,in2.fa,in3.fa \
		--in_gff=in1.gtf,in2.gtf,in3.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--position=0,5,-1
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing gtf")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)

	def test_case_16(self):
		"""
		Case 10:
		Use chromosome name from '--chrom' when chrom_id in in.gtf is different from ref.gtf
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/16/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--position=-1
		"""

		response = subprocess.getoutput(command)
		print(response)

		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
	
	def test_case_17(self):
		"""
		Case 17:
		The format of the insert files is different from that of the reference file.
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/17/')
		
		# ref.gff3 with in.gtf
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--in_fasta=in.fa \
		--in_gff=in.gtf \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gff3 \
		--position=5
		"""

		response = subprocess.getoutput(command)
		print(response)

		with open('gold.gff3', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gff3', 'r') as f:
			new_gff = f.read()
		print("Testing GFF3")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")

		# ref.gtf with in.gff3
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--in_fasta=in.fa \
		--in_gff=in.gff3 \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gtf \
		--position=5
		"""

		response = subprocess.getoutput(command)
		print(response)

		with open('gold.gtf', 'r') as f:
			gold_gff = f.read()
		with open('ref_reformed.gtf', 'r') as f:
			new_gff = f.read()
		print("Testing GTF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('ref_reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)

if __name__ == '__main__':
    unittest.main()
