import unittest
import filecmp
import subprocess
import os

class TestReform(unittest.TestCase):
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
		--in_gff=in.gff3 \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gff3 \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gff3', 'r') as f:
			gold_gff = f.read()
		with open('reformed.gff3', 'r') as f:
			new_gff = f.read()
		print("Testing GFF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)
		
	def test_case_2(self):
		"""
		Case 2:
		Indel resulting in truncating downstream side of feature.
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
		--in_gff=in.gff3 \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gff3 \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gff3', 'r') as f:
			gold_gff = f.read()
		with open('reformed.gff3', 'r') as f:
			new_gff = f.read()
		print("Testing GFF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('reformed.fa', 'r') as f:
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
		--in_gff=in.gff3 \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gff3 \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gff3', 'r') as f:
			gold_gff = f.read()
		with open('reformed.gff3', 'r') as f:
			new_gff = f.read()
		print("Testing GFF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('reformed.fa', 'r') as f:
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
		--in_gff=in.gff3 \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gff3 \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gff3', 'r') as f:
			gold_gff = f.read()
		with open('reformed.gff3', 'r') as f:
			new_gff = f.read()
		print("Testing GFF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('reformed.fa', 'r') as f:
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
		truncation upstream side of feature (first base), and the 
		remainder of the feature to be offset 9 bases downstream. 
		"""
		
		wd = os.getcwd()
		os.chdir('test_data/5/')
		
		command = """
		python3 ../../reform.py \
		--chrom="X" \
		--upstream_fasta=up.fa \
		--in_fasta=in.fa \
		--in_gff=in.gff3 \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gff3 \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gff3', 'r') as f:
			gold_gff = f.read()
		with open('reformed.gff3', 'r') as f:
			new_gff = f.read()
		print("Testing GFF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('reformed.fa', 'r') as f:
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
		--in_gff=in.gff3 \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gff3 \
		--downstream_fasta=down.fa
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gff3', 'r') as f:
			gold_gff = f.read()
		with open('reformed.gff3', 'r') as f:
			new_gff = f.read()
		print("Testing GFF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('reformed.fa', 'r') as f:
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
		--in_gff=in.gff3 \
		--ref_fasta=ref.fa \
		--ref_gff=ref.gff3 \
		--position=3
		"""

		response = subprocess.getoutput(command)
		print(response)
	
		with open('gold.gff3', 'r') as f:
			gold_gff = f.read()
		with open('reformed.gff3', 'r') as f:
			new_gff = f.read()
		print("Testing GFF")
		self.assertListEqual(list(gold_gff), list(new_gff))
		print("Done")
		
		with open('gold.fa', 'r') as f:
			gold_fa = f.read()
		with open('reformed.fa', 'r') as f:
			new_fa = f.read()
		print("Testing Fasta")
		self.assertListEqual(list(gold_fa), list(new_fa))
		print("Done")
		
		os.chdir(wd)


if __name__ == '__main__':
    unittest.main()
