"""
Filter Fastq Version

Author:  Lauren Kleine
 
This function parses through fastq files to determine if the version is correct for this
pipeline. If the version is inefficient, then the function 

"""

def trim_fastq(fastq_input_file, fastq_output_file):
	"""
	This function processes Fastq files by trimming line 2 (DNA sequence) and 
	line 4 (quality values) by specified numeric values on the 5' and 3' ends. 
	These numeric values are entered as the user arguements: 'first_base' and 
	'last_base'. This function requires a log file be entered as an arguement to 
	output statistics about the trimmed sequences, or to clarify errors.
	"""

	try:
		file_handle = open(fastq_input_file, 'r')
		fout = open(fastq_output_file, 'w')
	except:
		return -1
	with file_handle :
		count = 0
		for line in file_handle:
			count += 1
			if count % 4 == 1:
				seq1 = line
				fout.write(seq1)

		"""
		add if name = main stuff 
		"""

