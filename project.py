""" 
	RNA Alignment Assignment
	
	Implement each of the functions below using the algorithms covered in class.
	You can construct additional functions and data structures but you should not
	change the functions' APIs.

	You will be graded on the helper function implementations as well as the RNA alignment, although
	you do not have to use your helper function.
	
	*** Make sure to comment out any print statement so as not to interfere with the grading script
"""

import sys # DO NOT EDIT THIS
from shared import *
import textwrap

ALPHABET = [TERMINATOR] + BASES

def get_suffix_array(s):
	"""
	Naive implementation of suffix array generation (0-indexed). You do not have to implement the
	KS Algorithm. Make this code fast enough so you have enough time in Aligner.__init__ (see bottom).

	Input:
		s: a string of the alphabet ['A', 'C', 'G', 'T'] already terminated by a unique delimiter '$'
	
	Output: list of indices representing the suffix array

	>>> get_suffix_array('GATAGACA$')
	[8, 7, 5, 3, 1, 6, 4, 0, 2]
	"""

	bucket_keys = [a+b for a in BASES for b in ALPHABET] + [TERMINATOR]
	bucket_keys = sorted(bucket_keys)

	def add_to_bucket(d, bucket_key, value):
		if bucket_key in d:
			d[bucket_key].append(value)
		else:
			d[bucket_key] = [value]

	def generate_buckets(indices, offset):
		d = dict()
		for i in indices:
			j = i + offset
			if j == len(s) - 1:
				bucket_key = s[j]
			else:
				bucket_key = s[j] + s[j+1]
			add_to_bucket(d, bucket_key, i)
		return d

	def radix_sort(indices, offset):
		if len(indices) <= 1:
			return indices

		buckets = generate_buckets(indices, offset)

		result = list()

		for bucket_key in bucket_keys:
			if bucket_key in buckets:
				result += radix_sort(buckets[bucket_key], offset+2)

		return result

	return radix_sort(range(len(s)), 0)

def get_bwt(s, sa):
	"""
	Input:
		s: a string terminated by a unique delimiter '$'
		sa: the suffix array of s

	Output:
		L: BWT of s as a string
	"""
	L = [s[i-1] for i in sa ]
	return "".join(L)
	
def get_F(L):
	"""
	Input: L = get_bwt(s)

	Output: F, first column in Pi_sorted
	"""
	return "".join(sorted(L))

def get_M(F):
	"""
	Returns the helper data structure M (using the notation from class). M is a dictionary that maps character
	strings to start indices. i.e. M[c] is the first occurrence of "c" in F.

	If a character "c" does not exist in F, you may set M[c] = -1
	"""
	M = {}
	for i in range(len(F)):
		if F[i] in M.keys():
			continue
		else:
			M[F[i]] = i    
		M = {key:M[key] for key in sorted(M.keys())}
	return M 

def get_occ(L):
	"""
	Returns the helper data structure OCC (using the notation from class). OCC should be a dictionary that maps 
	string character to a list of integers. If c is a string character and i is an integer, then OCC[c][i] gives
	the number of occurrences of character "c" in the bwt string up to and including index i
	"""
	alphabet = list(set(L))
	Occ = {}
	for char in alphabet:
		Occ[char]=[]
 
	for i,c in enumerate(L):
		for char in alphabet:
			if len(Occ[char]) == 0:
				prev = 0
			else:
				prev = Occ[char][-1]
			if c == char:
				Occ[char].append(prev+1)
			else:
				Occ[char].append(prev)        
	return Occ
	

def exact_suffix_matches(p, M, occ):
	"""
	Find the positions within the suffix array sa of the longest possible suffix of p 
	that is a substring of s (the original string).
	
	Note that such positions must be consecutive, so we want the range of positions.

	Input:
		p: the pattern string
		M, occ: buckets and repeats information used by sp, ep

	Output: a tuple (range, length)
		range: a tuple (start inclusive, end exclusive) of the indices in sa that contains
			the longest suffix of p as a prefix. range=None if no indices matches any suffix of p
		length: length of the longest suffix of p found in s. length=0 if no indices matches any suffix of p

		An example return value would be ((2, 5), 7). This means that p[len(p) - 7 : len(p)] is
		found in s and matches positions 2, 3, and 4 in the suffix array.

	>>> s = 'ACGT' * 10 + '$'
	>>> sa = get_suffix_array(s)
	>>> sa
	[40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0, 37, 33, 29, 25, 21, 17, 13, 9, 5, 1, 38, 34, 30, 26, 22, 18, 14, 10, 6, 2, 39, 35, 31, 27, 23, 19, 15, 11, 7, 3]
	>>> L = get_bwt(s, sa)
	>>> L
	'TTTTTTTTTT$AAAAAAAAAACCCCCCCCCCGGGGGGGGGG'
	>>> F = get_F(L)
	>>> F
	'$AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT'
	>>> M = get_M(F)
	>>> sorted(M.items())
	[('$', 0), ('A', 1), ('C', 11), ('G', 21), ('T', 31)]
	>>> occ = get_occ(L)
	>>> type(occ) == dict, type(occ['$']) == list, type(occ['$'][0]) == int
	(True, True, True)
	>>> occ['$']
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
	>>> exact_suffix_matches('ACTGA', M, occ)
	((1, 11), 1)
	>>> exact_suffix_matches('$', M, occ)
	((0, 1), 1)
	>>> exact_suffix_matches('AA', M, occ)
	((1, 11), 1)
	"""
	
	# get the ordered list of keys of M, in order to find the "next" char 
	Order_char = list(sorted(list(M.keys())))
	# print(Order_char)
	
	#initialize:
	start = len(p)-1
	SP = M[p[start]]
	next_index = (1 + Order_char.index(p[start]))%len(Order_char)   
	EP=M[Order_char[next_index]]-1
	# print(SP,EP)
	
	for i in range(len(p)-1,-1,-1):

		start = i
		# Initialize search pointers:
		if(M[p[start]]+occ[p[start]][SP-1]<M[p[start]]+occ[p[start]][EP]-1):
			SP = M[p[start]]+occ[p[start]][SP-1]
			EP = M[p[start]]+occ[p[start]][EP]-1
		else:
			return (SP,EP+1),len(p)-i
	
	

def inexact_substring(str, substr, mis, first=True):
	
	results = []

	for i in range(0, len(str) - len(substr) + 1):
		
		mismatch = 0
		
		for j in range(0, len(substr)):
			
			if str[i+j] != substr[j]:
				mismatch += 1
				
				if mismatch > mis:
					break

		if mismatch <= mis:
			results.append((i, mismatch))
			if first:
				return results

	return results

MIN_INTRON_SIZE = 20
MAX_INTRON_SIZE = 10000

class Aligner:
	def __init__(self, genome_sequence, known_genes):
		"""
		Initializes the aligner. Do all time intensive set up here. i.e. build suffix array.

		genome_sequence: a string (NOT TERMINATED BY '$') representing the bases of the of the genome
		known_genes: a python set of Gene objects (see shared.py) that represent known genes. You can get the isoforms 
					 and exons from a Gene object

		Time limit: 500 seconds maximum on the provided data. Note that our server is probably faster than your machine, 
					so don't stress if you are close. Server is 1.25 times faster than the i7 CPU on my computer

		"""
		def exon_key(exon):
			return exon.start

		self.genome = genome_sequence
		self.transcriptome = dict()

		for gene in known_genes:
			for isoform in gene.isoforms:
				# append all exons, noting all junction points
				conc = ""
				index = 0
				isoform.exons = sorted(isoform.exons, key=exon_key)
				
				conc_map = list()

				for exon in isoform.exons:

					exon_str = genome_sequence[exon.start:exon.end]

					conc += exon_str

					conc_map.append((index, exon.start, len(exon_str)))

					index += len(exon_str)

				self.transcriptome[isoform.id] = (conc, conc_map)

		self.kmers_map = dict()

		for i in range(len(genome_sequence) - 6):
			key = genome_sequence[i:i+7]

			if key not in self.kmers_map:
				self.kmers_map[key] = [i]
			else:
				self.kmers_map[key].append(i)


	def align(self, read_sequence):
		"""
		Returns an alignment to the genome sequence. An alignment is a list of pieces. 
		Each piece consists of a start index in the read, a start index in the genome, and a length 
		indicating how many bases are aligned in this piece. Note that mismatches are count as "aligned".

		Note that <read_start_2> >= <read_start_1> + <length_1>. If your algorithm produces an alignment that 
		violates this, we will remove pieces from your alignment arbitrarily until consecutive pieces 
		satisfy <read_start_2> >= <read_start_1> + <length_1>

		Return value must be in the form (also see the project pdf):
		[(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

		If no good matches are found: return the best match you can find or return []

		Time limit: 0.5 seconds per read on average on the provided data.
		"""
		a = self.align_to_transcriptome(read_sequence)

		if not a:
			a = self.align_to_genome(read_sequence)

		return a

	def align_to_transcriptome(self, read_sequence):

		alignment = []

		for isoform_id, isoform in self.transcriptome.items():

			conc = isoform[0]
			conc_map = isoform[1]

			matches = inexact_substring(conc, read_sequence, 6)
			if len(matches) > 0:
				# match
				match, num_mismatches = matches[0]

				read_index = 0

				for mapping in conc_map:
					conc_start = mapping[0]
					genome_start = mapping[1]
					length = mapping[2]

					if match < conc_start + length:
						align_length = min(conc_start + length - match, len(read_sequence) - read_index)

						alignment.append([read_index, genome_start + (match - conc_start), align_length])

						read_index += align_length
						match += align_length

						if match >= matches[0][0] + len(read_sequence):
							break

				return alignment

		return []


	def align_to_genome(self, read_sequence):
		seed_len = 7
		seeds = textwrap.wrap(read_sequence, 7)

		align = []
		mismatches = 7

		for i, seed in enumerate(seeds):
			if seed not in self.kmers_map:
				continue

			kmers = self.kmers_map[seed]

			extend_left = i*seed_len
			
			for kmer in kmers:
				if kmer - extend_left + len(read_sequence) < len(self.genome) and kmer - extend_left >= 0:
					genome_substr = self.genome[kmer - extend_left : kmer - extend_left + len(read_sequence)]

					match = inexact_substring(read_sequence, genome_substr, 6)
					if len(match) > 0 and match[0][1] < mismatches:
						align = [(0, kmer - extend_left, len(read_sequence))]
						mismatches = match[0][1]

		return align
