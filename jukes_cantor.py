import re
import math

NUCLEOBASES = 'ACGT'

# Jukes Cantor distance formula: (-3/4)ln[1-p*(4/3)]
def distance (seq1, seq2):
    p = percentage_of_character_changes(seq1, seq2)
    return -0.75*math.log(1-(p*4/3))

# percentage of characters that are different between two sequences
def percentage_of_character_changes (seq1, seq2):
	
	# variable to hold number of nucleobase changes
	change_count = 0

	#variable to hold number of characters concidered in the calculation
	total_characters_count = 0.0 # this should be float to give correct percent at the end

	# loop through two sequences and count number of character changes and total characters considered in the calculation
	for a, b in zip(seq1, seq2):
		if (a in NUCLEOBASES) and (b in NUCLEOBASES):
			total_characters_count += 1
			if a != b:
				change_count += 1
	
	# return percentage of characters that are different between the two sequences
	return change_count / total_characters_count if total_characters_count else 0.0


# store previously read sequences
sequences_read = {}

# regex to find sequence id 
pattern = re.compile(r'(?<=\>)(.+?)(?=\|)')

# read input file
with open('Public_Eacles_Morpho_Aglia.fas', 'r') as input_file:
	while True:
		line = input_file.readline()
		if line.strip() == '': break
		
        # extract sequence id and read the sequence
		new_seq_id = pattern.search(line).group()
		new_seq = input_file.readline()
		
		# calculate distance with regard to the previously read sequences
		for seq_id, seq in sequences_read.items():
			d = distance(new_seq, seq)
            # discard zeros and distances greater than 4%
			if d and d <= 0.04: print ('{0}\t{1}\t{2}'.format(new_seq_id, seq_id, d))
		
		# add this sequence to previously read sequences
		sequences_read.update({new_seq_id: new_seq})
