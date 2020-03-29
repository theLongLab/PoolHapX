import sys

def convert2ms(timept,popsize):

	# Find out how many haplotypes there are. 
        lines = []
	with open(timept, "r") as file: 
		lines = file.readlines()
	numhaps = len(lines[0].strip().split("\t")) - 1

	# Based on the within-pool frequencies, find out how many to make of each kind of haplotype. 
	hapcts = [] # The number to generate for ms format
	hapfreqs = lines[1].strip().split("\t") # The actual frequency of the haplotypes
	sumcts = 0
	targetct = int(popsize) # Target number of individuals to make.
	for f in range(1,len(hapfreqs)): # For frquencies except the raw name and the last one...
		ct = round(float(hapfreqs[f]) * targetct)
		# print(str(f) + " " + str(hapfreqs[f]) + " " + str(ct))
		hapcts.append(ct)
		sumcts += ct
	difference = sumcts - targetct
	if difference > 0:
		reduction = float(1 - difference/sumcts)
		for h in range(0,numhaps):
			hapcts[h] = round(hapcts[h] * reduction)
			# print(str(hapcts[h]), end = " ")

	# Generate the allelic composition of each kind of haplotype. Position by position, not haplotype by haplotype. 
	tmpcomps = lines[2].strip().split("\t")[1:numhaps + 1] # Append each allele to the list. This creates a space for each haplotype.
	for l in range(3,len(lines)):
		alleles = lines[l].strip().split("\t")[1:numhaps + 1]
		for a in range(0,numhaps):
			tmpcomps[a] += alleles[a]
	# print(tmpcomps)

	# hapcomps is the master list of the within-pool haplotypes.  
	hapcomps = []
	for h in range(0,numhaps):
		for c in range(int(hapcts[h])):
			hapcomps.append(tmpcomps[h]) # For each haplotype in the list, write it count times.
			# print(hapcomps)
	return hapcomps

def get_snps(example,reflen):
	positions = []
	with open(example, "r") as file:
		lines = file.readlines() # Skip the first two lines describing the haplotypes and their within-sample frequencies.
		for l in range(2,len(lines)):
			positions.append(float(lines[l].split("\t")[0].split(";")[1])/reflen) # Store (position number)/(reference sequence length) i.e.: ms style
	return positions

def printms(workdir,h,hapcomps,positions):
	tot_hapnum = len(hapcomps)
	segsites = len(positions)
	with open(workdir + "/" + str(h) + ".ms.txt","w+") as msf:
		msf.write("./ms " + str(tot_hapnum) + " 1\n0 0 0\n\n//\nsegsites: " + str(segsites) + "\npositions:") # Write ./ms\n1 2 3\n\n//\nsegsites: s\npositions: \plist\n
		for p in range(0,segsites):
			msf.write(" " + str(positions[p]))
		msf.write("\n")
		for h in range(0,tot_hapnum):
			msf.write(hapcomps[h] + "\n") # For each haplotype in the list, write it in. Counts have been taken care of by convert2ms.

def main():
	workdir = sys.argv[1]
	positions = get_snps(workdir + "/p0.txt", float(sys.argv[2]))

	for h in range(int(sys.argv[3]),int(sys.argv[4])): 
		hapcomps = convert2ms(workdir + "/p" + str(h) + ".txt", sys.argv[5])
		printms(sys.argv[6],h,hapcomps,positions)

if __name__ == '__main__':
	main()
