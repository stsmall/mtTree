#A collection of python modules for the mtdna pipeline
#Author: Aaron Steele
import subprocess
import sys

#Shift a fasta file so the center of the sequence is the ends.  Used to reduce edge effects of circular genomes
def fasta_shift(fastaIn):
	infile = open(fastaIn,"r")
	fsplit = fastaIn.split("/")
	fsplit[-1] = "shift_" + fsplit[-1]
	fastaOut = "/".join(fsplit)
	outfile = open(fastaOut,"w")

	seq = ""
	lines = infile.readlines()
	outfile.write(lines[0])
	for line in lines[1:]:
		seq += line.rstrip()
	
	length = len(seq)
	mid = length/2

	for i in range(mid,mid+length):
		outfile.write(seq[i%length])

	infile.close()
	outfile.close()

	return fastaOut

	
	

#Random samples X reads from a given fastq file
def fq_se_sample(fqIn, fqOut, sampleSize):
	
	import random
	infile = open(fqIn,"r")
	outfile = open(fqOut,"w")

	#Count the number of lines in the fast file
	read_count = 0
	line = infile.readline()
	while line:
		read_count += 1
		line = infile.readline()

	read_count = read_count/4
	infile.seek(0)
	if sampleSize >= read_count:
		sys.stderr.write("Sample size exceeds # of reads in " + fqIn)
		sys.stderr.write("\nOutputting all reads\n")
		line = infile.readline()
		while line:
			outfile.write(line)
			line = infile.readline()
		infile.close()
		outfile.close()
		return 0
	
	sample = random.samples(xrange(read_count),sampleSize)
	
	#Output the selected sample
	line = infile.readline()
	readnum = 0
	while line:
		if readnum in sample:
			outfile.write(line)
			outfile.write(infile.readline())
			outfile.write(infile.readline())
			outfile.write(infile.readline())
		else:
			infile.readline()
			infile.readline()
			infile.readline()
		line = infile.readline()
		readNum += 1
	infile.close()
	outfile.close()

def sample_pe_fq(fq1In, fq2In, fq1Out, fq2Out, sampleSize):
	import random

	infile1 = open(fq1In,"r")
	infile2 = open(fq2In,"r")
	outfile1 = open(fq1Out,"w")
	outfile2 = open(fq2Out,"w")
	
	read_count1 = 0
	read_count2 = 0
	
	#Count the reads in fq1
	line = infile1.readline()
	while line:
		read_count1 += 1
		line = infile1.readline()
	read_count1 = read_count1 / 4
	infile1.seek(0)

	#Count the reads in fq2
	line = infile2.readline()
	while line:
		read_count2 += 1
		line = infile2.readline()
	read_count2 = read_count2 / 4
	infile2.seek(0)

	if read_count1 != read_count2:
		sys.stderr.write("Fastq files are not the same size...Exiting")
		sys.exit(1)

	#If sampleSize is larger than reads, output all reads
	read_count = read_count1
	if sampleSize > read_count:
		sys.stderr.write("Sample size exceeds # of reads in " + fq1In + " and " + fq2In)
		sys.stderr.write("\nOutputting all reads\n")
		line = infile1.readline()
		while line:
			outfile1.write(line)
			line = infile1.readline()
		line = infile2.readline()
		while line:
			outfile2.write(line)
			line = infile2.readline()
		infile1.close()
		infile2.close()
		outfile1.close()
		outfile2.close()
		return 0

	#Random sample and output
	sample = random.sample(xrange(read_count),sampleSize)
	
	#Output the forward reads
	line = infile1.readline()
	readNum = 0
	while line:
		if readNum in sample:
			outfile1.write(line)
			outfile1.write(infile1.readline())
			outfile1.write(infile1.readline())
			outfile1.write(infile1.readline())
		else:
			infile1.readline()
			infile1.readline()
			infile1.readline()
		line = infile1.readline()
		readNum += 1
	infile1.close()
	outfile1.close()
	
	line = infile2.readline()
	readNum = 0
	while line:
		if readNum in sample:
			outfile2.write(line)
			outfile2.write(infile2.readline())
			outfile2.write(infile2.readline())
			outfile2.write(infile2.readline())
		else:
			infile2.readline()
			infile2.readline()
			infile2.readline()
		line = infile2.readline()
		readNum += 1
	infile2.close()
	outfile2.close()

def sam_2_pe(samfile,pe1,pe2):
	
	#Extract forward reads
	command = "cat " + samfile + " | grep -v ^@ | awk 'NR%2==1 {print \"@\"$1\"_1\\n\"$10\"\\n+\\n\"$11}' > " + pe1
	proc = subprocess.Popen(command,shell=True)
	proc.wait()

	#Extract reverse reads
	command = "cat " + samfile + " | grep -v ^@ | awk 'NR%2==0 {print \"@\"$1\"_2\\n\"$10\"\\n+\\n\"$11}' > " + pe2
	proc = subprocess.Popen(command,shell=True)
	proc.wait()

def getRefLength(reference):
    infile = open(reference,"r")
    lines = infile.readlines()
    infile.close()

    size = 0    
    for line in lines:
        if not line.startswith(">"):
            size += len(line.rstrip())
    
    return size



