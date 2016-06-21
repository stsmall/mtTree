#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 16:17:20 2016
mtTree.py
this is a fork of the project by Aaron Steele to assemble mitochondrial genomes from WGS data.
The code follows the methods suggested by Prado-Martinez on great apes
@author: stsmall
"""

#dependencies: BWA, samtools, hapsembler, nucmer (from mummer3), Picard

import sys
import os
import subprocess
import getopt
import mtLib  #seperate module
from cSequenceBuilder import cSequenceBuilder  #seperate module

class mtTree:
    def __init__(self,argv):
        self.argv = argv
        self.bwa = "bwa"
        self.bowtie2 = ""
        self.nucmer = "nucmer"
        self.samtools = "samtools"
        self.picardtools ="SamToFastq"
        self.fq1 = ""
        self.fq2 = ""
        self.threads = 20
        self.reference = ""
        self.readLength = 250
        self.coverage = 300
        self.fq1Mapped = ""
        self.fq2Mapped = ""
        self.cwd = ""  #current working directory
        self.hap = "" #path to hapsembler working directory

    def usage(self):
         print("\nUsage: mtTree.py [OPTIONS]")
         print("\t -h\t Prints this message")
         print( "[REQUIRED]")
         print("\t -f\t First fastq file")
         print("\t -g\t Second fastq file")
         print("\t -r\t Reference assembly")
         print("\t -b\t Location of bwa executable (default=bwa)")
         print("\t -bt2\t Location of bowtie2 executable (default=bowtie2)")         
         print("\t -s\t Location of samtools executable (default=samtools)")
         print("\t -n\t Location of nucmer executable (default=nucmer)")
         print("\t -a\t Hapsembler installation directory")
         print("\t -p\t Location of SamToFastq executable (defualt=SamToFastq)")
         print("[OPTIONAL]")
         print("\t -c\t Downsample to coverage for assembly (default=300)")
         print("\t -l\t Average length of read (default = 250)")
         print("\t -t\t # of threads to use (default=20)")
    
    #TODO: Need checks to make sure required inputs are set
    def processArguments(self):
        try:
            opts, args = getopt.getopt(self.argv, "h:r:f:g:c:l:t:b:bt2:s:n:a:p:",["help","reference","fq1","fq2","coverage","readlength","threads","bwa","bowtie2","samtools","nucmer","hap","picardtools"])
        except getopt.GetoptError:
            print "Unrecongized option"
            self.usage()
            sys.exit(2)
    
        for opt,arg in opts:
            if opt in ("-h","--help"):
                self.usage()
                sys.exit(0)
            elif opt in ("-f","--fq1"):
                self.fq1 = os.path.realpath(arg)
                #Set the cwd
                path = os.path.realpath(self.fq1)
                psplit = path.split("/")
                self.cwd = "/".join(psplit[:-1]) + "/"
            elif opt in ("-g","--fq2"):
                self.fq2 = os.path.realpath(arg)
            elif opt in ("-r","--reference"):
                self.reference = os.path.realpath(arg)
            elif opt in ("-c","--coverage"):
                self.coverage = int(arg)    
            elif opt in ("-t","--threads"):
                self.threads = int(arg)
            elif opt in ("-b","--bwa"):
                self.bwa = os.path.realpath(arg)
            elif opt in ("-bt2","--bowtie2"):
                self.bowtie2 = os.path.realpath(str(arg))
            elif opt in ("-s","--samtools"):
                self.samtools = os.path.realpath(arg)
            elif opt in ("-n","--nucmer"):
                self.nucmer = os.path.realpath(arg)
            elif opt in ("-l","--readLength"):
                self.readLength = arg
            elif opt in ("-a", "--hap"):
                self.hap = os.path.realpath(arg)
            elif opt in ("-p", "--picardtools"):
                self.picardtools = os.path.realpath(arg)
    def alignMEM(self,outputSam,reference):      
#            command = self.bwa + " index " + reference
#            proc = subprocess.Popen(command, shell=True)
#            proc.wait()
#            
#            command = self.bwa + " mem -t " + str(self.threads) + " " + reference + " " + self.fq1 + " " + self.fq2 + " > tmp.sam"
#            proc = subprocess.Popen(command, shell=True)
#            proc.wait()
        command = self.bowtie2 + "-build -f " + reference
        print command
        proc = subprocess.Popen(command, shell=True)
        proc.wait()  

        command = self.bowtie2 + "bowtie2 -p " + str(self.threads) + " --no-unal -R 5 -N 1 -L 12 -D 25 -i S,2,.25 -x " + reference + " -1 " + self.fq1 + " -2 " + self.fq2 + " -S tmp.sam"
        proc = subprocess.Popen(command, shell=True)
        proc.wait()     
    
        #mappped reads only
        command = self.samtools + " view -q 15 -F 4 tmp.sam > " + outputSam
        proc = subprocess.Popen(command, shell=True)
        proc.wait()

        command = "rm tmp.sam"
        proc = subprocess.Popen(command, shell=True)
        proc.wait()

    def buildAssemblies(self,startCount,endCount,sam):

        #Convert the sam file to fastq
        command = self.picardtools + " I=" + sam + " F=mit_1.fq F2=mit_2.fq"
        proc = subprocess.Popen(command, shell=True)
        proc.wait()
        
        #Determine sample size using coverage and read length
        refLength = mtLib.getRefLength(self.reference)
        sampleSize = int((refLength * self.coverage/2)/self.readLength)

        for i in xrange(startCount,endCount+1):
            #Random sample
            mtLib.sample_pe_fq("mit_1.fq", "mit_2.fq", "mit_1.fq.tmp", "mit_2.fq.tmp", sampleSize)
            
            command = self.hap + "preprocr -p illumina -f mit_1.fq.tmp -x mit_2.fq.tmp -o mit.fq.tmp -d 33"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = self.hap + "overlappr -p illumina -f mit.fq.tmp -o mitK -g 15 -t " + str(self.threads)
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = self.hap + "hapsemblr -r mitK -c contigs.fa -g 15"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = self.hap + "consensr -p illumina -f mit.fq.tmp -c contigs.fa -o mit_contigs." + str(i) + ".fa -d 33"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = self.nucmer + " --mum -p mit_aln." + str(i) + " " + self.reference + " mit_contigs." + str(i) + ".fa"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = "rm mit_contigs." + str(i) + ".fa.tmp mit_1.fq.tmp mit_2.fq.tmp mit.fq.tmp"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()
    
        
        #Clean up
        command = "rm mit_1.fq mit_2.fq"
        proc = subprocess.Popen(command, shell=True)
        proc.wait()

    def run(self):
        self.processArguments()

        os.chdir(self.cwd)
        #Shift the reference
        shiftRef = mtLib.fasta_shift(self.reference)
        
        #Perform Regular ppieline
        sys.stderr.write("Performing regular Pipeline\n")
        self.alignMEM("mit_mapped_norm.sam", self.reference)
        self.buildAssemblies(1,5,"mit_mapped_norm.sam")

        #Perform Shifted pipeline
        sys.stderr.write("Performing Shifted Pipelinen\n")
        self.alignMEM("mit_mapped_shift.sam", shiftRef)
        self.buildAssemblies(6,10,"mit_mapped_shift.sam")        

        #Build consensus sequence
        sys.stderr.write("Builing a consensus sequence\n")
        c = cSequenceBuilder(self.cwd, self.reference)
        c.buildSequence("mit_consensus.fa")

        #Cleanup 
        sys.stderr.write("Cleaning up temp files\n")
        command = "rm mit_contigs.*.fa mit_aln.*.delta mitK.* " + shiftRef + "*"
        proc = subprocess.Popen(command,shell=True)
        proc.wait() 
    

if __name__ == "__main__":
    mtt = mtTree(sys.argv[1:])
    if len(sys.argv[1:]) < 2:
        mtt.usage()
        sys.exit(2)
    mtt.run()
