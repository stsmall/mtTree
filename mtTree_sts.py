# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 11:27:27 2016
usage: mtTree.py -h
this is a fork of the project by Aaron Steele https://bitbucket.org/steelea/ to assemble mitochondrial genomes from WGS data.
The code follows the methods suggested by Prado-Martinez Nature 2013: Great ape genetic diversity and population history.

Dependencies: bowtie2, hapsembler, nucmer (from mummer3), cSequenceBuilder

@author: stsmall
"""
import argparse
import logging
import os
import subprocess
import sys
import mtLibsts  #seperate module
from cSequenceBuilder import cSequenceBuilder  #seperate module


def get_args():
    parser = argparse.ArgumentParser(description='Assembles mitochondrial genomes from paired read info by mapping then assembly and consensus')  
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-bt2','--bowtie2', help='path to bowtie2 directory containing executable for bowtie2 and bowtie2-build')
    group.add_argument('-bwa','--bwa', help='path to bwa directory containing executable')    
    parser.add_argument('-f1','--fastq1', required=True, help='fastq containing read pair 1')      
    parser.add_argument('-f2','--fastq2', required=True, help='fastq containing read pair 2')
    parser.add_argument('-r','--reference', required=True, help='fasta containing reference')    
    parser.add_argument('-n','--nucmer',help='path to nucmer')
    parser.add_argument('-a','--hapsemblr',help='path to hapsemblr') 
    parser.add_argument('-c','--coverage',help='coverage for downsampling',type=int,default=300) 
    parser.add_argument('-l','--read_length',help='estiamted read length',type=int,default=250) 
    parser.add_argument('-t','--threads',help='number of threads for bowtie2',type=int,default=1)   
    parser.add_argument('-s','--samtools',help='path to samtools')
    parser.add_argument('-bdt','--bedtools',help='path to bedtools')
    parser.add_argument('-smb','--sambamba',help='path to sambamba')    
    parser.add_argument('--debug',action='store_true',help='increase output for code debugging')
    args = parser.parse_args()
    return args

class mtTree:
    def __init__(self,args):
        self.bowtie2 = os.path.join(args.bowtie2,"bowtie2")
        self.fastq1 = os.path.realpath(args.fastq1)
        self.fastq2 = os.path.realpath(args.fastq2)
        self.reference = os.path.realpath(args.reference)
        self.nucmer = os.path.join(args.nucmer,"nucmer")
        self.hapsemblr = os.path.realpath(args.hapsemblr) #executables in bin; path should end in bin
        self.coverage = args.coverage
        self.read_length = args.read_length
        self.threads = args.threads
        self.cwd = os.path.split(self.fastq1)[0]
        self.samtools = os.path.realpath(args.samtools)
        self.bedtools = os.path.realpath(args.bedtools) #executables in bin; path should end in bin
        self.sambamba = os.path.realpath(args.sambamba)
    def align(self,outputSam,reference):
        '''align reads from fastq files using bowtie2'''
        
        #make index        
        command = self.bowtie2 + "-build -f " + self.reference + " " + self.reference
        proc = subprocess.Popen(command, shell=True)
        proc.wait()
        
        #run bowtie2 alignment
        #as bam        
        #command = self.bowtie2 + " -p " + str(self.threads) + " --no-unal -R 5 -N 1 -L 12 -D 25 -i S,2,.25 -x " + self.reference + " -1 " + self.fastq1 + " -2 " + self.fastq2 + " | " + self.samtools + " view -bS > out.bam" 
        #as sam        
        command = self.bowtie2 + " -p " + str(self.threads) + " --no-unal -R 5 -N 1 -L 12 -D 25 -i S,2,.25 -x " + self.reference + " -1 " + self.fastq1 + " -2 " + self.fastq2 + " > out.sam" 
        print command        
        proc = subprocess.Popen(command, shell=True)
        proc.wait() 
        
        #sort bam        
        command = self.sambamba + " sort -n -t " + str(self.threads) + " out.sam -o " + outputSam
        proc = subprocess.Popen(command, shell=True)
        proc.wait()       
      
    def assemble(self,startCount,endCount,sam):    
        '''assemble reads using hapsemblr'''
        
        #bam to paired-end        
        #command = os.path.join(self.bedtools,"bamToFastq") + " -i " + sam + " -fq mt_1.fq -fq2 mt_2.fq"
        mtLibsts.sam_2_pe(sam,"mit_1.fq","mit_2.fq")
        
        
        #Determine sample size using coverage and read length
        refLength = mtLibsts.getRefLength(self.reference)
        sampleSize = int((refLength * self.coverage/2)/self.read_length)      
        
        #run hapsemblr
        for i in xrange(startCount,endCount+1):
            #random sample
            mtLibsts.sample_pe_fq("mit_1.fq", "mit_2.fq", "mit_1.fq.tmp", "mit_2.fq.tmp", sampleSize)
            
            command = os.path.join(self.hapsemblr,"preprocr") + " -p illumina -f mit_1.fq.tmp -x mit_2.fq.tmp -o mit.fq.tmp -d 33"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = os.path.join(self.hapsemblr,"overlappr") + " -p illumina -f mit.fq.tmp -o mitK -g 15 -t " + str(self.threads)
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = os.path.join(self.hapsemblr,"hapsemblr") + " -r mitK -c contigs.fa -g 15"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = os.path.join(self.hapsemblr,"consenr") + " -p illumina -f mit.fq.tmp -c contigs.fa -o mit_contigs." + str(i) + ".fa -d 33"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = self.nucmer + " --mum -p mit_aln." + str(i) + " " + self.reference + " mit_contigs." + str(i) + ".fa"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = "rm mit_contigs." + str(i) + ".fa.tmp mit_1.fq.tmp mit_2.fq.tmp mit.fq.tmp"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()
            
    def run(self):
        '''run everything at once'''
        os.chdir(self.cwd)
        #shift reference        
        shiftRef = mtLibsts.fasta_shift(self.reference)        
        #run align        
        sys.stderr.write("Performing regular Pipeline\n")        
        self.align("mit_mapped.bam",shiftRef)
        self.assemble(1,5,"mit_mapped.bam")
        #build consensus        
        sys.stderr.write("Builing a consensus sequence\n")
        c = cSequenceBuilder(self.cwd, self.reference)
        c.buildSequence("mit_consensus.fa")
        #Cleanup 
        #sys.stderr.write("Cleaning up temp files\n")
        #command = "rm mit_contigs.*.fa mit_aln.*.delta mitK.*"
        #proc = subprocess.Popen(command,shell=True)
        #proc.wait()
  
def main():
    args = get_args()
    if args.debug:
        logging.basicConfig(filename='mttree.debug.log',level=logging.DEBUG)
    else:
        logging.basicConfig(filename='mttree.log',level=logging.INFO)
    mtt = mtTree(args)
    mtt.run()

if __name__ == '__main__':
    main()