# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 11:27:27 2016
usage: mtTree.py -h
this is a fork of the project by Aaron Steele https://bitbucket.org/steelea/ to assemble mitochondrial genomes from WGS data.
The code follows the methods suggested by Prado-Martinez Nature 2013: Great ape genetic diversity and population history.

Dependencies: bowtie2 path, hapsembler path, samtools path

@author: stsmall
"""
import argparse
import logging
import os
import subprocess
import sys
import mtLibsts  #seperate module

def get_args():
    parser = argparse.ArgumentParser(description='Assembles mitochondrial genomes from paired read info by mapping then assembly and consensus')  
    parser.add_argument('-bt2','--bowtie2', help='path to bowtie2 directory containing executable for bowtie2 and bowtie2-build')
    parser.add_argument('-f1','--fastq1', required=True, help='fastq containing read pair 1')      
    parser.add_argument('-f2','--fastq2', required=True, help='fastq containing read pair 2')
    parser.add_argument('-r','--reference', required=True, help='fasta containing reference')    
    parser.add_argument('-a','--hapsemblr',help='path to hapsemblr') 
    parser.add_argument('-c','--coverage',help='coverage for downsampling',type=int,default=300) 
    parser.add_argument('-l','--read_length',help='estiamted read length',type=int,default=250) 
    parser.add_argument('-t','--threads',help='number of threads for bowtie2',type=int,default=1)   
    parser.add_argument('-s','--samtools',help='path to samtools')    
    parser.add_argument('--debug',action='store_true',help='increase output for code debugging')
    args = parser.parse_args()
    return args

class mtTree:
    def __init__(self,args):
        self.bowtie2 = os.path.join(args.bowtie2,"bowtie2") #~/bin,bowtie2 =/bin/bowtie2      
        self.fastq1 = os.path.realpath(args.fastq1)
        self.fastq2 = os.path.realpath(args.fastq2)
        self.reference = os.path.realpath(args.reference)
        self.hapsemblr = os.path.realpath(args.hapsemblr) #executables in bin; path should end in bin since join is later 
        self.coverage = args.coverage
        self.read_length = args.read_length
        self.threads = args.threads
        self.cwd = os.path.split(self.fastq1)[0]
        self.samtools = os.path.join(args.samtools,"samtools") #~/bin,samtools =/bin/samtools
    
    def align(self,outputSam,reference): #shift_ref contains complete path
        '''align reads from fastq files using bowtie2'''        
       
        #check if index exists
        if os.path.isfile(reference + ".1.bt2"):
            pass
        else:
            #make index        
            command = self.bowtie2 + "-build -f " + reference + " " + reference
            proc = subprocess.Popen(command, shell=True)
            proc.wait()
        
        #run bowtie2 alignment
        command = self.bowtie2 + " -p " + str(self.threads) + " --no-unal -R 5 -N 1 -L 12 -D 25 -i S,2,.25 -x " + reference + " -1 " + self.fastq1 + " -2 " + self.fastq2 + " > out.sam" 
        print command        
        proc = subprocess.Popen(command, shell=True)
        proc.wait()

        #samtools cull quality
        command = self.samtools + " view -q 30 out.sam > out.q30.sam"
        proc = subprocess.Popen(command, shell=True)
        proc.wait() 

        #sort with samtools
        command = self.samtools + " sort -n -@ " + str(self.threads) + " out.q30.sam -o " + outputSam
        #print command         
        proc = subprocess.Popen(command, shell=True)
        proc.wait()            
      
    def assemble(self,startCount,endCount,sam):    
        '''assemble reads using hapsemblr'''
        
        #bam to paired-end        
        mtLibsts.sam_2_pe(sam,"mit_1.fq","mit_2.fq")        
        
        #Determine sample size using coverage and read length
        refLength = mtLibsts.getRefLength(self.reference)
        sampleSize = int((refLength * self.coverage/2)/self.read_length)      
        #run hapsemblr
        for i in xrange(startCount,endCount+1):
            #random sample
            mtLibsts.write_random_records("mit_1.fq", "mit_2.fq", sampleSize)
            
            command = os.path.join(self.hapsemblr,"preprocr") + " -p illumina -f mit_1.fq.subset -x mit_2.fq.subset -o mit.fq.tmp -d 33"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = os.path.join(self.hapsemblr,"overlappr") + " -p illumina -f mit.fq.tmp -o mitK -g 15 -t " + str(self.threads)
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = os.path.join(self.hapsemblr,"hapsemblr") + " -r mitK -c contigs.fa -g 15"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = os.path.join(self.hapsemblr,"consensr") + " -p illumina -f mit.fq.tmp -c contigs.fa -o mit_contigs." + str(i) + ".fa -d 33"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()

            command = "rm mit_1.fq.subset mit_2.fq.subset mit.fq.tmp"
            proc = subprocess.Popen(command, shell=True)
            proc.wait()
        
        #cat all mit_contigs.*.fa to mit_contigs.f.fa and change names so unique
        filenames = ["mit_contigs.1.fa","mit_contigs.2.fa","mit_contigs.3.fa","mit_contigs.4.fa","mit_contigs.5.fa"]
        with open("mit_contigs.f2.fa",'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    outfile.write(infile.read())
        #rename the headers so no dups
        j = 0
        with open("mit_contigs.f.fa",'w') as outfile:
            with open('mit_contigs.f2.fa','r') as infile:
                for line in infile:
                    if line.startswith(">"):
                        outfile.write(">mit_%i\n" %j)
                        j += 1
                    else:
                        outfile.write(line)                
    
    def run(self):
        '''run everything at once'''
        print self.cwd        
        os.chdir(self.cwd)
        print self.reference

        #shift reference        
        shiftRef = mtLibsts.fasta_shift(self.reference)        

        #run align        
        sys.stderr.write("Performing regular Pipeline\n")        
        self.align("mit_mapped.sam",shiftRef)
        self.assemble(1,5,"mit_mapped.sam")

       #Cleanup 
        sys.stderr.write("Cleaning up temp files\n")
        command = "rm mit_contigs.{1,2,3,4,5}.fa mit_contigs.f2.fa out.sam out.q30.sam && " + self.samtools + " view -Sb mit_mapped.sam > mit_mapped.sort.bam"
        proc = subprocess.Popen(command,shell=True)
        proc.wait()
  
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