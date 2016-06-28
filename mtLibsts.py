# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 11:18:33 2016
remake of mtLib.py
@author: stsmall
"""
import subprocess
import sys
import os

def fasta_shift(fastaIn):
    infile = open(fastaIn,"r")
    path_ref = os.path.split(fastaIn)  
    fastaOut = os.path.join(path_ref[0], "shift_" + path_ref[1])
    outfile = open(fastaOut,"w")

    seq = ""
    lines = infile.readlines()
    outfile.write(lines[0])
    for line in lines[1:]:
        seq += line.rstrip()
    
    length = len(seq)
    mid = length/2

    outfile.write(seq+seq[0:mid])
  
    infile.close()
    outfile.close()

    return fastaOut

def getRefLength(reference):
    infile = open(reference,"r")
    lines = infile.readlines()
    infile.close()
    size = 0    
    for line in lines:
        if not line.startswith(">"):
            size += len(line.rstrip()) 
    return size

def sam_2_pe(samfile,pe1,pe2):
          
    #assumes a name sorted sam file; e.g., samtools sort -n or sambamba sort -n -t20
    r1 = subprocess.Popen(['head', '-n 2', samfile],stdout=subprocess.PIPE)
    (out,err) = r1.communicate()
    #print out
    if out.split("\n")[0].split("\t")[0] in out.split("\n")[1].split("\t")[0]:  
        
        #Extract forward reads
        command = "cat " + samfile + " | grep -v ^@ | awk 'NR%2==1 {print \"@\"$1\"_1\\n\"$10\"\\n+\\n\"$11}' > " + pe1 
        print command        
        proc = subprocess.Popen(command,shell=True)
        proc.wait()
        
        #extract reverse reads
        command = "cat " + samfile + " | grep -v ^@ | awk 'NR%2==0 {print \"@\"$1\"_2\\n\"$10\"\\n+\\n\"$11}' > " + pe2
        proc = subprocess.Popen(command,shell=True)
        proc.wait()
    else:
        print "file not sorted by query name"
        sys.exit(1)
        
def write_random_records(fqa, fqb, N):
    """ get N random headers from a fastq file without reading the
    whole thing into memory"""
    import random    
    records = sum(1 for _ in open(fqa)) / 4
    rand_records = sorted([random.randint(0, records - 1) for _ in xrange(N)])

    fha, fhb = open(fqa),  open(fqb)
    suba, subb = open(fqa + ".subset", "w"), open(fqb + ".subset", "w")
    rec_no = - 1
    for rr in rand_records:

        while rec_no < rr:
            rec_no += 1       
            for i in range(4): fha.readline()
            for i in range(4): fhb.readline()
        for i in range(4):
            suba.write(fha.readline())
            subb.write(fhb.readline())
        rec_no += 1 # (thanks @anderwo)

    print ("wrote to %s, %s" %(suba.name, subb.name))
