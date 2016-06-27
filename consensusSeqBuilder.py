#!/usr/bin/env python

#   This is a program constructs a consensus sequence (fasta) using 
#   multiple alignments (.delta) to a common reference. 

#   Copyright (C) <2014>  <Andres Martin>
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import sys
import operator
from os import listdir
from os.path import isfile, join
import os

class cSequenceBuilder:
	def __init__(self,startDir,ref):
		self.startDir = startDir

		#Build the consensus sequence dictionary and populate
		self.cSequence = {}	#{base position : {A : votes} {G : votes} {T : votes} {C : votes} {n : votes}}
		refFile = open(ref,"r")
		lines = refFile.readlines()
		refFile.close()
		size = 0
		for line in lines:
			if not line.startswith(">"):
				size += len(line.rstrip())

		for sequencePosition in range(size):
			self.cSequence[sequencePosition] = {"N" : 0, "A" : 0, "T" : 0, "C" : 0, "G" : 0}

		self.deltaArray = [] #filenames
		self.files = [f for f in listdir(startDir) if isfile(join(startDir+"/",f))]
		self.faArray = [] #filenames
		self.header = ""

	def addVote(self,sequencePosition,character):
		if sequencePosition in self.cSequence:
			self.cSequence[sequencePosition][character] += 1
		else:
			self.cSequence[sequencePosition] = {"n" : 0, "A" : 0, "T" : 0, "C" : 0, "G" : 0}
			self.cSequence[sequencePosition][character] += 1

	def complementBase(self,base):
		if base == "A":
			base = "T"
		elif base == "T":
			base = "A"
		elif base == "C":
			base = "G"
		elif base == "G":
			base = "C"
		return base

	def buildSequence(self,filename):
		outfile = open(filename,"w")
		self.header = ">" + os.path.realpath(filename).split("/")[-2]
		for i in range(1, 11):
			for f in self.files:
				if "." + str(i)+ "." in f:
					if "delta" in f:
						self.deltaArray.append(f)
					else:
						self.faArray.append(f)

		for i in range(len(self.deltaArray)):
			faStartDict = {} #{contig number : {instance : fa start position}}
			faEndDict = {} #{contig number : {instance : fa end position}}
			refStartDict = {} #{conting number : {instance : ref start position}}
			refEndDict = {} #{contig number : {instance : ref end position}}
			indelDict = {} #{contig number : {instance : [indel locations]}}
			coordinantsNext = False
			indelsNext = False
			for line in open(self.startDir + "/" + self.deltaArray[i], 'r'):
				if indelsNext is True:
					if int(line) is 0:
						indelsNext = False
						coordinantsNext = True
						instance += 1
					elif contigNumber in indelDict:
						if instance in indelDict[contigNumber]:
							if int(line) < 0 and indelDict[contigNumber][instance][-1] < 0:
								line = int(line) + indelDict[contigNumber][instance][-1]
							elif int(line) > 0 and indelDict[contigNumber][instance][-1] < 0:
								line = int(line) - indelDict[contigNumber][instance][-1]
							elif int(line) < 0 and indelDict[contigNumber][instance][-1] > 0:
								line = int(line) - indelDict[contigNumber][instance][-1]
							elif int(line) > 0 and indelDict[contigNumber][instance][-1] > 0:
								line = int(line) + indelDict[contigNumber][instance][-1]
							indelDict[contigNumber][instance].append(int(line))
						else:
							if int(line) < 0:
								if faStartDict[contigNumber][instance] < faEndDict[contigNumber][instance]:
									line = int(line) - int(faStartDict[contigNumber][instance])
								else:
									line = int(line) - int(faEndDict[contigNumber][instance])
							else:
								if faStartDict[contigNumber][instance] < faEndDict[contigNumber][instance]:
									line = int(line) + int(faStartDict[contigNumber][instance])
								else:
									line = int(line) + int(faEndDict[contigNumber][instance])
							indelDict[contigNumber].update({instance : [int(line)]})
					else:
						if int(line) < 0:
							if faStartDict[contigNumber][instance] < faEndDict[contigNumber][instance]:
								line = int(line) - int(faStartDict[contigNumber][instance])
							else:
								line = int(line) - int(faEndDict[contigNumber][instance])
						else:
							if faStartDict[contigNumber][instance] < faEndDict[contigNumber][instance]:
								line = int(line) + int(faStartDict[contigNumber][instance])
							else:
								line = int(line) + int(faEndDict[contigNumber][instance])
						indelDict[contigNumber] = {instance : [int(line) - 1]}
				elif ">" in line:
					instance = 0
					contigNumber = int(line.split(" ")[1].split("contig_")[1])
					coordinantsNext = True
				elif coordinantsNext is True:
					lineArray = line.split(" ")
					if contigNumber in faStartDict:
						refStartDict[contigNumber].update({instance : lineArray[0]})
						refEndDict[contigNumber].update({instance : lineArray[1]})
						faStartDict[contigNumber].update({instance : lineArray[2]})
						faEndDict[contigNumber].update({instance : lineArray[3]})
					else:
						refStartDict[contigNumber] = {instance : lineArray[0]}
						refEndDict[contigNumber] = {instance : lineArray[1]}
						faStartDict[contigNumber] = {instance : lineArray[2]}
						faEndDict[contigNumber] = {instance : lineArray[3]}
					coordinantsNext = False
					indelsNext = True
			for line in open(self.startDir + "/" + self.faArray[i], 'r'):
				if ">" in line:
					contigNumber = int(line.split("contig_")[1])
					basePosition = 0
				else:
					line = line.rstrip()
					for character in line:
						basePosition += 1
						for instance in faStartDict[contigNumber]:
							if basePosition >= int(faStartDict[contigNumber][instance]) and basePosition <= int(faEndDict[contigNumber][instance]):
								sequencePosition = basePosition - int(faStartDict[contigNumber][instance]) + int(refStartDict[contigNumber][instance])
								if contigNumber in indelDict:
									if instance in indelDict[contigNumber]:
										if basePosition in indelDict[contigNumber][instance]:
											if sequencePosition in self.cSequence:
												self.cSequence[sequencePosition]["N"] += 1
											else:
												self.cSequence[sequencePosition] = {"N" : 1, "A" : 0, "T" : 0, "C" : 0, "G" : 0}
										elif -basePosition in indelDict[contigNumber][instance]:
											indelDict[contigNumber][instance].remove(-basePosition)
											basePosition -= 1
											faEndDict[contigNumber][instance] = int(faEndDict[contigNumber][instance]) - 1 
										else:
											self.addVote(sequencePosition,character)
									else:
										self.addVote(sequencePosition,character)
								else:
									self.addVote(sequencePosition,character)
							if basePosition <= int(faStartDict[contigNumber][instance]) and basePosition >= int(faEndDict[contigNumber][instance]):
								sequencePosition = int(refEndDict[contigNumber][instance]) - (basePosition - int(faEndDict[contigNumber][instance]))
								if contigNumber in indelDict:
									if instance in indelDict[contigNumber]:
										if basePosition in indelDict[contigNumber][instance]:
											if sequencePosition in self.cSequence:
												self.cSequence[sequencePosition]["N"] += 1
											else:
												self.cSequence[sequencePosition] = {"N" : 1, "A" : 0, "T" : 0, "C" : 0, "G" : 0}
										elif -basePosition in indelDict[contigNumber][instance]:
											indelDict[contigNumber][instance].remove(-basePosition)
											basePosition -= 1
											faStartDict[contigNumber][instance] = int(faStartDict[contigNumber][instance]) - 1
										else:
											character = self.complementBase(character)
											self.addVote(sequencePosition,character)
									else:
										character = self.complementBase(character)
										self.addVote(sequencePosition, character)
								else:
									character = self.complementBase(character)
									self.addVote(sequencePosition, character)
		
		outfile.write(self.header + "\n")
		for i in range(1,max(self.cSequence.keys())+1):
			if i in self.cSequence.keys():
				keys = self.cSequence[i].keys()
				values = [self.cSequence[i][k] for k in keys]
				index,value = max(enumerate(values),key=operator.itemgetter(1))
				if value > 0:
					outfile.write(keys[index])
				else:
					outfile.write('N')

		outfile.close()

# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
