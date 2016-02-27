#!/usr/bin/env python

import sys
import time
from operator import itemgetter
import os
import subprocess

def main():

	python_dir = "/mnt/wigclust5/data/safe/kendall/Python-2.7.1"
	bowtie_dir = "/mnt/wigclust5/data/safe/kendall/bowtie-0.12.7"
	tophat2_dir = "/mnt/wigclust5/data/safe/kendall/tophat-2.0.4.Linux_x86_64"
	bowtie2_dir = "/mnt/wigclust5/data/safe/kendall/bowtie2-2.0.0-beta6"
	samtools_dir = "/mnt/wigclust5/data/safe/kendall/samtools-0.1.16"
	project_dir = "/mnt/wigclust5/data/safe/paula/anirban12"
	barcoded_dir = "/mnt/wigclust5/data/safe/paula/anirban12/barcode.split"
	infilename = project_dir + "/mm9/meta/anirban12.guide.20150522.txt"
	
	
	guide = fileToGuide(infilename)
	ncells = len(guide["seq.unit.id"])
	for i in range(ncells):
		#if i > 158:
		#	break
                if guide["process"][i] == "0":
                        continue
                        
                thisGenome = guide["genome"][i]
                thisDnaRna = guide["dna.rna"][i]
                thisRnaBarcodeGroup = guide["rna.barcode.group"][i]
                thisRnaBarcode = guide["rna.barcode"][i]
                thisSeqUnitId = guide["seq.unit.id"][i]
                thisFlowcell = guide["flowcell"][i]
                thisLane = guide["lane"][i]
                
                thisMappedDir = project_dir + "/" + thisGenome + "/mapped"
                thisProcessedDir = project_dir + "/" + thisGenome + "/processed"
                thisTophatWorkDir1 = thisMappedDir + "/" + thisSeqUnitId + ".r1"
                thisTophatWorkDir2 = thisMappedDir + "/" + thisSeqUnitId + ".r2"
                
                if thisRnaBarcode == "":
                        thisSeqfile1 = barcoded_dir + "/" + thisFlowcell + "/s_" + thisLane + "_1_sequence.0.r.1.txt.gz"
                        thisSeqfile2 = barcoded_dir + "/" + thisFlowcell + "/s_" + thisLane + "_2_sequence.0.r.2.txt.gz"
                elif thisRnaBarcode == "0":
                        thisSeqfile1 = barcoded_dir + "/" + thisFlowcell + "/s_" + thisLane + "_1_sequence.0.r.1.txt.gz"
                        thisSeqfile2 = barcoded_dir + "/" + thisFlowcell + "/s_" + thisLane + "_2_sequence.0.r.1.txt.gz"
                else:
                        thisSeqfile1 = barcoded_dir + "/" + thisFlowcell + "/s_" + thisLane + "_1_sequence." + thisRnaBarcodeGroup + thisRnaBarcode + ".r.1.txt.gz"
                        thisSeqfile2 = barcoded_dir + "/" + thisFlowcell + "/s_" + thisLane + "_2_sequence." + thisRnaBarcodeGroup + thisRnaBarcode + ".r.2.txt.gz"

		if os.path.isdir(thisTophatWorkDir1):
			pass
		else:
                        thisCommand = "mkdir " + thisTophatWorkDir1
			os.system(thisCommand)

		if os.path.isdir(thisTophatWorkDir2):
			pass
		else:
                        thisCommand = "mkdir " + thisTophatWorkDir2
			os.system(thisCommand)
		
		memoryNeeded = "16G"
		
		#if os.path.isfile(thisSeqfile):
		#	memoryNeeded = str( (os.path.getsize(thisSeqfile) / 400000000 ) + 5) + "G"
		#else:
		#	memoryNeeded = "4G"
		
		#thisCommand = "gunzip -c " + thisSeqfile + " | python solexa.quals01.py"
		#print thisCommand
		#p = subprocess.Popen(thisCommand, shell=True, stdout=subprocess.PIPE)		
		#thisSolexaQuals = p.communicate()[0].rstrip()
		
		#print "thisSolexaQuals ", thisSolexaQuals
		
                qsubFile = thisMappedDir + "/" + thisSeqUnitId + ".r.tophat2.qsub"
                bashFile = thisMappedDir + "/" + thisSeqUnitId + ".r.tophat2.bash"
                
                QSUB = open(qsubFile, "w")
                BASH = open(bashFile, "w")
                outline = '#$ -S /bin/bash\n'
                QSUB.write(outline)
                #outline = "#$ -l virtual_free=" + memoryNeeded + "\n"
                outline = "#$ -l virtual_free=4G\n"
                QSUB.write(outline)
                outline = '#$ -pe threads 7\n'
                QSUB.write(outline)
                outline = "export PATH=$PATH:" + bowtie2_dir + ":" + samtools_dir + "\n"
                BASH.write(outline)
                if thisGenome == "hg19":
                	outline = tophat2_dir +"/tophat2 -p 7 -o " + thisTophatWorkDir1 + " -G /mnt/wigclust5/data/safe/kendall/sequences/hg19.ucsc.refseq.gtf /mnt/wigclust5/data/safe/kendall/bowtie2-2.0.0-beta6/indexes/hg19 <(gunzip -c " + thisSeqfile1 + " )\n"
                elif thisGenome == "mm9":
                	outline = tophat2_dir +"/tophat2 -p 7 -o " + thisTophatWorkDir1 + " -G /mnt/wigclust5/data/safe/kendall/anirban05/mm9ercc92.refseq.gtf /mnt/wigclust5/data/safe/kendall/bowtie2-2.0.0-beta6/indexes/mm9ercc92 <(gunzip -c " + thisSeqfile1 + " )\n"                
                else:
                	outline = "## UNKNOWN GENOME"
                BASH.write(outline)
                if thisGenome == "hg19":
                	outline = tophat2_dir +"/tophat2 -p 7 -o " + thisTophatWorkDir2 + " -G /mnt/wigclust5/data/safe/kendall/sequences/hg19.ucsc.refseq.gtf /mnt/wigclust5/data/safe/kendall/bowtie2-2.0.0-beta6/indexes/hg19 <(gunzip -c " + thisSeqfile2 + " )\n"
                elif thisGenome == "mm9":
                	outline = tophat2_dir +"/tophat2 -p 7 -o " + thisTophatWorkDir2 + " -G  /mnt/wigclust5/data/safe/kendall/anirban05/mm9ercc92.refseq.gtf /mnt/wigclust5/data/safe/kendall/bowtie2-2.0.0-beta6/indexes/mm9ercc92 <(gunzip -c " + thisSeqfile2 + " )\n"
                else:
                	outline = "## UNKNOWN GENOME"
                BASH.write(outline)
                outline = "bash " + bashFile + "\n"
                QSUB.write(outline)
                QSUB.close()
                BASH.close()

                time.sleep(1)
                thisCommand = "chmod 755 " + bashFile
                os.system(thisCommand)
                thisCommand = "chmod 755 " + qsubFile
                os.system(thisCommand)
                #thisCommand = "qsub -q all.q@wigclust5 " + qsubFile
                thisCommand = "qsub " + qsubFile
                os.system(thisCommand)



def fileToGuide(infilename):
	INFILE = open(infilename, "r")
	
	guide = dict()
	x = INFILE.readline()
	colnames = x.rstrip().split("\t")
	for c in colnames:
		c = c.strip()
	ncols = len(colnames)

	a = []
	for x in INFILE:
		arow = x.rstrip().split("\t")
		for b in arow:
			b = b.strip()
		if len(arow) < ncols:
			n = len(arow)
			for i in range(n, ncols):
				arow.append("")
		a.append(arow)

	z = zip(*a)
	if len(z) == ncols:
		for i in range(ncols):
			if guide.has_key(colnames[i]):
				print "ERROR: Duplicate column name."
			else:
				guide[colnames[i]] = z[i]
	else:
		print "ERROR: Data and colnames length not equal."
		
	INFILE.close()
	return guide


if __name__ == "__main__":
	main()
