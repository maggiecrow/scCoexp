#!/usr/bin/env python

import sys
import time
from operator import itemgetter
import os
import subprocess

def main():

	###  input guide table sorted by sample id. library id, rna barcode group, rna barcode.
	
	"""
	EXAMPLE COMMAND LINE (2 COMMANDS!)
	TAB=`echo -e "\t"`
	python aggregate.varietal.tag.counts03.strand.py <(cat <(head -n 1 ../meta/anirban12.guide.20150522.txt) <(sort -t"$TAB" -k 3,3 -k 16,16 -k 10,10 -k 11,11 <(tail -n +2 ../meta/anirban12.guide.20150522.txt)))
	"""
	
	
	infilename = sys.argv[1]


	python_dir = "/mnt/wigclust5/data/safe/kendall/Python-2.7.1"
	bowtie_dir = "/mnt/wigclust5/data/safe/kendall/bowtie-0.12.7"
	tophat2_dir = "/mnt/wigclust5/data/safe/kendall/tophat-2.0.4.Linux_x86_64"
	bowtie2_dir = "/mnt/wigclust5/data/safe/kendall/bowtie2-2.0.0-beta6"
	samtools_dir = "/mnt/wigclust5/data/safe/kendall/samtools-0.1.16"
	project_dir = "/mnt/wigclust5/data/safe/paula/anirban12"
	barcoded_dir = "/mnt/wigclust5/data/safe/paula/anirban12/barcode.split"
	#infilename = project_dir + "/mm9/meta/anirban12.guide.20150415.txt"
	
	guide = fileToGuide(infilename)
	ncells = len(guide["seq.unit.id"])
	
	commandTemplate = 'python /mnt/wigclust5/data/safe/kendall/pythonpgms/coverage.tags.bins04.strand.py <(/mnt/wigclust5/data/safe/kendall/samtools-0.1.16/samtools merge - <MERGE_LIST> | /mnt/wigclust5/data/safe/kendall/samtools-0.1.16/samtools view -) /mnt/wigclust5/data/safe/kendall/ucscData/<GENOME>/refFlat.20120902.summary.sorted.txt /mnt/wigclust5/data/safe/paula/anirban06/<GENOME>/processed/<CELL_ID>.<GENOME>.refFlat.20120902.counts.04.strand.txt > /mnt/wigclust5/data/safe/paula/anirban06/<GENOME>/processed/<CELL_ID>.<GENOME>.refFlat.20120902.counts.04.strand.stdout'
	
	prevMergeList = ""
	prevCellId = ""
	prevGenome = ""
		
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
		thisLib = guide["lib"][i]
		thisSampleName = guide["sample.name"][i]

		thisCellId = thisSampleName + "_" + thisLib + "_" + thisRnaBarcodeGroup + thisRnaBarcode

		thisMappedDir = project_dir + "/" + thisGenome + "/mapped"
		thisProcessedDir = project_dir + "/" + thisGenome + "/processed"

		thisBamfile1 = thisProcessedDir + "/" + thisSeqUnitId + ".r1.sorted.bam"
		thisBamfile2 = thisProcessedDir + "/" + thisSeqUnitId + ".r2.sorted.bam"

        	if thisCellId == prevCellId:
        		if prevMergeList == "":
        			prevMergeList = thisBamfile1 + " " + thisBamfile2
        		else:
        			prevMergeList += (" " + thisBamfile1 + " " + thisBamfile2)
		else:
        		if prevCellId == "":
        			pass
        		else:
        			a = commandTemplate.replace("<MERGE_LIST>", prevMergeList)
        			a = a.replace("<CELL_ID>", prevCellId)
        			a = a.replace("<GENOME>", prevGenome)
        			
				memoryNeeded = "16G"

				qsubFile = thisProcessedDir + "/avtc_" + prevCellId + ".r.aggregate.varietal.tag.counts03.strand.qsub"
				bashFile = thisProcessedDir + "/avtc_" + prevCellId + ".r.aggregate.varietal.tag.counts03.strand.bash"

				QSUB = open(qsubFile, "w")
				BASH = open(bashFile, "w")
				outline = '#$ -S /bin/bash\n'
				QSUB.write(outline)
				outline = "#$ -l virtual_free=" + memoryNeeded + "\n"
				QSUB.write(outline)
				outline = a + "\n"
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
				print thisCommand
				os.system(thisCommand)

			prevMergeList = thisBamfile1 + " " + thisBamfile2

		prevCellId = thisCellId
		prevGenome = thisGenome
		

	if prevCellId == "":
		pass
	else:
		a = commandTemplate.replace("<MERGE_LIST>", prevMergeList)
		a = a.replace("<CELL_ID>", prevCellId)
		a = a.replace("<GENOME>", prevGenome)

		memoryNeeded = "16G"

		qsubFile = thisProcessedDir + "/avtc_" + prevCellId + ".r.aggregate.varietal.tag.counts03.strand.qsub"
		bashFile = thisProcessedDir + "/avtc_" + prevCellId + ".r.aggregate.varietal.tag.counts03.strand.bash"

		QSUB = open(qsubFile, "w")
		BASH = open(bashFile, "w")
		outline = '#$ -S /bin/bash\n'
		QSUB.write(outline)
		outline = "#$ -l virtual_free=" + memoryNeeded + "\n"
		QSUB.write(outline)
		outline = a + "\n"
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
		print thisCommand
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
