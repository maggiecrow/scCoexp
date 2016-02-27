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
                        thisTagfile1 = barcoded_dir + "/" + thisFlowcell + "/s_" + thisLane + "_2_sequence.0.r.1.txt.gz"
                        thisTagfile2 = barcoded_dir + "/" + thisFlowcell + "/s_" + thisLane + "_1_sequence.0.r.2.txt.gz"
                else:
                        thisSeqfile1 = barcoded_dir + "/" + thisFlowcell + "/s_" + thisLane + "_1_sequence." + thisRnaBarcodeGroup + thisRnaBarcode + ".r.1.txt.gz"
                        thisSeqfile2 = barcoded_dir + "/" + thisFlowcell + "/s_" + thisLane + "_2_sequence." + thisRnaBarcodeGroup + thisRnaBarcode + ".r.2.txt.gz"
                        thisTagfile1 = barcoded_dir + "/" + thisFlowcell + "/s_" + thisLane + "_2_sequence." + thisRnaBarcodeGroup + thisRnaBarcode + ".r.1.txt.gz"
                        thisTagfile2 = barcoded_dir + "/" + thisFlowcell + "/s_" + thisLane + "_1_sequence." + thisRnaBarcodeGroup + thisRnaBarcode + ".r.2.txt.gz"

		memoryNeeded = "16G"

				
                qsubFile = thisProcessedDir + "/" + thisSeqUnitId + ".r.add.varietal.tag01.qsub"
                bashFile = thisProcessedDir + "/" + thisSeqUnitId + ".r.add.varietal.tag01.bash"
                
                QSUB = open(qsubFile, "w")
                BASH = open(bashFile, "w")
                outline = '#$ -S /bin/bash\n'
                QSUB.write(outline)
                outline = "#$ -l virtual_free=" + memoryNeeded + "\n"
                QSUB.write(outline)
                outline = "/mnt/wigclust5/data/safe/kendall/cpppgms/tagmatch02.anirban01.exe <(gunzip -c " + thisTagfile1 + ") <(" + samtools_dir + "/samtools view -h " + thisMappedDir + "/" + thisSeqUnitId + ".r1/accepted_hits.bam) | " + samtools_dir + "/samtools view -Sb -o - - | " + samtools_dir + "/samtools sort - " + thisProcessedDir + "/" + thisSeqUnitId + ".r1.sorted\n"
                BASH.write(outline)
                outline = "python /mnt/wigclust5/data/safe/kendall/pythonpgms/coverage.tags.bins04.py <(" + samtools_dir + "/samtools view " + thisProcessedDir + "/" + thisSeqUnitId + ".r1.sorted.bam) /mnt/wigclust5/data/safe/kendall/ucscData/" + thisGenome + "/refFlat.20120902.summary.sorted.txt " + thisProcessedDir + "/" + thisSeqUnitId + ".r1." + thisGenome + ".refFlat.20120902.counts.04.txt > " + thisProcessedDir + "/" + thisSeqUnitId + ".r1." + thisGenome + ".refFlat.20120902.counts.04.stdout\n"
                BASH.write(outline)
                outline = "/mnt/wigclust5/data/safe/kendall/cpppgms/tagmatch02.anirban01.exe <(gunzip -c " + thisTagfile2 + ") <(" + samtools_dir + "/samtools view -h " + thisMappedDir + "/" + thisSeqUnitId + ".r2/accepted_hits.bam) | " + samtools_dir + "/samtools view -Sb -o - - | " + samtools_dir + "/samtools sort - " + thisProcessedDir + "/" + thisSeqUnitId + ".r2.sorted\n"
                BASH.write(outline)
                outline = "python /mnt/wigclust5/data/safe/kendall/pythonpgms/coverage.tags.bins04.py <(" + samtools_dir + "/samtools view " + thisProcessedDir + "/" + thisSeqUnitId + ".r2.sorted.bam) /mnt/wigclust5/data/safe/kendall/ucscData/" + thisGenome + "/refFlat.20120902.summary.sorted.txt " + thisProcessedDir + "/" + thisSeqUnitId + ".r2." + thisGenome + ".refFlat.20120902.counts.04.txt > " + thisProcessedDir + "/" + thisSeqUnitId + ".r2." + thisGenome + ".refFlat.20120902.counts.04.stdout\n"
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
