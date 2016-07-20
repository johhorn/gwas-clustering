
import sys
import string
"""
Inputfile contains a dump of rsids that has been extracted from the HyperBrowser.

Will create multiple files, each with a subset of the rsids, 
for expansion with snps in ld.
"""
inputfile = sys.argv[1]

diseaseFile = open(inputfile, 'r')
outputFile = open('rsquare_rsids/rsids1.txt', 'w')
i = 0
j = 0
rsids = []

for line in diseaseFile:
	outputFile.write(line)	
	if i == 10:
		j += i
		outputFile.close()
		filename = 'rsquare_rsids/rsids' + str(j) + '.txt'
		outputFile = open(filename, 'w')
		i = 0
	else:
		i += 1
outputFile.close()
diseaseFile.close()