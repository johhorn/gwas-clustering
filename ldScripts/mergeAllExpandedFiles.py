
import sys, os
import string
"""
Inputfile contains a dump of rsids that has been extracted from the HyperBrowser.
Inputfolder contains files of expanded SNPs from the inputfile.

Will find rsids that have not been expanded.
"""
path = sys.argv[1]
inputfile = sys.argv[2]

rsidsFile = open(inputfile, 'r')
alternateFile = open('unaccounted_new.txt', 'w')
unaccountedFile = open('unaccounted_old.txt', 'w')
finishedFile = open('expanded_with_ld.txt', 'w')

finishedFile.write('chrnum\tldPos\tldSNP\ttagPos\ttagSNP\tr2\n')

rsids_original= {}
rsids_expanded = {}

for line in rsidsFile:
	rsid = line.strip()
	rsids_original[rsid] = 1

not_expanded = set()

for root, dirs, files in os.walk(path, topdown=False):
	for name in files:

		if not name.startswith('out'):
			continue

		inputfile = os.path.join(root, name)
		SNPfile = open(inputfile, 'r')

		for snp in SNPfile:
			ld = snp.strip().split('\t')
			chrnum = ld[0]
			ldPos = ld[1]
			ldSNP = ld[2]
			tagPos = ld[3]
			tagSNP = ld[4]
			r2 = ld[5]
			
			if tagSNP.startswith('rs') and tagSNP not in rsids_original:
				not_expanded.add(tagSNP)

			if tagSNP.startswith('rs') and ldSNP.startswith('rs'):
				rsids_expanded[tagSNP] = 1
				finishedFile.write(chrnum + "\t" + ldPos + "\t" + ldSNP + "\t" + tagPos + "\t" + tagSNP + "\t" + r2 + "\n")

		SNPfile.close()

for tagSNP in rsids_original.keys():
	if tagSNP not in rsids_expanded:
		unaccountedFile.write(tagSNP + '\n')

for snp in not_expanded:
	alternateFile.write(snp + "\n")

alternateFile.close()
unaccountedFile.close()
finishedFile.close()