from quick.util.CommonFunctions import ensurePathExists
from quick.webtools.clustering.RsidMapper import RsidMapper


class GSuitePrimaryTrackModifier():
    """
    Functions for conversion of different primary tracks in a GSuite.
    """

    # Valued point track column headers
    RSID = 'id'
    SEQID = 'seqid'
    POS = 'start'
    VAL = 'value'

    # Sumstat column headers
    P = 'P'
    Z = 'Z'
    SNP = 'SNP'

    GTRACK_COLS = [SEQID, POS, RSID, VAL]

    @classmethod
    def _getAttributeNames(cls):
        return cls.GTRACK_COLS

    @classmethod
    def convertSumstat(cls, inFn, outFn, rsidDict, shouldLogTransform, valueFilter=None):
        """
        Converts a sumstat primary track to a valued point track.
        The sumstat track must be formatted as defined in the summary statistic file definition of Bulik-Sullivan et
        al., (2015). See specification at https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format

        :param inFn: path to original track
        :param outFn: path to new track
        :param rsidDict: dictionary for mapping of rsids to reference genome
        :param valueFilter: upper threshold for p/z values in original track
        """
        from math import log

        ensurePathExists(outFn)
        inFile = open(inFn, 'r')
        outFile = open(outFn, 'w')

        # Find columns
        inFileLines = [x.strip().split('\t') for x in inFile.readlines()]
        colNames = [col.upper() for col in inFileLines[0]]
        valueColNum = colNames.index(cls.P) if cls.P in colNames else colNames.index(cls.Z)
        idColNum = colNames.index(cls.SNP)

        # Track header information
        outFile.write('##track type: valued points\n')
        outFile.write("##1-indexed: False\n")
        outFile.write('###' + '\t'.join(cls.GTRACK_COLS) + '\n')

        # Convert each line
        for cols in inFileLines[1:]:

            rsid = cols[idColNum]
            value = cols[valueColNum]

            if not valueFilter or float(value) <= valueFilter:

                if shouldLogTransform:
                    try:
                        """Convert values to -log(pval), the values of GWAS Catalog SNPs"""
                        value = str(-log(float(value)))
                    except:
                        """For SNPs with reported p-value of 0.000, assume high significance"""
                        value = str(-log(0.0005))

                seq, pos = RsidMapper.getPosition(rsid, rsidDict)
                if seq and pos:
                    outFile.write('\t'.join([seq, pos, rsid, value]) + '\n')

        inFile.close()
        outFile.close()

    @classmethod
    def liftOverGTrack(cls, inFn, outFn, rsidDict):
        """
        Liftover for primary point tracks. The tracks must have a column 'id', with the rsid of the SNPs in each row.
        In addition, 'seqid' and 'start' is needed in the original tracks, as these columns will be the only ones
        modified for each track element.

        :param inFn: path to original track
        :param outFn: path to new track
        :param rsidDict: dictionary for mapping of rsids to reference genome
        """
        ensurePathExists(outFn)
        inFile = open(inFn, 'r')
        outFile = open(outFn, 'w')

        rsidCol = 0
        seqCol = 0
        startCol = 0

        # Lift over each line
        for line in inFile.readlines():
            if line.startswith('###'):
                cols = line[3:].strip().split('\t')
                rsidCol = cols.index(cls.RSID)
                seqCol = cols.index(cls.SEQID)
                startCol = cols.index(cls.POS)

            if line.startswith("##1-indexed:"):
                """
                The rsID-mapping is based on the dbSNP positions, which are 0-indexed.
                We need to make sure this attribute is correctly set in our tracks.
                """
                outFile.write("##1-indexed: False\n")

            elif line.startswith('#'):
                outFile.write(line)

            else:
                cols = line.strip().split('\t')
                rsid = cols[rsidCol]
                seq, pos = RsidMapper.getPosition(rsid, rsidDict)
                if seq and pos:
                    cols[seqCol] = str(seq)
                    cols[startCol] = pos
                    outFile.write('\t'.join(cols) + '\n')

        inFile.close()
        outFile.close()
