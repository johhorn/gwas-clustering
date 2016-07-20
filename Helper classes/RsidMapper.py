
class RsidMapper(object):
    """
    Uses static files stored in a folder /software/galaxy/galaxy_clustering/static/hyperbrowser/files/ + <genome>_bed
    to create a mapping between variant rsids and their positions.

    Returns a dictionary with rsids as keys, and a tuple of chromosome and start position for as corresponding value.
    This dictionary can further be used to uniformly map rsids of the same or different tracks (for instance within a
    GSuite) to the same reference genome.

    Is only defined for SNPs.
    """

    @classmethod
    def _getChromCount(cls):
        return 24

    @classmethod
    def _getChromFilenames(cls, genome):
        """
        Get all 24 chromosome file names for the given reference genome.
        Each file contain all dbSNP rsids for the chromosome it represents.

        :param genome: Chosen reference genome
        :return: List of dbSNP chromosome file names
        """
        from quick.util.StaticFile import StaticFile
        staticFileFolder = StaticFile(['files', genome])
        path = staticFileFolder.getDiskPath()
        chrIndexes = range(1, cls._getChromCount() - 1) + ['X', 'Y']
        return [path + '/bed_chr_' + str(x) + '.bed' for x in chrIndexes]

    @classmethod
    def createRsidMappingFromStaticFiles(cls, progressViewer, genome='hg38'):
        """

        :param progressViewer: HyperBrowser object for indicating progress in reading chromosome files
        :param genome: Chosen reference genome for rsid mapping
        :return: Dictionary with rsid keys and (chr, pos) values
        """
        fileNames = cls._getChromFilenames(genome)

        rsidMap = {}
        for fileName in fileNames:
            chrFile = open(fileName, 'r')
            chrFileLines = [x.strip().split('\t') for x in chrFile]

            for snp in chrFileLines[1:]:
                chr = snp[0]
                start = snp[1]
                stop = snp[2]
                rsid = snp[3]

                if int(stop) - int(start) > 1:  # Only keep SNPs
                    continue

                rsidMap[rsid] = (chr, start)

            progressViewer.update()

        return rsidMap

    @classmethod
    def getPosition(cls, rsid, rsidDict):
        if rsid in rsidDict:
            return rsidDict[rsid]
        else:
            return None, None






