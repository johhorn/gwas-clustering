
class CreateLDTrack():
    """
    Functions for creating a LD linked point track based on LD information from static HyperBrowser file.
    """

    @classmethod
    def getLDDict(cls, r2):
        """
        Read LD information from static file

        :param r2: Threshold of rsquare, lower limit of LD correlation
        :return: Dictionary with SNPs and their LD variants
        """
        from quick.util.StaticFile import StaticFile
        staticFileFolder = StaticFile(['files', 'linkage_disequilibrium'])
        path = staticFileFolder.getDiskPath()
        LDFile = open(path + '/significant_expanded_ld.txt', 'r')
        edgeDict = cls._addEdges(LDFile, float(r2))
        LDFile.close()
        return edgeDict

    @classmethod
    def _addEdges(cls, LDfile, r2_threshold):
        """
        Takes in a file with the following tab-separated columns: chrnum  pos_ldSNP  ldSNP  pos_tagSNP  tagSNP  r2

        Will generate a dictionary where each element have the following format, and is accessed by a rsid key
        for the given SNP:

        {
          'start': start,
          'seqid': seqid,
          'edges': {edgeRsid: r2, ...}
        }

        Only stores variants in LD with higher rsquare than the given r2_threshold.

        :param LDfile: File of LD information.
        :param r2_threshold: Lower limit of LD correlation
        :return:
        """
        edgeDict = {}
        for line in LDfile:

            snp = line.strip().split('\t')

            tagSNP = snp[4]
            ldSNP = snp[2]

            if not tagSNP.startswith('rs') and not ldSNP.startswith('rs'):
                continue
            if tagSNP == ldSNP:
                continue

            chrnum = snp[0]
            tagPos = snp[3]
            ldPos = snp[1]
            r2 = float(snp[5])

            if r2 >= float(r2_threshold):
                cls._addEdge(edgeDict, rsid=tagSNP, start=tagPos, seqid=chrnum, edgeRsid=ldSNP, r2=r2)
                cls._addEdge(edgeDict, rsid=ldSNP, start=ldPos, seqid=chrnum, edgeRsid=tagSNP, r2=r2)

        return edgeDict

    @classmethod
    def _addEdge(cls, edgeDict, rsid, start, seqid, edgeRsid, r2):
        """
        Adds a variant, given by its rsid, to the given LD dictionary.
        The variant given by edgesRsid is in LD with the main variant, and is stored as an edge.

        The rsid will be the key of the entry, and its corresponding value represents the variant in the following
        format:

         {
          'start': start,
          'seqid': seqid,
          'edges': {edgeRsid: r2, ...}
        }

        If the rsid variant is already present in the dictionary, only information of the edge will be added to the
        dictionary.

        :param edgeDict: LD dictionary
        :param rsid: Rsid of variant
        :param start: Position of variant
        :param seqid: Chromosome number of variant
        :param edgeRsid: Correlated LD variant
        :param r2: rsquare value of the LD correlation between rsid and edgeRsid
        :return:
        """
        if rsid in edgeDict:
            if edgeRsid not in edgeDict[rsid]['edges']:
                edgeDict[rsid]['edges'][edgeRsid] = float(r2)
        else:
            edgeDict[rsid] = {
                'start': start,
                'seqid': seqid,
                'edges': {edgeRsid: float(r2)}
            }

    @classmethod
    def _getLDSNP(cls, snpInfo, rsid):
        """
         Get a linked point track formatted line for a specific rsid.

        :param snpInfo: The LD dictionary entry of the rsid
        :param rsid: rsid of track element to format
        :return:
        """

        start = snpInfo['start']
        seqid = 'chr' + str(snpInfo['seqid'])
        edges = snpInfo['edges']

        text_edges = ""
        separator = ""
        for snp in edges:

            edge = separator + snp + '=' + str(edges[snp])
            text_edges += edge
            separator = ';'

        if text_edges != "":
            return seqid + '\t' + start + '\t' + rsid + '\t' + text_edges + '\n'
        else:
            return ''

    @classmethod
    def getExpansionDict(cls, rsids, ldDict):
        """
        Get subset of the master LD dictionary, as created in getLDDict.
        Takes in list of rsids, and returns a dictionary only containing LD information between the given rsids.

        :param rsids: Rsids of interest
        :param ldDict: Master LD dictionary
        :return:
        """
        expansions = {}
        for rsid in rsids:

            if rsid not in ldDict:
                print 'Not in LD:', rsid
                continue

            snpInfo = ldDict[rsid]
            start = snpInfo['start']
            seqid = snpInfo['seqid']
            edges = snpInfo['edges']

            for edgeRsid, r2 in edges.items():
                cls._addEdge(expansions, rsid, start, seqid, edgeRsid, r2)
                cls._addEdge(expansions, edgeRsid, ldDict[edgeRsid]['start'], ldDict[edgeRsid]['seqid'], rsid, r2)

        return expansions

    @classmethod
    def formatLinkedPointTrack(cls, expansionDict, isUndirected):
        """
        Format a linked point track, given an expansionDict.

        :param expansionDict: LD dictionary of all variants that should be formatted into linked point track
        :param isUndirected: Boolean parameter of whether or not the graph is undirected
        :return:
        """
        output = "##gtrack version: 1.0\n" \
                 "##track type: linked points\n" \
                 "##undirected edges: " + str(isUndirected) + "\n" \
                 "##edge weights: True\n" \
                 "##edge weight type: number\n" \
                 "##edge weight dimension: scalar\n" \
                 "##uninterrupted data lines: true\n" \
                 "##no overlapping elements: true\n" \
                 "##1-indexed: False\n" \
                 "###seqid\tstart\tid\tedges\n"

        for rsid, snpInfo in expansionDict.items():
            output += cls._getLDSNP(snpInfo, rsid)

        return output


    @classmethod
    def _getSNP(cls, rsid, rsidDict):
        from quick.webtools.clustering.RsidMapper import RsidMapper

        seqid, start = RsidMapper.getPosition(rsid, rsidDict)
        if seqid and start:
            return "\t".join([seqid, start, rsid]) + "\n"
        else:
            return ""


    @classmethod
    def formatPointTrack(cls, expansionDict, rsidDict, originalSNPs):
        """
        Format a point track, given an expansionDict.

        :param expansionDict: LD dictionary of all variants that should be formatted into linked point track
        """
        output = "##gtrack version: 1.0\n" \
                 "##track type: points\n" \
                 "##no overlapping elements: true\n" \
                 "###seqid\tstart\tid\n"

        for rsid in expansionDict.keys():
            output += cls._getSNP(rsid, rsidDict)

        # Make sure all original SNPs are present in the new track.
        for rsid in originalSNPs:
            if rsid not in expansionDict:
                output += cls._getSNP(rsid, rsidDict)

        return output


    @classmethod
    def parseFileIntoPointTrack(cls, inFn, outFn, ldDict, rsidDict):
        """
        Loops through a primary track and creates a new linked point track for the given track elements.
        The primary track must have the column header 'snps', whose column elements are rsids.

        :param inFn: Path to original track
        :param outFn: Path to new linked point track (LD graph)
        :param ldDict: Master LD dictionary
        :param edgeDir: Boolean parameter of whether or not the graph is undirected
        :return:
        """
        from quick.util.CommonFunctions import ensurePathExists
        ensurePathExists(outFn)
        inFile = open(inFn, 'r')
        outFile = open(outFn, 'w')
        rsids = cls.getUniqueRsids(inFile)
        expansionDict = CreateLDTrack.getExpansionDict(rsids, ldDict)
        outFile.write(CreateLDTrack.formatPointTrack(expansionDict, rsidDict, rsids))

        inFile.close()
        outFile.close()

    @classmethod
    def getUniqueRsids(cls, gtrackFile):
        """
        Takes in a matrix, where each row is a track line, each column a track attribute.
        Assumes a format where rsids is in the 'snps' column.
        Returns a list of unique rsids.
        """
        inFileLines = [x.strip().split('\t') for x in gtrackFile.readlines()]

        rsids = []
        colIndex = None
        for line in inFileLines:
            if line[0].startswith('###'):
                colIndex = line.index('snps')
            elif line[0].startswith('#'):
                continue
            else:
                rsids.append(line[colIndex])

        return list(set(rsids))
