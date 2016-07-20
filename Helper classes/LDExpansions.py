class LDExpansions(object):
    """
    Functions for creating and using an LD graph.

    In addition, it has functionality for finding positions from an LD graph track (linked point track).
    """

    @classmethod
    def createRSquareGraph(cls, ldGraphTrackName, r2_threshold):
        """
        Creates a dictionary of all pairs in a linked point track.
        Variants in LD must have rsquare >= the rsquare threshold passed to the function.

        :param ldGraphTrackName: linked point track, as chosen in tool (choices.ldtrack)
        :param r2_threshold: Lower limit of square value
        :return: Dictionary of all ld-pairs with sorted key = (rsid1, rsid2), value = rSquare
        """
        from quick.application.ExternalTrackManager import ExternalTrackManager
        from gold.origdata.GtrackGenomeElementSource import GtrackGenomeElementSource

        fileName = ExternalTrackManager.extractFnFromGalaxyTN(ldGraphTrackName)
        suffix = ExternalTrackManager.extractFileSuffixFromGalaxyTN(ldGraphTrackName)
        gtSource = GtrackGenomeElementSource(fileName, suffix=suffix)

        r2graph = {}

        for ge in gtSource:
            rsid = ge.id
            edges = ge.edges
            weights = ge.weights

            for i in range(0, len(edges)):
                ldRsid = edges[i]
                r2 = weights[i]

                if r2 >= float(r2_threshold):
                    cls.addEdge(r2graph, rsid, ldRsid, r2)

        return r2graph

    @classmethod
    def sortRsidTuple(cls, rsid1, rsid2):
        id1 = int(rsid1[2:])
        id2 = int(rsid2[2:])
        return (rsid1, rsid2) if id1 < id2 else (rsid2, rsid1)

    @classmethod
    def addEdge(cls, r2graph, rsid1, rsid2, r2):
        key = cls.sortRsidTuple(rsid1, rsid2)
        if key not in r2graph:
            r2graph[key] = r2

    @classmethod
    def getEdge(cls, rsid1, rsid2, graph):
        """
        Sorts the node ids, i.e. rsids, to get uniform key and finds the edge value between the nodes.
        :param rsid1: Node id 1
        :param rsid2: Node id 2
        :param graph: Graph to find node pair in
        :return: rsquare value between nodes, 1 if completely similar, 0 if no correlation (not in graph)
        """
        if rsid1 == rsid2:
            return 1

        key = cls.sortRsidTuple(rsid1, rsid2)
        if key in graph:
            return graph[key]
        else:
            return 0

    @classmethod
    def generateTracksAndLabels(cls, gSuite, analysisBins):
        """
        Used in bipartite matching along with rSquare graph.
        For each track, a list of its rsids is generated. A list of these lists, along with a list of track labels,
        is returned.
        """
        from gold.application.HBAPI import doAnalysis
        from gold.description.AnalysisDefHandler import AnalysisSpec
        from quick.statistic.UniquePointTrackStat import UniquePointTrackStat
        from gold.track.Track import Track

        analysisSpec = AnalysisSpec(UniquePointTrackStat)

        tracks = []
        labels = []

        for gSuiteTrack in gSuite.allTracks():
            track = Track(gSuiteTrack.trackName)
            result = doAnalysis(analysisSpec, analysisBins, [track]).getGlobalResult()
            if 'Result' in result:
                tracks.append(result['Result'])
                labels.append(gSuiteTrack.title)

        return tracks, labels

    @classmethod
    def createPositionDict(cls, ldGraphTrackName):
        """
        Creates position dictionary from linked point track. To be used for empiric exploration of positions,
        based on LD correlation (rsquare values)
        :param ldGraphTrackName: linked point track, as chosen in tool (choices.ldtrack)
        :return: Dictionary of all nodes in track with key = rsid, value = position

        """
        from quick.application.ExternalTrackManager import ExternalTrackManager
        from gold.origdata.GtrackGenomeElementSource import GtrackGenomeElementSource

        fileName = ExternalTrackManager.extractFnFromGalaxyTN(ldGraphTrackName)
        suffix = ExternalTrackManager.extractFileSuffixFromGalaxyTN(ldGraphTrackName)
        gtSource = GtrackGenomeElementSource(fileName, suffix=suffix)

        positionDict = {}

        for ge in gtSource:
            rsid = ge.id
            position = ge.start
            cls.addPosition(positionDict, rsid, position)

        return positionDict

    @classmethod
    def addPosition(cls, positionDict, rsid, position):
        if rsid not in positionDict:
            positionDict[rsid] = position

    @classmethod
    def getPosition(cls, positionDict, rsid):
        if rsid in positionDict:
            return positionDict[rsid]
        else:
            return -1
