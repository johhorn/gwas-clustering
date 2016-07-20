# Author: Johanne
from gold.statistic.GraphStat import GraphStat
from gold.statistic.MagicStatFactory import MagicStatFactory
from gold.statistic.RawDataStat import RawDataStat
from gold.statistic.Statistic import Statistic
from gold.track.TrackFormat import TrackFormatReq


class ExpandWithLDVariantStat(MagicStatFactory):
    """
    For a track, return a dictionary of mapping between rsids, where the value tagSNP represent the LD cluster
    the key SNP is part of.

    Used mainly in createChildren of other statistics.
    """
    pass


class ExpandWithLDVariantStatUnsplittable(Statistic):

    def _init(self, rsquareLimit):
        self._rsquareLimit = rsquareLimit

    def _compute(self):
        track = self._children[0].getResult()
        graph = self._children[1].getResult()

        rsids = track.extrasAsNumpyArray('snps')
        positions = track.startsAsNumpyArray()

        return self.createLDClusterDict(graph, positions, rsids)

    def createLDClusterDict(self, graph, positions, rsids):
        """
        Create dictionary where SNPs in LD within the track is mapped to the same tagSNP that represent the cluster
        """
        expansionDict = {}
        keys = {}
        for index in range(0, len(positions)):
            rsid = rsids[index]
            pos = positions[index]
            if rsid in keys:
                rsid = keys[rsid]
            else:
                keys[rsid] = rsid

            if pos not in expansionDict:
                expansionDict[pos] = rsid

            if not graph.hasNode(rsid):
                continue

            node = graph.getNode(rsid)

            for neighborEdge in node.getNeighborIter():

                if neighborEdge.weight < self._rsquareLimit:
                    continue
                neighbor = neighborEdge.toNode
                neighborPos = neighbor.start()
                neighborRsid = neighbor.id()

                if neighborPos not in expansionDict:
                    expansionDict[neighborPos] = rsid

                if neighborRsid not in keys:
                    keys[neighborRsid] = rsid
        return expansionDict

    def _createChildren(self):
        self._addChild(RawDataStat(self._region, self._track, TrackFormatReq(allowOverlaps=True)))
        self._addChild(GraphStat(self._region, self._track2))
