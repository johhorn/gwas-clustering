# Author: Johanne
from gold.statistic.MagicStatFactory import MagicStatFactory
from gold.statistic.Statistic import StatisticSplittable, MultipleTrackStatistic
from quick.statistic.ExpandWithLDVariantStat import ExpandWithLDVariantStat


class ExpandTrackAndMatchStat(MagicStatFactory):
    """
    Statistic which calculates positive matches (1, 1) and mismatches (0, 1) or (1, 0)
    between two tracks. Is made for valued point tracks, and return a genome-wide count.

    Each binary feature represents an LD cluster.

    Can be used for distance computations as described in the article "A survey of Binary Similarity and Distance
    Measures" by Choi, Cha and Tappert.

    Returns a dictionary containing scores for all these cases. The dictionary have the following keys:
    a: Overlap between LD clusters
    b: Count of unmatched LD clusters in track1
    c: Count of unmatched LD clusters in track2
    d: Negative matches, undefined (-1)

    See example of usage in quick/webtools/clustering/LDExpansionClusteringTool.py
    """
    pass


class ExpandTrackAndMatchStatSplittable(StatisticSplittable):
    def _combineResults(self):
        countDict = {'a': 0, 'b': 0, 'c': 0, 'd': -1}
        for child in self._childResults:
            countDict['a'] += child[0]
            countDict['b'] += child[1]
            countDict['c'] += child[2]

        return countDict


class ExpandTrackAndMatchStatUnsplittable(MultipleTrackStatistic):
    """
    The expandedTracks are dictionaries which map LD structure between rsids.
    A positive match is counted as unique mappings between two ldClusters. This mapping is represented by the pair of
    LD cluster tagSNPs, in sorted order.

    For the region, a, b, c and d is computed and returned.
    """

    def _getTagSNPCount(self, expandedTrack):
        rsids = set(expandedTrack.values())
        return len(rsids)

    def _getRsid(self, dict, key):
        if key in dict:
            return dict[key]
        else:
            return None

    def _sortRsidTuple(self, rsid1, rsid2):
        id1 = int(rsid1[2:])
        id2 = int(rsid2[2:])
        return (rsid1, rsid2) if id1 < id2 else (rsid2, rsid1)

    def _compute(self):
        # Get LD cluster dictionaries
        expandedTrack = self._children[0].getResult()
        expandedTrack2 = self._children[1].getResult()

        # Get count of unmatched LD clusters
        filteredCount = self._getTagSNPCount(expandedTrack)
        filteredCount2 = self._getTagSNPCount(expandedTrack2)

        a = self.findLDClusterMatches(expandedTrack, expandedTrack2)
        b = filteredCount - a
        c = filteredCount2 - a
        d = self._region.getTotalBpSpan() - a - b - c
        return a, b, c, d

    def findLDClusterMatches(self, expandedTrack, expandedTrack2):
        """
        Traverses through the smallest track and finds all positive matches between LD clusters
        """
        loopingDict = expandedTrack if len(expandedTrack) <= len(expandedTrack2) else expandedTrack2

        matches = set()
        for position in loopingDict:

            rsid1 = self._getRsid(expandedTrack, position)
            rsid2 = self._getRsid(expandedTrack2, position)

            if rsid1 is None or rsid2 is None:
                continue

            match = self._sortRsidTuple(rsid1, rsid2)
            matches.add(match)

        return len(matches)

    def _createChildren(self):
        rsquareLimit = 0.8  # Default value if no limit is set (consensus lower limit of rsquare)

        if 'rsquareLimit' in self._kwArgs:
            rsquareLimit = self._kwArgs['rsquareLimit']

        linkedPointTrack = self._tracks[2]

        self._addChild(ExpandWithLDVariantStat(
            region=self._region,
            track=self._track,
            track2=linkedPointTrack,
            rsquareLimit=rsquareLimit
        ))

        self._addChild(ExpandWithLDVariantStat(
            region=self._region,
            track=self._track2,
            track2=linkedPointTrack,
            rsquareLimit=rsquareLimit
        ))
