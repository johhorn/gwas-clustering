# Author: Johanne
from gold.statistic.MagicStatFactory import MagicStatFactory
from gold.statistic.RawDataStat import RawDataStat
from gold.statistic.Statistic import Statistic, StatisticSplittable
from gold.track.TrackFormat import TrackFormatReq
from gold.track.TrackView import TrackView
from quick.statistic.FilterSNPStat import FilterSNPStat


class DistanceMetricsFoundationStat(MagicStatFactory):
    """
    Statistic which calculates positive matches (1, 1), negative matches (0, 0) and mismatches (0, 1) or (1, 0)
    between two tracks. Is made for valued point tracks, and return a genome-wide count.

    Each count represents a base pair in the tracks. Negative matches is calculated as chromosome size, minus
    a, b and c.

    Can be used for distance computations as described in the article "A survey of Binary Similarity and Distance
    Measures" by Choi, Cha and Tappert.

    Returns a dictionary containing scores for all these cases. The dictionary have the following keys:
    a: Positive matches
    b: Mismatches (1, 0)
    c: Mismatches (0, 1)
    d: Negative matches

    See example of usage in quick/webtools/clustering/BinaryClusteringTool.py

    Two subclasses exist, which build on the implementation below:
    - DistanceMetricsBlockFoundationStat
    - DistanceMetricsFuzzyFoundationStat
    """
    pass


class DistanceMetricsFoundationStatSplittable(StatisticSplittable):

    def _combineResults(self):
        """Combines the matches found for each chromosome in _compute of the Unsplittable class."""

        count = {'a': 0, 'b': 0, 'c': 0, 'd': 0}
        for child in self._childResults:
            if child:
                count['a'] += child[0]
                count['b'] += child[1]
                count['c'] += child[2]
                count['d'] += child[3]

        return count


class DistanceMetricsFoundationStatUnsplittable(Statistic):
    @classmethod
    def _getDistances(cls, snps1, snps2, region):
        """
        Takes in two lists of SNPs, and compute matches and mismatches for the lists, where each SNP
        represents its own binary feature.
        """
        import numpy as np

        intersectSNPs = np.intersect1d(snps1, snps2)
        a = intersectSNPs.size
        b = snps1.size - intersectSNPs.size
        c = snps2.size - intersectSNPs.size
        d = region.getTotalBpSpan() - a - b - c
        return a, b, c, d

    def _compute(self):
        """
        For each bin (chromosome), get child data for each track, and count matches for each list.
        Return these matches as a, b, c and d, defined above.
        Pass these lists to the Splittable class.
        """
        snps1 = self._children[0].getResult()
        snps2 = self._children[1].getResult()
        chromosome = self._region

        if isinstance(snps1, TrackView) and isinstance(snps2, TrackView):
            snps1 = snps1.startsAsNumpyArray()
            snps2 = snps2.startsAsNumpyArray()

        return self._getDistances(snps1, snps2, chromosome)

    def _createChildren(self):
        """
        Create children corresponding to each track. If a filter threshold is set (as in the definitions of
        genetic loci), each child is a filtered list of SNPs.
        If no threshold is set, the children are TrackViews with all SNPs reported for each track.
        """

        if 'filterThreshold' in self._kwArgs:
            self._addChild(FilterSNPStat(self._region, self._track, filterThreshold=self._kwArgs['filterThreshold']))
            self._addChild(FilterSNPStat(self._region, self._track2, filterThreshold=self._kwArgs['filterThreshold']))

        else:
            self._addChild(RawDataStat(self._region, self._track, TrackFormatReq(allowOverlaps=False)))
            self._addChild(RawDataStat(self._region, self._track2, TrackFormatReq(allowOverlaps=False)))
