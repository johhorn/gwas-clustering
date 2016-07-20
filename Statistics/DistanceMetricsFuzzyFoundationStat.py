# Author: Johanne
from gold.statistic.MagicStatFactory import MagicStatFactory
from quick.statistic.DistanceMetricsBlockFoundationStat import DistanceMetricsBlockFoundationStatUnsplittable, \
    DistanceMetricsBlockFoundationStatSplittable


class DistanceMetricsFuzzyFoundationStat(MagicStatFactory):
    """
    Statistic which calculates positive matches (1, 1) and mismatches (0, 1) or (1, 0)
    between two tracks. Is made for valued point tracks, and return a genome-wide count.

    Each count represents a match for a genetic locus. This genetic locus is defined as a gaussian curve, where SNPs
    within the same block, count as a match, but rather than adding a full 1, the gaussian value is used.
    Positive matches with this gaussian definition is more fuzzy, as matches not at the exact same position in both
    track, will add a score to a, and divide 1 - a between b and c.

    A filter threshold is set, which filter out SNPs that map to the same locus (see DistnaceMetricsFoundationStat).
    An overlap threshold is also set, typically threshold / 4, which define the boundaries for a match between SNPs
    of different tracks.

    Negative matches are not well defined, and -1 is returned in stead.

    Can be used for distance computations as described in the article "A survey of Binary Similarity and Distance
    Measures" by Choi, Cha and Tappert.

    Returns a dictionary containing scores for all these cases. The dictionary have the following keys:
    a: Positive matches
    b: Mismatches (1, 0)
    c: Mismatches (0, 1)
    d: -1

    See example of usage in quick/webtools/clustering/LociClusteringTool.py
    """
    pass


class DistanceMetricsFuzzyFoundationStatSplittable(DistanceMetricsBlockFoundationStatSplittable):
    pass


class DistanceMetricsFuzzyFoundationStatUnsplittable(DistanceMetricsBlockFoundationStatUnsplittable):

    def _computeMatch(self, snp1, snp2, threshold):
        """
        The threshold is half the genetic loci used in filtering. It behaves much like the one in the block definition,
        but for computing a fuzzy positive match, it is split in half again.
        The curve fits somewhat a model of LD decay, with respect to physical distance, when the sigma of the
        function (below: bell_width) is set to a quarter of the full genetic loci range.
        """

        from math import exp
        bell_width = threshold / 2.0
        snp1 = float(snp1)
        snp2 = float(snp2)
        numerator = (snp2 - snp1) ** 2
        denominator = 2 * bell_width ** 2
        a = exp(- (numerator / denominator))
        b = (1 - a) / 2
        c = (1 - a) / 2
        return a, b, c
