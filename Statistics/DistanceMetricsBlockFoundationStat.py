# Author: Johanne
from gold.statistic.MagicStatFactory import MagicStatFactory
from quick.statistic.DistanceMetricsFoundationStat import DistanceMetricsFoundationStatUnsplittable, \
    DistanceMetricsFoundationStatSplittable


class DistanceMetricsBlockFoundationStat(MagicStatFactory):
    """
    Statistic which calculates positive matches (1, 1), and mismatches (0, 1) or (1, 0)
    between two tracks. Is made for valued point tracks, and return a genome-wide count.

    Each count represents a match for a genetic locus. This genetic locus is defined as a block, where SNPs within the
    same block count as a match.
    A filter threshold is set, which filter out SNPs that map to the same locus (see DistnaceMetricsFoundationStat).
    An overlap threshold is also set, typically threshold / 2, which define the boundaries for a match between SNPs
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


class DistanceMetricsBlockFoundationStatSplittable(DistanceMetricsFoundationStatSplittable):
    pass


class DistanceMetricsBlockFoundationStatUnsplittable(DistanceMetricsFoundationStatUnsplittable):

    def _computeMatch(self, snp1, snp2, threshold):
        a = 1
        b = 0
        c = 0
        return a, b, c

    def _getDistances(self, snps1, snps2, region):

        threshold = self._getThreshold()

        a = 0
        b = 0
        c = 0
        snp1len = len(snps1)
        snp2len = len(snps2)

        i, j = (0, 0)
        while i < snp1len and j < snp2len:

            snp1 = snps1[i]
            snp2 = snps2[j]

            if abs(snp1 - snp2) <= threshold:
                a_val, b_val, c_val = self._computeMatch(snp1, snp2, threshold)
                a += a_val
                b += b_val
                c += c_val

                i += 1
                j += 1

            elif snp1 < snp2:
                b += 1
                i += 1
            elif snp1 > snp2:
                c += 1
                j += 1

        if i < snp1len and snp2len == 0:
            b += snp1len
        elif snp1len == 0 and j < snp2len:
            c += snp2len
        elif i < snp1len and snp2len != 0:
            snp = snps2[j - 1]
            while i < snp1len and abs(snp - snps1[i]) <= threshold:
                a_val, b_val, c_val = self._computeMatch(snps1[i], snp, threshold)
                a += a_val
                b += b_val
                c += c_val
                i += 1
            b += snp1len - i
        elif j < snp2len and snp1len != 0:
            snp = snps1[i - 1]
            while j < snp2len and abs(snp - snps2[j]) <= threshold:
                a_val, b_val, c_val = self._computeMatch(snps2[j], snp, threshold)
                a += a_val
                b += b_val
                c += c_val
                j += 1
            c += snp2len - j
        return a, b, c, -1

    def _getThreshold(self):
        """
        Set thresholds for overlap of SNPs at the same genetic loci.
        Use half of the value used in the filtering.
        """
        threshold = 0
        if 'filterThreshold' in self._kwArgs:
            threshold = int(self._kwArgs['filterThreshold']) / 2
        return threshold

