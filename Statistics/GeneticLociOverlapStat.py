# Author: Johanne
from gold.statistic.MagicStatFactory import MagicStatFactory
from gold.statistic.Statistic import Statistic, StatisticSumResSplittable
from quick.statistic.FilterSNPStat import FilterSNPStat


class GeneticLociOverlapStat(MagicStatFactory):
    """
    Calculate sum of overlapping genetic loci between two tracks.

    A genetic locus is defined by a threshold set as parameter for the analysis. The tracks are filtered to
    only contain the most significant SNP if multiple SNPs map to the same locus.
    The overlap is then computed for half the threshold, so that all SNPs can only overlap with one locus at the other
    track.
    Count of overlaps, summed over all chromosomes, are returned.
    """
    pass


class GeneticLociOverlapStatSplittable(StatisticSumResSplittable):
    pass


class GeneticLociOverlapStatUnsplittable(Statistic):

    def _countOverlappingLoci(self, snps1, snps2):

        if len(snps1) == 0 or len(snps2) == 0:
            return 0

        threshold = self._getGeneticLociThreshold()
        count = 0
        i, j = (0, 0)
        while i < len(snps1) and j < len(snps2):

            snp1 = snps1[i]
            snp2 = snps2[j]

            if abs(snp1 - snp2) <= threshold:
                count += 1
                i += 1
                j += 1

            elif snp1 < snp2:
                i += 1
            elif snp1 > snp2:
                j += 1
        return count

    def _getGeneticLociThreshold(self):
        threshold = 0
        if 'filterThreshold' in self._kwArgs:
            threshold = int(self._kwArgs['filterThreshold']) / 2

        return threshold

    def _compute(self):
        snps1 = self._children[0].getResult()
        snps2 = self._children[1].getResult()
        return self._countOverlappingLoci(snps1, snps2)

    def _createChildren(self):

        self._addChild(FilterSNPStat(self._region, self._track, filterThreshold=self._getGeneticLociThreshold() * 2))
        self._addChild(FilterSNPStat(self._region, self._track2, filterThreshold=self._getGeneticLociThreshold() * 2))

