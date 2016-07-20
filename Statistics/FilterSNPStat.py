from gold.statistic.MagicStatFactory import MagicStatFactory
from gold.statistic.Statistic import Statistic
from gold.statistic.RawDataStat import RawDataStat
from gold.track.TrackFormat import TrackFormatReq
# Author: Johanne


class FilterSNPStat(MagicStatFactory):
    """
    Filter SNPs that map to the same genetic locus.
    Work for valued point tracks.

    A genetic locus is defined by the filterThreshold parameter passed to the Unsplittable class.
    This class will traverse the track and only keep the most significant point if multiple points
    lies within the same filterThreshold range. The most significant point is the one with the lowest value.

    Used mainly in createChildren of other statistics.
    """


class FilterSNPStatUnsplittable(Statistic):

    """
    Filters each track bin (chromosome)
    A SNP is defined as a valued point in the track
    1) For each SNP, check if there are other SNPs within the genetic locus threshold passed to class.
    2) If there are several within genetic locus threshold, only keep most significant, i.e SNP with lowest p-value
    """

    def _init(self, filterThreshold=0):
        self._filterThreshold = int(filterThreshold)

    def _keepSNP(self, snps, pval, notNanIndexes, indexList):

        filterIndexes = notNanIndexes[indexList]
        notNanArray = pval[indexList][filterIndexes]
        minPVal = min(notNanArray)
        index = list(pval[indexList]).index(minPVal)
        return snps[indexList][index]

    def filterSNP(self, snps, pval):

        import numpy as np
        filterThreshold = self._filterThreshold

        tmp = []
        filtered = []
        SNPcount = len(snps)

        maxWindowFlagSet = False
        currSNP = snps[0]
        nextSNP = snps[1]
        notNanIndexes = ~np.isnan(pval)

        i = 0
        while i < SNPcount - 1:

            currSNP = snps[i]
            nextSNP = snps[i+1]

            if nextSNP - currSNP <= filterThreshold and not maxWindowFlagSet:
                if len(tmp) > 0 and nextSNP - snps[tmp[0]] > filterThreshold:
                    maxWindowFlagSet = True
                    i -= 1
                    continue
                else:
                    tmp.append(i)

            elif len(tmp) > 0:
                if nextSNP - currSNP > filterThreshold:
                    tmp.append(i)

                filtered.append(self._keepSNP(snps, pval, notNanIndexes, np.asarray(tmp)))
                tmp = []
                maxWindowFlagSet = False
            else:
                filtered.append(currSNP)
            i += 1

        if nextSNP - currSNP <= filterThreshold:
            tmp.append(SNPcount - 1)
            filtered.append(self._keepSNP(snps, pval, notNanIndexes, tmp))
        else:
            filtered.append(nextSNP)

        return filtered

    def _compute(self):
        """
        Does computation on each bin/chromosome. Only need to define the computation for one chromosome, as the
        underlying structure of the HyperBrowser will apply this to each bin/chromosome

        Returns a filtered list of SNPs for the given chromosome region
        """
        rawData = self._children[0].getResult()
        snps = rawData.startsAsNumpyArray()
        pval = rawData.valsAsNumpyArray()

        if len(snps) > 1:
            return self.filterSNP(snps, pval)
        else:
            return list(snps)

    def _createChildren(self):
        self._addChild(RawDataStat(self._region, self._track, TrackFormatReq(allowOverlaps=True)))
