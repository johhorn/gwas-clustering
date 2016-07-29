# Author: Johanne
from gold.statistic.MagicStatFactory import MagicStatFactory
from gold.statistic.Statistic import Statistic, StatisticSplittable
from gold.statistic.RawDataStat import RawDataStat
from gold.track.TrackFormat import TrackFormatReq
from math import isinf


class OverlappingValsListStat(MagicStatFactory):
    """
    Return two lists containing values, one for each track.
    The features of the vectors represent positions where both tracks have a SNP.
    """
    pass


class OverlappingValsListStatSplittable(StatisticSplittable):
    """
    Combines the vectors created for each bin (chromosome) into a genome-wide vector of values.
    """

    def _combineResults(self):

        pvalsX = []
        pvalsY = []

        for child in self._childResults:
            pvalsX += child[0]
            pvalsY += child[1]

        return {'X': pvalsX, 'Y': pvalsY}


class OverlappingValsListStatUnsplittable(Statistic):
    """
    For each region, finds all positions present in both tracks, and store the lowest values of the track position,
    one for each track. Returns two vectors with values for all such overlapping positions in the corresponding
    track.
    """

    def _compute(self):

        track1 = self._children[0].getResult()
        track2 = self._children[1].getResult()

        positions1 = track1.startsAsNumpyArray()
        positions2 = track2.startsAsNumpyArray()
        track1Length = track1.getNumElements()
        track2Length = track2.getNumElements()

        if track1Length == 0 or track2Length == 0:
            return [], []

        pvals1 = track1.valsAsNumpyArray()
        pvals2 = track2.valsAsNumpyArray()

        pvalsX = {}
        pvalsY = {}

        i, j = (0, 0)

        while i < track1Length and j < track2Length:

            pos1 = positions1[i]
            pos2 = positions2[j]
            pval1 = pvals1[i]
            pval2 = pvals2[j]

            if self._isSamePosition(pos1, pos2) and self._notInfinity(pval1, pval2):
                self._appendValue(pos1, pval1, pvalsX)
                self._appendValue(pos2, pval2, pvalsY)
                i += 1
                j += 1
            elif pos1 < pos2:
                self._updateIfPresent(pos1, pval1, pvalsX)
                i += 1
            elif pos2 < pos1:
                self._updateIfPresent(pos2, pval2, pvalsY)
                j += 1

        # If any values remain in the longer list that have the same position as the last match
        if i == track1Length and j < track2Length:
            pos1 = positions1[i - 1]
            while j < track2Length:
                pos2 = positions2[j]
                if pos2 != pos1:
                    break

                pval2 = pvals2[j]
                self._updateIfPresent(pos2, pval2, pvalsY)
                j += 1
        elif j == track2Length and i < track1Length:
            pos2 = positions2[j - 1]
            while i < track1Length:
                pos1 = positions1[i]
                if pos2 != pos1:
                    break

                pval1 = pvals1[i]
                self._updateIfPresent(pos1, pval1, pvalsX)
                i += 1

        return self._dictToVector(pvalsX), self._dictToVector(pvalsY)

    def _notInfinity(self, pval1, pval2):
        return not isinf(pval1) or not isinf(pval2)

    def _isSamePosition(self, pos1, pos2):
        return pos1 == pos2

    def _appendValue(self, pos, value, valueDict):
        """
        I the children are created with allowOverlaps=False, SNPs that have several P-values (reported in different
        studies) will get a value of nan. This is common for tracks generated from the GWAS Catalog via the
        HyperBrowser.

        To in stead get the most significant value of these positions, while at the same time not creating multiple
        features for the same position, a dictionary store all vector features for each track.
        """
        if pos in valueDict:
            self._updateValue(pos, value, valueDict)
        else:
            valueDict[pos] = value

    def _updateValue(self, pos, val, valueDict):
        """
        Called to update value, when position is key in dictionary.
        Values of -log(pval) are assumed, and the most significant SNP of the two is defined as the one with
        the highest value.
        """
        value = valueDict[pos]

        if val > value:
            valueDict[pos] = val

    def _updateIfPresent(self, pos, val, valueDict):
        """
        Check to see if the value maps to a feature that is stored for the track vector.
        If it is, and the new value has higher significance (lower value), the new value replace the old.
        """
        if isinf(val):
            return

        if pos in valueDict:
            self._updateValue(pos, val, valueDict)

    def _dictToVector(self, valueDict):
        """
        Return vector of unique positions, in sorted order
        """
        sortedPositions = sorted(valueDict.keys())
        sortedValues = [valueDict[key] for key in sortedPositions]
        return sortedValues

    def _createChildren(self):
        self._addChild(RawDataStat(self._region, self._track, TrackFormatReq(allowOverlaps=True)))
        self._addChild(RawDataStat(self._region, self._track2, TrackFormatReq(allowOverlaps=True)))



