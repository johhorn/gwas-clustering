# Author: Johanne
from gold.statistic.MagicStatFactory import MagicStatFactory
from gold.statistic.Statistic import Statistic, StatisticSplittable
from gold.statistic.RawDataStat import RawDataStat
from gold.track.TrackFormat import TrackFormatReq


class UniquePointTrackStat(MagicStatFactory):
    """
    Returns list of unique rsIDs in track
    """
    pass


class UniquePointTrackStatSplittable(StatisticSplittable):

    def _combineResults(self):
        rsids = set()
        for track in self._childResults:
            if track.hasExtra('snps'):
                rsids.update(track.extrasAsNumpyArray('snps'))
        return list(rsids)


class UniquePointTrackStatUnsplittable(Statistic):

    def _compute(self):
        return self._children[0].getResult()

    def _createChildren(self):
        self._addChild(RawDataStat(self._region, self._track, TrackFormatReq(allowOverlaps=True)))
