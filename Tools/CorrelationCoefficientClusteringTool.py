from gold.application.HBAPI import doAnalysis
from gold.description.AnalysisDefHandler import AnalysisSpec
from gold.track.Track import Track
from quick.application.UserBinSource import GlobalBinSource
from quick.multitrack.MultiTrackCommon import getGSuiteFromGalaxyTN
from quick.statistic.OverlappingValsListStat import OverlappingValsListStat
from quick.webtools.clustering.CommonCorrelationFunctions import CommonCorrelationFunctions
from quick.webtools.mixin.DebugMixin import DebugMixin


class CorrelationCoefficientClusteringTool(DebugMixin, CommonCorrelationFunctions):

    @staticmethod
    def getToolName():
        return "Find correlation between vectors of values from pairs of tracks"

    @staticmethod
    def isPublic():
        return True

    @classmethod
    def getInputBoxNames(cls):
        return cls.getCommonCorrelationInputBoxes() + cls.getInputBoxNamesForDebug()

    @staticmethod
    def getInputBoxOrder():
        return [
            'gSuite',
            'corrStat',
            'linkageCriterion'
        ]

    @classmethod
    def execute(cls, choices, galaxyFn=None, username=''):
        import time
        start = time.clock()

        # HTML settings
        from gold.result.HtmlCore import HtmlCore
        htmlCore = HtmlCore()
        htmlCore.divBegin(style=cls.HTML_STYLE)

        # Set debug environment
        cls._setDebugModeIfSelected(choices)

        # Analysis environment
        gSuite = getGSuiteFromGalaxyTN(choices.gSuite)
        analysisBins = GlobalBinSource(gSuite.genome)

        # Print info about tool
        cls.htmlClusterTitle(cls.getToolName(), htmlCore)
        cls.htmlClusterSubtext(choices.corrStat, [cls.CORR_PEARSON, cls.CORR_SPEARMAN],
                               choices.linkageCriterion, htmlCore)
        cls.htmlVectorExplanation(htmlCore)

        corrDict, labels = cls.trackOverlapValuesCorrelation(analysisBins, gSuite)
        if corrDict and labels:
            cls.printCorrPlots(corrDict, labels, choices.corrStat, choices.linkageCriterion, galaxyFn, htmlCore)

        cls.htmlClusterTime(str(time.clock() - start), htmlCore)
        htmlCore.divEnd()
        print htmlCore

    @classmethod
    def trackOverlapValuesCorrelation(cls, analysisBins, gSuite):
        """
        Represent each track as a vector with values at positions that are present in both tracks
        """
        corrDict = cls.createDistDict([cls.CORR_PEARSON, cls.CORR_SPEARMAN])
        analysisSpec = AnalysisSpec(OverlappingValsListStat)
        gSuiteSize = gSuite.numTracks()
        labels = []

        for i in range(0, gSuiteSize):

            gSuiteTrack = gSuite.getTrackFromIndex(i)
            labels.append(gSuiteTrack.title)

            for j in range(i + 1, gSuiteSize):

                gSuiteTrack2 = gSuite.getTrackFromIndex(j)
                track1 = Track(gSuiteTrack.trackName)
                track2 = Track(gSuiteTrack2.trackName)

                track1List = []
                track2List = []

                result = doAnalysis(analysisSpec, analysisBins, [track1, track2]).getGlobalResult()

                if 'X' in result and 'Y' in result:
                    track1List = result['X']
                    track2List = result['Y']

                cls.updateCorrDict(corrDict, track1List, track2List)

        return corrDict, labels

    @staticmethod
    def validateAndReturnErrors(choices):
        """
        Should validate the selected input parameters. If the parameters are
        not valid, an error text explaining the problem should be returned. The
        GUI then shows this text to the user (if not empty) and greys out the
        execute button (even if the text is empty). If all parameters are
        valid, the method should return None, which enables the execute button.
        """
        errorString = CommonCorrelationFunctions.validateGSuite(choices)
        if errorString:
            return errorString

        if choices.corrStat == CommonCorrelationFunctions.DEFAULT_SELECT:
            return 'Please select a correlation coefficient'

        if choices.linkageCriterion == CommonCorrelationFunctions.DEFAULT_SELECT:
            return 'Please select a linkage criterion'

    @staticmethod
    def getToolDescription():
        from gold.result.HtmlCore import HtmlCore

        core = HtmlCore()
        core.header("Pearson/Spearman correlation of values")
        core.paragraph(""
                       "This tool finds correlation for all pairs of valued point tracks within the given GSuite. "
                       "Each track is represented as a vector. The vector features are the values in the respective "
                       "tracks, for SNPs that are present in both tracks."
                       "A correlation coefficient is then computed for the vector, using the Pearson or Spearman "
                       "r statistics. The result will be a value between -1 and 1, "
                       "where 0 denotes no correlation, -1 negative correlation, and 1 positive correlation. "
                       )
        core.paragraph("Note these methods require a large number of samples to work as expected. "
                       "A vector size of at least 100, i.e. at least 100 overlapping SNPs between two tracks to be "
                       "compared, is recommended.")
        core.divider()
        core.smallHeader("Conversion to distance")
        core.paragraph('To cluster the tracks, the correlation coefficients must be converted into standardized '
                       'distance measures. <br>'
                       'In the tool, this conversion is done by computing (corr + 1) / 2. '
                       'High positive correlation gives small distance, while high '
                       'negative correlation gives large distance. '
                       )
        return str(core)

    @staticmethod
    def getOutputFormat(choices=None):
        return 'html'

    @staticmethod
    def isDebugMode():
        return False

    @classmethod
    def htmlVectorExplanation(cls, htmlCore):
        htmlCore.line("The vectors of this computation contains the point values of the different tracks. ")
        htmlCore.line("Each vector consists of the values of the corresponding track, for SNPs present in both tracks.")