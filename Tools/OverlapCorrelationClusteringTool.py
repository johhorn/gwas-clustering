from gold.application.HBAPI import doAnalysis
from gold.description.AnalysisDefHandler import AnalysisSpec
from gold.track.Track import Track
from quick.application.UserBinSource import GlobalBinSource
from quick.multitrack.MultiTrackCommon import getGSuiteFromGalaxyTN
from quick.statistic.GeneticLociOverlapStat import GeneticLociOverlapStat
from quick.webtools.clustering.CommonCorrelationFunctions import CommonCorrelationFunctions
from quick.webtools.mixin.DebugMixin import DebugMixin


class OverlapCorrelationClusteringTool(DebugMixin, CommonCorrelationFunctions):
    @staticmethod
    def getToolName():
        return "Find correlation between vectors of genetic overlap to other tracks"

    @staticmethod
    def isPublic():
        return True

    @classmethod
    def getInputBoxNames(cls):
        return cls.getCommonCorrelationInputBoxes() + [
            ('Select size of genetic locus', 'geneticLocus')
        ] + cls.getInputBoxNamesForDebug()

    @staticmethod
    def getOptionsBoxGeneticLocus(choices):
        return '500000'

    @staticmethod
    def getInputBoxOrder():
        return [
            'gSuite',
            'geneticLocus',
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
        analysisSpec = AnalysisSpec(GeneticLociOverlapStat)
        analysisSpec.addParameter('filterThreshold', int(choices.geneticLocus))

        # Print tool information:
        cls.htmlClusterTitle(cls.getToolName(), htmlCore)
        cls.htmlClusterSubtext(choices.corrStat, [cls.CORR_PEARSON, cls.CORR_SPEARMAN],
                               choices.linkageCriterion, htmlCore)
        cls.htmlVectorHandling(htmlCore)

        # Get correlations
        overlapMatrix, labels = cls.getOverlapMatrix(analysisBins, analysisSpec, gSuite)
        corrDict = cls.getTriangularCorrMatrix(overlapMatrix)
        cls.printCorrPlots(corrDict, labels, choices.corrStat, choices.linkageCriterion, galaxyFn, htmlCore)

        cls.htmlClusterTime(str(time.clock() - start), htmlCore)
        htmlCore.divEnd()
        print htmlCore

    @classmethod
    def getOverlapMatrix(cls, analysisBins, analysisSpec, gSuite):
        overlapMatrix = []
        labels = []

        i = 0
        for gSuiteTrack in gSuite.allTracks():

            labels.append(gSuiteTrack.title)
            overlapMatrix.append([])
            track = Track(gSuiteTrack.trackName)

            for gSuiteTrack2 in gSuite.allTracks():

                if gSuiteTrack == gSuiteTrack2:
                    overlapMatrix[i].append('.')  # This value must be removed before correlating the values.
                    continue

                track2 = Track(gSuiteTrack2.trackName)
                result = doAnalysis(analysisSpec, analysisBins, [track, track2]).getGlobalResult()
                overlapMatrix[i].append(result['Result'])

            i += 1

        return overlapMatrix, labels

    @classmethod
    def getTriangularCorrMatrix(cls, overlapMatrix):
        corrDict = cls.createDistDict([cls.CORR_PEARSON, cls.CORR_SPEARMAN])

        size = len(overlapMatrix)
        for i in range(0, size):
            for j in range(i + 1, size):
                vector1 = cls.modifyVector(overlapMatrix[i], i, j)
                vector2 = cls.modifyVector(overlapMatrix[j], j, i)
                cls.updateCorrDict(corrDict, vector1, vector2)

        return corrDict

    @classmethod
    def modifyVector(cls, track, trackIndex, otherIndex):

        lower = trackIndex
        higher = otherIndex

        if otherIndex < trackIndex:
            lower = otherIndex
            higher = trackIndex

        return track[:lower] + track[lower + 1:higher] + track[higher + 1:]

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

        try:
            value = float(choices.geneticLocus)
            if not value.is_integer():
                return 'Please define size of genetic loci as a positive integer'
        except:
            return 'Please define size of genetic loci as a positive integer'

        if choices.corrStat == CommonCorrelationFunctions.DEFAULT_SELECT:
            return 'Please select a correlation coefficient'
        if choices.linkageCriterion == CommonCorrelationFunctions.DEFAULT_SELECT:
            return 'Please select a linkage criterion'

    @staticmethod
    def getToolDescription():
        from gold.result.HtmlCore import HtmlCore

        core = HtmlCore()
        core.header("Pearson/Spearman correlation of genetic overlap to other tracks")
        core.paragraph(""
                       "This tool finds correlation for all pairs of valued point tracks within the given GSuite. "
                       "Each track is represented as a vector of overlap with other tracks. "
                       "A correlation coefficient is then computed for the vector, using the Pearson or Spearman "
                       "r statistics. The result will be a value between -1 and 1, "
                       "where 0 denotes no correlation, -1 negative correlation, and 1 positive correlation. ")
        core.paragraph("Note these methods require a large number of samples to work as expected. "
                       "A vector size of at least 100, i.e. at least 100 overlapping SNPs between two tracks to be "
                       "compared, is recommended.")
        core.divider()
        core.smallHeader("Vector representation")
        core.paragraph(
            "The vectors of the different tracks are constructed using overlap to other tracks. "
            "Overlap is determined using a genetic locus definition."
            'A genetic loci is defined as an '
            'area of a particular size surrounding a SNP, for instance a block of 500kb. A preliminary '
            'filtering step done on all tracks, only keeps the most significant SNP within a genetic loci. '
            'Overlap is then computed for all pairs of tracks, and if two SNPs map to the same loci, their overlap '
            'is counted. '
            "To avoid bias of tracks overlapping with themselves, these features are removed before every pairwise "
            "comparison."
        )
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
        return True

    @classmethod
    def htmlVectorHandling(cls, htmlCore):
        htmlCore.line("For the vectors, the features where the tracks overlap with themselves, are removed.")
