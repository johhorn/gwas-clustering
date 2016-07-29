from gold.application.HBAPI import doAnalysis
from gold.description.AnalysisDefHandler import AnalysisSpec
from gold.track.Track import Track
from quick.application.ExternalTrackManager import ExternalTrackManager
from quick.application.UserBinSource import GlobalBinSource
from quick.multitrack.MultiTrackCommon import getGSuiteFromGalaxyTN
from quick.statistic.ExpandTrackAndMatchStat import ExpandTrackAndMatchStat
from quick.webtools.GeneralGuiTool import GeneralGuiTool
from quick.webtools.clustering.CommonClusteringFunctions import CommonClusteringFunctions
from quick.webtools.mixin.DebugMixin import DebugMixin


class LDExpansionClusteringTool(CommonClusteringFunctions, DebugMixin):

    # List of available measures
    CLUSTER_LIST = CommonClusteringFunctions.DISTLIST

    @staticmethod
    def getToolName():
        return "Clustering using definition of haplotype blocks"

    @staticmethod
    def isPublic():
        return True

    @classmethod
    def getInputBoxNames(cls):
        return cls.getCommonClusteringInputBoxNames() + [
            ('Select LD graph track', 'ldTrack'),
            ('Set r<sup>2</sup> threshold', 'rSquare')
        ]

    @staticmethod
    def getInputBoxOrder():
        return [
            'gSuite',
            'ldTrack',
            'rSquare',
            'distanceMeasure',
            'linkageCriterion',
            'debugMode'
        ]

    @staticmethod
    def getOptionsBoxLdTrack(choices):
        return GeneralGuiTool.getHistorySelectionElement('gtrack')

    @staticmethod
    def getInfoForOptionsBoxLdTrack(choices):
        return 'Linked point track that describe LD relations between track elements of the GSuite. Can be ' \
               'generated with the "LD track generator" tool.'

    @staticmethod
    def getOptionsBoxRSquare(choices):
        return [CommonClusteringFunctions.DEFAULT_SELECT, '1.0', '0.9', '0.8', '0.7']

    @staticmethod
    def getInfoForOptionsBoxRSquare(choices):
        return 'Lower limit of r<sup>2</sup> for LD correlation between two variants. NB: Will only filter if the ' \
               'given LD graph track was created with a lower r<sup>2</sup> limit.'

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

        # Print tool information
        cls.htmlClusterTitle(cls.getToolName(), htmlCore)
        cls.htmlClusterSubtext(choices.distanceMeasure, cls.CLUSTER_LIST, choices.linkageCriterion, htmlCore)
        htmlCore.line('Threshold of r<sup>2</sup>: ' + choices.rSquare)

        # Analysis environment
        gSuite = getGSuiteFromGalaxyTN(choices.gSuite)
        analysisBins = GlobalBinSource(gSuite.genome)
        analysisSpec = AnalysisSpec(ExpandTrackAndMatchStat)

        splitName = choices.ldTrack.split(":")
        trackName = ExternalTrackManager.getPreProcessedTrackFromGalaxyTN(gSuite.genome, splitName)
        linkedPointTrack = Track(trackName)

        # Find distance/correlation matrix
        labels = []
        distDict = cls.createDistDict(cls.CLUSTER_LIST)
        size = gSuite.numTracks()
        for i in range(0, size):
            gSuiteTrack1 = gSuite.getTrackFromIndex(i)
            labels.append(gSuiteTrack1.title)
            for j in range(i + 1, size):
                gSuiteTrack2 = gSuite.getTrackFromIndex(j)
                track1 = Track(gSuiteTrack1.trackName)
                track2 = Track(gSuiteTrack2.trackName)
                count = doAnalysis(analysisSpec, analysisBins, [track1, track2, linkedPointTrack]).getGlobalResult()
                cls.updateDistDict(distDict, count)

        # Cluster and print plots
        cls.printDistPlots(distDict, labels, choices.distanceMeasure, choices.linkageCriterion, galaxyFn, htmlCore)

        cls.htmlClusterTime(str(time.clock() - start), htmlCore)
        htmlCore.divEnd()
        print htmlCore

    @staticmethod
    def isDebugMode():
        return False

    @staticmethod
    def getOutputFormat(choices=None):
        return 'html'

    @staticmethod
    def validateAndReturnErrors(choices):
        errorString = CommonClusteringFunctions.validateGSuite(choices)
        if errorString:
            return errorString

        if not choices.ldTrack:
            return 'Please select an LD graph track (linked point track)'

        if choices.rSquare == CommonClusteringFunctions.DEFAULT_SELECT:
            return 'Please select a threshold for r<sup>2</sup>'

        errorString = CommonClusteringFunctions.checkClusterOptions(
            choices.distanceMeasure,
            choices.linkageCriterion
        )
        if errorString:
            return errorString

        return None

    @staticmethod
    def getToolDescription():
        from gold.result.HtmlCore import HtmlCore

        core = HtmlCore()
        core.smallHeader('Clustering of haplotype blocks')
        core.paragraph('This tool makes use of actual LD information when it compares genetic loci.'
                       'The LD is given in the form of a linked point track. Base pairs in LD will be seen as '
                       'a correlated LD cluster within the same track. Overlap between these LD clusters on '
                       'different is then counted as a match. ')
        core.paragraph('An LD track for the GSuite can be created with the "LD track generator" tool.')
        core.divider()
        core.smallHeader('Matching definitions')
        core.paragraph('We compute the following matching values for each pair of tracks:'
                       '<ul>'
                       '<li>a: Positive matches: Overlap between LD clusters </li>'
                       '<li>b: Count of unmatched LD clusters in first track</li>'
                       '<li>c: Count of unmatched LD clusters in second track</li>'
                       '<li>d: Not defined</li>'
                       '</ul>')
        core.divider()
        core.smallHeader('Similarity definitions')
        core.paragraph(''
                       '<ul>'
                       '<li>Jaccard: a / (a + b + c) </li>'
                       '<li>Cosine: a / (sqrt(a + b) * sqrt(a + c))</li>'
                       '<li>Simpson: a / min(a + b, a + c) </li>'
                       '<li>Otsuka: a / sqrt((a + b) * (a + c)) </li>'
                       '<li>Sorgenfrei: a ** 2 / ((a + b) * (a + c)) </li>'
                       '<li>Kulczynski: ((a / 2) * (2 * a + b + c)) / ((a + b) * (a + c))</li>'
                       '<li>McConnaughey: ((a ** 2) - (b * c)) / (a + b) * (a + c) </li>'
                       '</ul>'
                       'All measures are defined from binary operational taxonomic units (OTUs).'
                       'They evaluate how similar two tracks are, and the different measures focus on different '
                       'properties.')

        core.divider()
        core.smallHeader('Converting similarity and correlation into distance')
        core.paragraph('To cluster the tracks, given the different definitions of correlation and similarity, the '
                       'measures must be converted into standardized distance measures. <br>'
                       'The six first similarity measures is converted into distances by subtracting them from 1.')
        core.paragraph('The last measure, McConnaughey, output a value between -1 and 1, and is normalized to a '
                       'distance '
                       'measure by (corr + 1) / 2. High positive correlation gives small distance, while high '
                       'negative correlation gives large distance.')
        core.divider()
        core.smallHeader('Output')
        core.paragraph('For the measure chosen for the clustering, the output is a heatmap of all pairwise '
                       'similarities '
                       'between tracks in the given GSuite. In addition, a clustering dendrogram show '
                       'the relations between the tracks, computed from the standardized distance of the measure. '
                       'Raw text values for similarity and distance matrixes, as well as the linkage matrix, are '
                       'included. Lastly, a rank matrix is printed, to show the order of the standardized pairwise '
                       'distances. The smallest rank show which tracks have the smallest distance between them.')
        return str(core)
