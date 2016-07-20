from quick.multitrack.MultiTrackCommon import getGSuiteFromGalaxyTN
from gold.application.HBAPI import Track, doAnalysis, AnalysisSpec, GlobalBinSource
from quick.statistic.DistanceMetricsBlockFoundationStat import DistanceMetricsBlockFoundationStat
from quick.statistic.DistanceMetricsFuzzyFoundationStat import DistanceMetricsFuzzyFoundationStat
from quick.webtools.clustering.CommonClusteringFunctions import CommonClusteringFunctions
from quick.webtools.mixin.DebugMixin import DebugMixin


class LociClusteringTool(DebugMixin, CommonClusteringFunctions):

    # Type of comparison
    COMP_GAUSS = 'Gaussian definition of genetic loci'
    COMP_BLOCK = 'Block definition of genetic loci'

    # List of available measures
    CLUSTER_LIST = CommonClusteringFunctions.DISTLIST

    @staticmethod
    def isPublic():
        return True

    @staticmethod
    def isDebugMode():
        return True

    @staticmethod
    def getToolName():
        return "Clustering using definition of genetic loci"

    @classmethod
    def getInputBoxNames(cls):
        return \
            CommonClusteringFunctions.getCommonClusteringInputBoxNames() + [
                ('Select definition for genetic loci', 'similarityCase'),
                ('Select range of genetic loci', 'geneticLocus')
            ]

    @staticmethod
    def getInputBoxOrder():
        return [
            'gSuite',
            'similarityCase',
            'geneticLocus',
            'distanceMeasure',
            'linkageCriterion',
            'debugMode'
        ]

    @staticmethod
    def getOptionsBoxSimilarityCase(choices):
        return [
            CommonClusteringFunctions.DEFAULT_SELECT,
            LociClusteringTool.COMP_BLOCK,
            LociClusteringTool.COMP_GAUSS
        ]

    @staticmethod
    def getOptionsBoxGeneticLocus(choices):
        return '500000'

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

        # Print tool info
        cls.htmlClusterTitle(choices.similarityCase, htmlCore)
        cls.htmlClusterSubtext(choices.distanceMeasure, cls.CLUSTER_LIST, choices.linkageCriterion, htmlCore)
        htmlCore.line('Range of genetic loci: ' + choices.geneticLocus)

        # Analysis environment
        gSuite = getGSuiteFromGalaxyTN(choices.gSuite)
        analysisBins = GlobalBinSource(gSuite.genome)
        distDict = None
        labels = None

        if choices.similarityCase == LociClusteringTool.COMP_GAUSS:
            analysisSpec = AnalysisSpec(DistanceMetricsFuzzyFoundationStat)
            analysisSpec.addParameter('filterThreshold', int(choices.geneticLocus))

            distDict, labels = cls.computeDistance(gSuite, analysisSpec, analysisBins)

        elif choices.similarityCase == LociClusteringTool.COMP_BLOCK:
            analysisSpec = AnalysisSpec(DistanceMetricsBlockFoundationStat)
            analysisSpec.addParameter('filterThreshold', int(choices.geneticLocus))

            distDict, labels = cls.computeDistance(gSuite, analysisSpec, analysisBins)

        # Cluster and print plots
        if distDict and labels:
            cls.printDistPlots(distDict, labels, choices.distanceMeasure, choices.linkageCriterion, galaxyFn, htmlCore)

        cls.htmlClusterTime(str(time.clock() - start), htmlCore)
        htmlCore.divEnd()
        print htmlCore

    @staticmethod
    def getOutputFormat(choices=None):
        return 'html'

    @classmethod
    def computeDistance(cls, gSuite, differentTracksAnalysis, analysisBins):

        distDict = cls.createDistDict(cls.DISTLIST)
        labels = []

        trackCount = gSuite.numTracks()
        for i in range(0, trackCount):
            gSuiteTrack = gSuite.getTrackFromIndex(i)
            track = Track(gSuiteTrack.trackName)
            labels.append(gSuiteTrack.title)

            for j in range(i + 1, trackCount):
                track2 = Track(gSuite.getTrackFromIndex(j).trackName)
                tracks = [track, track2]
                count = doAnalysis(differentTracksAnalysis, analysisBins, tracks).getGlobalResult()

                cls.updateDistDict(distDict, count)

        return distDict, labels

    @staticmethod
    def validateAndReturnErrors(choices):
        """
        Should validate the selected input parameters. If the parameters are
        not valid, an error text explaining the problem should be returned. The
        GUI then shows this text to the user (if not empty) and greys out the
        execute button (even if the text is empty). If all parameters are
        valid, the method should return None, which enables the execute button.
        """
        errorString = CommonClusteringFunctions.validateGSuite(choices)
        if errorString:
            return errorString

        if choices.similarityCase == CommonClusteringFunctions.DEFAULT_SELECT:
            return 'Please select a genetic loci definition'

        try:
            value = float(choices.geneticLocus)
            if not value.is_integer():
                return 'Please define size of genetic loci as a positive integer'
        except:
            return 'Please define size of genetic loci as a positive integer'

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
        core.smallHeader('Clustering using predefined genetic loci')
        core.paragraph('This clustering tool uses a notion of genetic loci to represent binary features.'
                       'A genetic loci is defined as an '
                       'area of a particular size surrounding a SNP, for instance a block of 500kb. A preliminary '
                       'filtering step done on all tracks, only keeps the most significant SNP within a genetic loci. '
                       'Then, SNPs on different tracks that map to the same genetic loci, is counted. '
                       'The block version counts 1 for each pair of SNPs at the same genetic loci. '
                       'A gaussian version is also available, where positive matches are updated with a '
                       'gaussian score. The further '
                       'away two SNPs are from each other, the smaller the overlap score is. The gaussian overlap '
                       'definition does not continue indefinitely, but stops counting overlap at the border of the '
                       'genetic loci')
        core.divider()
        core.smallHeader('Matching definitions')
        core.paragraph('We compute the following matching values for each pair of tracks:'
                       '<ul>'
                       '<li>a: SNPs map to same genetic loci</li>'
                       '<li>b: SNP in first track has no overlap</li>'
                       '<li>c: SNP in second track has no overlap</li>'
                       '<li>d: Undefined</li>'
                       '</ul>'
                       'With the gaussian definition, overlap is updated for <i>a</i> with a gaussian score, based on '
                       'the distance between the SNPs that map to the same loci, while b and c is updated with '
                       '(1 - gaussian) / 2 each. If there is no overlap, <i>b</i> and <i>c</i> '
                       'is updated similar to the block definition. <br>'
                       'Negative matches, <i>d</i>, and consequently, <i>n</i>, is not defined with this feature '
                       'representation.')
        core.divider()
        core.smallHeader('Similarity/correlation definitions')
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
                       'The six first elements of the list are measures of similarity. The last measure is a'
                       'correlation coefficient, and output a value between -1 and 1.')

        core.divider()
        core.smallHeader('Converting similarity and correlation into distance')
        core.paragraph('To cluster the tracks, given the different definitions of correlation and similarity, the '
                       'measures must be converted into standardized distance measures. <br>'
                       'The six first similarity measures is converted into distances by subtracting them from 1.'
                       'McConnaughey output a correlation coefficient between -1 and 1. '
                       'A value of -1 indicate the two tracks are highly independent with respect to location, '
                       'while 1 means high correlation and co-occurrence of SNPs. It is normalized to a distance '
                       'measure by (corr + 1) / 2. High positive correlation gives small distance, while high '
                       'negative correlation gives large distance.')
        core.divider()
        core.smallHeader('Output')
        core.paragraph('For the measure chosen for the clustering, the output is a heatmap of all pairwise distances '
                       'between tracks in the given GSuite. For McConnaughey and Pearson, '
                       'the heatmap output the correlation coefficients. In addition, a clustering dendrogram show '
                       'the relations between the tracks, computed from the standardized distance of the measure. '
                       'Raw text values for correlation or distance matrixes, as well as the linkage matrix, are '
                       'included. Lastly, a rank matrix is printed, to show the order of the standardized pairwise '
                       'distances. The smallest rank show which tracks have the smallest distance between them.')
        return str(core)
