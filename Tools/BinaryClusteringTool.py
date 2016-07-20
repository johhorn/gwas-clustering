from numpy import nonzero
from quick.multitrack.MultiTrackCommon import getGSuiteFromGalaxyTN
from gold.application.HBAPI import Track, doAnalysis, AnalysisSpec, GlobalBinSource
from quick.statistic.DistanceMetricsFoundationStat import DistanceMetricsFoundationStat
from quick.statistic.PointCountPerMicroBinV2Stat import PointCountPerMicroBinV2Stat
from quick.webtools.clustering.CommonClusteringFunctions import CommonClusteringFunctions
from quick.webtools.mixin.DebugMixin import DebugMixin


class BinaryClusteringTool(DebugMixin, CommonClusteringFunctions):
    # Type of comparison
    COMP_DIRECT = 'Direct base pair occurrences'
    COMP_BINS = 'Occurrences within user defined bins'

    @staticmethod
    def isPublic():
        return True

    @staticmethod
    def isDebugMode():
        return True

    @staticmethod
    def getToolName():
        return "Clustering using binary vector definition"

    @classmethod
    def getInputBoxNames(cls):
        return \
            CommonClusteringFunctions.getCommonClusteringInputBoxNames() + [
                ('Select vector representation', 'vectorDef'),
                ('Select microbin size', 'microBin')
            ]

    @staticmethod
    def getInputBoxOrder():
        return [
            'gSuite',
            'vectorDef',
            'microBin',
            'distanceCorrMeasure',
            'linkageCriterion',
            'debugMode'
        ]

    @staticmethod
    def getOptionsBoxVectorDef(choices):
        return [
            CommonClusteringFunctions.DEFAULT_SELECT,
            BinaryClusteringTool.COMP_DIRECT,
            BinaryClusteringTool.COMP_BINS
        ]

    @staticmethod
    def getOptionsBoxMicroBin(choices):
        if choices.vectorDef == BinaryClusteringTool.COMP_BINS:
            return '500000'
        else:
            return None

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
        distDict = None
        labels = None

        # Print tool information:
        cls.htmlClusterTitle(choices.vectorDef, htmlCore)
        cls.htmlClusterSubtext(choices.distanceCorrMeasure, cls.CORR_DISTLIST, choices.linkageCriterion, htmlCore)

        # Get distance / correlation matrixes
        if choices.vectorDef == BinaryClusteringTool.COMP_DIRECT:
            distDict, labels = cls.directVectorDistance(gSuite, analysisBins)
        elif choices.vectorDef == BinaryClusteringTool.COMP_BINS:
            distDict, labels = cls.microBinDistance(gSuite, analysisBins, choices)

        # Print plots
        if distDict and labels:
            cls.printDistPlots(distDict, labels, choices.distanceCorrMeasure, choices.linkageCriterion, galaxyFn,
                               htmlCore)

        cls.htmlClusterTime(str(time.clock() - start), htmlCore)
        htmlCore.divEnd()
        print htmlCore

    @staticmethod
    def getOutputFormat(choices=None):
        return 'html'

    @classmethod
    def directVectorDistance(cls, gSuite, analysisBins):
        """
        Each base pair represents its own feature.
        """
        analysisSpec = AnalysisSpec(DistanceMetricsFoundationStat)
        labels = []
        size = gSuite.numTracks()
        distDict = cls.createDistDict(cls.CORR_DISTLIST)

        for i in range(0, size):
            gSuiteTrack1 = gSuite.getTrackFromIndex(i)
            labels.append(gSuiteTrack1.title)

            for j in range(i + 1, size):
                gSuiteTrack2 = gSuite.getTrackFromIndex(j)
                track1 = Track(gSuiteTrack1.trackName)
                track2 = Track(gSuiteTrack2.trackName)

                count = doAnalysis(analysisSpec, analysisBins, [track1, track2]).getGlobalResult()
                cls.updateDistDict(distDict, count)

        return distDict, labels

    @classmethod
    def microBinDistance(cls, gSuite, analysisBins, choices):
        """
        Each bin represents a feature.
        """
        bins = []
        labels = []
        size = gSuite.numTracks()
        distDict = cls.createDistDict(cls.CORR_DISTLIST)

        analysisSpec = AnalysisSpec(PointCountPerMicroBinV2Stat)
        analysisSpec.addParameter('microBin', int(choices.microBin))

        # Get bins:
        for gSuiteTrack in gSuite.allTracks():
            labels.append(gSuiteTrack.title)
            track = [Track(gSuiteTrack.trackName)]
            res = doAnalysis(analysisSpec, analysisBins, track).getGlobalResult()
            if 'Result' in res:
                result = res['Result']
                bins.append(result)

        for i in range(0, size):
            for j in range(i + 1, size):
                bin1 = bins[i]
                bin2 = bins[j]
                count = {'a': 0, 'b': 0, 'c': 0, 'd': 0}

                # Create numpy masks
                nonzero_intersect = nonzero((bin1 != 0) & (bin2 != 0))
                nonzero_snps1 = nonzero(bin1)
                nonzero_snps2 = nonzero(bin2)

                # Use masks to find a, b, c, d
                count['a'] = len(bin1[nonzero_intersect])
                count['b'] = len(bin1[nonzero_snps1]) - count['a']
                count['c'] = len(bin2[nonzero_snps2]) - count['a']
                count['d'] = len(bin1) - count['a'] - count['b'] - count['c']
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

        if choices.vectorDef == CommonClusteringFunctions.DEFAULT_SELECT:
            return 'Please select a vector representation'

        if choices.vectorDef == BinaryClusteringTool.COMP_BINS:
            try:
                value = float(choices.microBin)
                if not value.is_integer():
                    return 'Please define bin size as a positive integer'
            except:
                return 'Please define bins size as a positive integer'

        errorString = CommonClusteringFunctions.checkClusterOptions(
            choices.distanceCorrMeasure,
            choices.linkageCriterion
        )
        if errorString:
            return errorString

        return None

    @staticmethod
    def getToolDescription():
        from gold.result.HtmlCore import HtmlCore

        core = HtmlCore()
        core.header('Binary vector clustering')
        core.paragraph('This tool compares all pairs of tracks in a GSuite, by using a standard binary vector '
                       'definition to represent each track. The features of the vectors is either the direct base '
                       'pairs of the entire genome, or a standard bin spanning a smaller region of each chromosome. '
                       'Each feature is a 0 or 1, denoting absence or presence of a track element at that feature. '
                       'For each pair of tracks, a binary matching representation is computed, as defined below. '
                       )
        core.divider()
        core.smallHeader('Matching definitions')
        core.paragraph('We compute the following matching values for each pair of tracks. Each pair of features '
                       'in the track vectors are either 0 (absence) or 1 (presence). The total number of features are '
                       'defined as n = a + b + c + d.'
                       'The match counts are further used to compute measures of distance or correlation.'
                       '<ul>'
                       '<li>a: Positive matches: (1, 1) </li>'
                       '<li>b: Mismatch: (1, 0) in first track</li>'
                       '<li>c: Mismatch: (0, 1) in second track</li>'
                       '<li>d: Negative matches: (0, 0) </li>'
                       '</ul>')
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
                       '<li>Forbes: (n * a) / ((a + c) * (a + b))</li>'
                       '<li>McConnaughey: ((a ** 2) - (b * c)) / (a + b) * (a + c) </li>'
                       '<li>Pearson: ((a * d) - (b * c)) / (sqrt((a + b) * (d + c) * (a + c) * (d + b)))'
                       '</ul>'
                       'The six first elements of the list are measures of similarity. The next, Forbes, measures '
                       'expected overlap (A value of 10 means it overlaps 10 times more than what is expected, given '
                       'independence. A value of 0.01 means it overlaps 10 times less). The last two measures are'
                       'correlation coefficients, and output a value between -1 and 1.')

        core.divider()
        core.smallHeader('Converting similarity and correlation into distance')
        core.paragraph('To cluster the tracks, given the different definitions of correlation and similarity, the '
                       'measures must be converted into standardized distance measures. <br>'
                       'The six first similarity measures is converted into distances by subtracting them from 1.'
                       'The distance value for Forbes is computed by 1 / max(1, log(forbes) + 1), as we are interested '
                       'in the cases where overlap happens more than expected.')
        core.paragraph('The last two measures, McConnaughey and Pearson, output coefficients between -1 and 1. '
                       'A value of -1 indicate the two tracks are highly independent with respect to location, '
                       'while 1 means high correlation and co-occurrence of SNPs. They are normalized to a distance '
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
