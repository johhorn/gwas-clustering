from quick.webtools.GeneralGuiTool import GeneralGuiTool
from quick.multitrack.MultiTrackCommon import getGSuiteFromGalaxyTN
from gold.application.HBAPI import GlobalBinSource, doAnalysis
from quick.webtools.clustering.BipartiteMatching import BipartiteMatching
from quick.webtools.clustering.CommonClusteringFunctions import CommonClusteringFunctions
from quick.webtools.clustering.LDExpansions import LDExpansions
from quick.webtools.mixin.DebugMixin import DebugMixin


class LDBipartiteMatchingTool(DebugMixin, CommonClusteringFunctions):

    # List of available measures
    CLUSTER_LIST = CommonClusteringFunctions.DISTLIST

    LD_GREEDY = 'Greedy bipartite graph match'
    LD_LAPJV = 'Jonker-Volgenant bipartite graph match'

    @staticmethod
    def isPublic():
        return True

    @staticmethod
    def getToolName():
        return "Clustering using LD bipartite matching score"

    @classmethod
    def getInputBoxNames(cls):
        return \
            CommonClusteringFunctions.getCommonClusteringInputBoxNames() + [
                ('Select LD graph track', 'ldGraphTrack'),
                ('Select algorithm  for bipartite  matching', 'ldGraphMatching'),
                ('Select r<sup>2</sup> threshold', 'rSquare')
            ]

    @staticmethod
    def getInputBoxOrder():
        return [
            'gSuite',
            'ldGraphTrack',
            'rSquare',
            'ldGraphMatching',
            'distanceMeasure',
            'linkageCriterion',
            'debugMode'
        ]

    @staticmethod
    def getOptionsBoxLdGraphMatching(choices):
        return [
            CommonClusteringFunctions.DEFAULT_SELECT,
            LDBipartiteMatchingTool.LD_GREEDY,
            LDBipartiteMatchingTool.LD_LAPJV
        ]

    @staticmethod
    def getOptionsBoxLdGraphTrack(choices):
        return GeneralGuiTool.getHistorySelectionElement('gtrack')

    @staticmethod
    def getInfoForOptionsBoxLdGraphTrack(choices):
        return 'Linked point track that describe LD relations between the track elements of the GSuite. Can be ' \
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
        cls.htmlClusterTitle(choices.ldGraphMatching, htmlCore)
        cls.htmlClusterSubtext(choices.distanceMeasure, cls.CLUSTER_LIST, choices.linkageCriterion, htmlCore)
        htmlCore.line('Threshold of r<sup>2</sup>: ' + choices.rSquare)

        # Analysis environment
        gSuite = getGSuiteFromGalaxyTN(choices.gSuite)
        analysisBins = GlobalBinSource(gSuite.genome)

        graph = LDExpansions.createRSquareGraph(choices.ldGraphTrack, float(choices.rSquare))
        tracks, labels = LDExpansions.generateTracksAndLabels(gSuite, analysisBins)

        # Find distance/correlation matrices
        size = gSuite.numTracks()
        distDict = cls.createDistDict(cls.CLUSTER_LIST)
        for i in range(0, size):
            for j in range(i + 1, size):

                track1 = tracks[i]
                track2 = tracks[j]

                cost_matrix = BipartiteMatching.generateCostMatrix(track1, track2, graph)

                if choices.ldGraphMatching == LDBipartiteMatchingTool.LD_GREEDY:
                    count = BipartiteMatching.greedyBipartite(cost_matrix)
                else:
                    count = BipartiteMatching.lapjvBipartite(cost_matrix)

                cls.updateDistDict(distDict, count)

        # Cluster and print plots
        cls.printDistPlots(distDict, labels, choices.distanceMeasure, choices.linkageCriterion, galaxyFn, htmlCore)

        cls.htmlClusterTime(str(time.clock() - start), htmlCore)
        htmlCore.divEnd()
        print htmlCore

    @staticmethod
    def getOutputFormat(choices=None):
        return 'html'

    @staticmethod
    def isDebugMode():
        """
        Specifies whether the debug mode is turned on.
        """
        return True

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

        if not choices.ldGraphTrack:
            return 'Please select an LD graph track (linked point track)'

        if choices.rSquare == CommonClusteringFunctions.DEFAULT_SELECT:
            return 'Please select a threshold for r<sup>2</sup>'

        if choices.ldGraphMatching == CommonClusteringFunctions.DEFAULT_SELECT:
            return 'Please select algorithm for graph matching'

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
        core.header('Bipartite matching')
        core.paragraph('This tool compares all pairs of tracks in a GSuite using a bipartite matching technique. '
                       'Each track represent a set of nodes, with weighted edges as given in the LD graph track. '
                       'For each pair of tracks, a binary matching representation is computed, as defined below. '
                       )
        core.paragraph('An LD track for the GSuite can be created with the "LD track generator" tool.')
        core.divider()
        core.smallHeader('Matching definitions')
        core.paragraph('For each match that contribute to the overall matching score, the following is updated: '
                       'a += Bipartite matching score b += (1 - a) / 2 c += (1 - a) / 2 '
                       'Count of SNPs with no edges, i.e no chance of a match score, in each track is added to b and c.')
        core.divider()
        core.smallHeader('Algorithms')
        core.paragraph('The greedy bipartite matching algorithm is fast, but susceptible to local optima. The optimal '
                       'matcher, based on the Jonker-Volgenant algorithm, is much slower, with worst time complexity '
                       'O(n^3). '
                       'It is however guaranteed to find the highest rsquare matching between two tracks.')
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
