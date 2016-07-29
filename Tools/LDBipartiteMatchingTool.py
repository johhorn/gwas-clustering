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
        return False

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
        core.paragraph('This tool uses a bipartite matching technique to find the best matching, in terms of maximal '
                       'weights overall. A SNP can only be matched once with a SNP in the other track, and there '
                       'is no intra-track matching.'
                       'For each match that contribute to the overall matching score, the following is updated: '
                       '<ul>'
                       '<li>a: Overall bipartite matching score </li>'
                       '<li>b: Number of SNPs without edges in first track</li>'
                       '<li>c: Number of SNPs without edges in second track</li>'
                       '<li>d: Not defined</li>'
                       '</ul>')
        core.divider()
        core.smallHeader('Algorithms')
        core.paragraph('Two algorithms are provided to compute the matching scores. The first is a greedy heuristic, '
                       'which is fast, but susceptible to local optima. The optimal '
                       'matcher, based on the Jonker-Volgenant algorithm, is much slower, with worst time complexity '
                       'O(n<sup>3</sup>), '
                       'but is guaranteed to find the highest <sup>2</sup> matching between two tracks.')
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
