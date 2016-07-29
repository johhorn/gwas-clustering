from collections import OrderedDict
from gold.result.HtmlCore import HtmlCore
from quick.webtools.GeneralGuiTool import GeneralGuiTool
from quick.webtools.clustering.CommonClusteringFunctions import CommonClusteringFunctions
from quick.webtools.clustering.LDExpansions import LDExpansions
from quick.webtools.clustering.MatplotlibPlots import MatplotlibPlots
from quick.util.StaticFile import GalaxyRunSpecificFile
from quick.webtools.mixin.DebugMixin import DebugMixin


class EmpiricLDDataTool(GeneralGuiTool, DebugMixin):

    # Empirical cases related to LD
    DEFAULT_SELECT = '--- Select ---'
    LD_CORRELATION = 'Distances between SNPs in LD'
    LD_RSQUARE = 'Distance curves, relative to different r<sup>2</sup>'

    @staticmethod
    def getToolName():
        return "Empirical exploration of linked point tracks"

    @staticmethod
    def isDebugMode():
        return False

    @staticmethod
    def isPublic():
        return True

    @classmethod
    def getInputBoxNames(cls):
        return [
                   ('Select ld graph track', 'ldGraphTrack'),
                   ('Select empirical case', 'empiricalCase'),
                   ('Select r<sup>2</sup> threshold (value between 0 and 1)', 'rSquare'),
                   ('Select multiple r<sup>2</sup> thresholds for comparison', 'multipleRSquares')
               ] + cls.getInputBoxNamesForDebug()

    @staticmethod
    def getOptionsBoxLdGraphTrack():
        return GeneralGuiTool.getHistorySelectionElement('gtrack')

    @staticmethod
    def getOptionsBoxEmpiricalCase(choices):
        return [
            EmpiricLDDataTool.DEFAULT_SELECT,
            EmpiricLDDataTool.LD_CORRELATION,
            EmpiricLDDataTool.LD_RSQUARE
        ]

    @staticmethod
    def getOptionsBoxRSquare(choices):
        if choices.empiricalCase == EmpiricLDDataTool.LD_CORRELATION:
            return '0.8'

    @staticmethod
    def getOptionsBoxMultipleRSquares(choices):
        """
        The rSquare values between 0.7 and 1 can be seen as biologically significant and is therefore checked by
        default. RSquare below 0.7 are seen as noise, but for empirical exploration of how these values behave,
        one can choose to check all values down to 0.
        """
        if choices.empiricalCase == EmpiricLDDataTool.LD_RSQUARE:
            return OrderedDict([
                ('1.0', True),
                ('0.9', True),
                ('0.8', True),
                ('0.7', True),
                ('0.6', False),
                ('0.5', False),
                ('0.4', False),
                ('0.3', False),
                ('0.2', False),
                ('0.1', False),
                ('0.0', False)
            ])

    @classmethod
    def execute(cls, choices, galaxyFn=None, username=''):

        cls._setDebugModeIfSelected(choices)
        htmlCore = HtmlCore()
        htmlCore.divBegin(style=CommonClusteringFunctions.HTML_STYLE)

        if choices.empiricalCase == EmpiricLDDataTool.LD_CORRELATION:
            htmlCore.bigHeader(EmpiricLDDataTool.LD_CORRELATION)
            htmlCore.divider(True)
            cls.getLDDistances(choices.ldGraphTrack, float(choices.rSquare), galaxyFn, htmlCore)
        elif choices.empiricalCase == EmpiricLDDataTool.LD_RSQUARE:
            htmlCore.bigHeader(EmpiricLDDataTool.LD_RSQUARE)
            htmlCore.divider(True)
            cls.getLDDistancesOfMultipleRsquares(choices.ldGraphTrack, choices.multipleRSquares, galaxyFn, htmlCore)

        htmlCore.divEnd()
        print htmlCore

    @classmethod
    def getLDDistances(cls, ldGraphTrack, rSquare, galaxyFn, htmlCore):
        graph = LDExpansions.createRSquareGraph(ldGraphTrack, rSquare)
        positions = LDExpansions.createPositionDict(ldGraphTrack)

        distances = cls.findAllDistancesInLD(graph, positions, rSquare, htmlCore)
        bins = range(0, 525000, 25000)
        cls.plotDistances(distances, galaxyFn, bins, rSquare, htmlCore)

        htmlCore.divider(True)
        htmlCore.header('Exact LD pair count within different distance blocks')
        htmlCore.divBegin(style=CommonClusteringFunctions.TABLE_STYLE)
        htmlCore.tableHeader(['Distance start', 'Distance end', 'Count'])

        from numpy import histogram
        binData = histogram(distances, bins)
        for binNumber in range(0, len(binData[0])):

            binCount = binData[0][binNumber]
            binStart = binData[1][binNumber]
            binEnd = binData[1][binNumber + 1]

            htmlCore.tableRowBegin()
            htmlCore.tableCell(str(binStart))
            htmlCore.tableCell(str(binEnd))
            htmlCore.tableCell(str(binCount))
            htmlCore.tableRowEnd()
        htmlCore.divEnd()

    @classmethod
    def findAllDistancesInLD(cls, graph, positionDict, r2Filter, htmlCore):
        distances = []
        for rsidPair, r2 in graph.items():
            if r2 < r2Filter:
                continue

            rsid1 = rsidPair[0]
            rsid2 = rsidPair[1]

            if rsid1 == rsid2:
                continue

            pos1 = LDExpansions.getPosition(positionDict, rsid1)
            pos2 = LDExpansions.getPosition(positionDict, rsid2)

            if pos1 == -1 or pos2 == -1:
                continue

            distance = abs(pos1 - pos2)
            distances.append(distance)

        cls.printSummary(distances, htmlCore, r2Filter)
        return distances

    @classmethod
    def printSummary(cls, distances, htmlCore, r2Filter):
        from numpy import median, mean
        htmlCore.header('Summary of physical distances between LD pairs with r<sup>2</sup> >= ' + str(r2Filter))
        htmlCore.line('Max distance: ' + str(max(distances)))
        htmlCore.line('Min distance: ' + str(min(distances)))
        htmlCore.line('Average distance: ' + str(mean(distances)))
        htmlCore.line('Median distance: ' + str(median(distances)))
        htmlCore.line('Number of distances: ' + str(len(distances)))

    @classmethod
    def plotDistances(cls, distances, galaxyFn, bins, r2, htmlCore):
        distanceBins, occurrences = cls.standardizeLineGraph(distances)
        distFile = GalaxyRunSpecificFile(['distancegraph.pdf'], galaxyFn)
        MatplotlibPlots.pointGraph(x=distanceBins[1:], y=occurrences, fileLocation=distFile,
                                   xlabel='Distances between SNP pairs in LD', ylabel='LD-pair count')

        dist = sorted(distances)
        histFile = GalaxyRunSpecificFile(['histogram.pdf'], galaxyFn)
        MatplotlibPlots.histogramPlot(dist, bins, histFile, 'Distances between SNP pairs in LD')

        htmlCore.divider(True)
        htmlCore.header('Line plot with distribution of distances between LD pairs, r<sup>2</sup>  >= ' + str(r2))
        htmlCore.line(distFile.getEmbeddedImage())
        htmlCore.link('PDF of distances between tracks here', distFile.getURL())

        htmlCore.divider(True)
        htmlCore.header('Histogram with distribution of distances between LD pairs, r<sup>2</sup> >= ' + str(r2))
        htmlCore.line(histFile.getEmbeddedImage())
        htmlCore.link('PDF histogram of distances here', histFile.getURL())

    @classmethod
    def getLDDistancesOfMultipleRsquares(cls, ldGraphTrack, rSquareThresholds, galaxyFn, htmlCore):

        graph = LDExpansions.createRSquareGraph(ldGraphTrack, 0)
        positions = LDExpansions.createPositionDict(ldGraphTrack)
        ldDistances = []
        rSquareLabels = []
        bins = []

        for rSquare, isSet in rSquareThresholds.items():
            if not isSet:
                continue

            rSquareLabels.append(rSquare)
            rSquare = float(rSquare)
            distances = cls.findAllDistancesInLD(graph, positions, rSquare, htmlCore)
            bins, ldPairCount = cls.standardizeLineGraph(distances)
            ldDistances.append(ldPairCount)

        graphFile = GalaxyRunSpecificFile(['multipleLines.pdf'], galaxyFn)

        MatplotlibPlots.multipleLineGraph(
            bins[1:],
            ldDistances,
            rSquareLabels,
            graphFile,
            'Distance',
            'LD-pair count'
        )

        if len(bins) == 0:
            return

        htmlCore.divider(True)
        htmlCore.header('Distribution of distances between LD pairs with different thresholds of r<sup>2</sup>')
        htmlCore.line(graphFile.getEmbeddedImage())
        htmlCore.link('PDF of distances between tracks here', graphFile.getURL())

    @classmethod
    def standardizeLineGraph(cls, distances):
        from numpy import histogram
        ldPairCount = []
        step = 5000
        bins = range(0, 500000 + step, step)
        binData = histogram(distances, bins)

        for binNumber in range(0, len(binData[0])):
            binCount = binData[0][binNumber]
            ldPairCount.append(binCount)

        return bins, ldPairCount

    @staticmethod
    def getOutputFormat(choices=None):
        return 'html'

    @staticmethod
    def validateAndReturnErrors(choices):

        if not choices.ldGraphTrack:
            return 'Please select an LD graph track (linked point track)'
        if choices.empiricalCase == EmpiricLDDataTool.DEFAULT_SELECT:
            return 'Please select an empirical case'

        if choices.empiricalCase == EmpiricLDDataTool.LD_CORRELATION:
            try:
                r2 = float(choices.rSquare)
                if r2 < 0 or r2 > 1:
                    return 'r<sup>2</sup> must be a value between 0 and 1'
            except:
                return 'r<sup>2</sup> must be numeric'

    @staticmethod
    def getToolDescription():
        from gold.result.HtmlCore import HtmlCore

        core = HtmlCore()
        core.header('Empiric exploration of LD tracks')
        core.paragraph('This tool explores the distribution of physical distance between SNPs in LD, given a threshold '
                       'of r<sup>2</sup>. It takes as input an LD graph as a linked point track, where edges between '
                       'the nodes of the linked point track connect SNPs in LD, weighted by r<sup>2</sup>. The LD '
                       'track generator tool can be used to create an LD track for a specific GSuite.')
        core.divider()
        core.smallHeader(EmpiricLDDataTool.LD_CORRELATION)
        core.paragraph('With this tool option, one value of r<sup>2</sup> is investigated for the given LD track. A '
                       'histogram and a line graph is plotted, showing distribution of the LD-pairs relative to '
                       'the physical distance between them.')
        core.divider()
        core.smallHeader(EmpiricLDDataTool.LD_RSQUARE)
        core.paragraph('This tool option allows investigation of several thresholds of r<sup>2</sup>. Can be used to '
                       'argue effect of r<sup>2</sup> thresholds on the data set.')
        core.divider()
        core.smallHeader('NB! Preset thresholds of r<sup>2</sup> in linked point track')
        core.paragraph('Note that the linked point track might have a threshold of r<sup>2</sup> already set, and '
                       'no edges with values below this threshold. This tool cannot detect which thresholds have been'
                       'applied previously. If values are unexpected, check that the tool threshold is not lower than'
                       'the actual r<sup>2</sup> of the linked point edges. A way to discover a preset threshold '
                       'within the linked point track, is to use the ' + EmpiricLDDataTool.LD_RSQUARE + ' option with '
                       'all values of r<sup>2</sup>. If the graph have fewer lines than the number of r<sup>2</sup> '
                       'labels, it will indicate a preset threshold.')

        return core
