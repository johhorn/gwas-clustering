from gold.result.HtmlCore import HtmlCore
from gold.statistic.NearestPointDistsStat import NearestPointDistsStat
from quick.statistic.PointGapsStat import PointGapsStat
from quick.statistic.UniquePointTrackStat import UniquePointTrackStat
from quick.webtools.GeneralGuiTool import GeneralGuiTool
from quick.multitrack.MultiTrackCommon import getGSuiteFromGalaxyTN
from gold.application.HBAPI import doAnalysis, AnalysisSpec, GlobalBinSource
from quick.webtools.clustering.CommonClusteringFunctions import CommonClusteringFunctions
from quick.webtools.clustering.MatplotlibPlots import MatplotlibPlots
from quick.util.StaticFile import GalaxyRunSpecificFile
from numpy import asarray, log
from gold.application.HBAPI import Track
from quick.webtools.mixin.DebugMixin import DebugMixin
from quick.webtools.restricted.visualization.visualizationGraphs import visualizationGraphs


class SNPStatsTool(GeneralGuiTool, DebugMixin):

    # Empirical cases related to SNP stats of gSuite
    SNP_COUNT = 'SNP frequencies in GSuite'
    SNP_WITHIN = 'Distances within tracks'
    SNP_BETWEEN = 'Distances between tracks'

    @staticmethod
    def getToolName():
        return "Tool for exploration of GSuite of valued point tracks"

    @staticmethod
    def isPublic():
        return True

    @staticmethod
    def isDebugMode():
        return False

    @classmethod
    def getInputBoxNames(cls):
        return [
                   ('Select GSuite from history', 'gSuite'),
                   ('Select empirical case', 'getEmpiricalCase')
               ] + cls.getInputBoxNamesForDebug()

    @staticmethod
    def getOptionsBoxGSuite():
        return GeneralGuiTool.getHistorySelectionElement('gsuite')

    @staticmethod
    def getOptionsBoxGetEmpiricalCase(choices):
        return [
                    CommonClusteringFunctions.DEFAULT_SELECT,
                    SNPStatsTool.SNP_COUNT,
                    SNPStatsTool.SNP_WITHIN,
                    SNPStatsTool.SNP_BETWEEN
            ]

    @classmethod
    def execute(cls, choices, galaxyFn=None, username=''):

        cls._setDebugModeIfSelected(choices)

        gSuite = getGSuiteFromGalaxyTN(choices.gSuite)
        analysisBins = GlobalBinSource(gSuite.genome)

        htmlCore = HtmlCore()
        htmlCore.divBegin(style=CommonClusteringFunctions.HTML_STYLE)

        if choices.getEmpiricalCase == SNPStatsTool.SNP_COUNT:
            htmlCore.bigHeader(SNPStatsTool.SNP_COUNT)
            htmlCore.divider(True)
            cls.getSNPFrequencyStats(analysisBins, gSuite, galaxyFn, htmlCore)
        elif choices.getEmpiricalCase == SNPStatsTool.SNP_BETWEEN:
            htmlCore.bigHeader(SNPStatsTool.SNP_BETWEEN)
            htmlCore.divider(True)
            cls.getDistancesBetween(analysisBins, gSuite, galaxyFn, htmlCore)
        elif choices.getEmpiricalCase == SNPStatsTool.SNP_WITHIN:
            htmlCore.bigHeader(SNPStatsTool.SNP_WITHIN)
            htmlCore.divider(True)
            cls.getDistancesWithin(analysisBins, gSuite, galaxyFn, htmlCore)

        print htmlCore

    @classmethod
    def getDistancesWithin(cls, analysisBins, gSuite, galaxyFn, htmlCore):
        """
        Finds all distances between points on the same tracks, for all tracks in a GSuite
        """
        distances = []
        analysisSpec = AnalysisSpec(PointGapsStat)

        for gSuiteTrack in gSuite.allTracks():
            track = Track(gSuiteTrack.trackName)
            result = doAnalysis(analysisSpec, analysisBins, [track])
            cls.addDistances(distances, result.getGlobalResult()['Result'])

        cls.printStats(distances, 'distance', htmlCore)
        cls.plotDistances(distances, galaxyFn, 'within', htmlCore)

    @classmethod
    def getDistancesBetween(cls, analysisBins, gSuite, galaxyFn, htmlCore):
        """
        Computes distances for one pair of track in a GSuite at a time.
        Finds the smallest distance to a point in the other track for all points in each track.
        These distances might be asymmetric, i.e. smallest distance between point in track1 -> closest point in track2
        could be different from smallest distance between the same point in track2 -> closest point in track1
        """

        distances = []
        analysisSpec = AnalysisSpec(NearestPointDistsStat)
        allTracks = gSuite.allTracks()

        for gSuiteTrack in allTracks:
            track = Track(gSuiteTrack.trackName)
            for gSuiteTrack2 in allTracks:
                if gSuiteTrack.trackName != gSuiteTrack2.trackName:
                    track2 = Track(gSuiteTrack2.trackName)
                    tracks = [track, track2]
                    result = doAnalysis(analysisSpec, analysisBins, tracks)
                    cls.addDistances(distances, result.getGlobalResult()['Result'])

        cls.printStats(distances, 'distance', htmlCore)
        htmlCore.line('Asymmetries related to alternating shortest distance: ' + str(cls.countAsymmetries(distances)))
        cls.plotDistances(distances, galaxyFn, 'between', htmlCore)

    @classmethod
    def printStats(cls, objectList, type, htmlCore):
        from numpy import mean, median

        htmlCore.header('Summary statistics of ' + type + 's:')
        htmlCore.line('Average ' + type + ' length: ' + str(mean(asarray(objectList))))
        htmlCore.line('Median ' + type + ' length: ' + str(median(asarray(objectList))))
        htmlCore.line('Max ' + type + ' length: ' + str(max(objectList)))
        htmlCore.line('Min ' + type + ' length: ' + str(min(objectList)))
        htmlCore.line('Number of ' + type + 's: ' + str(len(objectList)))

    @staticmethod
    def getOutputFormat(choices=None):
        return 'html'

    @classmethod
    def getSNPFrequencyStats(cls, bins, gSuite, galaxyFn, htmlCore):
        rsIDs = set()
        snpCount = []
        analysisSpec = AnalysisSpec(UniquePointTrackStat)
        trackLabels = []

        for gSuiteTrack in gSuite.allTracks():
            track = Track(gSuiteTrack.trackName)
            trackLabels.append(gSuiteTrack.title)
            result = doAnalysis(analysisSpec, bins, [track])
            if 'Result' in result.getGlobalResult():
                observations = result.getGlobalResult()['Result']
                snpCount.append(len(observations))
                rsIDs.update(observations)

        snpcountFile = GalaxyRunSpecificFile(['snpfrequencies.pdf'], galaxyFn)
        MatplotlibPlots.pointGraphY(snpCount, snpcountFile, ylabel='SNP counts',
                                    xticks=trackLabels)

        snpdistributionFile = GalaxyRunSpecificFile(['snpfreqhistogram.pdf'], galaxyFn)
        MatplotlibPlots.histogramRugPlot(snpCount, 10, snpdistributionFile, 'SNP counts')

        totalSNPCount = sum(snpCount)
        cls.printStats(snpCount, 'track', htmlCore)
        htmlCore.line('Total number of SNPs: ' + str(totalSNPCount))
        htmlCore.line('Unique SNPs: ' + str(len(rsIDs)))
        htmlCore.line('Overlapping rsIDs: ' + str(totalSNPCount - len(rsIDs)))

        htmlCore.divider(True)
        htmlCore.header('Graph of SNP frequencies in GSuite tracks')
        htmlCore.line(snpcountFile.getEmbeddedImage())
        htmlCore.link('PDF of SNP frequency graph', snpcountFile.getURL())

        htmlCore.divider(True)
        htmlCore.header('Histogram of SNP frequencies in GSuite tracks')
        htmlCore.line(snpdistributionFile.getEmbeddedImage())
        htmlCore.link('PDF of SNP frequency histogram', snpdistributionFile.getURL())

        cls.getInteractiveColumnChartWithLabels(snpCount, trackLabels, htmlCore)

    @classmethod
    def plotDistances(cls, distances, galaxyFn, distCase, htmlCore):

        # Plot distance graph
        xdata, ydata = cls.standardizeLineGraph(distances)
        distFile = GalaxyRunSpecificFile(['distancegraph.pdf'], galaxyFn)
        MatplotlibPlots.pointGraph(xdata[1:], ydata, distFile, 'Smallest distance for each point',
                                   'Distance point count')
        # Write distance graph
        htmlCore.divider(True)
        htmlCore.header('Graph of smallest distances for all points ' + distCase + ' tracks in GSuite')
        htmlCore.line(distFile.getEmbeddedImage())
        htmlCore.link('PDF of distance graph', distFile.getURL())

        # Plot distance histograms
        dist = sorted(distances)
        bins = 20
        histFile = GalaxyRunSpecificFile(['histogram.pdf'], galaxyFn)
        loghistFile = GalaxyRunSpecificFile(['loghistogram.pdf'], galaxyFn)
        MatplotlibPlots.histogramRugPlot(dist, bins, histFile, 'Distances')
        MatplotlibPlots.histogramRugPlot(log(dist), bins, loghistFile, 'Log of distances')
        helperText = 'The rugs/vertical lines at the bottom show the distribution of point distances.'

        # Write distance histograms
        htmlCore.divider(True)
        htmlCore.header('Histogram of smallest distances for all points ' + distCase +
                        ' tracks in GSuite')
        htmlCore.line(helperText)
        htmlCore.line(histFile.getEmbeddedImage())
        htmlCore.link('PDF of distance histogram', histFile.getURL())

        htmlCore.header('Histogram of log of smallest distances for all points ' + distCase +
                        ' tracks in GSuite')
        htmlCore.line(helperText)
        htmlCore.line(loghistFile.getEmbeddedImage())
        htmlCore.link('PDF of log distance histogram', loghistFile.getURL())

        # Plot and write interactive bar chart
        cls.getInteractiveColumnChart(dist, distCase, htmlCore)

    @classmethod
    def standardizeLineGraph(cls, distances):
        from numpy import histogram
        distanceCount = []
        step = 1000000

        bins = range(0, int(1e8)+step, step)
        binData = histogram(distances, bins)

        for binNumber in range(0, len(binData[0])):
            binCount = binData[0][binNumber]
            distanceCount.append(binCount)

        return bins, distanceCount

    @classmethod
    def countAsymmetries(cls, distances):
        count = 0
        distDict = {}
        for distance in distances:
            if distance in distDict:
                distDict[distance] += 1
            else:
                distDict[distance] = 1

        for key, value in distDict.items():
            if value % 2:
                count += 1

        return count

    @classmethod
    def getInteractiveColumnChart(cls, distances, distCase, htmlCore):
        htmlCore.divider()
        htmlCore.header('Interactive distance plot')
        htmlCore.line('To see the exact distance for the different points ' + distCase +
                      ' tracks, hover over the bars of interest.')
        htmlCore.divEnd()  # End div set in execute to remove font settings

        vg = visualizationGraphs()
        lines = vg.drawColumnChart(
            dataY=distances,
            titleText='Smallest distances for all points ' + distCase + ' tracks',
            yAxisTitle='Distances'
        )

        htmlCore.line(lines)

    @classmethod
    def getInteractiveColumnChartWithLabels(cls, data, labels, htmlCore):
        htmlCore.divider()
        htmlCore.header('Interactive SNP frequency plot')
        htmlCore.line('To see exact SNP counts for the different tracks, hover over the bars of interest.')
        htmlCore.divEnd()  # End div set in execute to remove font settings

        vg = visualizationGraphs()
        lines = vg.drawColumnChart(
            dataY=data,
            titleText='SNP frequencies in tracks',
            yAxisTitle='SNP count',
            categories=labels,
            xAxisRotation=90
        )

        htmlCore.line(lines)

    @classmethod
    def addDistances(cls, distances, result):
        for distance in result:
            if distance:
                distances.append(distance)

    @staticmethod
    def validateAndReturnErrors(choices):
        errorString = CommonClusteringFunctions.validateGSuite(choices)
        if errorString:
            return errorString

        if choices.getEmpiricalCase == CommonClusteringFunctions.DEFAULT_SELECT:
            return 'Please select an empirical case'

    @staticmethod
    def getToolDescription():
        from gold.result.HtmlCore import HtmlCore

        core = HtmlCore()
        core.header('Empiric exploration of valued point tracks')
        core.paragraph('This tool prints out SNP information of a given GSuite.  Its purpose is to give a general '
                       'overview of track summary statistics, and show internal and external relationships of '
                       'distances between points within and between tracks.')
        core.divider()
        core.smallHeader(SNPStatsTool.SNP_COUNT)
        core.paragraph('This option generate summary statistics for SNP count across all tracks. Different '
                       'plots illustrate the SNP frequencies for each individual track.')
        core.divider()
        core.smallHeader('Distances between and within tracks')
        core.paragraph('These options generate summary statistics of distances of all points, either within or'
                       'between tracks. In addition, several plots are generated to illustrate the properties of '
                       'distances between points in the given GSuite.')

        return str(core)