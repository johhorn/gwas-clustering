from gold.gsuite import GSuiteConstants
from quick.multitrack.MultiTrackCommon import getGSuiteFromGalaxyTN
from quick.util.StaticFile import GalaxyRunSpecificFile
from quick.webtools.GeneralGuiTool import GeneralGuiTool
from quick.webtools.clustering.MatplotlibPlots import MatplotlibPlots
from quick.webtools.mixin.DebugMixin import DebugMixin


class CommonClusteringFunctions(GeneralGuiTool):
    """
    Suite of functions to use in clustering tools.

    If getCommonClusteringInputBoxNames are used, it is recommended to define getInputBoxOrder, in order to get the
    needed input fields in the correct order. Any key not in getInputBoxOrder, although defined in this class,
    will not appear in the actual tool, so that one can customize the common list.
    """

    # HTML style settings
    HTML_STYLE = 'font-family:serif; font-size:12pt'
    TABLE_STYLE = 'font-family:sans-serif; font-size:10pt'

    # Similarity measures
    SIM_JACCARD = 'Jaccard'
    SIM_COSINE = 'Cosine'
    SIM_SIMPSON = 'Simpson'
    SIM_SORGENFREI = 'Sorgenfrei'
    SIM_OTSUKA = 'Otsuka'
    SIM_KULCZYNSKI = 'Kulczynski'
    SIM_FORBES = 'Forbes'
    SIM_MCCONNAUGHEY = 'McConnaughey'

    # Correlation measure:
    CORR_PEARSON = 'Pearson'

    # For comparison of methods:
    ALL_MEASURES = 'All measures in list'

    # Lists of measures:
    DISTLIST = [
        SIM_JACCARD,
        SIM_COSINE,
        SIM_SORGENFREI,
        SIM_OTSUKA,
        SIM_SIMPSON,
        SIM_KULCZYNSKI,
        SIM_MCCONNAUGHEY
     ]

    CORR_DISTLIST = DISTLIST + [
        SIM_FORBES,
        CORR_PEARSON
    ]

    # Linkage criterion
    SINGLE = 'single'
    AVG = 'average'
    COMPLETE = 'complete'
    WEIGHTED = 'weighted'

    # Error handling
    GSUITE_ALLOWED_FILE_FORMATS = [GSuiteConstants.PREPROCESSED]
    GSUITE_ALLOWED_LOCATIONS = [GSuiteConstants.LOCAL]
    GSUITE_ALLOWED_TRACK_TYPES = [GSuiteConstants.VALUED_POINTS]
    GRAPH_ALLOWED_FILE_FORMATS = [GSuiteConstants.LINKED_POINTS]

    DEFAULT_SELECT = '--- Select ---'

    @classmethod
    def getCommonClusteringInputBoxNames(cls):
        return [
            ('Select GSuite from history', 'gSuite'),
            ('Select distance measure', 'distanceMeasure'),
            ('Select distance measure', 'distanceCorrMeasure'),
            ('Select linkage criterion', 'linkageCriterion')
        ] + DebugMixin.getInputBoxNamesForDebug()

    @staticmethod
    def getOptionsBoxGSuite():
        return GeneralGuiTool.getHistorySelectionElement('gsuite')

    @staticmethod
    def getOptionsBoxLinkageCriterion(choices):
        return [
            CommonClusteringFunctions.DEFAULT_SELECT,
            CommonClusteringFunctions.SINGLE,
            CommonClusteringFunctions.COMPLETE,
            CommonClusteringFunctions.AVG,
            CommonClusteringFunctions.WEIGHTED
        ]

    @staticmethod
    def getOptionsBoxDistanceMeasure(prevChoices):
        return [
            CommonClusteringFunctions.DEFAULT_SELECT,
            CommonClusteringFunctions.ALL_MEASURES,
        ] + CommonClusteringFunctions.DISTLIST

    @staticmethod
    def getOptionsBoxDistanceCorrMeasure(prevChoices):
        return [
                   CommonClusteringFunctions.DEFAULT_SELECT,
                   CommonClusteringFunctions.ALL_MEASURES,
               ] + CommonClusteringFunctions.CORR_DISTLIST


    @classmethod
    def printClusterPlots(cls, correlationMatrix, linkageMatrix, galaxyFn, distanceMeasure, labels, htmlCore):
        from numpy import amax, amin, isnan
        maxVal = amax(correlationMatrix)
        minVal = amin(correlationMatrix)

        seabornFile = GalaxyRunSpecificFile(['Image', distanceMeasure + 'seabornHeatmap.pdf'], galaxyFn)
        dendrogramFile = GalaxyRunSpecificFile(['Image', distanceMeasure + 'dendrogram.pdf'], galaxyFn)

        if minVal < 0 or isnan(minVal):
            MatplotlibPlots.seabornHeatmapPlot(
                correlationMatrix,
                labels,
                max=maxVal if maxVal >= 1 else 1,
                min=minVal if minVal <= -1 else -1,
                fileName=seabornFile,
                cmap="RdBu_r"
            )
        else:
            MatplotlibPlots.seabornHeatmapPlot(
                correlationMatrix,
                labels,
                max=maxVal if maxVal >= 1 else 1,
                min=minVal if minVal <= 0 else 0,
                fileName=seabornFile
            )

        MatplotlibPlots.dendrogramClusteringPlot(linkageMatrix, labels, dendrogramFile)
        htmlCore.line(seabornFile.getEmbeddedImage())
        htmlCore.link('PDF of similarity matrix', seabornFile.getURL())
        htmlCore.line(dendrogramFile.getEmbeddedImage())
        htmlCore.link('PDF of dendrogram', dendrogramFile.getURL())

    @classmethod
    def printTextMatrixes(cls, correlationMatrix, linkageMatrix, distanceMatrix, galaxyFn, filename, htmlCore):

        # Print correlation matrix
        corrMatrixFile = GalaxyRunSpecificFile(['corr_matrix_result_' + filename + '.txt'], galaxyFn)
        corrMatrixPath = corrMatrixFile.getDiskPath(True)
        open(corrMatrixPath, 'w').write(str(correlationMatrix))
        htmlCore.link('<br><br>View the raw text similarity/correlation matrix for this analysis',
                      corrMatrixFile.getURL())

        # Print distance matrix
        distMatrixFile = GalaxyRunSpecificFile(['dist_matrix_result_' + filename + '.txt'], galaxyFn)
        distMatrixPath = distMatrixFile.getDiskPath(True)
        open(distMatrixPath, 'w').write(str(distanceMatrix))
        htmlCore.link('<br><br>View the raw text triangular distance matrix for this analysis', distMatrixFile.getURL())

        # Print linkage matrix
        linkMatrixFile = GalaxyRunSpecificFile(['linkage_matrix_result_' + filename + '.txt'], galaxyFn)
        linkMatrixPath = linkMatrixFile.getDiskPath(True)
        open(linkMatrixPath, 'w').write(str(linkageMatrix))
        htmlCore.link('<br><br>View the raw text linkage matrix for this analysis', linkMatrixFile.getURL())

    @classmethod
    def validateGSuite(cls, choices):
        errorString = cls._checkGSuiteFile(choices.gSuite)
        if errorString:
            return errorString

        gSuite = getGSuiteFromGalaxyTN(choices.gSuite)
        errorString = GeneralGuiTool._checkGSuiteRequirements(
            gSuite,
            allowedLocations=CommonClusteringFunctions.GSUITE_ALLOWED_LOCATIONS,
            allowedFileFormats=CommonClusteringFunctions.GSUITE_ALLOWED_FILE_FORMATS,
            allowedTrackTypes=CommonClusteringFunctions.GSUITE_ALLOWED_TRACK_TYPES)
        if errorString:
            return errorString

        errorString = GeneralGuiTool._checkGSuiteTrackListSize(gSuite, 2)
        if errorString:
            return errorString


    @classmethod
    def roundFloat(cls, distance):
        """
        Rounds of a float to avoid dissimilarity due to floating points
        """
        return int((distance * 1000000) + 0.000001) / 1000000.0

    @classmethod
    def getDistance(cls, res, similarityMetric):
        from math import sqrt

        a = float(res['a'])
        b = float(res['b'])
        c = float(res['c'])

        denominator = 1
        nominator = 0

        if similarityMetric == CommonClusteringFunctions.SIM_JACCARD:
            denominator = (a + b + c)
            nominator = a
        elif similarityMetric == CommonClusteringFunctions.SIM_COSINE:
            denominator = (sqrt(a + b) * sqrt(a + c))
            nominator = a
        elif similarityMetric == CommonClusteringFunctions.SIM_SIMPSON:
            denominator = min(a + b, a + c)
            nominator = a
        elif similarityMetric == CommonClusteringFunctions.SIM_OTSUKA:
            denominator = sqrt((a + b) * (a + c))
            nominator = a
        elif similarityMetric == CommonClusteringFunctions.SIM_SORGENFREI:
            nominator = a ** 2
            denominator = (a + b) * (a + c)
        elif similarityMetric == CommonClusteringFunctions.SIM_KULCZYNSKI:
            nominator = (a / 2) * (2 * a + b + c)
            denominator = (a + b) * (a + c)
        elif similarityMetric == CommonClusteringFunctions.SIM_FORBES:
            """
            Originally a measure for expected overlap. Normalized to a measure of distance
            NB: Only to be used on binary count dictionaries that have negative matches (d) defined.
            """
            from math import log
            d = float(res['d'])
            n = a + b + c + d
            nominator = n * a
            denominator = (a + c) * (a + b)
            if denominator > 0:
                forbes = nominator / denominator
                return 1 / max(1, log(forbes) + 1) if forbes != 0 else 1
            return 1

        elif similarityMetric == CommonClusteringFunctions.SIM_MCCONNAUGHEY:
            """
            A correlation coefficient, the measure will be between -1 and 1
            """
            nominator = (a ** 2) - (b * c)
            denominator = (a + b) * (a + c)
            return nominator / denominator if denominator != 0 else 0

        elif similarityMetric == CommonClusteringFunctions.CORR_PEARSON:
            """
            A correlation coefficient, the measure will be between -1 and 1
            NB: Only to be used on binary count dictionaries that have negative matches (d) defined.
            """
            d = float(res['d'])
            nominator = (a * d) - (b * c)
            denominator = sqrt((a + b) * (d + c) * (a + c) * (d + b))
            return nominator / denominator

        similarity = nominator / denominator if denominator != 0 else 0
        return 1 - cls.roundFloat(similarity)

    @classmethod
    def createDistDict(cls, allMeasures):
        """
        Create dictionary for a range of different distance/correlation measures.
        Made to contain triangular distance/correlation matrices of different measures.
        The specified measures are stored inside the dictionary.
        """
        distDict = {'availableMetrics': allMeasures}
        for dist in allMeasures:
            distDict[dist] = []

        return distDict

    @classmethod
    def updateDistDict(cls, distDict, matches):
        """
        Update all distance/correlation measures of the dictionary with the matches passed to the function.
        Will add the results to the correct triangular distance/correlation matrix.
        """
        for dist in cls.getDistDictKeys(distDict):
            distance = cls.getDistance(matches, dist)
            distDict[dist].append(distance)

    @classmethod
    def getDistDictKeys(cls, distDict):
        """Get the distance/correlation measures the dictionary contain."""
        return distDict['availableMetrics']

    @classmethod
    def getDistMatrixes(cls, distDict, distMeasure, linkageCriterion):
        """
        Find and return the correlation matrix, linkage matrix and distance matrix for the distance/correlation
        measure given with distMeasure parameter.
        """
        from scipy.spatial.distance import squareform
        from numpy import ones, fill_diagonal
        from scipy.cluster.hierarchy import linkage

        if distMeasure == cls.CORR_PEARSON or distMeasure == cls.SIM_MCCONNAUGHEY:
            '''As these measures generate values between -1 and 1, need special handling'''

            # Cluster distances, i.e. convert correlation into distance between 0 and 1
            triangularCorrMatrix = distDict[distMeasure]
            triangularDistMatrix = ones(len(triangularCorrMatrix)) - [(x + 1) / 2 for x in triangularCorrMatrix]
            linkageMatrix = linkage(cls.removeNanDistances(triangularDistMatrix), linkageCriterion)

            # Make correlation matrix square
            correlationMatrix = squareform(triangularCorrMatrix)
            fill_diagonal(correlationMatrix, 1)
        else:

            # Cluster distances
            triangularDistMatrix = distDict[distMeasure]
            linkageMatrix = linkage(cls.removeNanDistances(triangularDistMatrix), linkageCriterion)

            # Convert triangular distances into square correlation matrix
            squareDistMatrix = squareform(triangularDistMatrix)
            squareSize = len(squareDistMatrix)
            correlationMatrix = ones((squareSize, squareSize)) - squareDistMatrix

        return correlationMatrix, linkageMatrix, triangularDistMatrix

    @classmethod
    def removeNanDistances(cls, triangularDistanceMatrix):
        from numpy import isnan
        for i in range(0, len(triangularDistanceMatrix)):
            if isnan(triangularDistanceMatrix[i]):
                triangularDistanceMatrix[i] = 1

        return triangularDistanceMatrix

    @classmethod
    def printDistPlots(cls, distDict, labels, distanceMeasure, linkageCriterion, galaxyFn, htmlCore):

        if distanceMeasure == cls.ALL_MEASURES:
            for measure in cls.getDistDictKeys(distDict):
                cls.printForOneDistanceMeasure(distDict, galaxyFn, htmlCore, labels, linkageCriterion, measure)
        else:
            cls.printForOneDistanceMeasure(distDict, galaxyFn, htmlCore, labels, linkageCriterion, distanceMeasure)

    @classmethod
    def printForOneDistanceMeasure(cls, distDict, galaxyFn, htmlCore, labels, linkageCriterion, measure):
        htmlCore.divider(True)
        measureType = 'Correlation' if measure == cls.CORR_PEARSON else 'Similarity'
        htmlCore.smallHeader(measureType + ' matrix and clustering of distances with ' + measure)
        htmlCore.line('<br>')
        corr, linkage, distance = cls.getDistMatrixes(distDict, measure, linkageCriterion)
        cls.printClusterPlots(corr, linkage, galaxyFn, measure, labels, htmlCore)
        cls.printTextMatrixes(corr, linkage, distance, galaxyFn, measure, htmlCore)
        cls.findRanking(distance, labels, measure, htmlCore)


    @classmethod
    def htmlClusterTitle(cls, clusterCase, htmlCore):
        htmlCore.header('Clustering of ' + clusterCase.lower())

    @classmethod
    def htmlClusterSubtext(cls, distanceMeasure, distList, linkageCriterion, htmlCore):

        if distanceMeasure == cls.ALL_MEASURES:
            htmlCore.line('Clustering is done for all available distance and correlation measures, to facilitate '
                          'comparison of the different measures on binary data representations: '
                          '<ul>' + ' '.join(['<li>' + dist + '</li>' for dist in distList]) +
                          '</ul> ')

        else:
             htmlCore.line('Chosen distance/correlation measure for clustering on binary data representations: ' + distanceMeasure)

        htmlCore.line('Linkage criterion for clustering dendrograms is ' + linkageCriterion)

    @classmethod
    def htmlClusterTime(cls, time, htmlCore):
        htmlCore.divider()
        htmlCore.line('Total execute time: ' + time)

    @classmethod
    def checkClusterOptions(cls, dist, link):
        if dist == cls.DEFAULT_SELECT:
            return 'Please select a distance measure'

        if link == cls.DEFAULT_SELECT:
            return 'Please select a linkage criterion'

    @classmethod
    def findRanking(cls, triangularDistance, labels, measure, htmlCore):
        """
        Prints order of smallest distance. It is an alternative to the clustering dendrograms and linkage matrix.
        The rank is 1 for the smallest distance between a pair of tracks, and increases with distance until all the
        pairs have a rank.
        """

        rankMatrix = cls.createRankMatrix(triangularDistance)

        htmlCore.smallHeader('<br><br><br>Rank matrix of ' + measure)
        htmlCore.paragraph(
            'Matrix containing order, or rank, of distances between tracks. A rank of 1 denotes the pair '
            'of tracks with the smallest distance between them, a rank of ' + str(len(triangularDistance)) +
            ' means the pair has the largest distance between them.')
        htmlCore.divEnd()
        htmlCore.divBegin(style=CommonClusteringFunctions.TABLE_STYLE)
        htmlCore.tableHeader([''] + labels)

        labelIndex = 0
        for row in rankMatrix:
            htmlCore.tableRowBegin()

            htmlCore.tableCell(labels[labelIndex])
            for cell in row:
                if cell == 0:
                    htmlCore.tableCell('X')
                else:
                    htmlCore.tableCell(str(int(cell)))

            htmlCore.tableRowEnd()
            labelIndex += 1

        htmlCore.tableFooter()
        htmlCore.divEnd()
        htmlCore.divBegin(style=CommonClusteringFunctions.HTML_STYLE)

    @classmethod
    def createRankMatrix(cls, distances):
        """
        Finds the rank order of the triangular distance matrix passed to the function.
        Returns a square (symmetric) matrix of the ranks, where the diagonal (track distance with itself) is set to 0.
        """
        from numpy import argmin, zeros
        from scipy.spatial.distance import squareform

        distanceCount = len(distances)
        rankings = zeros(distanceCount)
        rank = 1
        while rank <= distanceCount:
            pos = argmin(distances)
            rankings[pos] = rank
            distances[pos] = 10
            rank += 1

        return squareform(rankings)




