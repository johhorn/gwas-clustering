from quick.webtools.clustering.CommonClusteringFunctions import CommonClusteringFunctions


class CommonCorrelationFunctions(CommonClusteringFunctions):
    """
    Suite of functions to use in correlation clustering tools.

    Similar to the CommonClusteringFunctions, which it inherits. Some of the functions in the super class have been
    declared for correlation specific purposes. Others are used directly.
    """

    # Correlation measures
    CORR_PEARSON = 'Pearson'
    CORR_SPEARMAN = 'Spearman'

    @classmethod
    def getCommonCorrelationInputBoxes(cls):
        return [
                    ('Select GSuite from history', 'gSuite'),
                    ('Select correlation coefficient', 'corrStat'),
                    ('Select linkage criterion', 'linkageCriterion')
                ]

    @staticmethod
    def getOptionsBoxCorrStat(choices):
        return [
            CommonClusteringFunctions.DEFAULT_SELECT,
            CommonClusteringFunctions.ALL_MEASURES,
            CommonCorrelationFunctions.CORR_PEARSON,
            CommonCorrelationFunctions.CORR_SPEARMAN
        ]

    @classmethod
    def updateCorrDict(cls, corrDict, vector1, vector2):
        """Similar to the updateDistDict of CommonClusteringFunctions"""

        for corr in cls.getDistDictKeys(corrDict):
            corrDict[corr].append(cls.getCorrelation(vector1, vector2, corr))

    @classmethod
    def getCorrelation(cls, vector1, vector2, corrStat):
        from scipy.stats import pearsonr, spearmanr

        if corrStat == cls.CORR_PEARSON:
            return pearsonr(vector1, vector2)[0]
        elif corrStat == cls.CORR_SPEARMAN:
            return spearmanr(vector1, vector2)[0]

    @classmethod
    def _toSquareMatrix(cls, triangularMatrix):
        """
        Takes in correlation matrix with values from -1 to 1. Makes it square and adds 1 for the diagonals,
        where tracks are perfectly correlated (with themselves).
        """
        from scipy.spatial.distance import squareform
        square = squareform(triangularMatrix)
        size = len(square)
        for i in range(0, size):
            square[i, i] = 1
        return square

    @classmethod
    def printCorrPlots(cls, corrDict, labels, corrStat, linkageCriterion, galaxyFn, htmlCore):
        if corrStat == cls.ALL_MEASURES or corrStat == cls.CORR_PEARSON:
            cls.printForOneCorrCoeff(corrDict, galaxyFn, htmlCore, labels, linkageCriterion, cls.CORR_PEARSON)
        if corrStat == cls.ALL_MEASURES or corrStat == cls.CORR_SPEARMAN:
            cls.printForOneCorrCoeff(corrDict, galaxyFn, htmlCore, labels, linkageCriterion, cls.CORR_SPEARMAN)

    @classmethod
    def printForOneCorrCoeff(cls, corrDict, galaxyFn, htmlCore, labels, linkageCriterion, corrStat):
        """
        Prints clustering, raw text matrixes and ranking for the correlation matrixes. Uses functions defined
        in CommonClusteringFunctions.
        """

        htmlCore.divider(True)
        htmlCore.smallHeader('Correlation with ' + corrStat)
        corr, linkage, distance = cls.getDistMatrixes(corrDict, corrStat, linkageCriterion)
        cls.printClusterPlots(corr, linkage, galaxyFn, corrStat, labels, htmlCore)
        cls.printTextMatrixes(corr, linkage, distance, galaxyFn, corrStat, htmlCore)
        cls.findRanking(distance, labels, corrStat, htmlCore)

    @classmethod
    def getDistMatrixes(cls, corrDict, corrCoeff, linkageCriterion):
        """
        Overrides function in class CommonClusteringFunctions.
        Is defined for Pearson and Spearman only, and not the distance measures defined in the super class.
        Needed for the functions called in printForOneCorrCoeff.
        """
        from scipy.spatial.distance import squareform
        from numpy import ones, fill_diagonal
        from scipy.cluster.hierarchy import linkage

        # Cluster distances, i.e. convert correlation into distance between 0 and 1
        triangularCorrMatrix = corrDict[corrCoeff]
        triangularDistMatrix = ones(len(triangularCorrMatrix)) - [(x + 1) / 2 for x in triangularCorrMatrix]
        linkageMatrix = linkage(cls.removeNanDistances(triangularDistMatrix), linkageCriterion)

        # Make correlation matrix square
        correlationMatrix = squareform(triangularCorrMatrix)
        fill_diagonal(correlationMatrix, 1)

        return correlationMatrix, linkageMatrix, triangularDistMatrix

    @classmethod
    def htmlClusterTitle(cls, clusterCase, htmlCore):
        htmlCore.header(clusterCase)

    @classmethod
    def htmlClusterSubtext(cls, distanceMeasure, distList, linkageCriterion, htmlCore):

        if distanceMeasure == cls.ALL_MEASURES:
            htmlCore.line('Clustering is done for all available correlation measures, to facilitate '
                          'comparison of the different measures on continuous data representations: '
                          '<ul>' + ' '.join(['<li>' + dist + '</li>' for dist in distList]) +
                          '</ul> ')

        else:
            htmlCore.line('Chosen correlation measure for '
                          'clustering on binary data representations: ' + distanceMeasure)

        htmlCore.line('Linkage criterion for clustering dendrograms is ' + linkageCriterion)