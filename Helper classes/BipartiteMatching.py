from quick.webtools.clustering.LDExpansions import LDExpansions


class BipartiteMatching(object):
    """
    Functions for bipartite matching of tracks. The cost matrix, is a len(track1) x len(track2) matrix, where the
    cells are weights between the nodes in the tracks.

    Both greedyBipartite and lapjvBipartite returns a dictionary containing scores for the bipartite matching scores.
    For each match that contribute to the overall matching score, the following is updated:
    a += Bipartite matching score
    b += (1 - a) / 2
    c += (1 - a) / 2
    d = undefined (-1)

    A dictionary of the a, b, c values is returned. Count of SNPs with no edges, i.e no chance of a match score,
    in each track is added to b and c.

    The greedy algorithm is susceptible to local optima, and cannot guarantee that the overall score is the maximal
    score possible.

    One can use generateCostMatrix and the functions in LDExpansions.py to generate a cost matrix based on a linked
    point track. In LDExpansions, a graph is created from a linked point track, and the function assume that
    the weights are rsquare values of LD correlation. However, any graph representation in a linked point track,
    which is symmetric, with edges and weights, can use these functions for cost matrix generation.

    For example usage, see the LDBipartiteMatchingTool.
    """

    @classmethod
    def greedyBipartite(cls, cost_matrix):
        """
        Greedy algorithm for bipartite matching of two tracks.
        Takes in a len(track1) x len(track2) matrix, where the cells are weights between the nodes in
        the tracks.

        Finds a, b and c for the tracks, as defined above.
        """
        from numpy import zeros, shape, max, argmax, unravel_index, matrix

        dimensions = shape(cost_matrix)
        cost = {'a': 0.0, 'b': 0.0, 'c': 0.0, 'd': -1}

        cost['b'] += cls.getZeroOccurrences(cost_matrix)
        cost['c'] += cls.getZeroOccurrences(matrix.transpose(cost_matrix))

        while max(cost_matrix) > 0:
            pos = argmax(cost_matrix)
            row, col = unravel_index([pos], dimensions)
            dimensions = list(dimensions)

            val = max(cost_matrix)
            cost_matrix[row[0], :] = zeros(dimensions[1])
            cost_matrix[:, col[0]] = zeros(dimensions[0])
            cost['a'] += val

        return cost

    @classmethod
    def lapjvBipartite(cls, cost_matrix):
        """
        Optimal algorithm for bipartite matching of two tracks.
        Takes in len(track1) x len(track2) matrix, where the cells are weights between the nodes in
        the tracks.

        Finds a, b and c for the tracks, as defined above.
        """
        from quick.webtools.clustering.JonkerVolgenant import JonkerVolgenant
        from numpy import matrix
        matches = JonkerVolgenant.findJonkerVolgenant(cost_matrix)

        cost = {'a': 0, 'b': 0, 'c': 0, 'd': -1}

        track1Length = len(cost_matrix)
        track2Length = len(cost_matrix[0])

        cost['b'] += cls.getZeroOccurrences(cost_matrix)
        cost['c'] += cls.getZeroOccurrences(matrix.transpose(cost_matrix))

        if track1Length < track2Length:
            for matchID in range(0, track1Length):
                val = cost_matrix[matchID, matches[matchID]]
                cost['a'] += val
        else:
            for matchID in range(0, track2Length):
                val = cost_matrix[matches[matchID], matchID]
                cost['a'] += val

        return cost

    @classmethod
    def generateCostMatrix(cls, track1, track2, graph):
        """
        Take in an LD graph, as is generated from the LDExpansions.createRSquareGraph function, and two tracks,
        represented as lists of the track SNP rsids.
        Return a cost matrix for the two tracks, where cells of value mark edges with specific weights between
        the nodes. Weights are rsquare values that weights edges between nodes (SNPs).
        """
        from numpy import zeros

        dim1 = len(track1)
        dim2 = len(track2)

        cost_matrix = zeros((dim1, dim2))
        i, j = (0, 0)
        while i < dim1:
            while j < dim2:
                r2 = LDExpansions.getEdge(track1[i], track2[j], graph)
                cost_matrix[i, j] = r2
                j += 1
            i += 1
            j = 0

        return cost_matrix

    @classmethod
    def getZeroOccurrences(cls, matrix):
        """
        Takes in cost matrix, and returns the number of rows with only zeros in them.
        """
        zeroCount = 0
        for row in matrix:
            if all([x == 0 for x in row]):
                zeroCount += 1

        return zeroCount
