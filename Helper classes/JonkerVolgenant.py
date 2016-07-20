
class JonkerVolgenant(object):
    """
    A python implementation of the Hungarian algorithm as described in Jonker and Volgenants article "A Shortest
    Augmenting Path Algorithm for Dense and Sparse Linear Assignment Problems" (1986)

    It is an optimal bipartite matcher, which finds the highest total weight between two sets of nodes X and Y.
    Nodes in X and Y can be linked together, as in a graph, with a weight between them. This is represented as a
    cost matrix, where the nodes represent rows and columns, and the cell contains the weight of the edge between them.

    One area of usage is to match tracks of SNPs in LD.

    The Java coded version of the Jonker Volgenant algorithm presented in the git repo below was used to a great extent
    for inspiration and help when translating the original code from Pascal:
    https://github.com/dberm22/Bipartite-Solver
    """

    @staticmethod
    def findJonkerVolgenant(cost_matrix):

        from numpy import zeros, min, argmin

        rdim = len(cost_matrix)
        cdim = len(cost_matrix[0])

        # Pad matrix to ensure it is square
        cost_matrix = JonkerVolgenant.pad_matrix(cost_matrix)
        dim = len(cost_matrix)
        v = zeros(dim)

        row_solution = zeros(dim)
        col_solution = zeros(dim)
        free = zeros(dim)
        numfree = 0
        minimum = 0
        last = 0
        index = 0

        matches = zeros(dim)     # Init how many times a row will be assigned in the column reduction

        # Column reduction:
        for j in range(dim - 1, -1, -1):    # Reverse order gives better results
            minimum = min(cost_matrix[:, j])
            i_minimum = argmin(cost_matrix[:, j])
            v[j] = minimum

            matches[i_minimum] += 1
            if matches[i_minimum] == 1:  # Init assignment if min row assigned for first time
                row_solution[i_minimum] = j
                col_solution[j] = i_minimum
            else:
                col_solution[j] = -1    # Row already assigned, column not assigned

        # Reduction transfer
        for index in range(0, dim):
            if matches[index] == 0:     # Fill list of unassigned "free" rows
                free[numfree] = index
                numfree += 1
            elif matches[index] == 1:   # Transfer reduction from rows that are assigned once
                j1 = row_solution[index]
                minimum = 2        # higher than max possible cost

                for j in range(0, dim):
                    tmp_cost = cost_matrix[index, j] - v[j]
                    if j != j1 and tmp_cost < minimum:
                        minimum = tmp_cost
                v[j1] -= minimum

        # Augmenting row reduction
        for loopcount in range(0, 2):
            k = 0
            prev_numfree = numfree
            numfree = 0     # start list of rows still free after augmenting row reduction

            while k < prev_numfree:
                index = free[k]
                k += 1

                # find minimum and second minimum reduced cost over columns
                umin = cost_matrix[index][0] - v[0]
                j1 = 0
                usubmin = 2

                for j in range(1, dim):
                    h = cost_matrix[index, j] - v[j]
                    if usubmin > h >= umin:
                        usubmin = h
                        j2 = j
                    elif h < usubmin and h < umin:
                        usubmin = umin
                        umin = h
                        j2 = j1
                        j1 = j

                i0 = col_solution[j1]
                if umin < usubmin:
                    v[j1] -= (usubmin - umin)
                elif i0 >= 0:
                    j1 = j2     # Swap columns j1 and j2, as j2 may be unassigned
                    i0 = col_solution[j2]

                # (Re)assign i to j1, possible de-assigning an i0
                row_solution[index] = j1
                col_solution[j1] = index

                # Put in current k and go back to that k
                # Continue augmenting path i - j1 with i0
                if i0 >= 0 and umin < usubmin:
                    k -= 1
                    free[k] = i0
                # No further augmenting reduction possible
                # Store i0 in list of free rows for next phase
                elif i0 >= 0:
                    free[numfree] = i0
                    numfree += 1

        # Augment solution for each free row
        d = zeros(dim)
        pred = zeros(dim)
        col_list = zeros(dim)
        endofpath = 0
        for f in range(0, numfree):
            freerow = free[f]

            # Dijkstra shortest path algorithm
            # Runs until unassigned column added to shortest path tree
            for j in range(0, dim):
                d[j] = cost_matrix[freerow, j] - v[j]
                pred[j] = freerow
                col_list[j] = j

            low = 0     # Columns in [0 ... low - 1] which are ready
            up = 0      # Columns in [low ... up - 1] which are to be scanned for current minimum
            unassignedfound = False

            while not unassignedfound:
                if up == low:
                    last = low - 1

                    # Scan columns for [up ... dim - 1] to find all indices for which new minimum occurs
                    # Store these indices between [low .. up - 1] (increasing up)
                    minimum = d[col_list[up]]
                    up += 1

                    for k in range(up, dim):
                        j = col_list[k]
                        h = d[j]
                        if h <= minimum:

                            if h < min:     # New minimum
                                up = low
                                minimum = h

                            # New index with same minimum, put on index up and extend list
                            col_list[k] = col_list[up]
                            col_list[up] = j
                            up += 1

                    # Check if any of the minimum columns happens to be unassigned
                    # If so, we have an augmenting path
                    for k in range(low, up):
                        if col_solution[col_list[k]] < 0:
                            endofpath = col_list[k]
                            unassignedfound = True
                            break

                if not unassignedfound:
                    j1 = col_list[low]
                    low += 1
                    index = col_solution[j1]
                    h = cost_matrix[index, j1] - v[j1] - minimum

                    for k in range(up, dim):
                        j = col_list[k]
                        v2 = cost_matrix[index, j] - v[j] - h
                        if v2 < d[j]:

                            pred[j] = index
                            if v2 == minimum:
                                if col_solution[j] < 0:     # If unassigned, shortest augmenting path is complete
                                    endofpath = j
                                    unassignedfound = True
                                    break
                                else:                       # Add to list to be scanned right away
                                    col_list[k] = col_list[up]
                                    col_list[up] = j
                                    up += 1
                            d[j] = v2

            # Update column prices
            for k in range(0, last + 1):
                j1 = col_list[k]
                v[j1] += d[j1] - minimum

            # Reset row and column assignments along the alternating path
            while index != freerow:
                index = pred[endofpath]
                col_solution[endofpath] = index
                j1 = endofpath
                endofpath = row_solution[index]
                row_solution[index] = j1

        # Trim and return
        if rdim < cdim:
            return row_solution[0:rdim]
        else:
            return col_solution[0:cdim]

    @staticmethod
    def pad_matrix(cost_matrix):
        from numpy import ones, where

        rdim = len(cost_matrix)
        cdim = len(cost_matrix[0])
        dim = max(rdim, cdim)

        new_matrix = ones((dim, dim))
        indices = where(cost_matrix > 0)

        for index in indices[0]:
            for j in indices[1]:
                new_matrix[index, j] = 1 - cost_matrix[index, j]

        return new_matrix
