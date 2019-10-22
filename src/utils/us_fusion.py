"""
This file was created by the CMU DELPHI Epi-forecasting group.
https://github.com/cmu-delphi/

===============
=== Purpose ===
===============

Prepares the variables necessary for performing sensor fusion of signals from
the various regions and states within the US. This includes, in particular, the
matrices H and W.

H maps from state space (columns) to input space (rows). W maps from state
space (columns) to output space (rows).

In general, this file takes as input the locations for which sensor readings
are available and returns as output H, W, and a list of locations which make up
output space.
"""

# standard library
from fractions import Fraction
import functools

# third party
import numpy as np

# first party
from src.utils.geo.locations import Locations
from src.utils.geo.populations import get_population


class UsFusion:
    """Prepares for sensor fusion of signals based on US regions and states."""

    __known_statespace = {}

    @staticmethod
    def get_weight_row(location, season, atoms):
        """
        Return a list of the population weights of all atoms, with respect to the
        given location. Atoms not within the location will have a weight of zero.
        The returned weights will sum to one.
        """

        # check each atom individually
        total_population = 0
        atom_populations = []
        for atom in atoms:
            if atom in Locations.region_map[location]:
                if season:
                    population = get_population(atom, season)
                else:
                    population = get_population(atom)
            else:
                population = 0
            total_population += population
            atom_populations.append(population)

        # sanity check
        if total_population == 0:
            raise Exception(('location has no constituent atoms', location))

        # return list of fractional populations
        get_fraction = lambda pop: Fraction(pop, total_population)
        return list(map(get_fraction, atom_populations))

    @staticmethod
    def get_weight_matrix(locations, season, atoms):
        """
        Return a matrix of weights, where rows correspond to the given locations
        and columns correspond to the given atomic locations.
        """

        # stack rows for each location and return the matrix
        get_row = lambda loc: UsFusion.get_weight_row(loc, season, atoms)
        return np.array(list(map(get_row, locations)))

    @staticmethod
    def eliminate(X):
        """
        Compute the canonical reduced row echelon form of the given matrix. The
        Gauss-Jordan algorithm is used to compute the elimination. The matrix is
        modified in-place.

        For numerical stability, it is strongly suggested that the elements of the
        input matrix be Fractions. Although discouraged, matrices of floats are also
        supported.

        input:
          X: the input matrix

        output:
          the matrix in reduced row echelon form
        """

        # dimensions
        num_r, num_c = X.shape

        # forward elimination
        r, c = 0, 0
        while r < num_r and c < num_c:
            values = [float(x) for x in X[r:, c]]
            i = r + np.argmax(np.abs(values))
            if X[i, c] != 0:
                if i != r:
                    temp = X[i, :].copy()
                    X[i, :] = X[r, :]
                    X[r, :] = temp
                X[r, c:] /= X[r, c]
                for i in range(r + 1, num_r):
                    X[i, c:] -= X[i, c] * X[r, c:]
                r += 1
            c += 1

        # backward substitution
        for r in range(num_r - 1, -1, -1):
            for c in range(num_c):
                if X[r, c] != 0:
                    for i in range(r - 1, -1, -1):
                        X[i, c:] -= X[i, c] * X[r, c:]
                    break

        # return the result
        return X

    @staticmethod
    def matmul(*matrices):
        """
        Compute the product of the given matrices. The matrices must all have
        elements of type Fraction or float. The type of the output will be the same
        as the type of the input.

        This function is not particularly efficient -- O(n^3) -- and is intended only
        for computing the product of matrices of fractions. The product of matrices
        of floats can be computed more efficiently by numpy or scipy.

        input:
          *matrices: the input matrices

        output:
          the product of inputs matrices
        """

        if len(matrices) == 1:
            return matrices[0]
        elif len(matrices) == 2:
            A, B = matrices
            (rows, size), (temp, cols) = A.shape, B.shape
            if size != temp:
                raise Exception('matrix dimensions do not match')
            dot = lambda U, V: sum(u * v for (u, v) in zip(U, V))
            vals = [[dot(A[r, :], B[:, c]) for c in range(cols)] for r in
                    range(rows)]
            return np.array(vals)
        else:
            return UsFusion.matmul(matrices[0], UsFusion.matmul(*matrices[1:]))

    @staticmethod
    def _determine_statespace(H0, W0):
        """
        Return matrices mapping from latent statespace to input space and output
        space. These are the matrices H and W, respectively, used in the sensor
        fusion kernel. Since some outputs may be indeterminate, the indices of the
        fully determined rows are returned. This may be used, for example, to find
        the set of outputs which make up the rows of the returned W matrix.

        inputs:
          H0: map from full statespace to inputs (I x S)
          W0: map from full statespace to outputs (O x S)

        outputs:
          - the matrix H, mapping subspace to inputs (I x S')
          - the matrix W, mapping subspace to outputs (O' x S')
          - list of row indices of W0 that make up W (O')

        notes:
          - S' <= S and O' <= O
          - for numerical stability, inputs should be matrices of Fractions
        """

        # helper function to convert a float matrix into a fraction matrix
        fractions = lambda X: np.array(
            [[Fraction(x) for x in row] for row in X])

        # Find a set of basis vectors that span the same subspace (of the full
        # statespace) that is spanned by the input vectors in H0. The result is a
        # minimal set of elements from which all inputs can be unambiguously
        # determined.
        B = UsFusion.eliminate(H0.copy())

        # the dimensions of full statespace (number of columns)
        size = B.shape[1]

        # the dimensions of the subspace (number of non-empty rows)
        rank = np.sum(np.sum(np.abs(B), axis=1) > 0)

        # B should be a square matrix with rows of zeros below rows of basis vectors
        num_rows = B.shape[0]
        if num_rows < size:
            Z = fractions(np.zeros((size - num_rows, size)))
            B = np.vstack((B, Z))
        elif num_rows > size:
            B = B[:size, :]

        # Attempt to build each input and output vector as a linear combination of
        # the subspace basis vectors. Since B may not be full rank, it may not be
        # invertible. Instead, solve by eliminating the augmented matrix of B
        # (transposed) with the identity matrix. After elimination, the (transposed)
        # inverse of B is contained within the augmented matrix.
        I = fractions(np.eye(size))
        BtI = np.hstack((B.T, I))
        IBit = UsFusion.eliminate(BtI)
        Bi = IBit[:, size:].T

        # possible, or "actual", solutions are in the leftmost columns
        # impossible, or "pseudo", solutions are in the rightmost columns
        Bi_actual, Bi_pseudo = Bi[:, :rank], Bi[:, rank:]

        # compute H, the map from statespace B to inputs
        # all inputs are within the span of statespace B
        H = UsFusion.matmul(H0, Bi_actual)

        # compute W, the map from statespace B to outputs
        # outputs not within the span of statespace B must be excluded
        W_actual = UsFusion.matmul(W0, Bi_actual)
        W_pseudo = UsFusion.matmul(W0, Bi_pseudo)

        # only keep rows where the coeficient of all pseudo basis vectors is zero
        actual_rows = np.flatnonzero(np.sum(np.abs(W_pseudo), axis=1) == 0)
        W = W_actual[actual_rows, :]

        # return H, W, and the indices of the rows of W0 that make up W
        return H, W, actual_rows

    @staticmethod
    @functools.lru_cache(maxsize=16)
    def determine_statespace(
            input_locations,
            season=None,
            exclude_locations=()):
        """
        Return matrices mapping from latent statespace to input space and output
        space. These are the matrices H and W, respectively, used in the sensor
        fusion kernel. A list of output locations corresponding to the rows of W is
        also returned.

        Results are cached for better performance.

        inputs:
          input_locations: a tuple of sensor locations
          season (optional): The season (year) in which the nowcast is being made.
            This is generally only helpful for retrospective nowcasts where
            historical population weights should be used. By default, the most
            recent population weights are used. (See populations.py)
          exclude_locations (optional): A tuple of atoms to exclude from
            statespace. This is generally only helpful for retrospective nowcasts
            where it is known that some state or territory was not reporting and,
            therefore, was not included in regional or national wILI.

        outputs:
          - the matrix H, mapping subspace to inputs
          - the matrix W, mapping subspace to outputs
          - tuple of output locations, corresponding to rows of W
        """

        # quick sanity check
        if set(exclude_locations) & set(input_locations):
            raise Exception('input contains excluded locations')

        # function to filter out excluded atoms
        atom_filter = lambda a: a not in exclude_locations

        # list of all locations, including nat, hhs, cen, and atoms
        all_locations = list(filter(atom_filter, Locations.region_list))

        # list of atomic locations only
        atoms = list(filter(atom_filter, Locations.atom_list))

        # precursors of the H and W matrices, assuming that statespace is US atoms
        get_matrix = lambda locs: UsFusion.get_weight_matrix(locs, season,
                                                             atoms)
        H0 = get_matrix(input_locations)
        W0 = get_matrix(all_locations)

        # optimization for the typical case where all US atoms are represented
        if set(input_locations) >= set(atoms):
            # statespace is all US atoms, so H and W are already correct
            H, W, output_locations = H0, W0, all_locations
        else:
            # determine optimal H and W matrices
            H, W, selected_rows = UsFusion._determine_statespace(H0, W0)
            # select the output locations
            output_locations = [all_locations[i] for i in selected_rows]

        # convert fractions to floats and return the result
        return H.astype(np.float), W.astype(np.float), output_locations
