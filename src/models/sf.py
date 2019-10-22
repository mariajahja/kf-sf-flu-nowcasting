"""
Contains sensor fusion methods.

We use the python interface to the Gurobi solver. To speed up computation time
for a grid of parameters, we reset the Gurobi model to an unsolved state (rather
than reconstruct a new model).
"""

# standard
import logging

import numpy as np
# third party
from gurobipy import *


def sf_l2(X, Z, H, lams):
    """Fit ridge regression constrained to be equivalent to sensor fusion.

    Args:
        X: matrix of historical wILI data (weeks x states)
        Z: matrix of sensor data (weeks x sensors)
        H: matrix of population weights (sensors x states)
        lams: Array of positive penalty parameters

    Returns:
        Beta: matrix of solutions to the minimization problem
    """
    t, k = X.shape  # weeks x states
    d = Z.shape[1]  # num sensors
    Ht = H.T
    Beta = np.zeros((len(lams), d, k)) * np.nan

    # Set up Gurobi model
    m = Model()
    m.setParam('OutputFlag', False)
    Bj = [m.addVar(lb=-GRB.INFINITY) for _ in range(d)]  # Beta_j
    G = [m.addVar(lb=-GRB.INFINITY) for _ in range(t)]  # G = X_j-ZBeta_j
    m.update()

    # Add constraints
    # HtBj = ej
    m.addConstrs((quicksum([Ht[r, j] * Bj[j] for j in range(d)]) == 0 \
                  for r in range(k)), name='HtB')
    G_constrs = None  # Updated in loop
    m.update()

    for region in range(k):
        Xj = X[:, region]  # n_weeks by 1

        # reset HtB constraint to 1 for that region
        m.remove(m.getConstrByName(f'HtB[{region}]'))
        m.addConstr((quicksum([Ht[region, j] * Bj[j] for j in range(d)]) == 1),
                    name=f'HtB[{region}]')
        m.update()

        # reset G eq
        if G_constrs:
            m.remove(G_constrs)
            m.update()

        G_constrs = m.addConstrs(
            (Xj[i] - quicksum(Z[i, j] * Bj[j] for j in range(d)) == G[i] for i
             in range(t)))
        m.update()

        if not len(m.getConstrs()) == t + k:
            logging.error("Error in removing constraints for ridge regression")
            raise Exception

        for q, lam in enumerate(lams):
            # Set objective
            m.setObjective(quicksum([G[i] * G[i] for i in range(t)]) + \
                           lam * quicksum([Bj[j] * Bj[j] for j in range(d)]),
                           GRB.MINIMIZE)
            m.update()
            m.optimize()

            Beta[q, :, region] = [Bj[i].X for i in range(d)]
            m.reset()  # reset model to unsolved state

        # reset HtB constraint to 0
        m.remove(m.getConstrByName(f'HtB[{region}]'))
        m.addConstr((quicksum([Ht[region, j] * Bj[j] for j in range(d)]) == 0),
                    name=f'HtB[{region}]')
        m.update()

    for q, lam in enumerate(lams):
        if not np.allclose(np.dot(Ht, Beta[q, :, :]), np.eye(k)):  # identity
            logging.error(f'Const. Ridge constraints not met.'
                          f' Solution is bad for lam {lam}.')
        else:
            logging.info(f'Solved const. ridge problem for lam {lam}.')

    return Beta


def sf_l1(X, Z, H, lams):
    """Fit lasso regression with sensor fusion constraints.

        Args:
            X: matrix of historical wILI data (weeks x states)
            Z: matrix of sensor data (weeks x sensors)
            H: matrix of population weights (sensors x states)
            lams: Array of positive penalty parameters

        Returns:
            Beta: matrix of solutions to the minimization problem
        """
    t, k = X.shape  # weeks x states
    d = Z.shape[1]  # num sensors
    Ht = H.T
    Beta = np.zeros((len(lams), d, k)) * np.nan

    # Set up Gurobi model
    m = Model()
    m.setParam('OutputFlag', False)
    Bj = [m.addVar(lb=-GRB.INFINITY) for _ in range(d)]  # Beta_j
    aBj = [m.addVar(lb=0) for _ in range(d)]  # lasso variable
    G = [m.addVar(lb=-GRB.INFINITY) for _ in range(t)]  # G = X_j-ZBeta_j
    m.update()

    # HtB = e constraint
    m.addConstrs((quicksum([Ht[r, j] * Bj[j] for j in range(d)]) == 0 \
                  for r in range(k)), name='HtB')
    G_constrs = None  # updated in loop
    m.update()

    # Set lasso objective
    for j in range(d): m.addConstr(aBj[j] == abs_(Bj[j]))
    m.update()

    for region in range(k):
        Xj = X[:, region]  # n_weeks by 1

        # reset HtB constraint to 1 for that region
        m.remove(m.getConstrByName(f'HtB[{region}]'))
        m.addConstr((quicksum([Ht[region, j] * Bj[j] for j in range(d)]) == 1),
                    name=f'HtB[{region}]')
        m.update()

        # reset G eq
        if G_constrs:
            m.remove(G_constrs)
            m.update()

        G_constrs = m.addConstrs(
            (Xj[i] - quicksum(Z[i, j] * Bj[j] for j in range(d)) == G[i] for i
             in range(t)))
        m.update()

        for q, lam in enumerate(lams):
            m.setObjective((1 / (2 * t)) * quicksum(
                [G[i] * G[i] for i in range(t)]) + lam * quicksum(
                [aBj[j] for j in range(d)]), GRB.MINIMIZE)
            m.update()
            m.optimize()

            Beta[q, :, region] = [Bj[i].X for i in range(d)]
            m.reset()  # reset model to unsolved state

        # reset HtB constraint to 0
        m.remove(m.getConstrByName(f'HtB[{region}]'))
        m.addConstr(
            (quicksum([Ht[region, j] * Bj[j] for j in range(d)]) == 0),
            name=f'HtB[{region}]')
        m.update()

    for q, lam in enumerate(lams):
        if not np.allclose(np.dot(Ht, Beta[q, :, :]), np.eye(k)):  # identity
            logging.error(f'Const. lasso constraints not met.'
                          f' Solution is bad for lam {lam}.')
        else:
            logging.info(f'Solved const. lasso problem for lam {lam}.')

    return Beta


if __name__ == '__main__':
    H = np.random.randn(100, 5)
    X = np.random.randn(1000, 5)
    Z = X @ H.T + np.random.randn(1000, 100)

    ridge_Beta = sf_l2(X, Z, H, [0, 1, 10, 20])
    lasso_Beta = sf_l1(X, Z, H, [0, 1, 10, 20])

    G = Z - X @ H.T
    cov_G = G.T @ G / X.shape[0]
    Ri = np.linalg.inv(cov_G)
    kf_Beta = np.linalg.inv(H.T @ Ri @ H) @ H.T @ Ri

    assert np.allclose(ridge_Beta[0], lasso_Beta[0])
    assert np.allclose(kf_Beta.T, ridge_Beta[0])

    # test lasso - will need to remove H.TB=I constraint in func for equality
    # from sklearn.linear_model import Lasso
    #
    # lasso_Beta = sf_l1(X[:, 1].reshape(-1, 1), Z, H, [0.25])
    # lasso_Beta = lasso_Beta[0]
    # lasso_Beta[np.abs(lasso_Beta) < 1e-6] = 0.0
    #
    # mod = Lasso(alpha=0.25, fit_intercept=False)
    # mod.fit(Z, X[:, 1])
    # print(np.sum(np.abs(mod.coef_ - lasso_Beta.T)))
