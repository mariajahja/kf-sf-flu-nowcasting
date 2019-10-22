"""
Contains (unconstrained) regression methods.

We use the python interface to the Gurobi solver. To speed up computation time
for a grid of parameters, we reset the Gurobi model to an unsolved state (rather
than reconstruct a new model).
"""

# third party
import numpy as np
from gurobipy import *


def ridge(X, Z, H, lams):
    """Fit ridge regression.

    Args:
        X: matrix of historical wILI data (weeks x states)
        Z: matrix of sensor data (weeks x sensors)
        H: unused (this method is unconstrained)
        lams: Array of positive penalty parameters

    Returns:
        Beta: matrix of solutions to the minimization problem
    """
    t, k = X.shape  # weeks x states
    d = Z.shape[1]  # num sensors
    Beta = np.zeros((len(lams), d, k)) * np.nan

    # Set up Gurobi model
    m = Model()
    m.setParam('OutputFlag', False)
    Bj = [m.addVar(lb=-GRB.INFINITY) for _ in range(d)]  # Beta_j
    G = [m.addVar(lb=-GRB.INFINITY) for _ in range(t)]  # G = X_j-ZBeta_j
    m.update()

    G_constrs = None  # Updated in loop

    for region in range(k):
        Xj = X[:, region]  # n_weeks by 1

        # reset G eq
        if G_constrs:
            m.remove(G_constrs)
            m.update()

        G_constrs = m.addConstrs(
            (Xj[i] - quicksum(Z[i, j] * Bj[j] for j in range(d)) == G[i] for i
             in range(t)))
        m.update()

        for q, lam in enumerate(lams):
            m.setObjective(quicksum([G[i] * G[i] for i in range(t)]) + \
                           lam * quicksum([Bj[j] * Bj[j] for j in range(d)]),
                           GRB.MINIMIZE)
            m.update()
            m.optimize()

            Beta[q, :, region] = [Bj[i].X for i in range(d)]
            m.reset()  # reset model to unsolved state

    return Beta


def lasso(X, Z, H, lams):
    """Fit lasso regression.

        Args:
            X: matrix of historical wILI data (weeks x states)
            Z: matrix of sensor data (weeks x sensors)
            H: unused (this method is unconstrained)
            lams: Array of positive penalty parameters

        Returns:
            Beta: matrix of solutions to the minimization problem
        """
    t, k = X.shape  # weeks x states
    d = Z.shape[1]  # num sensors
    Beta = np.zeros((len(lams), d, k)) * np.nan

    # Set up Gurobi model
    m = Model()
    m.setParam('OutputFlag', False)
    Bj = [m.addVar(lb=-GRB.INFINITY) for _ in range(d)]  # Beta_j
    aBj = [m.addVar(lb=0) for _ in range(d)]  # lasso variable
    G = [m.addVar(lb=-GRB.INFINITY) for _ in range(t)]  # G = X_j-ZBeta_j
    m.update()

    G_constrs = None  # updated in loop
    for j in range(d): m.addConstr(aBj[j] == abs_(Bj[j]))
    m.update()

    for region in range(k):
        Xj = X[:, region]  # n_weeks by 1

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

    return Beta


if __name__ == '__main__':
    H = np.random.randn(100, 5)
    X = np.random.randn(1000, 5)
    Z = X @ H.T + np.random.randn(1000, 100)

    ridge_Beta = ridge(X[:, 1].reshape(-1, 1), Z, H, [0, 1, 10])
    lasso_Beta = lasso(X[:, 1].reshape(-1, 1), Z, H, [0, 0.5, 1])

    # with 0 regularization, solutions should be exactly the same
    assert np.allclose(ridge_Beta[0], lasso_Beta[0])

    # check against sklearn implementation
    from sklearn.linear_model import Ridge, Lasso

    mod = Ridge(alpha=1, fit_intercept=False)
    mod.fit(Z, X[:, 1])
    assert np.allclose(mod.coef_, ridge_Beta[1].T, atol=1e-4)

    mod = Lasso(alpha=0.5, fit_intercept=False)
    mod.fit(Z, X[:, 1])
    assert np.allclose(mod.coef_, lasso_Beta[1].T, atol=1e-4)
