from scipy.io import loadmat
from dataclasses import dataclass
import numpy as np
import cvxpy as cp
from numpy.linalg import matrix_power
import matplotlib.pyplot as plt

from plotting import *


global n_agents
global nx
global nu
# global ny
n_agents = 4
nx = 2
nu = 1
# ny = 2


@dataclass
class Agent:
    A: np.ndarray
    B: np.ndarray
    x0: np.ndarray
    M: np.ndarray
    b: np.ndarray
    h: np.ndarray
    idx: int


@dataclass
class Subproblem:
    pr: cp.Problem
    xfinal: cp.Expression
    f: cp.Expression
    constraint: cp.Expression
    augment: cp.Expression
    u: cp.Variable
    lambdas: cp.Parameter
    xf: cp.Parameter
    rho: cp.Parameter

    def dual_update(self, xf_new, rho):
        self.xf.value = xf_new.reshape((2, 1))
        self.lambdas.value = (
            self.lambdas.value + rho * self.xfinal.value - rho * xf_new.reshape((2, 1))
        )
        return self.lambdas.value

    def primal_update(self):
        self.pr.solve()
        return self.u.value, self.xfinal.value, self.lambdas.value


def gen_M_b(A, B, x0):
    nt = int(Tfinal)
    M = np.zeros((2 * nt, nt))
    b = np.zeros((2 * nt))

    powers = np.zeros((nt + 1, 2, 2))
    for i in range(nt + 1):
        powers[i, :, :] = matrix_power(A, i)

    # b matrix
    for i in range(nt):
        idx = 2 * i
        b[idx : idx + 2] = (powers[i + 1, :, :] @ x0).flatten()

    # A matrix
    for row in range(nt):
        for col in range(row + 1):
            rowidx = 2 * row
            M[rowidx : rowidx + 2, col] = (powers[row - col, :, :] @ B).flatten()
    return M, b


def h(idx):
    # h for hardcoding
    I = np.eye(2)
    zero = np.zeros([2, 2])
    match idx:
        case 0:
            # return np.array([I, zero, zero, -I])
            return np.array([I, zero, zero, -I]).reshape(8, 2)
        case 1:
            return np.array([-I, I, zero, zero]).reshape(8, 2)
        case 2:
            return np.array([zero, -I, I, zero]).reshape(8, 2)
        case 3:
            return np.array([zero, zero, -I, I]).reshape(8, 2)


def load_variables():
    a = loadmat("agents.mat")
    global Tfinal
    Tfinal = a["Tfinal"][0, 0]

    M, b = gen_M_b(a["A1"], a["B1"], a["x01"])
    agent_1 = Agent(a["A1"], a["B1"], a["x01"], M, b, h(0), 0)

    M, b = gen_M_b(a["A2"], a["B2"], a["x02"])
    agent_2 = Agent(a["A2"], a["B2"], a["x02"], M, b, h(1), 1)

    M, b = gen_M_b(a["A3"], a["B3"], a["x03"])
    agent_3 = Agent(a["A3"], a["B3"], a["x03"], M, b, h(2), 2)

    M, b = gen_M_b(a["A4"], a["B4"], a["x04"])
    agent_4 = Agent(a["A4"], a["B4"], a["x04"], M, b, h(3), 3)

    return [agent_1, agent_2, agent_3, agent_4], a["umax"][0, 0], a["Tfinal"][0, 0]


def sub_problem(agent: Agent, rho) -> cp.Problem:
    nt = int(Tfinal)
    u = cp.Variable((Tfinal, 1))
    lambdas = cp.Parameter((2, 1))
    xf = cp.Parameter((2, 1))
    # rho = cp.Parameter(pos=True, value=1)
    # xf_var = cp.Variable((2, 1))

    xfinal = agent.M[2 * nt - 2 :, :] @ u + agent.b[2 * nt - 2 :].reshape(2, 1)
    Mf = agent.M[2 * nt - 2 :, :]
    bf = agent.b[2 * nt - 2 :].reshape(2, 1)
    f = (
        cp.quad_form(u, (agent.M.T @ agent.M + np.eye(20)))
        + 2 * (agent.b @ agent.M) @ u
    )
    # f = (
    #     cp.quad_form(
    #         u,
    #         (agent.M.T @ agent.M + np.eye(20) + (rho) * 0.25 * Mf.T @ Mf),
    #         assume_PSD=True,
    #     )
    #     + 2 * (agent.b @ agent.M + rho / 2 * (bf - xf).T @ Mf) @ u
    # )
    augment = rho / 2 * cp.sum_squares(xfinal - xf)

    # augment = rho / 2 * (cp.quad_form(u, Mf.T @ Mf) + 2 * (bf - xf).T @ Mf @ u)

    constraint = xfinal - xf
    cost = f + lambdas.T @ constraint + augment
    # augment = 1
    return Subproblem(
        cp.Problem(cp.Minimize(cost), [-u <= u_max, u <= u_max]),
        xfinal,
        f,
        constraint,
        augment,
        u,
        lambdas,
        xf,
        rho,
    )


def admm_optimization(
    subproblems,
    agents,
    u_cent: np.ndarray,
    tol=1e-16,
    n_iter=400,
    rho=1.0,
    lambdas=None,
):
    if lambdas is None:
        lambdas = np.zeros([8, 1])
    nt = int(Tfinal)

    lambdas = lambdas
    iterations = []
    xfinal_hist = np.zeros([n_iter, 2, 4])
    e_hist = np.zeros([n_iter])

    # initialize with optimization parameters
    for problem in subproblems:
        # problem.rho.value = rho
        problem.lambdas.value = np.array([[0], [0]], np.float64)
        problem.xf.value = np.array([[0], [0]], np.float64)

    for i in range(n_iter):
        # primal update
        u_calc = []
        lambdas = []
        xfinals = []

        for subproblem in subproblems:
            u_sub, x_sub, lambdas_i = subproblem.primal_update()
            u_calc.append(u_sub)
            xfinals.append(x_sub)
            lambdas.append(lambdas_i)

        # consensus update
        xfinal_new = 0.25 * sum(xfinals) + 0.25 * rho * sum(lambdas_i)
        iterations.append(i)
        xfinal_hist[i, :, :] = np.array(xfinals).T

        # dual update
        for problem in subproblems:
            problem.dual_update(xfinal_new, rho)

        error_norm = np.linalg.norm(np.array(u_calc) - u_cent)
        e_hist[i] = error_norm
        # if error_norm < tol:
        #     break

    print(f"optimization finished, {iterations[-1]+1} update steps taken")
    u_conv = [problem.u.value for problem in subproblems]
    x_conv = [
        ((agent.M @ u_res)[:, 0] + agent.b).reshape(nt, 2)
        for agent, u_res in zip(agents, u_conv)
    ]
    return (
        iterations,
        xfinal_hist[: i + 1, :, :],
        e_hist[: i + 1],
        u_conv,
        x_conv,
    )


def get_centralized_variables(agent):
    nt = int(Tfinal)
    u = cp.Variable((Tfinal, 1))
    xfinal = agent.M[2 * nt - 2 :, :] @ u + agent.b[2 * nt - 2 :].reshape(2, 1)
    f = cp.quad_form(u, (agent.M.T @ agent.M + np.eye(20))) + 2 * agent.b @ agent.M @ u
    return u, xfinal, f


def centralized_solution(agents):
    nt = int(Tfinal)

    u1, xfinal1, f1 = get_centralized_variables(agents[0])
    u2, xfinal2, f2 = get_centralized_variables(agents[1])
    u3, xfinal3, f3 = get_centralized_variables(agents[2])
    u4, xfinal4, f4 = get_centralized_variables(agents[3])

    cost = f1 + f2 + f3 + f4

    constraints = [
        -u1 <= u_max,
        u1 <= u_max,
        -u2 <= u_max,
        u2 <= u_max,
        -u3 <= u_max,
        u3 <= u_max,
        -u4 <= u_max,
        u4 <= u_max,
        xfinal1 == xfinal2,
        xfinal2 == xfinal3,
        xfinal3 == xfinal4,
        xfinal4 == xfinal1,
    ]

    pr = cp.Problem(cp.Minimize(cost), constraints)
    pr.solve()

    u_result = [u1.value, u2.value, u3.value, u4.value]
    x_result = [
        ((agent.M @ u_res)[:, 0] + agent.b).reshape(nt, 2)
        for agent, u_res in zip(agents, u_result)
    ]

    return u_result, x_result


def admm_sol(agents, rho):
    subproblems = [sub_problem(agent, rho) for agent in agents]
    lambdas = np.zeros([8, 1])

    u_cent, x_cent = centralized_solution(agents)
    return (
        admm_optimization(
            subproblems, agents, np.array(u_cent), tol=1e-10, n_iter=150, rho=rho
        ),
        u_cent,
        x_cent,
    )


def rho_effect(agents):
    rhos = [0.5, 1, 1.5, 1.9]
    results = [
        admm_sol(
            agents,
            rho=rho,
        )
        for rho in rhos
    ]

    # results_normal = subgradient_optimization(
    #     subproblems,
    #     agents,
    #     np.array(u_cent),
    #     tol=1e-10,
    #     n_iter=100,
    #     alpha=5.0,
    #     steps="constant",
    # )

    iterations, e_hist, x_conv, x_cent = [
        [res[i] for res in results] for i in (0, 2, 4, 6)
    ]

    admm_rho_effect(
        iterations,
        e_hist,
        rhos,
    )


if __name__ == "__main__":
    global u_max
    agents, u_max, Tfinal = load_variables()

    rho = 0.5
    # subproblems = [sub_problem(agent) for agent in agents]
    # lambdas = np.zeros([8, 1])

    u_cent, x_cent = centralized_solution(agents)
    # iterations, xfinal_hist, e_hist, u_conv, x_conv = admm_optimization(
    #     subproblems, agents, np.array(u_cent), tol=1e-10, n_iter=300, rho=0.5
    # )

    [iterations, xfinal_hist, e_hist, u_conv, x_conv], u_cent, x_cent = admm_sol(
        agents, rho=1
    )
    plt.figure()
    convergence_error(iterations, e_hist)
    plt.figure()
    compare_cent_decentr(u_cent, x_cent, u_conv, x_conv)
    plt.figure()
    final_state_convergence(iterations, xfinal_hist)
    # plt.figure()
    # rho_effect(agents)
    plt.figure()
    plot_agent_subplots(x_cent, x_conv)
    plt.show()
