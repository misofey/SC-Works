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
    cost: cp.Expression
    u: cp.Variable
    lambdas: cp.Parameter


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


def sub_problem(agent: Agent) -> cp.Problem:
    nt = int(Tfinal)
    u = cp.Variable((Tfinal, 1))
    lambdas = cp.Parameter((8, 1))

    xfinal = agent.M[2 * nt - 2 :, :] @ u + agent.b[2 * nt - 2 :].reshape(2, 1)
    f = cp.quad_form(u, (agent.M.T @ agent.M + np.eye(20))) + 2 * agent.b @ agent.M @ u
    constraint = agent.h @ xfinal
    cost = f + lambdas.T @ constraint

    return Subproblem(
        cp.Problem(cp.Minimize(cost), [-u <= u_max, u <= u_max]),
        xfinal,
        f,
        constraint,
        cost,
        u,
        lambdas,
    )


def get_subgradient(lambdas, subproblems, agents):
    gradient = np.zeros([8, 1])
    u_results = np.zeros([20, 4])
    xfinals = np.zeros([2, 4])
    for spr, agent, i in zip(subproblems, agents, range(4)):
        spr.lambdas.value = lambdas
        spr.pr.solve()
        u_calc = spr.u.value
        gradient += spr.constraint.value
        xfinals[
            :,
            i,
        ] = spr.xfinal.value.flatten()
    return gradient, xfinals


def subgradient_optimization(
    subproblems,
    agents,
    u_cent: np.ndarray,
    tol=1e-16,
    n_iter=400,
    alpha=1.0,
    lambdas=None,
    steps="constant",
):
    if lambdas is None:
        lambdas = np.zeros([8, 1])
    nt = int(Tfinal)

    lambdas = lambdas
    grad_hist = []
    iterations = []
    xfinal_hist = np.zeros([n_iter, 2, 4])
    e_hist = np.zeros([n_iter])

    for i in range(0, n_iter):
        grad_new, xfinals = get_subgradient(lambdas, subproblems, agents)
        norm = np.linalg.norm(grad_new)
        iterations.append(i)
        grad_hist.append(norm)
        xfinal_hist[i, :, :] = xfinals

        u_calc = np.array([problem.u.value for problem in subproblems])
        error_norm = np.linalg.norm(u_calc - u_cent)
        e_hist[i] = error_norm
        # if error_norm < tol:
        #     break
        match steps:
            case "constant":
                size = 1
            case "variable":
                size = 10 / (i + 1)

        lambdas = lambdas + alpha * size * (grad_new)

    print(f"optimization finished, {iterations[-1]} subgradient steps taken")
    u_conv = [problem.u.value for problem in subproblems]
    x_conv = [
        ((agent.M @ u_res)[:, 0] + agent.b).reshape(nt, 2)
        for agent, u_res in zip(agents, u_conv)
    ]
    return (
        iterations,
        grad_hist,
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


def accelerated_subgradient_optimization(
    subproblems,
    agents,
    u_cent: np.ndarray,
    tol=1e-16,
    n_iter=400,
    alpha=1.0,
    lambdas=None,
    steps="constant",
):
    if lambdas is None:
        lambdas = np.zeros([8, 1])
    nt = int(Tfinal)

    lambdas = lambdas
    lambda_prev = lambdas
    grad_hist = []
    iterations = []
    xfinal_hist = np.zeros([n_iter, 2, 4])
    e_hist = np.zeros([n_iter])
    beta = 0

    for i in range(0, n_iter):
        stuff = (1 + np.sqrt(1 + 4 * beta * beta)) / 2
        gammas = lambdas + (beta - 1) / stuff * (lambdas - lambda_prev)
        beta = stuff
        grad_new, xfinals = get_subgradient(gammas, subproblems, agents)
        norm = np.linalg.norm(grad_new)
        iterations.append(i)
        grad_hist.append(norm)
        xfinal_hist[i, :, :] = xfinals

        u_calc = np.array([problem.u.value for problem in subproblems])
        error_norm = np.linalg.norm(u_calc - u_cent)
        e_hist[i] = error_norm
        # if error_norm < tol:
        #     break
        match steps:
            case "constant":
                size = 1
            case "variable":
                size = 10 / (i + 1)

        lambdas = lambdas + alpha * size * (grad_new)

    print(f"optimization finished, {iterations[-1]} subgradient steps taken")
    u_conv = [problem.u.value for problem in subproblems]
    x_conv = [
        ((agent.M @ u_res)[:, 0] + agent.b).reshape(nt, 2)
        for agent, u_res in zip(agents, u_conv)
    ]
    return (
        iterations,
        grad_hist,
        xfinal_hist[: i + 1, :, :],
        e_hist[: i + 1],
        u_conv,
        x_conv,
    )


# def convergence_for_alpha(agents, subproblems, alpha):
#     iterations, grad_hist, xfinal_hist, e_hist, u_conv, x_conv = (
#         subgradient_optimization(
#             subproblems, agents, np.array(u_cent), tol=1e-10, n_iter=100, alpha=alpha
#         )
#     )
#     return iterations, e_hist


def problem_a(agent, subproblems, u_cent, x_cent):
    iterations, grad_hist, xfinal_hist, e_hist, u_conv, x_conv = (
        subgradient_optimization(
            subproblems, agents, np.array(u_cent), tol=1e-10, n_iter=100, alpha=1.0
        )
    )
    plt.figure()
    convergence_error(iterations, e_hist)
    plt.figure()
    plot_agent_subplots(x_cent, x_conv)
    plt.figure()
    final_state_convergence(iterations, xfinal_hist)


def problem_b(agents, subproblems, u_cent, x_cent):
    stepsizes = [7.5, 5, 2, 1]
    # stepsizes = [1.0]
    results_variable = [
        subgradient_optimization(
            subproblems,
            agents,
            np.array(u_cent),
            tol=1e-10,
            n_iter=100,
            alpha=alpha,
            steps="variable",
        )
        for alpha in stepsizes
    ]

    results_constant = [
        subgradient_optimization(
            subproblems,
            agents,
            np.array(u_cent),
            tol=1e-10,
            n_iter=100,
            alpha=alpha,
            steps="constant",
        )
        for alpha in stepsizes
    ]

    iterations_variable, e_hist_variable = [
        [res[i] for res in results_variable] for i in (0, 3)
    ]

    iterations_constant, e_hist_constant = [
        [res[i] for res in results_constant] for i in (0, 3)
    ]
    plt.figure()
    plot_variable_convergence(
        iterations_constant,
        e_hist_constant,
        iterations_variable,
        e_hist_variable,
        stepsizes,
    )


def problem_c(agents, subproblems, u_cent, x_cent):
    # betas = [0.5, 1, 2, 3.5, 5]
    result_accelerated = accelerated_subgradient_optimization(
        subproblems,
        agents,
        np.array(u_cent),
        tol=1e-10,
        n_iter=100,
        alpha=2.5,
        steps="constant",
    )

    results_normal = subgradient_optimization(
        subproblems,
        agents,
        np.array(u_cent),
        tol=1e-10,
        n_iter=100,
        alpha=2.5,
        steps="constant",
    )

    # iterations_accelerated, e_hist_accelerated = [
    #     [res[i] for res in result_accelerated] for i in (0, 3)
    # ]

    accelerated_method(
        results_normal[0],
        results_normal[3],
        result_accelerated[0],
        result_accelerated[3],
    )


if __name__ == "__main__":
    global u_max
    agents, u_max, Tfinal = load_variables()

    subproblems = [sub_problem(agent) for agent in agents]
    lambdas = np.zeros([8, 1])

    u_cent, x_cent = centralized_solution(agents)

    # problem_a(agents, subproblems, u_cent, x_cent)
    # problem_b(agents, subproblems, u_cent, x_cent)
    problem_c(agents, subproblems, u_cent, x_cent)
    plt.show()
