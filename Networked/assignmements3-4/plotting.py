import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["text.usetex"] = True


def convergence_error(iterations, error):
    plt.semilogy(iterations, error)
    plt.grid(which="both")
    plt.xlabel("Iterations")
    plt.ylabel("Error")


def final_state_convergence(iterations, xfinal_hist):
    plt.plot(iterations, xfinal_hist[:, 0, :], "-")
    plt.plot(xfinal_hist[:, 1, :], "--")
    plt.grid()
    plt.xlabel("Iterations")
    plt.ylabel("Final tates")
    proxy_central = plt.Line2D([0], [0], color="black", linestyle="-", label="State 1")
    proxy_decentr = plt.Line2D([0], [0], color="black", linestyle="--", label="State 2")
    plt.legend(handles=[proxy_central, proxy_decentr], loc="best")


def compare_cent_decentr(u_cent: [], x_cent: [], u_dec: [], x_dec: []):
    nt = len(u_cent[0])

    t = range(1, nt + 1)

    plt.plot(t, x_dec[0][:, 0], label="Decentralized", linestyle="-")
    plt.plot(t, x_cent[0][:, 0], label="Centralized", linestyle="--")
    # plt.plot(t, x_cent[0][:, 1])
    # plt.plot(t, x_dec[0][:, 1])
    plt.xlabel("t")
    plt.ylabel("$x_1(t)[0]$")
    plt.legend()
    plt.grid()


def plot_agent_subplots(x_cent, x_conv):
    """
    Creates a 4x2 grid of subplots for 4 agents and 2 state dimensions,
    comparing centralized (solid) vs decentralized (dashed) trajectories.

    Parameters:
    - t       : 1D array of time steps (length NT)
    - x_cent  : list of 4 arrays, each of shape (NT, 2)
    - x_conv  : list of 4 arrays, each of shape (NT, 2)
    """
    nt = 20
    t = range(1, nt + 1)

    fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(12, 16), sharex=True)
    for i in range(4):
        for j in range(2):
            ax = axes[i, j]
            ax.plot(t, x_cent[i][:, j], label="Centralized")
            ax.plot(t, x_conv[i][:, j], "--", label="Decentralized")
            ax.set_title(f"Agent {i+1}, State dim {j+1}")
            ax.set_ylabel("Value")
            ax.grid(True)
    # only add legend once in top-left subplot
    axes[0, 0].legend(loc="best")
    plt.xlabel("Time step")
    plt.tight_layout()
    plt.show()


def plot_2d_convergence(iterations, e_hist, xfinal_hist, x_dec, x_cent):
    # plt.figure()
    # plt.figure()
    T, _, n_agents = xfinal_hist.shape
    for i in range(n_agents):
        traj = xfinal_hist[:, :, i]
        plt.scatter(
            traj[:, 0], traj[:, 1], marker="o", linewidths=0.1, label=f"Agent {i+1}"
        )
    # Overlay centralized optimal states
    plt.scatter(
        x_cent[0][-1, 0],
        x_cent[0][-1, 1],
        marker="x",
        s=100,
        label="Centralized Optimum",
    )
    # plt.title("Agent State Trajectories in State-Space")
    plt.xlabel("State Dimension 1")
    plt.ylabel("State Dimension 2")
    plt.legend(loc="best")
    plt.grid(True)

    plt.show()


def plot_variable_convergence(
    iterations_constant, e_hist_constant, iterations_variable, e_hist_variable, alphas
):
    for i in range(len(alphas)):
        plt.semilogy(
            iterations_variable[i],
            e_hist_variable[i],
            label="alpha = " + str(alphas[i]),
        )
    proxy_constant = plt.Line2D(
        [0], [0], color="black", linestyle="--", label="Constant step size"
    )
    proxy_variable = plt.Line2D([0], [0], color="black", label="Diminishing step size")

    a = plt.gcf().legend(loc="outside right upper")
    plt.gcf().add_artist(a)

    plt.gcf().legend(
        handles=[proxy_constant, proxy_variable], loc="outside right lower"
    )
    plt.gca().set_prop_cycle(None)

    for i in range(len(alphas)):
        plt.semilogy(
            iterations_constant[i],
            e_hist_constant[i],
            # label="$alpha=$" + str(alphas[i]),
            linestyle="--",
        )
    plt.grid()
    plt.xlabel("Iterations")
    plt.ylabel("Error")


def accelerated_method(
    iterations_normal, e_hist_normal, iterations_accelerated, e_hist_accelerated
):
    plt.semilogy(
        iterations_accelerated,
        e_hist_accelerated,
        label="alpha = 2.5, Nesterov's method",
    )
    # proxy_constant = plt.Line2D(
    #     [0], [0], color="black", linestyle="--", label="Constant step size"
    # )
    # proxy_variable = plt.Line2D([0], [0], color="black", label="Accelerated method")
    # a = plt.gcf().legend(loc="outside right upper")
    # plt.gcf().add_artist(a)

    # plt.gcf().legend(
    #     handles=[proxy_constant, proxy_variable], loc="outside right lower"
    # )
    plt.gca().set_prop_cycle(None)
    plt.semilogy(
        iterations_normal,
        e_hist_normal,
        label="alpha=2.5, constant step size",
        linestyle="--",
    )
    plt.legend(loc="best")

    plt.xlabel("Iterations")
    plt.ylabel("Error")
    plt.grid()


def admm_rho_effect(iterations, e_hist, rhos):
    for i in range(len(rhos)):
        plt.semilogy(
            iterations[i],
            e_hist[i],
            label="beta = " + str(rhos[i]),
        )
    # proxy_constant = plt.Line2D(
    #     [0], [0], color="black", linestyle="--", label="Constant step size"
    # )
    # proxy_variable = plt.Line2D([0], [0], color="black", label="Accelerated method")
    # a = plt.gcf().legend(loc="outside right upper")
    # plt.gcf().add_artist(a)

    # plt.gcf().legend(
    #     handles=[proxy_constant, proxy_variable], loc="outside right lower"
    # )
    # plt.gca().set_prop_cycle(None)
    # plt.semilogy(
    #     iterations_normal,
    #     e_hist_normal,
    #     label="alpha=5, constant step size",
    #     linestyle="--",
    # )
    plt.legend(loc="best")

    plt.xlabel("Iterations")
    plt.ylabel("Error")
    plt.grid(which="both")
