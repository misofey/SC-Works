import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def plot_environment(safe_set, obstacle, target, reference_trajectory_for_position=None):
    def plot_rect(ax, rect, facecolor, edgecolor='black', label=None, textcolor='black'):
        lower, upper = rect.lower[0:2], rect.upper[0:2]
        width, height = (upper - lower).tolist()
        r = Rectangle(lower.tolist(), width, height,
                      facecolor=facecolor, edgecolor=edgecolor, linewidth=2)
        ax.add_patch(r)
        if label:
            cx = (lower[0] + upper[0]) / 2
            cy = (lower[1] + upper[1]) / 2
            ax.text(cx, cy, label, color=textcolor, ha='center', va='center',
                    fontsize=12, fontweight='bold')

    fig, ax = plt.subplots(figsize=(5, 5))

    # Plot sets
    plot_rect(ax, safe_set, facecolor='white', edgecolor='black')
    plot_rect(ax, obstacle, facecolor='black', edgecolor='black', label='obs', textcolor='white')
    plot_rect(ax, target, facecolor='gray', edgecolor='black', label='tgt', textcolor='black')

    # Plot reference trajectory if provided
    if reference_trajectory_for_position is not None:
        ref = reference_trajectory_for_position.detach().cpu().numpy()
        ax.plot(ref[:, 0], ref[:, 1], 'ro-', label='Reference trajectory')  # red dots + line

    lower = safe_set.lower[0:2]
    upper = safe_set.upper[0:2]
    ax.set_xlim(lower[0] - 0.5, upper[0] + 0.5)
    ax.set_ylim(lower[1] - 0.5, upper[1] + 0.5)
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('$x_1$')
    ax.set_ylabel('$x_2$')
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.legend()
    plt.show()


def plot_position_trajectories(traj_nn, traj_oracle, labels=("Neural Network", "Oracle")):
    # Extract x, y for both trajectories
    x1, y1 = traj_nn[:, 0].numpy(), traj_nn[:, 1].numpy()
    x2, y2 = traj_oracle[:, 0].numpy(), traj_oracle[:, 1].numpy()

    # Create scatter plot
    plt.figure(figsize=(5, 5))
    plt.scatter(x1, y1, s=10, alpha=0.7, label=labels[0], color="tab:blue")
    plt.scatter(x2, y2, s=10, alpha=0.7, label=labels[1], color="tab:orange")

    plt.xlabel("x1")
    plt.ylabel("x2")
    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.show()