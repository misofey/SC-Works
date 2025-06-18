from scipy.io import loadmat
from dataclasses import dataclass
import numpy as np


@dataclass
class Agent:
    A: np.ndarray
    B: np.ndarray
    x0: np.ndarray


def load_variables():
    a = loadmat("agents.mat")

    agent_1 = Agent(a["A1"], a["B1"], a["x01"])
    agent_2 = Agent(a["A2"], a["B2"], a["x02"])
    agent_3 = Agent(a["A3"], a["B3"], a["x03"])
    agent_4 = Agent(a["A4"], a["B4"], a["x04"])

    return agent_1, agent_2, agent_3, agent_4, a["umax"], a["Tfinal"]


if __name__ == "__main__":
    print(load_variables())
