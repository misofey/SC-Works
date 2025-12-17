import torch

class Data:
    def __init__(self, x, u, next_x):
        self.x = x
        self.u = u
        self.extended_input = torch.cat((x, u), dim=-1)
        self.output = next_x


class Dataset:
    def __init__(self, dynamics, state_support, control_support, dataset_size):
        self.dynamics = dynamics
        self.state_support = state_support
        self.control_support = control_support
        self.dataset_size = dataset_size

        self.data = self._generate()

    def _generate(self):

        x = self.state_support.lower + (self.state_support.upper - self.state_support.lower) * torch.rand(
            (self.dataset_size, self.state_support.lower.shape[0]))
        u = self.control_support.lower + (self.control_support.upper - self.control_support.lower) * torch.rand(
            (self.dataset_size, self.control_support.lower.shape[0]))

        x_next = self.dynamics.forward(x, u)

        return Data(x, u, x_next)
