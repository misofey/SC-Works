import torch.nn as nn
import torch.nn.functional as F

class NeuralNetworkDynamics(nn.Module):  # inherit the nn.Module class for backpropagation and training functionalities
    def __init__(self):
        super(NeuralNetworkDynamics, self).__init__()
        # TODO: Implement a NN with at least 3 layers, and ReLU activation
        # TODO: Feel free to use any kind of structure (e.g. dense feedforward, residual, etc). This is an engineering choice
        # TODO: You can consult, for instance: https://docs.pytorch.org/tutorials/beginner/basics/buildmodel_tutorial.html