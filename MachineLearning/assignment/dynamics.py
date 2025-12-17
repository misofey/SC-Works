import math
import torch

class DubinsCar:
    def __init__(self, tau=0.1, **kwargs):
        self.tau = tau

    def forward(self, x: torch.Tensor, u: torch.Tensor):
        steering_angle = x[..., -1]
        heading_velocity, omega = u[..., 0], u[..., 1]

        # Dubins car update
        dx = self.tau * heading_velocity * torch.cos(steering_angle)
        dy = self.tau * heading_velocity * torch.sin(steering_angle)
        dtheta = self.tau * omega

        next_x = x.clone()
        next_x[..., 0] += dx
        next_x[..., 1] += dy
        next_x[..., 2] += dtheta

        # Keep in angle range
        next_x[..., 2] = next_x[..., 2] % (2 * math.pi)

        return next_x