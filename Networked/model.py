import numpy as np
import scipy as sp
import control as ct

student_number = 5476747
_a, _b, _c = 5, 7, 7
_h = 1

A = np.array([[0, 0.5 - _c], [0.2 + _a - _b]])
B = np.array([[1.0], [0.0]])
C = np.eye(2)
D = np.array([[0.0], [0.0]])

base_sys = ct.ss(A, B, C, D)
# def lsim(sys, u, x0 = None):
#     A = sys.A
#     B = sys.B

#     N = len(u)
#     x = np.zeros((2, N+1))
#     if x0:
#         x[:, 0] = 0
#     for i in range(1, N):
#         x[:, i] = A @ x[:, i-1] + B @ u[i-1]

#     return x
