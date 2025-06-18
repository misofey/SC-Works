import numpy as np
from scipy import optimize

H = np.array([[2, 1, 1], [1, 1, 1], [1, 1, 2]])
c = np.array([0, 1, 1])
u_og = np.array([1, 0, 1])

# update step 1
one = np.array([1, 0, 0])
update1_a = one @ H @ one
update1_b = one @ H @ ((1 - one) * u_og) + 0.0 * (u_og * (1 - one)) @ H @ one + c @ one
update1_u1 = -update1_b / update1_a
print(update1_a, update1_b)
print(update1_u1)

# update step 2
two = np.array([0, 1, 0])
update2_a = two @ H @ two
update2_b = two @ H @ ((1 - two) * u_og) + 0.0 * (u_og * (1 - two)) @ H @ two + c @ two
update2_u2 = -update2_b / update2_a
print(update2_a, update2_b)
print(update2_u2)

# update step 3
three = np.array([0, 0, 1])
update3_a = three @ H @ three
update3_b = (
    three @ H @ ((1 - three) * u_og)
    # + 0.5 * (u_og * (1 - three)) @ H @ three
    + c @ three
)
update3_u3 = -update3_b /3 update3_a
print(update3_a, update3_b)
print(update3_u3)


u_new = np.array([update1_u1, update2_u2, update3_u3])
cost_new = 0.5 * u_new @ H @ u_new + c @ u_new
cost_og = 0.5 * u_og @ H @ u_og + c @ u_og
u_weighted = u_og / 2 + u_new / 2
cost_weighted = 0.5 * u_weighted @ H @ u_weighted + c @ u_weighted

print(u_new)
print(cost_og)
print(cost_new)
print(cost_weighted)
