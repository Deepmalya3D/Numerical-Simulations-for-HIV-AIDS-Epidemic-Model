import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

'''
Epidemic Model for HIV/AIDS under consideration
S˙ = A − βIS − μS − dS,
I˙ = βIS + r1T − dI − k1I − k2I,
A˙ = k1I − (δ1 + d)A + r2T,
T˙ = k2I − r1T − (d + δ2 + r2)T,
R˙ = μS − dR.
'''


def f(l, t):
    a = 0.61
    b = 0.028
    u = 0.055
    d = 0.0195
    k1 = 0.15
    k2 = 0.35
    r1 = 0.08      # variable
    r2 = 0.03      # variable
    d1 = 0.0807
    d2 = 0.0657
    s = l[0]
    i = l[1]
    a_ = l[2]
    z = l[3]

    ds_dt = a - b*i*s - u*s - d*s
    di_dt = b*i*s + r1*z - d*i - k1*i - k2*i
    da_dt = k1*i - (d1+d)*a_ + r2*z
    dz_dt = k2*i - r1*z - (d+d2+r2)*z
    return [ds_dt, di_dt, da_dt, dz_dt]


t = np.linspace(0, 200)
s0 = [37.7, 25, 15, 10]

s = odeint(f, s0, t)

plt.plot(t, s[:, 0], "b-", lw=0.5)
plt.plot(t, s[:, 1], "g-", lw=0.5)
plt.plot(t, s[:, 2], "r-", lw=0.5)
plt.plot(t, s[:, 3], "m-", lw=0.5)
plt.xlabel("Time")
plt.ylabel("Population Size")
plt.margins(0)
plt.legend(["S", "I", "A", "T"])
plt.show()
