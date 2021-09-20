# CONSTANT PRESSURE, TEMPERATURE, VOLUME

from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt

eps_ = 1
k = 5.0  # dm^3/mol/s
Ca_0 = 0.2  # mol/dm3
v_0 = 1  # dm3/s


def Xfun(X): return (1.0 - eps_ * X)*(1.0 - eps_ * X) / (1.0 - X) / (1.0 - X)


# Fa = v*Ca
#Fa_0 = v_0 * Ca_0
V = v_0/k/Ca_0*quad(Xfun, 0, 1.0)[0]
print(V)
