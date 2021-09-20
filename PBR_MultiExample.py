import numpy as np
from scipy.integrate import solve_ivp, odeint
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec

DelH1 = -25000  # cal/mol
DelH2 = -125000  # cal/mol
DelH3 = -100000  # cal/mol

Ua = 1.49  # cal/(K s kg_cat)
mcp = 1000  # cal/(s K)
sum_ = 45  # sum of Cp cal/mol-K

P0 = 5  # bar
Fa0 = 10  # mol/s
T0 = 850  # K

Ct0 = P0/.08314/T0  # mol/L
Ca0 = Ct0 * 0.33  # mol/L
Ft0 = 30  # mol/s
v0 = Ft0/Ct0


def k1(T):
    return 2e5 * np.exp(-17000/1.99/T)  # L^1.5 / (mol^0.5 kg-cat sec)


def k2(T):
    return 8e6 * np.exp(-23000/1.99/T)  # L^1.5 / (mol^0.5 kg-cat sec)


def k3(T):
    return k2(T)


def C(F_comp, F_t, T):
    return F_comp/v0 * F_t/Ft0 * T0/T


IniVal = [10, 20, 0, 0, 0, 850, 300]

d_dw = IniVal


def Reactor(VEC, W):
    FA = VEC[0]
    FB = VEC[1]
    FC = VEC[2]
    FD = VEC[3]
    FE = VEC[4]
    T = VEC[5]
    Ta = VEC[6]

    FT = np.sum(VEC[0:5])

    r1 = - k1(T) * C(FA, FT, T)**0.5 * C(FB, FT, T)
    r2 = - k2(T) * C(FA, FT, T)**0.5 * C(FB, FT, T)
    r3 = - k3(T) * C(FC, FT, T)**0.5 * C(FB, FT, T)

    d_dw[0] = r1+r2  # mol/(kg-cat s)
    d_dw[1] = r1/2.0 + 3*r2 + 5.0*r3/2.0  # mol/(kg-cat s)
    d_dw[2] = -r1 + r3  # mol/(kg-cat s)
    d_dw[3] = -2.0*r2 - 2.0*r3  # mol/(kg-cat s)
    d_dw[4] = d_dw[3]  # mol/(kg-cat s)

    d_dw[5] = (r1*DelH1 + r2*DelH2 - r3*DelH3 - Ua *
               (T - Ta)) / (FT*sum_)  # K/ kg-cat
    d_dw[6] = Ua/mcp * (T-Ta)  # K/ kg-cat

    return d_dw


w = np.linspace(0, 70, 500)

fig = plt.figure(constrained_layout=True)
gs = gridspec.GridSpec(2, 1, figure=fig)

# x = solve_ivp(Reactor, [0, 10], IniVal, method="RK45")
# ax.plot(x.t, x.y)

x = odeint(Reactor, IniVal, w)

ax = fig.add_subplot(gs[0, 0])
ax.plot(w, x[:, 0:5], label=["FA", "FB", "FC", "FD", "FE"])
ax.grid(True)
# plt.gca().legend()
ax.set_xlabel('Massa de catalisador (kg)')
ax.set_ylabel('F (mol/s)')
ax.set_title("")

ax2 = fig.add_subplot(gs[1, 0])

ax2.plot(w, x[:, 5:])
ax2.grid(True)
ax2.set_xlabel('Massa de catalisador (kg)')
ax2.set_ylabel('Temperatura (K)')
ax2.set_title("")
# plot results


# text = '\n'.join((
#     r'$\mathrm{P_{entrada}}=%.2f\,kPa$' % (P, ),
#     r'$\mathrm{y_{H_{2}S}}=%.2e$' % (y_A_0, ),
#     r'$\mathrm{P_{H_{2}S}}=%.2f\,kPa$' % (P_H2S, ),
#     r'$\mathrm{Q_{H_{2}S}}=%.2e\,mol/s$' % (F_A, ),
#     r'$\mathrm{T}=%.2f\,K$' % (T, )))

# props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# ax.text(0.5, 0.45, text, transform=ax.transAxes, fontsize=14,
#         verticalalignment="top", bbox=props)


plt.show()
