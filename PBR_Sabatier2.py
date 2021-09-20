import numpy as np
from scipy.integrate import solve_ivp, odeint
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
# import matplotlib.gridspec as gridspec


# Program based on:
#  Miguel, C. V., Mendes, A., & Madeira, L. M. (2018).
#  Intrinsic kinetics of CO2 methanation over an industrial nickel-based catalyst.
#  Journal of CO2 Utilization, 25, 128–136. doi:10.1016/j.jcou.2018.03.011 

def k(T):
    return 0.8930/100 * np.exp(14275.0 / T)  # mol / (g-cat h)


def K_eq(T):
    return np.exp((1.0/1.987)*(56000.0/T/T + 34633.0/T - 16.4*np.log(T) + 0.00557*T) + 33.165)


def K_OH(T):
    return 0.4326 * np.exp(7414.0 / T)  # 1/kPa^0.5


P0 = 1000  # kPa
T = 100+273.15  # K
F_H2_ini = 2207.96
Ratio_H2_CO2 = 4.0


IniVal = [F_H2_ini/Ratio_H2_CO2, F_H2_ini, 0, 0]
d_dw = IniVal
y = IniVal


def Reactor(W, VEC):

    F_CO2 = VEC[0]
    F_H2 = VEC[1]
    F_CH4 = VEC[2]
    F_H2O = VEC[3]

    FT = np.sum(VEC)
    for i in range(len(VEC)):
        if (VEC[i] != 0):
            y[i] = np.abs(VEC[i]/FT)
        else:
            y[i] = 0

    # Ideal gas partial pressure calculation (kPa)
    P_CO2 = y[0]*P0
    P_H2 = y[1]*P0
    P_CH4 = y[2]*P0
    P_H2O = y[3]*P0

    # print(P_CO2, P_H2, P_CH4, P_H2O)
#
    r_ch4 = k(T)*pow(P_CO2, 0.5)*pow(P_H2, 0.5)/pow(1.0 + K_OH(T)*P_H2O /
                                                    pow(P_H2, 0.5), 2.0) * (1.0 - P_CH4*P_H2O*P_H2O/P_CO2/pow(P_H2, 4.0)/K_eq(T))
    # print(r_ch4)
    d_dw[0] = - r_ch4  # mol/(kg-cat s)
    d_dw[1] = - 4.0 * r_ch4  # mol/(kg-cat s)
    d_dw[2] = r_ch4  # mol/(kg-cat s)
    d_dw[3] = 2 * r_ch4  # mol/(kg-cat s)

    return d_dw


# w = np.linspace(0, 70, 1000)

fig = plt.figure(constrained_layout=True)
ax = fig.add_subplot()
# gs = gridspec.GridSpec(2, 1, figure=fig)

x = solve_ivp(Reactor, [0, 500e3], IniVal, method="RK45")
ax.plot(x.t[0:]/1000, x.y[0, 0:], '^-r', label=r'$\mathrm{F_{CO_2}}$')
ax.plot(x.t[1:]/1000, x.y[1, 1:], '--', label=r'$\mathrm{F_{H_2}}$')
ax.plot(x.t/1000, x.y[2, :], 'k-', label=r'$\mathrm{F_{CH_4}}$')
ax.plot(x.t/1000, x.y[3, :], 'b.-', label=r'$\mathrm{F_{H_{2}O}}$')


ax.grid(True)
ax.legend()
plt.gca().set_title(r" $\bf{Kinetics\,\, data\,\, from:}$" + "\n" + r"Miguel, C. V., Mendes, A., & Madeira, L. M. (2018)." + "\n" +
                    r"Intrinsic kinetics of $CO_2$ methanation over an industrial nickel-based catalyst." + "\n" + r"Journal of $CO_2$ Utilization, 25, 128–136. DOI: 10.1016/j.jcou.2018.03.011", fontsize=9)
ax.set_xlabel('Massa de catalisador (kg)')
ax.set_ylabel('Fluxo (mol/h)')

text = '\n'.join((
    r'$\mathrm{P}=%.2f\,kPa$' % (P0, ),
    r'$\mathrm{F_{H_{2}}}=%.2f\,mol/h$' % (F_H2_ini, ),
    r'$\mathrm{F_{H_2}/F_{CO_2}}=%.2f$' % (Ratio_H2_CO2, ),
    r'$\mathrm{T}=%.2f\,K$' % (T, )))

props = dict(boxstyle='round', Facecolor='white', alpha=1)
ax.text(0.1, 0.97, text, transform=ax.transAxes, fontsize=9,
        verticalalignment="top", bbox=props)


plt.show()
