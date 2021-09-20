import numpy as np
from scipy.integrate import solve_ivp, odeint
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
# import matplotlib.gridspec as gridspec


# def k_co(T):
#     return 3.35e3 * np.exp(-74102/8.314/T)  # mol / (kg-cat sec)


# def k_wg(T):
#     return 3.0421e4 * np.exp(-161740/8.314/T)  # mol / (kg-cat sec Pa^1.5)


# def KA(T):
#     return np.exp(-11.72-7361/T)


# def KB(T):
#     return np.exp(-23.95+8738/T)


# def KC(T):
#     return np.exp(-2.38-778/T)


# def KWG(T):
#     return np.exp(-4.4556+4683/T)

def K_FUN(K0, Q, T):
    return K0 * np.exp(-Q*1000/8.314/T)


def K(COMP, T):
    if COMP == "CO2":
        return K_FUN(4.7e-6, 56.8, T)
    elif COMP == "CO":
        return K_FUN(2.9e-9, 121.4, T)
    elif COMP == "H2":
        return K_FUN(4.3e-10, 105.1, T)
    elif COMP == "H2O":
        return K_FUN(7.6e-7, 94.4, T)
    else:
        print("Component not chosen")


def k0(RXN, T):
    if RXN == "RWGS":
        return K_FUN(1.1e-21, 242.4, T)
    elif RXN == "CO_M":
        return K_FUN(1.9e12, 144.0, T)
    elif RXN == "CO2_M":
        return K_FUN(1.6e7, 94.3, T)
    else:
        print("Reaction not chosen")


def main():
    F_H2 = 34.31

    IniVal = [F_H2/4.0, F_H2, 0, 0, 0]

    d_dw = IniVal

    y = np.zeros(5)


    P0 = 50.0  # bar
    T = 425+273.15
    # F_CO20 = 10  # mol/s
    # T0 = 850  # K

    # Ct0 = P0/.08314/T0  # mol/L
    # Ca0 = Ct0 * 0.33  # mol/L
    # Ft0 = 30  # mol/s
    # v0 = Ft0/Ct0


    def Reactor(VEC, W):

        F_CO2 = VEC[0]
        F_H2 = VEC[1]
        F_CO = VEC[2]
        F_H2O = VEC[3]
        F_CH4 = VEC[4]

        FT = np.sum(VEC)
        for i in range(len(VEC)):
            if (VEC[i] != 0):
                y[i] = np.abs(VEC[i]/FT)
            else:
                y[i] = 0

        P_CO2 = y[0]*P0
        P_H2 = y[1]*P0
        P_CO = y[2]*P0
        P_H2O = y[3]*P0
        print(P_CO2, P_H2, P_CO, P_H2O)

        # print(P_COO2, P_H2, P_CO, P_H2O)
        # r_co = k_co(T) * (KA(T)*P_CO**(0.5)*P_H2**(0.5)) / \
        #     ((1.0 + KA(T)*P_CO + KB(T)*P_H2O / (P_H2**(0.5)))**2.0)

        # r_w = k_wg(T) * (KC(T)*P_CO / (P_H2**(0.5))
        #                  * P_H2O - P_CO2*P_H2**(0.5)/KWG(T))/((1.0 + KA(T)*P_CO + KB(T)*P_H2O/(P_H2**(0.5)))**2.0)
        det = (1.0 + K("CO2", T)*P_CO2 + K("CO", T)*P_CO +
            K("H2", T)*P_H2 + K("H2O", T)**0.5 * P_H2O**0.5)

        r_co = k0("RWGS", T)*P_CO2 / det

        r_ch4_1 = k0("CO2_M", T) * K("CO2", T) * \
            K("H2", T) * P_CO2 * P_H2 / det/det

        r_ch4_2 = k0("CO_M", T) * K("CO", T) * \
            K("H2", T) * P_CO * P_H2 / det/det

        d_dw[0] = -r_ch4_1 - r_co  # mol/(kg-cat s)
        d_dw[1] = -4.0 * r_ch4_1 - r_co - 3.0*r_ch4_2  # mol/(kg-cat s)
        d_dw[2] = r_co - r_ch4_2  # mol/(kg-cat s)
        d_dw[3] = r_ch4_1 + r_ch4_2  # mol/(kg-cat s)
        d_dw[4] = 2.0*r_ch4_1 + r_co + r_ch4_2  # mol/(kg-cat s)

        return d_dw


    w = np.linspace(0, 70)

    fig = plt.figure(constrained_layout=True)
    # gs = gridspec.GridSpec(2, 1, figure=fig)

    # x = solve_ivp(Reactor, [0, 10], IniVal, method="RK45")
    # ax.plot(x.t, x.y)

    x = odeint(Reactor, IniVal, w)
    # gs[0, 0]
    ax = fig.add_subplot()
    ax.plot(w, x)
    # ax.grid(True)
    # plt.gca().legend()
    ax.set_xlabel('Massa de catalisador (g)')
    ax.set_ylabel('F (mol/min)')
    ax.set_title("")

    # ax2 = fig.add_subplot(gs[1, 0])

    # ax2.plot(w, x[:, 5:])
    # ax2.grid(True)
    # ax2.set_xlabel('Massa de catalisador (kg)')
    # ax2.set_ylabel('Temperatura (K)')
    # ax2.set_title("")
    # plot results


    # text = '\n'.join((
    #     r'$\mathrm{P_{entrada}}=%.2f\,kPa$' % (P, ),
    #     r'$\mathrm{y_{H_{2}S}}=%.2e$' % (y_A_0, ),
    #     r'$\mathrm{P_{H_{2}S}}=%.2f\,kPa$' % (P_H2S, ),
    #     r'$\mathrm{Q_{H_{2}S}}=%.2e\,mol/s$' % (F_A, ),
    #     r'$\mathrm{T}=%.2f\,K$' % (T, )))

    # props = dict(boxstyle='round', F_CO2cecolor='wheat', alpha=0.5)
    # ax.text(0.5, 0.45, text, transform=ax.transAxes, fontsize=14,
    #         verticalalignment="top", bbox=props)


    plt.show()  


if (__name__ == "__main__"):
    main()