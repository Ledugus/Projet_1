# Simulation d'une grue sur une barge flottante
import math
import matplotlib.pyplot as plt
import numpy as np

### Constantes

g = 9.81  # gravitation [m/s**2]

### Paramètres du système

m_barge = 10  # masse de la barge [kg]
m_charge = 0.7  # mass de la charge déplacée [kg]
l = 1  # longueur de la barge [m]
h_barge = 0.08  # hauteur de la barge [m]
d = 2  # distance de la charge
h_charge = 0.1
D = 0.3 # coefficient de frottement visqueux
def calculate_moment_inertie():
    moment_barge = (m_barge * (l ** 2 + h_barge ** 2)) / 12
    moment_charge = m_charge * (d**2+(h_charge-h_im)**2)
    moment_inertie = moment_barge + moment_charge
    return moment_inertie

# Grandeurs constantes calculées
m_tot = m_charge + m_barge 
f_poussee_gravite = m_tot * g
h_im = m_tot / ((l ** 2) * 1000)
z_g = (m_barge * (h_barge / 2) + m_charge * h_charge / 2) / m_tot
moment_inertie =  calculate_moment_inertie()
angle_stable = math.atan((m_charge*d)/(m_tot * ((l**2)/(12*h_im)-h_im/2-(z_g-h_im))))
print(moment_inertie)
print(z_g-h_im)
### Paramètres de la simulation

step = 0.001  # pas (dt) [s]
end = 20.0  # durée [s]
theta_0 = 0  # angle d'inclinaison du départ

t = np.arange(0, end, step)
x_g = np.empty_like(t)
theta = np.empty_like(t)
v_theta = np.empty_like(t)
a_theta = np.empty_like(t)
x_c = np.empty_like(t)


def get_x_c(theta):
    h1 = h_im + math.tan(theta)*l/2
    h2 = h_im - math.tan(theta)*l/2
    l_C = l * (h1 + 2*h2)/(3*(h1+h2))
    h_C = (h1**2+h1*h2+h2**2)/(3*(h1+h2))
    Y_C = -l/2 + l_C
    Z_C = -h_im + h_C
    # print("Theta", theta%(2*math.pi), "Y_C", Y_C, "Z_C", Z_C)
    y_C = Y_C * math.cos(theta) - Z_C * math.sin(theta)
    # print("Différence : ", y_C - (l**2 / (12 * h_im) - (h_im / 2)) * math.sin(theta))
    return y_C

def get_x_g(theta):
    return (math.sin(theta) * (z_g - h_im))


def simulation():
    """
    pre: 
    post: exécute une simulation jusqu'à t=end par pas de dt=step.
          Remplit les listes x, v, a des positions, vitesses et accélérations.
    """
    # conditions initiales
    theta[0] = theta_0

    for i in range(len(t) - 1):
        dt = step
        # calcul des centres de force
        x_g[i] = get_x_g(theta[i])
        x_c[i] = get_x_c(theta[i])
        # if i % 100 == 0:
        #     print(theta[i], "poussee", x_c[i], "gravité", x_g[i])
        # calcul des couples
        couple_destabilisateur = m_charge * g * d
        couple_redressement = -f_poussee_gravite * (x_c[i]-x_g[i])
        couple_tot = couple_destabilisateur + couple_redressement
        # calcul accélération, vitesse, position
        a_theta[i] = (couple_tot / moment_inertie) - D * v_theta[i]
        v_theta[i + 1] = v_theta[i] + a_theta[i] * dt
        theta[i + 1] = theta[i] + v_theta[i] * dt
        a_theta[i + 1] = a_theta[i]


def graphiques():
    plt.figure(1)
    plt.subplot(3, 1, 1)
    plt.plot(t, theta, label="Angle de la plateforme")
    plt.axhline(y=angle_stable)
    plt.legend()
    plt.subplot(3, 1, 2)
    plt.plot(t, v_theta, label="Vitesse angulaire de la plateforme")
    plt.legend()
    plt.subplot(3, 1, 3)
    plt.plot(t, a_theta, label="Accélération angulaire de la plateforme")
    plt.legend()
    plt.show()


"""
def graphiques_energie():
    e_ressort = np.empty_like(t)
    e_cin = np.empty_like(t)
    e_tot = np.empty_like(t)
    for i in range(len(t)):
        e_ressort[i] = k * x[i] ** 2 / 2
        e_cin[i] = m * v[i] ** 2 / 2
        e_tot[i] = e_ressort[i] + e_cin[i]

    plt.figure(2)
    plt.plot(t, e_ressort, label="ressort")
    plt.plot(t, e_cin, label="cinétique")
    plt.plot(t, e_tot, label="total")
    plt.legend()
    plt.show()

"""
### programme principal

simulation()
graphiques()
# graphiques_energie()

# mu = 0.1
# x_0 = +10.0
# simulation()
# graphiques()
# graphiques_energie()