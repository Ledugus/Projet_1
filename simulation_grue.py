# Simulation d'une grue sur une barge flottante
import math
import matplotlib.pyplot as plt
import numpy as np

### Constantes

g = 9.81  # gravitation [m/s**2]

### Paramètres du système
m_charge = 0.7  # mass de la charge déplacée [kg]
h_charge = 0.1 # hauteur de la charge au bas de la barge
d = 2  # distance de la charge au centre de la barge [m]

m_caisse = 0
h_caisse = 0.1

l_barge = 0.6  # longueur de la barge [m]
h_barge = 0.08  # hauteur de la barge [m]
vol_barge = h_barge * l_barge ** 2 # volume de la barge [m³]
rho_barge = 1000/8 # masse volumique de la barge [kg/m³]
m_barge = rho_barge * vol_barge # masse de la barge [kg]

D = 0.3 # coefficient de frottement visqueux

def calculate_moment_inertie():
    moment_barge = (m_barge * (l_barge ** 2 + h_barge ** 2)) / 12
    moment_caisse = m_caisse * (h_caisse-h_c)**2
    moment_charge = m_charge * (d**2+(h_charge-h_c)**2)
    moment_inertie = moment_barge + moment_caisse + moment_charge
    return moment_inertie

def calculate_stable_angle():
    return math.atan((m_charge*d)/(m_tot * ((l_barge**2)/(12*h_c)-h_c/2-(z_g-h_c))))


# Grandeurs constantes calculées
m_tot = m_charge + m_barge + m_caisse
f_poussee_gravite = m_tot * g
h_c = m_tot / ((l_barge ** 2) * 1000)
z_g = (m_barge * (h_barge / 2) + m_charge * h_charge + m_caisse * h_caisse / 2) / m_tot
moment_inertie = calculate_moment_inertie()
angle_stable = calculate_stable_angle()
angle_submersion = 0
angle_soulevement = 0

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
    return (l_barge**2 / (12 * h_c) - (h_c / 2)) * math.sin(theta)


def get_x_g(theta):
    return (math.sin(theta) * (z_g - h_c)) + ((math.cos(theta) * m_charge * d) / m_tot)


def simulation():
    """
    pre: 
    post: exécute une simulation jusqu'à t=end par pas de dt=step.
          Remplit les listes theta, v_theta, a_theta des angles, vitesses et accélérations angulaires.
    """
    # conditions initiales
    theta[0] = theta_0

    for i in range(len(t) - 1):
        dt = step
        # calcul des centres de force
        x_g[i] = get_x_g(theta[i])
        x_c[i] = get_x_c(theta[i])

        # calcul des couples
        couple_gravite = f_poussee_gravite * x_g[i]
        couple_archimede = -f_poussee_gravite * x_c[i]
        couple_tot = couple_gravite + couple_archimede

        # calcul accélération, vitesse, angle
        a_theta[i] = (couple_tot / moment_inertie) - D * v_theta[i]
        v_theta[i + 1] = v_theta[i] + a_theta[i] * dt
        theta[i + 1] = theta[i] + v_theta[i] * dt
        a_theta[i + 1] = a_theta[i]


def graphiques():
    plt.figure(1)
    plt.subplot(3, 1, 1)
    plt.plot(t, theta, label="Angle de la plateforme")
    # ligne de stabilité
    plt.axhline(y=angle_stable, label="Stabilité théorique", color="blue") 
    # ligne de submersion
    plt.axhline(y=angle_submersion, color="red") 
    plt.axhline(y=-angle_submersion, color="red")
    # ligne de soulevement
    plt.axhline(y=angle_soulevement, color="orange")
    plt.axhline(y=-angle_soulevement, color="orange")
    plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(t, v_theta, label="Vitesse angulaire de la plateforme")
    plt.legend()
    plt.subplot(3, 1, 3)
    plt.plot(t, a_theta, label="Accélération angulaire de la plateforme")
    plt.legend()
    plt.show()

### programme principal
if __name__ == "__main__":
    simulation()
    graphiques()