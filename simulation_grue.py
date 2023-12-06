# Simulation d'une grue sur une barge flottante
import math
import matplotlib.pyplot as plt
import numpy as np

# Constantes

g = 9.81  # gravitation [m/s**2]

# Paramètres du système

m_charge = 0.2  # mass de la charge déplacée [kg]
h_charge = 0.3  # hauteur de la charge au bas de la barge
# distance de la charge au centre de la barge [m]

m_caisse = 0  # masse du chargement sur la barge
h_caisse = 0.1  # hauteur du centre de gravité du chargement sur la barge

l_barge = 0.6  # longueur de la barge [m]
h_barge = 0.08  # hauteur de la barge [m]
vol_barge = h_barge * l_barge ** 2  # volume de la barge [m³]
rho_barge = 10 / vol_barge  # masse volumique de la barge [kg/m³]
m_barge = rho_barge * vol_barge  # masse de la barge [kg]
d = 2
D = 0.6  # coefficient de frottement visqueux


def calculate_moment_inertie():
    moment_barge = (m_barge * (l_barge ** 2 + h_barge ** 2)) / 12
    moment_caisse = m_caisse * (h_caisse - h_c) ** 2
    moment_charge = m_charge * (d ** 2 + (h_charge - h_c) ** 2)
    return moment_barge + moment_caisse + moment_charge


def calculate_stable_angle(distance: float, masse_charge: float) -> float:
    return math.atan((masse_charge * distance) / (m_tot * ((l_barge ** 2) / (12 * h_c) + h_c / 2 - (z_g))))


def d_by_time(time: float, time_max: float) -> float:
    return min(d*time/time_max + 0.3, d)


# Grandeurs constantes calculées
m_tot = m_charge + m_barge + m_caisse
f_poussee_gravite = m_tot * g
h_c = m_tot / ((l_barge ** 2) * 1000)
z_g = (m_barge * (h_barge / 2) + m_charge *
       h_charge + m_caisse * h_caisse / 2) / m_tot
moment_inertie = calculate_moment_inertie()
angle_stable = calculate_stable_angle(d, m_charge)
angle_submersion = math.atan(2 * (h_barge - h_c) / l_barge)
angle_soulevement = math.atan(2 * h_c / l_barge)

# Paramètres de la simulation

step = 0.001  # pas (dt) [s]
end = 20.0  # durée [s]
theta_0 = 0  # angle d'inclinaison du départ [rad]

t = np.arange(0, end, step)
x_g = np.empty_like(t)
theta = np.empty_like(t)
v_theta = np.empty_like(t)
a_theta = np.empty_like(t)
x_c = np.empty_like(t)


def get_x_c(theta: float) -> float:
    return (l_barge ** 2 / (12 * h_c) - (h_c / 2)) * math.sin(theta)


def get_x_g(theta: float, d: float) -> float:
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
        x_g[i] = get_x_g(theta[i], d_by_time(i*dt, end/4))
        x_c[i] = get_x_c(theta[i])

        # calcul des couples
        couple_gravite = f_poussee_gravite * x_g[i]
        couple_archimede = -f_poussee_gravite * x_c[i]
        couple_tot = couple_gravite + couple_archimede

        # calcul accélération et vitesse angulaire, angle
        a_theta[i] = (couple_tot / moment_inertie) - D * v_theta[i]
        v_theta[i + 1] = v_theta[i] + a_theta[i] * dt
        theta[i + 1] = theta[i] + v_theta[i] * dt
        a_theta[i + 1] = a_theta[i]


def graphique_angle_temps(draw_lines=True):
    plt.figure(1)
    plt.subplot(3, 1, 1)
    plt.plot(t, theta, label="Angle de la plateforme")
    plt.xlabel("Temps (s)")
    plt.ylabel("Angle (rad)")
    # ligne de stabilité
    plt.axhline(y=angle_stable, label="Stabilité théorique",
                color="blue", linestyle="dotted")
    if draw_lines:
        # ligne de submersion
        plt.axhline(y=angle_submersion, label="Submersion", xmin=0,
                    xmax=end, color="black", linestyle="dotted")
        plt.axhline(y=-angle_submersion, xmin=0, xmax=end,
                    color="black", linestyle="dotted")
        # ligne de soulevement
        plt.axhline(y=angle_soulevement, label="Soulèvement",
                    xmin=0, xmax=end, color="red", linestyle="dotted")
        plt.axhline(y=-angle_soulevement, xmin=0, xmax=end,
                    color="red", linestyle="dotted")
    plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(t, v_theta, label="Vitesse angulaire de la plateforme")
    plt.legend()
    plt.subplot(3, 1, 3)
    plt.plot(t, a_theta, label="Accélération angulaire de la plateforme")
    plt.legend()

    plt.savefig("angle_temps_charge_variable.png")
    plt.show()


def diagramme_des_phases():
    plt.figure(1)
    plt.plot(theta, v_theta, label="Phase du système")
    plt.xlabel("Angle (rad)")
    plt.ylabel("Vitesse (rad/s)")
    plt.legend()
    plt.plot(angle_stable, 0, 'bo')
    plt.savefig("diagramme_phases.png")
    plt.show()


def inclinaison_distance_masse():
    plt.figure(1)
    masse_max = 1
    step_masse = 0.2
    d_max = 2  # [m]
    d_step = 0.1
    m = np.arange(0, masse_max, step_masse)
    distance_array = np.arange(0, d_max, d_step)
    angle_stable_array = np.empty_like(distance_array)
    for masse in m:
        for index, distance in enumerate(distance_array):
            angle_stable_array[index] = calculate_stable_angle(distance, masse)
        plt.plot(distance_array, angle_stable_array, label=f'm = {masse}')
    plt.xlabel("Distance (m)")
    plt.ylabel("Angle stable (rad)")
    plt.legend()
    # ligne de submersion
    plt.axhline(y=angle_submersion, label="Submersion", xmin=0,
                xmax=end, color="black", linestyle="dotted")
    # ligne de soulevement
    plt.axhline(y=angle_soulevement, label="Soulèvement", xmin=0,
                xmax=end, color="red", linestyle="dotted")
    plt.legend()
    plt.savefig("inclinaison_distance_masse.png")
    plt.show()


# Programme principal
if __name__ == "__main__":
    simulation()
    graphique_angle_temps(draw_lines=False)
