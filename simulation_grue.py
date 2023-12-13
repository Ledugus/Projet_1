# Simulation d'une grue sur une barge flottante
import math
import matplotlib.pyplot as plt
import numpy as np

# Constantes
g = 9.81  # gravitation [m/s**2]

# Paramètres du système
m_charge = 0.2  # mass de la charge déplacée [kg]
h_charge = 0.3  # hauteur de la charge au bas de la barge [m]
d = 2  # distance de la charge au centre de la barge [m]

m_grue = 4  # masse du chargement sur la barge [kg]
h_grue = 0.5  # hauteur du centre de gravité du chargement sur la barge [m]

l = 0.6  # longueur de la barge [m]
h1 = 0.08  # hauteur de la barge [m]
m1 = 5  # masse de la barge [kg]

D = 0.6  # coefficient de frottement visqueux


def calculate_moment_inertie():
    """Retourne le moment d'inertie de la plateforme"""
    # moment d'un parallélipipède (plateforme)
    moment_barge = (m1 * (l ** 2 + h1 ** 2)) / 12
    # moment d'une masse ponctuelle (grue)
    moment_caisse = m_grue * (h_grue - h_c) ** 2
    # moment d'une masse ponctuelle (charge)
    moment_charge = m_charge * (d ** 2 + (h_charge - h_c) ** 2)

    return moment_barge + moment_caisse + moment_charge  # somme des moments


def calculate_stable_angle(distance: float, masse_charge: float) -> float:
    """Retourne l'angle de stabilité théorique"""
    return math.atan((masse_charge * distance) / (m_tot * ((l ** 2) / (12 * h_c) + h_c / 2 - (z_g))))


def d_by_time(time: float, time_max: float) -> float:
    """Retourne la distance de la charge en fonction du temps (charge variable)"""
    return min(d*time/time_max + 0.3, d)


# Autres paramètres calculés

m_tot = m_charge + m1 + m_grue  # masse totale [kg]
# Intensité des forces de poussée et de gravité [N]
f_poussee_gravite = m_tot * g
h_c = m_tot / ((l ** 2) * 1000)  # enfoncement de la plateforme [m]
z_g = (m1 * ((h1 / 2)-h_c) + m_charge *
       h_charge + m_grue * h_grue / 2) / m_tot  # hauteur du centre de gravité du système [m]
moment_inertie = calculate_moment_inertie()
angle_stable = calculate_stable_angle(d, m_charge)
angle_submersion = math.atan(2 * (h1 - h_c) / l)  # angle de submersion
angle_soulevement = math.atan(2 * h_c / l)  # angle de soulèvement

# Paramètres de la simulation

step = 0.001  # pas (dt) [s]
end = 20.0  # durée [s]
theta_0 = 0  # angle d'inclinaison du départ [rad]

# Initialisation des listes
t = np.arange(0, end, step)
x_g = np.empty_like(t)
theta = np.empty_like(t)
v_theta = np.empty_like(t)
a_theta = np.empty_like(t)
x_c = np.empty_like(t)


def get_x_c(theta: float) -> float:
    """Retourne la position du centre de poussée"""
    return (l ** 2 / (12 * h_c) - (h_c / 2)) * math.sin(theta)


def get_x_g(theta: float, d: float) -> float:
    """Retourne la position du centre de gravité"""
    return (math.sin(theta) * (z_g)) + ((math.cos(theta) * m_charge * d) / m_tot)


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

        # calcul des points d'application des forces
        x_g[i] = get_x_g(theta[i], d_by_time(dt*i, 10))
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
    plt.savefig("diagramme_phases_charge_variable.png")
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
    diagramme_des_phases()
