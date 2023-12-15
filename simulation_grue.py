# Simulation d'une grue sur une barge flottante
import math
import matplotlib.pyplot as plt
import numpy as np

# Constantes
g = 9.81  # gravitation [m/s**2]

# Paramètres du système
m_charge = 0.2  # mass de la charge déplacée [kg]
h_charge = 0.3  # hauteur de la charge au bas de la barge [m]
d_max = 2  # distance de la charge au centre de la barge [m]

m_grue = 4  # masse du chargement sur la barge [kg]
h_grue = 0.5  # hauteur du centre de gravité du chargement sur la barge [m]

l = 0.6  # longueur de la barge [m]
h1 = 0.08  # hauteur de la barge [m]
m1 = 5  # masse de la barge [kg]

D = 0.4  # coefficient de frottement visqueux


def calculate_moment_inertie():
    """Retourne le moment d'inertie de la plateforme"""
    # moment d'un parallélipipède (plateforme)
    moment_barge = (m1 * (l ** 2 + h1 ** 2)) / 12
    # moment d'une masse ponctuelle (grue)
    moment_caisse = m_grue * (h_grue - h_c) ** 2
    # moment d'une masse ponctuelle (charge)
    moment_charge = m_charge * (d_max ** 2 + (h_charge - h_c) ** 2)

    return moment_barge + moment_caisse + moment_charge  # somme des moments


def calculate_stable_angle(distance: float, masse_charge: float) -> float:
    """Retourne l'angle de stabilité théorique"""
    return math.atan((masse_charge * distance) / (m_tot * ((l ** 2) / (12 * h_c) + h_c / 2 - (z_g))))


def d_by_time(time: float, time_max: float) -> float:
    """Retourne la distance de la charge en fonction du temps (charge variable)"""
    return min(d_max*time/time_max + 0.2, d_max)


# Autres paramètres calculés

m_tot = m_charge + m1 + m_grue  # masse totale [kg]
# Intensité des forces de poussée et de gravité [N]
f_poussee_gravite = m_tot * g
h_c = m_tot / ((l ** 2) * 1000)  # enfoncement de la plateforme [m]
z_g = (m1 * ((h1 / 2)-h_c) + m_charge *
       h_charge + m_grue * h_grue / 2) / m_tot  # hauteur du centre de gravité du système [m]
moment_inertie = calculate_moment_inertie()
angle_stable = calculate_stable_angle(d_max, m_charge)
angle_submersion = math.atan(2 * (h1 - h_c) / l)  # angle de submersion
angle_soulevement = math.atan(2 * h_c / l)  # angle de soulèvement


def get_x_c(theta: float) -> float:
    """Retourne la position du centre de poussée"""
    return (l ** 2 / (12 * h_c) - (h_c / 2)) * math.sin(theta)


def get_x_g(theta: float, d: float) -> float:
    """Retourne la position du centre de gravité"""
    return (math.sin(theta) * (z_g)) + ((math.cos(theta) * m_charge * d) / m_tot)

# Paramètres de la simulation


step = 0.001  # pas (dt) [s]
end = 20.0  # durée [s]
theta_0 = 0  # angle d'inclinaison du départ [rad]


def simulation(charge_variable=False):
    """
    pre: 
    post: exécute une simulation jusqu'à t=end par pas de dt=step.
                                      Remplit les listes theta, v_theta, a_theta des angles, vitesses et accélérations angulaires.
    """
    # conditions initiales
    t = np.arange(0, end, step)
    x_g = np.empty_like(t)
    theta = np.empty_like(t)
    v_theta = np.empty_like(t)
    a_theta = np.empty_like(t)
    x_c = np.empty_like(t)
    theta[0] = theta_0
    for i in range(len(t) - 1):
        dt = step

        # calcul des points d'application des forces
        d = d_by_time(dt*i, 10) if charge_variable else d_max
        x_g[i] = get_x_g(theta[i], d)
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
    return t, theta, v_theta, a_theta


def get_donnees_tracker(filename):
    """Retourne les données Tracker du fichier spécifié"""
    f = open(filename)
    data = f.read()

    # Extraire les données en listes
    temps, position_y = [], []
    for ligne in data.split('\n'):
        if ligne:
            t, y = map(float, ligne.replace(',', '.').split())
            temps.append(t)
            position_y.append(y)
    # Convertir les listes en tableaux NumPy pour une manipulation plus aisée
    temps = np.array(temps)
    position_y = np.array(position_y)

    # Effectuez une régression polynomiale d'ordre 2 (ajustez cet ordre selon vos besoins)
    coefficients = np.polyfit(temps, position_y, len(temps)-2)

    # Créez une fonction polynomiale à partir des coefficients
    polynome = np.poly1d(coefficients)

    # Générez des valeurs de temps pour la courbe d'approximation
    temps_approx = np.linspace(temps.min(), temps.max(), 100)

    # Calculez les valeurs de position y correspondantes à la courbe d'approximation
    position_y_approx = polynome(temps_approx)-0.422
    cos_angle = position_y_approx/0.3
    angle_approx = np.arccos(np.clip(cos_angle, -1.0, 1.0))-1
    angles_en_degres = (np.degrees(angle_approx) - 90)
    return temps_approx, angle_approx


def get_periode(array):
    array = array[400:]
    min_array = min(array)
    max_array = max(array)
    return abs(np.where(array == max_array)[0][0]-np.where(array == min_array)[0][0])


def get_amp(array):
    min_array = min(array)
    max_array = max(array)
    return max_array-min_array


def graphique_angle_temps(title="angle_temps_charge_variable_simple.png", draw_lines=True, donnees_tracker=(), charge_variable=False):
    t, theta, v_theta, a_theta = simulation(charge_variable=charge_variable)
    plt.figure(1)
    if not draw_lines:
        plt.subplot(3, 1, 1)
    plt.plot(t, theta, label="Angle de la plateforme")
    plt.xlabel("Temps (s)")
    plt.ylabel("Angle (rad)")
    # ligne de stabilité
    plt.axhline(y=angle_stable, label="Stabilité théorique",
                color="blue", linestyle="dotted")
    if donnees_tracker and draw_lines:
        plt.plot(donnees_tracker[0], donnees_tracker[1],
                 label='Angle approximé en rad')
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
    if not draw_lines:
        plt.subplot(3, 1, 2)
        plt.plot(t, v_theta, label="Vitesse angulaire de la plateforme")
        plt.legend()
        plt.subplot(3, 1, 3)
        plt.plot(t, a_theta, label="Accélération angulaire de la plateforme")
        plt.legend()
    plt.savefig(title)
    plt.show()


def diagramme_des_phases(charge_variable=False):
    t, theta, v_theta, a_theta = simulation(charge_variable=charge_variable)
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
    temps_approx, angle_approx = get_donnees_tracker("tracker.txt")
    print(get_amp(angle_approx))
    for x in range(13, 14):
        D = x/10 + 0.1
        t, theta, v_theta, a_theta = simulation(charge_variable=False)
        print(D, get_periode(theta)*step, get_amp(theta))

        graphique_angle_temps(title=str(x),
                              draw_lines=True, donnees_tracker=(temps_approx, angle_approx-angle_approx[0]))
