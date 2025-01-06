import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
from Messdaten_Zeit import *
from Messdaten_Winkel_Experiment_Magnet import *

# Parameter
länge = 0.6255  # Kleine Kugel: 0.6255
winkel_anfang = 33.8
zeit = 40
zeitschritt = 0.0001
g = 9.806
rho = 1.225
c_w = 0.47
radius = 0.0155 # Kleine Kugel: 0.0155.
masse = 0.1312 # Kleine Kugel: 0.1312.

# Magnetische Parameter
x_magnet = 63.8  # Position des Magneten entlang der x-Achse
y_magnet = -51.3  # Position des Magneten entlang der y-Achse

# Anfangswerte
winkel = math.radians(winkel_anfang)
omega = 0.0
omega_werte = []
winkel_werte = []
zeiten = np.arange(0, zeit, zeitschritt)
fläche = math.pi * radius**3

# Hilfsfunktion: Berechnung der Winkelbeschleunigung (omega)
def berechne_omega(winkel, omega, k_magnet):
        
    F_grav = -g / länge * math.sin(winkel)
    v = omega * länge
    F_luft = -0.5 * rho * c_w * fläche * v * abs(v) / masse

    x_pendel = länge * math.sin(winkel)
    y_pendel = -länge * math.cos(winkel)

    dx = x_pendel - x_magnet
    dy = y_pendel - y_magnet
    abstand_magnet_pendel = math.sqrt(dx**2 + dy**2)

    F_mag = k_magnet / abstand_magnet_pendel**2

    return F_grav + F_luft + F_mag

def simuliere_pendel(beta, experiment_zeiten, winkel_anfang):
    
    winkel = math.radians(winkel_anfang)
    omega = 0.0  # Anfangswinkelgeschwindigkeit (rad/s)
    omega_werte = []
    winkel_werte = []

    zeiten = np.arange(0, experiment_zeiten[-1] + zeitschritt, zeitschritt)
    for t in zeiten:
        # Runge-Kutta 4. Ordnung
        k1_omega = zeitschritt * berechne_omega(winkel, omega, beta)
        k1_winkel = zeitschritt * omega

        k2_omega = zeitschritt * berechne_omega(winkel + k1_winkel / 2, omega + k1_omega / 2, beta)
        k2_winkel = zeitschritt * (omega + k1_omega / 2)

        k3_omega = zeitschritt * berechne_omega(winkel + k2_winkel / 2, omega + k2_omega / 2, beta)
        k3_winkel = zeitschritt * (omega + k2_omega / 2)

        k4_omega = zeitschritt * berechne_omega(winkel + k3_winkel, omega + k3_omega, beta)
        k4_winkel = zeitschritt * (omega + k3_omega)

        omega += (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega) / 6
        winkel += (k1_winkel + 2 * k2_winkel + 2 * k3_winkel + k4_winkel) / 6

        omega_werte.append(omega)
        winkel_werte.append(math.degrees(winkel))

    return np.interp(experiment_zeiten, zeiten, winkel_werte)

# Fehlerfunktion für curve_fit
def simulation_model(zeiten, k_magnet):
    return simuliere_pendel(k_magnet, zeiten, winkel_anfang)

# Optimierung der magnetischen Kosntante mit curve_fit
erstes_k_magnet = 0.1 # Startwert für die Optimierung
optimales_k_magnet, _ = curve_fit(simulation_model, zeit_messung, winkel_magnet, p0=[erstes_k_magnet])

# Finales Ergebnis simulieren
finales_k_magnet = round(optimales_k_magnet[0], 2)
simulation_mit_magnet = simuliere_pendel(finales_k_magnet, zeiten, winkel_anfang)



# Diagramm erstellen
plt.plot(zeit_messung, winkel_magnet, label="Experimentelle Daten")
plt.plot(zeiten, simulation_mit_magnet, label=fr"$k_{{\mathit{{mag}}}} = {finales_k_magnet}$")
plt.xlabel('Zeit (s)')
plt.ylabel('Winkel (Grad)')
plt.legend(loc="upper left")
plt.title('Simulation Pendelbewegung mit Magnet')
plt.grid(True)
plt.show()