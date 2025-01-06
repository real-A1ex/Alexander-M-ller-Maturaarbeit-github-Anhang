import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
from Messdaten_Zeit import *
from Messdaten_Winkel_1 import *
from Messdaten_Winkel_2 import *
from Messdaten_Winkel_3 import *
from Messdaten_Winkel_4 import *
from Messdaten_Winkel_5 import *

# Parameter
länge = 0.6255  # Kleine Kugel: 0.6255. Grosse Kugel: 0.636.
winkel_anfang = 30  # Kleine Kugel: 29°, 25°, 20°. Grosse Kugel: 30°, 27.5°.
zeit = 40
zeitschritt = 0.0001
g = 9.806
rho = 1.225
c_w = 0.47
radius = 0.02 # Kleine Kugel: 0.0155. Grosse Kugel: 0.02.
masse = 0.2862  # Kleine Kugel: 0.1312. Grosse Kugel: 0.2862.

# Anfangswerte
winkel = math.radians(winkel_anfang)
omega = 0.0
omega_werte = []
winkel_werte = []
zeiten = np.arange(0, zeit, zeitschritt)
fläche = math.pi * radius**2

# Funktion zur Berechnung der Beschleunigung
def berechne_omega(winkel, omega, beta):

    F_grav = -g / länge * math.sin(winkel)
    v = omega * länge
    F_luft = -0.5 * rho * c_w * fläche * v * abs(v) / masse
    F_reibung = -beta * omega  # Dämpfung durch Klemmstelle

    return F_grav + F_luft + F_reibung

# Runge-Kutta 4er Ordnung
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
def simulation_model(zeiten, beta):
    return simuliere_pendel(beta, zeiten, winkel_anfang)

# Optimierung der Dämpfungskonstante mit curve_fit
erstes_beta = 0.01  # Startwert für die Optimierung
optimales_beta, _ = curve_fit(simulation_model, zeit_messung, winkel_4, p0=[erstes_beta])

# Finales Ergebnis simulieren
finales_beta = optimales_beta[0]
simulation_mit_dämpfung = simuliere_pendel(finales_beta, zeiten, winkel_anfang)


# Diagramm erstellen
plt.plot(zeit_messung, winkel_4, label="Experimentelle Daten")
plt.plot(zeiten, simulation_mit_dämpfung, label=f"$\beta = {finales_beta}$")
plt.xlabel('Zeit (s)')
plt.ylabel('Winkel (Grad)')
plt.legend(loc="upper left")
plt.title('Simulation Pendelbewegung 4 mit Dämpfung')
plt.grid(True)
plt.show()