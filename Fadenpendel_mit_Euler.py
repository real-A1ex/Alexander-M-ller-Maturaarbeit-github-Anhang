import math
import numpy
import matplotlib.pyplot as plt
from Messdaten_Zeit import *
from Messdaten_Winkel_1 import *
from Messdaten_Winkel_2 import *
from Messdaten_Winkel_3 import *
from Messdaten_Winkel_4 import *
from Messdaten_Winkel_5 import *


# Parameter
länge = 0.636  # Kleine Kugel: 0.6255. Grosse Kugel: 0.636.
winkel_anfang = 27.5  # Kleine Kugel: 29°, 25°, 20°. Grosse Kugel: 30°, 27.5°.
zeit = 40 # Angegeben in Sekunden
zeitschritt = 0.0001 # Je kleiner, desto höher der Rechenaufwand und die Genauigkeit
g = 9.806 # Erdbeschleunigung am Ort der Experimentdurchführung
rho = 1.225 # Luftdichte am Ort der Experimentdurchführung
c_w = 0.47 # Strömungswiderstandskoeffizient einer Kugel
radius = 0.02 # Kleine Kugel: 0.0155. Grosse Kugel: 0.02.
masse = 0.2862  # Kleine Kugel: 0.1312. Grosse Kugel: 0.2862.

# Anfangswerte
winkel = math.radians(winkel_anfang)
omega = 0.0
omega_werte = []
winkel_werte = []
zeiten = numpy.arange(0, zeit, zeitschritt)
fläche = math.pi * radius**2

# Euler-Verfahren
for t in zeiten:
    F_grav = -g / länge * math.sin(winkel)
    
    v = omega * länge
    F_luft = -0.5 * rho * c_w * fläche * v * abs(v) / masse
    
    alpha = F_grav + F_luft
       
    winkel += omega * zeitschritt
    omega += alpha * zeitschritt 
    
    winkel_werte.append(math.degrees(winkel))
    omega_werte.append(omega)

# Diagramm erstellen
plt.plot(zeiten, winkel_werte, label="Simulation")
plt.plot(zeit_messung, winkel_5, label="Messungen Experiment") # Kleine Kugel: winkel_1: 29°, winkel_2: 25°, winkel_3: 20°.     Grosse Kugel: winkel_4: 30°, winkel_5: 27.5°.
plt.xlabel("Zeit (s)")
plt.ylabel("Winkel (Grad)")
plt.legend(loc="upper left")
plt.title("Simulation Pendelbewegung 4: Grosse Kugel, 27.5°")
plt.grid(True)
plt.show()
