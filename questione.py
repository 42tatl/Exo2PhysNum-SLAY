import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os

def read_in_file(filename):
    variables = {}
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith("#"):  
                key, value = line.split("=")
                key = key.strip()
                value = value.split("//")[0].strip()  
                
                try:
                    if "." in value:
                        variables[key] = float(value)
                    else:
                        variables[key] = int(value)
                except ValueError:
                    print(f"Warning: Could not convert '{value}' to number. Check {filename}.")
    return variables

executable = './Exe2'  # Remove .exe for Mac
repertoire = r"/Users/lilimouelle/Desktop/PHYSNUM/exo_2/Exo2PhysNum"  # Modify for correct directory
os.chdir(repertoire)

input_filename = 'confige.in'  # Modify according to configuration

params = read_in_file(input_filename)

theta0 = params.get("theta0", 0.0)
thetadot0 = params.get("thetadot0", 0.0)  
mu = params.get("mu", 0.0) 
B0 = params.get("B0", 0.0) 
B1 = params.get("B1", 0.0)
m = params.get("m", 0.0)  
L = params.get("L", 0.0) 
nsteps = params.get("nsteps", 0)  # Number of steps
N_excit = params.get("N_excit", 0)  # Number of excitations

Omega = 2 * np.sqrt(12 * mu * B0 / (m * L**2))  
T = 2 * np.pi / Omega  # Period of excitation

#Parameters of the simulation


#theta0 = np.array([2, -1, 5, -3.2, 0.8, 6.5, -7.3, 4.1, -2.9, 1.7, 3.6, -5.4, 0.2, -8.9, 7.8, -6.1, 9.4, -0.5, 1.2, -4.7], dtype=float)
#theta_dot0 = np.array([12, -8, 0.2, -15, 5.5, 18.3, -3.1, 10.7, -20.5, 7.9, -6.2, 14.8, -11.4, 22.1, -9.8, 3.3, -13.6, 25.2, -19.1, 8.6], dtype=float)
theta0 = np.array([2, -1], dtype=float) #1
theta_dot0 = np.array([12, -8], dtype=float)
nsimul = len(theta0)



#Creation des fichiers de simulations automatiquement pour diff nsteps
outputs = []  
for i in range(nsimul):
   # Création du nom du fichier de sortie incluant les deux paramètres
    output_file = f"theta0={theta0[i]:.5f}_thetadot0={theta_dot0[i]:.5f}.out"
    outputs.append(output_file)
    # Construction de la commande avec les deux paramètres
    cmd = f"{executable} {input_filename} theta0={theta0[i]:.15g} thetadot0={theta_dot0[i]:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

t_all = []
theta_all = []
theta_dot_all = []
for i in range(nsimul):
    data = np.loadtxt(outputs[i])
    t_all.append(data[:, 0])  
    #theta_all.append(np.mod(data[:, 1], 2*np.pi)) # Modulo 2pi
    theta_all.append((data[:, 1] + np.pi) % (2*np.pi) - np.pi)
    theta_dot_all.append(data[:, 2])

'''
total_points = sum(len(theta) for theta in theta_all)
print(f"Nombre total de points sélectionnés: {total_points}")
'''

#plt.figure(figsize=(8, 6))
#plt.scatter(theta_all, theta_dot_all, s=10, color='b', alpha=0.7)
colors = plt.cm.viridis(np.linspace(0, 1, nsimul))  # Palette de couleurs

for i in range(nsimul):
    plt.scatter(theta_all[i], theta_dot_all[i], s=1, color=colors[i], alpha=0.7)
plt.xlabel(r'$\theta$ [rad]')
plt.ylabel(r'$\dot{\theta}$ [rad $\cdot$ s$^{-1}$]')
plt.grid()
plt.savefig("poincare_e.png")
plt.show()
