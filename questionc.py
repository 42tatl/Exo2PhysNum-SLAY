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
repertoire = r"/Users/lilimouelle/Desktop/PHYSNUM/Exo2PhysNum-SLAY"  # Modify for correct directory
os.chdir(repertoire)

input_filename = 'configc.in'  # Modify according to configuration

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

nsimul = 50
theta0 = np.linspace(-np.pi, np.pi, nsimul) 
theta_dot0 = np.linspace(-10, 10, nsimul)
theta0 = np.append(theta0,1e-3)
theta_dot0 = np.append(theta_dot0,0)
nsimul += 1 

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
    theta_all.append((data[:, 1] + np.pi) % (2*np.pi) - np.pi)
    theta_dot_all.append(data[:, 2])


# Poincaré section
colors = plt.cm.rainbow(np.linspace(0, 1, nsimul))  

for i in range(nsimul):
    plt.scatter(theta_all[i], theta_dot_all[i], s=0.1,)
plt.xlabel(r'$\theta$ [rad]',fontsize=16)
plt.ylabel(r'$\dot{\theta}$ [rad $\cdot$ s$^{-1}$]',fontsize=16)
plt.xticks(fontsize=12)  
plt.yticks(fontsize=12)
plt.grid()
plt.tight_layout()
plt.savefig("poincare_c3.png",dpi=300)
plt.show()
