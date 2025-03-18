import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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

executable = './Exe'  # Remove .exe for Mac
repertoire = r"/Users/lilimouelle/Desktop/PHYSNUM/exo_2/Exo2PhysNum"  # Modify for correct directory
os.chdir(repertoire)

input_filename = 'configb.in'  # Modify according to configuration

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

#Parametres de simulation 
#nsteps = np.array([1000], dtype=int)  #first part of question b
#nsteps = np.array([30,35,40,50,60,70,80,90,100,300,500,1000,3000], dtype=int)
nsteps = np.array([100,150,200,300,500,1000,2500], dtype=int)
nsimul = len(nsteps) 
paramstr = 'nsteps'
param = nsteps

#Creation des fichiers de simulations automatiquement pour diff nsteps
outputs = []  
convergence_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    cmd = f"{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

T = 2 * np.pi / Omega  # Period of excitation
#Second part of question b (error of theta_fin)
tfins = []
theta_fins = []
delta_ts = []
delta_ts_squared = []
delta_ts_third = []
for i in range(nsimul):
    with open(outputs[i], 'r') as file:
         data = np.loadtxt(file)
    t = data[:, 0]
    tfins.append(t[-1])
    theta_fins.append(data[-1, 1])
    delta_ts.append(T/nsteps[i])
    delta_ts_squared.append((T/nsteps[i])**2)
    delta_ts_third.append((T/nsteps[i])**3)


plt.plot(delta_ts_squared, theta_fins, marker='o', linestyle='-', color='b', label='Error vs. Steps')
plt.xlabel(r'$(\Delta t)^2$')
plt.ylabel(r'$\theta_{fin}$')
plt.grid(True)
plt.legend()
plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
plt.ticklabel_format(style='sci', axis='x', scilimits=(-4, 4)) 
plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
plt.ticklabel_format(style='sci', axis='y', scilimits=(-4, 4)) 
plt.tight_layout()
plt.savefig("error_theta_fin.pdf")
plt.show()


'''
#First part of question b
try:
    data = np.loadtxt(output_file)
except ValueError:
    print(f"Erreur : le fichier '{output_file}' contient des données non numériques.")
    with open(output_file, "r") as file:
        print(file.readlines()[:10])  # Affiche les 10 premières lignes pour debug
    exit()

t = data[:, 0]
theta = data[:, 1]
theta_dot = data[:, 2]
emec = data[:, 3]
power = data[:, 4]

#theta as a function on time
plt.plot(t, theta, color='b')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$\theta$ [rad]')
plt.grid()
plt.savefig("theta_b.pdf")
plt.show()

#theta_dot as a function on time
plt.plot(theta, theta_dot, color='b')
plt.xlabel(r'$\theta$ [rad]')
plt.ylabel(r'$\dot{\theta}$ [rad $\cdot$ s$^{-1}$]')
plt.grid()
plt.savefig("phase_space_b.pdf")
plt.show()

#Emec as a function on time
plt.plot(t, emec, color='g')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$E_{mec}$ [J]')
plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
plt.ticklabel_format(style='sci', axis='y', scilimits=(-3, 3)) 
plt.grid()
plt.tight_layout()
plt.savefig("emec_b.pdf")
plt.show()

#Time derivative of mechanical energy
theta_dotdot = -12*mu/(m*pow(L,2))*(B0+B1*np.sin(Omega*t))*np.sin(theta)
dEmec_dt = 1.0/12.0*m*pow(L,2)*theta_dotdot*theta_dot+mu*B0*theta_dot*np.sin(theta)

# dEmec_dt as a function on time
plt.plot(t, dEmec_dt, color='r')
#plt.plot(t, power, label=r'$P_{nc}$', color='b')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$\frac{dE_{mec}}{dt}$ [W]')
plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
plt.ticklabel_format(style='sci', axis='y', scilimits=(-3, 3)) 
plt.grid()
plt.tight_layout()
plt.savefig("dEmec_dt.pdf")
plt.show()

# Power as a function on time
plt.plot(t, power, color='r')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$P_{nc}$ [W]')
plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
plt.ticklabel_format(style='sci', axis='y', scilimits=(-3, 3)) 
plt.grid()
plt.tight_layout()
plt.savefig("Power.pdf")
plt.show()

diff = np.abs(dEmec_dt - power)
#diff = dEmec_dt - power

# Difference between dEmec_dt and power as a function on time
plt.plot(t, diff, color='r')
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$\left|E_{mec} - P_{nc}\right|$ [W]')
plt.grid()
plt.savefig("diff.pdf")
plt.show()
'''
