import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from scipy.stats import linregress

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

repertoire = r"C:\Users\Avril\Desktop\Exo2PhysNum-SLAY"  
executable = os.path.join(repertoire, "Exe.exe")  

os.chdir(repertoire)

input_filename = "config.a"

params = read_in_file(input_filename)

theta0 = params.get("theta0", 0.0)
thetadot0 = params.get("thetadot0", 0.0)
mu = params.get("mu", 0.0)
B0 = params.get("B0", 0.0)
m = params.get("m", 0.0)
L = params.get("L", 0.0)

nsteps = np.array([1000, 2500, 5000, 10000, 15000, 20000, 3e4, 35000, 4e4, 5e4], dtype=int)  
nsimul = len(nsteps)
paramstr = "nsteps"
param = nsteps

outputs = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)

    cmd = f'"{executable}" {input_filename} {paramstr}={param[i]:.15g} output={output_file}'

    print(f"\n Running command: {cmd}")

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    print(" Command Output:", result.stdout)
    print(" Command Error:", result.stderr)

    if os.path.exists(output_file):
        print(f"✅ SUCCESS: Output file '{output_file}' was created!")
    else:
        print(f"❌ ERROR: The output file '{output_file}' was NOT created!")

print("\n🔍 Checking if output files exist before reading...\n")

tfins = []
theta_fins = []
theta_dot_fins = []
for i in range(nsimul):
    if os.path.exists(outputs[i]):  
        print(f"✅ Reading file: {outputs[i]}")
        data = np.loadtxt(outputs[i])
        t = data[:, 0]
        tfins.append(t[-1])
        theta_fins.append(data[-1, 1])
        theta_dot_fins.append(data[-1, 2])
    else:
        print(f"❌ Skipping missing file: {outputs[i]}")

deltas = []
delta_t = []
w0 = np.sqrt(12 * mu * B0 / (m * L**2))
A = theta0
B = thetadot0 / w0

def theta(t):
    return A * np.cos(w0 * t) + B * np.sin(w0 * t)

def theta_dot(t):
    return w0 * (B * np.cos(w0 * t) - A * np.sin(w0 * t))

for simul in range(nsimul):
    theta_tfin = theta_fins[simul]
    theta_dot_tfin = theta_dot_fins[simul]
    delta_t.append(tfins[simul] / nsteps[simul])
    delta = np.sqrt(w0**2 * (theta_tfin - theta(tfins[simul]))**2 + (theta_dot_tfin - theta_dot(tfins[simul]))**2)
    deltas.append(delta)
    print(f"Final error: {delta:.15g} for nsteps = {param[simul]:.15g}")

ref_x = np.linspace(min(delta_t), max(delta_t), 100)
ref_y = ref_x**2  

plt.figure(figsize=(8, 5))
plt.loglog(delta_t, deltas, marker='o', linestyle='-', color='b', label=r"Error on $\theta$ and $\dot\theta$")
plt.loglog(ref_x, ref_y, linestyle='--', color='r', label=r"$\sim \Delta t^2$")

plt.xlabel(r"$\Delta t$ [s]")
plt.ylabel(r"$\delta$ [rad/s]")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend()

plt.savefig("questiona2.png", dpi=300)
plt.show()
