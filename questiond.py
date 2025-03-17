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

repertoire = r"C:\Users\Avril\Desktop\Exo2PhysNum"  # Windows path format
executable = os.path.join(repertoire, "Exe.exe")
os.chdir(repertoire)

input_filename = "configb.in"

params = read_in_file(input_filename)

mu = params.get("mu", 0.0)
B0 = params.get("B0", 0.0)
m = params.get("m", 0.0)
L = params.get("L", 0.0)
w0 = np.sqrt(12 * mu * B0 / (m * L**2))

nsteps = np.array([2000], dtype=int)
ci = np.array([0, 1e-6, np.pi / 2 , np.pi / 2 + 1e-6], dtype=float)
nsimul = len(nsteps) * len(ci)

paramstr1 = "nsteps"
paramstr2 = "theta0"

outputs = []

for n in nsteps:
    for theta0 in ci:
        output_file = f"{paramstr1}={n}_theta0={theta0:.7f}.out"
        outputs.append(output_file)

        cmd = f'"{executable}" {input_filename} {paramstr1}={n} {paramstr2}={theta0:.15g} output={output_file}'

        print(f"\n Running command: {cmd}")

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if os.path.exists(output_file):
            print(f"SUCCESS: Output file '{output_file}' was created!")
        else:
            print(f"ERROR: The output file '{output_file}' was NOT created!")

plt.figure(figsize=(8, 5))
thetas = []
thetas_dot = []

for i in range(nsimul):
    if os.path.exists(outputs[i]):  
        print(f"Reading file: {outputs[i]}")
        data = np.loadtxt(outputs[i])
        print(data.shape)
        t = data[:, 0]
        thetas.append(data[:, 1])
        thetas_dot.append(data[:, 2])

thetas = np.column_stack(thetas)  # Shape: (n_time_steps, n_simul)
thetas_dot = np.column_stack(thetas_dot)  # Shape: (n_time_steps, n_simul)
for i in range(2):
    delta = np.sqrt(w0**2*(thetas[:,i+1] - thetas[:,i])**2 + (thetas_dot[:,i+1] - thetas_dot[:,i])**2)
    print(f"Max error for simulation {i+1}: {np.max(delta)}")
    plt.plot(t, delta)
    plt.xlabel('t [s]')
    plt.ylabel(r'$\delta_{ab}$')
    plt.legend()
    plt.grid()
    plt.savefig(f'error_{i+1}.png')
    plt.show()


