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

repertoire = r"C:\Users\Avril\Desktop\Exo2PhysNum-SLAY" 
executable = os.path.join(repertoire, "Exe.exe")
os.chdir(repertoire)

input_filename = "configb.in"

params = read_in_file(input_filename)

mu = params.get("mu")
B0 = params.get("B0")
m = params.get("m")
L = params.get("L")
w0 = np.sqrt(12 * mu * B0 / (m * L**2))

nsteps = np.array([1000], dtype=int)
ci = np.array([
    [1.808796, 0],
    [1.808797, 0],
    [-2.697327, 0],
    [-2.697326,0]
], dtype=float)

nsimul = len(nsteps) * len(ci)

paramstr1 = "nsteps"
paramstr2 = "theta0"
paramstr3 = "thetadot0"  

outputs = []

for n in nsteps:
    for theta0, thetadot0 in ci:
        output_file = f"{paramstr1}={n}_theta0={theta0:.10f}_thetadot0={thetadot0:.10f}.out"
        outputs.append(output_file)

        cmd = (
            f'"{executable}" {input_filename} {paramstr1}={n} {paramstr2}={theta0:.15g} {paramstr3}={thetadot0:.15g} output={output_file}'
        )

        print(f"\n Running command: {cmd}")

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if os.path.exists(output_file):
            print(f"SUCCESS: Output file '{output_file}' was created!")
        else:
            print(f"ERROR: The output file '{output_file}' was NOT created!")

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

thetas = np.column_stack(thetas)  
thetas_dot = np.column_stack(thetas_dot)  

plt.figure(figsize=(8, 5))

delta1 = np.sqrt(w0**2 * (thetas[:, 1] - thetas[:, 0])**2 + (thetas_dot[:, 1] - thetas_dot[:, 0])**2)
delta2 = np.sqrt(w0**2 * (thetas[:, 3] - thetas[:, 2])**2 + (thetas_dot[:, 3] - thetas_dot[:, 2])**2)

slope1, intercept1 = np.polyfit(t, np.log(delta1), 1)  
slope2, intercept2 = np.polyfit(t, np.log(delta2), 1)

plt.plot(t, delta1, label=rf'$\delta_{{\theta_1}}, \ \lambda = {slope1:.3f}$')
plt.plot(t, delta2, label=rf'$\delta_{{\theta_2}}, \ \lambda = {slope2:.3f}$')
plt.xlabel('t [s]')
plt.ylabel(r'$\delta_{ab}$')
plt.legend()
plt.grid()
plt.savefig('error.png')
plt.show()

plt.figure(figsize=(8, 5))
plt.semilogy(t, delta2, label=r'$\delta_{\theta_2}$')
plt.semilogy(t, np.exp(slope2 * t + intercept2), '--', label=rf'Fit : $\delta_{{\theta_2}}, \ \lambda = {slope2:.3f}$')
plt.xlabel('t [s]')
plt.ylabel(r'$ln(\delta_{ab})$')
plt.legend()
plt.grid()
plt.savefig('error1.png')
plt.show()