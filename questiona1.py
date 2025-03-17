import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os

# Function to read config file
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

# Define executable path (make sure it exists)
repertoire = r"C:\Users\Avril\Desktop\Exo2PhysNum"  # Windows path format
executable = os.path.join(repertoire, "Exe.exe")

# Change directory to where the executable is located
os.chdir(repertoire)

# Define input config file
input_filename = "config.a"

# Read parameters from config file
params = read_in_file(input_filename)

theta0 = params.get("theta0", 0.0)
thetadot0 = params.get("thetadot0", 0.0)
mu = params.get("mu", 0.0)
B0 = params.get("B0", 0.0)
m = params.get("m", 0.0)
L = params.get("L", 0.0)

w0 = np.sqrt(12 * mu * B0 / (m * L**2))
A = theta0
B = thetadot0 / w0

def theta_a(t):
    return A * np.cos(w0 * t) + B * np.sin(w0 * t)

def theta_dot(t):
    return w0 * (B * np.cos(w0 * t) - A * np.sin(w0 * t))

# Simulation parameters
nsteps = np.array([10, 20, 50, 100], dtype=int)
nsimul = len(nsteps)
paramstr = "nsteps"
param = nsteps

# Create simulation output files
outputs = []
convergence_list = []

for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)

    # Corrected command format for Windows
    cmd = f'"{executable}" {input_filename} {paramstr}={param[i]:.15g} output={output_file}'

    # Debugging print
    print(f"\n📢 Running command: {cmd}")

    # Run the subprocess and capture output
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Print debugging information
    print("🟢 Command Output:", result.stdout)
    print("🔴 Command Error:", result.stderr)

    # Check if output file was created
    if os.path.exists(output_file):
        print(f"✅ SUCCESS: Output file '{output_file}' was created!")
    else:
        print(f"❌ ERROR: The output file '{output_file}' was NOT created!")

print("\n🔍 Checking if output files exist before reading...\n")

# Initialize the figure (✅ FIXED: Only one figure)
plt.figure(figsize=(8, 5))

# Read and plot results on the same figure
for i in range(nsimul):
    if os.path.exists(outputs[i]):  # Ensure the file exists before reading
        print(f"✅ Reading file: {outputs[i]}")
        data = np.loadtxt(outputs[i])
        t = data[:, 0]
        theta = data[:, 1]
        # Plot numerical solutions for each `nsteps`
        plt.plot(t, theta, label=rf"$\theta(t)$ for {param[i]} steps", linewidth=1)

# Add the analytical solution θ_a(t) as a dashed red line
t_a = np.linspace(0, max(t), 1000)
plt.plot(t_a, theta_a(t_a), label=r"$\theta_a(t)$", color="black", linestyle="--", linewidth=1)

# Labels and legend
plt.xlabel("t [s]")
plt.ylabel(r"$\theta$ [rad]")
plt.grid(True)
plt.legend()
plt.savefig("questiona1.png", dpi=300)
plt.show()
