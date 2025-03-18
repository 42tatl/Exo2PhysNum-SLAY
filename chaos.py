import numpy as np
import os
import fonctions as fct  # Assuming this module contains `read_in_file` and `run_simulation`
import matplotlib.pyplot as plt

# Define paths
repertoire = r"C:\Users\Avril\Desktop\Exo2PhysNum-SLAY"
executable = os.path.join(repertoire, "Exe.exe")
os.chdir(repertoire)

# Define input config file
input_filename = "configb.in"

# Read parameters from config file
params = fct.read_in_file(input_filename)

# Extract parameters
theta0_a, thetadot0, mu, B0, m, L, w0 = fct.get_params(params)

# Function to run simulations for theta0 and theta0 + 1e-6
def run_sim(nsteps, ci_theta0, thetadot0):
    outputs = {}
    for n in nsteps:
        for theta0 in ci_theta0:
            perturbed_theta0 = theta0 + 1e-6  # Perturbation

            # Run simulation for theta0
            output_file_a, _ = fct.run_simulation(
                executable=executable,
                input_filename=input_filename,
                output_template="nsteps={nsteps}_theta0={theta0:.6f}.out",
                nsteps=n,
                theta0=theta0,
                thetadot0=thetadot0
            )

            # Run simulation for perturbed theta0
            output_file_b, _ = fct.run_simulation(
                executable=executable,
                input_filename=input_filename,
                output_template="nsteps={nsteps}_theta0={theta0:.6f}_perturbed.out",
                nsteps=n,
                theta0=perturbed_theta0,
                thetadot0=thetadot0
            )

            outputs[theta0] = (output_file_a, output_file_b)
    return outputs

# Function to compute the Lyapunov exponent for perturbed theta0 values
def compute_lyapunov_exponent(t, theta_a, theta_dot_a, theta_b, theta_dot_b, w0):
    delta = np.sqrt(w0**2 * (theta_b - theta_a)**2 + (theta_dot_b - theta_dot_a)**2)
    slope, intercept = np.polyfit(t, np.log(delta), 1)
    return slope

# Main function that runs simulations and computes Lyapunov exponents
def exp(nsteps, ci_theta0, thetadot0):

    output_files = run_sim(nsteps, ci_theta0, thetadot0)

    lyapunov_exponents = {}
    for theta0, (file_a, file_b) in output_files.items():
        if os.path.exists(file_a) and os.path.exists(file_b):
            print(f"✅ Reading files: {file_a} and {file_b}")
            data_a = np.loadtxt(file_a)
            data_b = np.loadtxt(file_b)

            t = data_a[:, 0]  # Time array
            theta_a = data_a[:, 1]  # Theta values for original
            theta_dot_a = data_a[:, 2]  # Theta dot values for original

            theta_b = data_b[:, 1]  # Theta values for perturbed
            theta_dot_b = data_b[:, 2]  # Theta dot values for perturbed

            lyapunov_exp = compute_lyapunov_exponent(
                t, theta_a, theta_dot_a, theta_b, theta_dot_b, w0
            )

            key = f"theta0={theta0:.6f}"
            lyapunov_exponents[key] = lyapunov_exp

    return lyapunov_exponents

nsteps = np.array([1000], dtype=int)
ci_theta0 = np.linspace(-np.pi, np.pi, 100)

lyapunov_results = exp(nsteps, ci_theta0, thetadot0)
values = list(lyapunov_results.values())
min_key = min(lyapunov_results, key=lyapunov_results.get)
max_key = max(lyapunov_results, key=lyapunov_results.get)

lyap_min = lyapunov_results[min_key]
lyap_max = lyapunov_results[max_key]

print(f"🔹 Minimum Lyapunov Exponent: {lyap_min:.6f} at {min_key}")
print(f"🔺 Maximum Lyapunov Exponent: {lyap_max:.6f} at {max_key}")

plt.plot(ci_theta0, values, '+')
plt.xlabel('theta0')
plt.ylabel('Lyapunov Exponent')
plt.grid()
plt.show()

