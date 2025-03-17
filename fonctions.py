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

def get_params(params):
    theta0 = params.get("theta0", 0.0)
    thetadot0 = params.get("thetadot0", 0.0)
    mu = params.get("mu", 0.0)
    B0 = params.get("B0", 0.0)
    m = params.get("m", 0.0)
    L = params.get("L", 0.0)
    w0 = np.sqrt(12 * mu * B0 / (m * L**2))
    return theta0, thetadot0, mu, B0, m, L, w0

def run_simulation(executable, input_filename, output_template, **params):
    """
    Runs a simulation with an arbitrary number of parameters.

    Parameters:
    - executable (str): Path to the executable file.
    - input_filename (str): Input configuration file name.
    - output_template (str): Template string for the output filename.
    - **params: Arbitrary number of named parameters for the command.

    Returns:
    - output_file (str): The generated output filename.
    - result (subprocess.CompletedProcess): The result of the subprocess execution.
    """
    # Create output filename dynamically based on parameters
    output_filename = output_template.format(**params)
    
    # Construct command dynamically
    param_str = " ".join(f"{key}={value:.15g}" for key, value in params.items())
    cmd = f'"{executable}" {input_filename} {param_str} output={output_filename}'

    # Debugging print
    print(f"\n📢 Running command: {cmd}")

    # Run the command
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Print command outputs
    print("🟢 Command Output:", result.stdout)
    print("🔴 Command Error:", result.stderr)

    # Check if the output file was created
    if os.path.exists(output_filename):
        print(f"✅ SUCCESS: Output file '{output_filename}' was created!")
    else:
        print(f"❌ ERROR: The output file '{output_filename}' was NOT created!")

    return output_filename, result