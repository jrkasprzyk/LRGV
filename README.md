# LRGV Project

This repository contains the source code and resources for the LRGV (Lower Rio Grande Valley) simulation and optimization project. The focus is on a command-line interface version of the model that accepts decision variables and writes outputs to std-out.

The original reference for the model is Kasprzyk et al. (2009) "Managing population and drought risks using many-objective water portfolio planning under uncertainty", Water Resources Research, https://doi.org/10.1029/2009WR008121 but this model also reflects additions added later.

## Structure

- `lrgv_src/` — C++ source code for the simulation model
- `lrgv_bin/` — Input data, Makefiles, and Python runner script

## Python Environment

Create your own virtual environment and install dependencies as follows:

1. Ensure you have Python 3.13.12 or newer installed.
2. (Optional but recommended) Create a virtual environment:
  ```sh
  python -m venv venv
  ```
3. Activate the virtual environment:
  - On Windows (cmd):
    ```sh
    venv\Scripts\activate.bat
    ```
  - On PowerShell:
    ```sh
    .\venv\Scripts\Activate.ps1
    ```
4. Install dependencies (none required except standard library):
  ```sh
  pip install -r requirements.txt
  ```

For the current Python demo, no external Python packages are required; only the standard library is used. The script was developed and tested with Python 3.13.12.

## C++ Build Instructions

To build the C++ code, two Makefiles are provided:

1. **On Windows (latest testing):**
  - Open a terminal in the repository root.
  - Run:
    ```sh
    make -f lrgv_bin/Makefile.msvc
    ```
  - This will build the simulation executable using the MSVC-compatible Makefile. This is the version most recently tested.

2. **On Linux/macOS (or with GNU Make):**
  - Open a terminal in the repository root.
  - Run:
    ```sh
    make -f lrgv_bin/MakefileSerial
    ```
  - This will build the simulation executable using the standard Makefile. This was originally tested when the repo was created in 2015.

## Running the Python Demo

- To run the simulation from Python, use:
  ```sh
  python lrgv_bin/run_from_python.py
  ```
- This demonstration calls the C++ executable and shows how input decision variables and output objective functions are communicated through console input and output.
- The executable should be compatible with [MOEAFramework](https://moeaframework.org/) and other tools that perform simulation-based optimization using standard input and output.

## Notes
- Input data files are in `lrgv_bin/`.
- C++ source code is in `lrgv_src/`.
- The project uses both C++ and Python for different components.