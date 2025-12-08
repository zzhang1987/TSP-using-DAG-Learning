# NeurIPS 2025 - Trace-Guided Cost Augmentation for ATSP

This repository contains the implementation of the NeurIPS 2025 paper **“Solving the Asymmetric Traveling Salesman Problem via Trace-Guided Cost Augmentation”** by **Zhen Zhang, Javen Qinfeng Shi, Wee Sun Lee**. It provides a Python interface built with **pybind11** and utilizes the **LKH‑3** heuristic as a baseline.

## Prerequisites

- **C++ compiler** (supports C++14 or newer)  
- **CMake** (≥ 3.15)  
- **Python** 3.12 (the project was built with this version)  
- **Git** (for submodule handling)  

> **Note:** The `third_party/pybind11` submodule is required for building the Python extension.  

## Build Steps

1. **Initialize submodules** – fetch the `pybind11` source:

   ```bash
   git submodule init
   git submodule update
   ```

2. **Compile the Python extension** – this builds the C++ code and creates the `partial_dag_model` module in‑place:

   ```bash
   python setup.py --build_ext --inplace
   ```

3. **Build and install LKH‑3** – the LKH‑3 solver is required for solving ATSP instances. Run the provided script:

   ```bash
   ./third_party/build_isntall_lkh3.sh
   ```

   The script configures, compiles, and installs LKH‑3 into the `third_party` directory.

## Usage Example

```python
python solve_atsp.py
```

## Project Structure

```
├─ csrc/                     # C++ source files
├─ third_party/
│   ├─ pybind11/             # Submodule – pybind11 headers
│   └─ build_isntall_lkh3.sh# Script to compile LKH‑3
├─ setup.py                  # Build script for the Python extension
├─ solve_atsp.py             # High‑level Python wrapper
└─ README.md                 # This file
```

## Troubleshooting

- **Missing compiler / CMake** – install via your package manager (e.g., `brew install cmake gcc` on macOS).
- **Submodule not found** – ensure you ran both `git submodule init` and `git submodule update`.
- **Build failures** – verify that the Python version matches the one used for `setup.py` and that the required development tools are in your `PATH`.

## Reference

If you use this code in your research, please cite the original paper:

```
@inproceedings{zhang2025trace,
  title={Solving the Asymmetric Traveling Salesman Problem via Trace-Guided Cost Augmentation},
  author={Zhang, Zhen and Shi, Javen Qinfeng and Lee, Wee Sun},
  booktitle={Advances in Neural Information Processing Systems},
  year={2025}
}
```

