# G4_NaI - NaI Detector Background Simulation

[![Geant4](https://img.shields.io/badge/Simulation-GEANT4-009688?logo=gnu)](https://geant4.web.cern.ch)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Geant4 simulation for characterizing internal background in Scionix NaI(Tl) detectors used in nuclear astrophysics experiments at the [LUNA (Laboratory for Underground Nuclear Astrophysics)](https://luna.lngs.infn.it) laboratory.

## Project Overview

This simulation package was developed to:
- Characterize internal background radiation in Scionix NaI(Tl) detectors
- Identify and quantify intrinsic radioimpurities in detector components
- Optimize experimental setup for underground measurements at LUNA

Key features:
- Detailed detector geometry implementation
- Radioactive decay simulations of common contaminants (K-40, Th/U chains, etc.)
- Particle tracking and energy deposition calculations

## Prerequisites

- [GEANT4](https://geant4.web.cern.ch) (version 10.7 or newer)
- CMake (version 3.16+)
- C++17 compatible compiler
- [ROOT Data Analysis Framework](https://root.cern) (recommended for analysis)

## Installation & Setup

1. Clone repository:
```bash
git clone https://github.com/skowrons94/G4_NaI.git
cd G4_NaI
```

2. Build with CMake:
```bash
mkdir build
cd build
cmake -DGeant4_DIR=/path/to/geant4/install ..
make -j4
```

## Running the Simulation

Execute from build directory:
```bash
./G4_NaI mac/[macro-file.mac]
```

In alternative, run the ```run.py``` script to simulate all the different radioisotope contributions inside the NaI detector.

Example macros:
- `vis.mac` - Visualization mode
