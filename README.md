# TDE-diploma-thesis
Package for numerical modelling of tidal disruption events.

## Repository structure

- **`TDESim/`**: The core Python package containing modularized functions.
  - `postnewtonian/`: Engine for simulation of TDEs in PN approximation (alias *pn*)
    - `pn_discrete_engine.py`: Simulations of trajectory represented by discrete ellipses with precessions
  - `utils/`: Other utilities used in the package
    - `statistics.py`: Statistical distributions and other specific functionalities used in MC simulations
    - `units.py`: Basic astrophysical constants, functions for unit conversion to/from geometrized units
  - `config.py`: Dataclass *SimulationParams* used for organizing parameters in simulations
  - `criteria.py`: Criteria used in simulations
  - `graph_3d.py`: Prepared 3D trajectory plotting utility (this is optional and can be removed, so that plotly package isn't needed
- **`01_Monte-Carlo-Guillochon.ipynb`**: Demo jupyter notebook containing first **Monte Carlo simulation following Guillochon et al. (2015)**

## Dependencies

To run the simulation notebook, you will need standard Python specific libraries: numpy, matplotlib, scipy and plotly (optional, see above).
