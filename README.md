# Gas Adsorption Model using Bayesian Analysis

This repository contains MATLAB code for modeling gas adsorption processes using Bayesian analysis and various isotherm models. The code includes functions for setting up and executing Bayesian analysis, forward modeling based on inferred parameters, and supporting functions for solving linear systems and calculating isotherm properties.

## Description

The code is divided into several key components:

### 1. Experimental Data Preparation and Bayesian Analysis Setup (`PI_GasAdsorption.m`)
This module processes experimental data and sets up the Bayesian analysis framework for gas adsorption models. It includes data loading and processing, Bayesian analysis setup, defining prior distributions, and configuring the Bayesian inversion solver and sampler.

### 2. IUQ Procedure Module (`uq_PhysicalAdsorption`)
This function executes the Inverse Uncertainty Quantification (IUQ) procedure for parameter estimation in physical adsorption models. It processes multiple data sets to estimate adsorption model parameters using Bayesian analysis.

### 3. Forward Modeling for Gas Physical Adsorption (`Adsorption_Model`)
This function performs forward modeling for gas physical adsorption. It simulates the adsorption process using various isotherm models and solves the system using an implicit solver with adaptive time stepping.

### 4. Supporting Functions
- **Thomas Algorithm (`thomas`)**: Solves tridiagonal linear systems using the Thomas algorithm.
- **Langmuir Isotherm (`Isotherm_Langmuir`)**: Calculates the adsorption quantity based on the Langmuir isotherm model.
- **Derivative of Langmuir Isotherm (`Deriv_Langmuir`)**: Computes the derivative of the Langmuir isotherm model.

## Authors

- **Yesol Hyun** - School of Mathematics and Computing (Computational Science and Engineering), Yonsei University - yesol2@yonsei.ac.kr
- **Geunwoo Oh** - School of Mathematics and Computing (Computational Science and Engineering), Yonsei University - gwoh@yonsei.ac.kr
- **Jung-Il Choi** - School of Mathematics and Computing (Computational Science and Engineering), Yonsei University - jic@yonsei.ac.kr


## Installation

Clone the repository to your local machine using:

```bash
git clone https://github.com/MPMC-Lab/1dGasAdsorption_PI_repo.git

## Contributing

Contributions to this project are welcome. Please fork the repository and submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

Special thanks to the contributors and researchers who have provided insights and feedback on this project.

## References

For more information, please refer to the reference paper(s) associated with this project, and for further academic context, visit the School of Mathematics and Computing (Computational Science and Engineering) at Yonsei University.

