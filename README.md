# Bayesian predictive model for gas adsorption

This repository contains MATLAB code dedicated to modeling gas adsorption processes. The codebase leverages Bayesian analysis and various isotherm models to provide a comprehensive tool for understanding and predicting gas adsorption behavior. Key functionalities include the setup and execution of Bayesian analysis, forward modeling based on inferred parameters, and essential functions for solving linear systems and calculating isotherm properties.

## Detailed Description

### 1. Comprehensive Bayesian Analysis and Results Visualization for Gas Adsorption Modeling (`BPM_GasAdsorption_Main.m`)
The `BPM_GasAdsorption_Main.m` module is central to the Bayesian analysis of gas adsorption processes. It is designed to process experimental data meticulously, laying the groundwork for sophisticated analysis. Key steps within this module include:
- **Data Loading and Processing**: Efficiently manages the ingestion and preliminary treatment of experimental data, ensuring readiness for analysis.
- **Bayesian Analysis Framework**: Establishes a robust Bayesian inference framework. This involves setting up algorithmic structures for Bayesian analysis, including defining prior distributions and configuring options for Bayesian solvers.
- **Inversion Solver Configuration**: Carefully configures the Bayesian inversion solver and sampler, optimizing them for the specific needs of gas adsorption modeling.
- **Results Reporting and Visualization**:Finalizes the analysis process by reporting and visually representing results, focusing on the interpretation and understanding of the Bayesian analysis outcomes.

### 2. IUQ Procedure Module (`uq_GasAdsorption.m`)
This function embodies the essence of the Inverse Uncertainty Quantification (IUQ) method, a critical aspect of modern parameter estimation techniques in physical adsorption models. It is designed to:
- **Handle Multiple Data Sets**: Processes various datasets simultaneously, ensuring a comprehensive approach to parameter estimation.
- **Estimate Model Parameters**: Utilizes Bayesian analysis to precisely estimate parameters within the gas adsorption model, enhancing both accuracy and reliability.

### 3. Forward Modeling for Gas Physical Adsorption (`Model_GasAdsorption.m`)
The `Adsorption_Model` function is a testament to the detailed simulation capabilities of the tool. It undertakes forward modeling for gas physical adsorption by:
- **Employing Various Isotherm Models**: Adapts to different isotherm models, providing versatility in modeling the adsorption process.
- **Solving with Implicit Solver**: Utilizes an advanced implicit solver with iterative linearization using the Newton method, ensuring rapid and precise simulation results even in complex 
    scenarios.

### 4. Supporting Functions
These functions form the backbone of the computational processes within the tool:
- **Thomas Algorithm (`thomas`)**: A robust solution for tridiagonal linear systems, pivotal in numerical methods and simulations.
- **Langmuir Isotherm (`Isotherm_Langmuir`)**: Calculates adsorption quantities, implementing the renowned Langmuir isotherm model for single-layer adsorption.
- **Derivative of Langmuir Isotherm (`Deriv_Langmuir`)**: Provides the derivative calculations of the Langmuir isotherm model, essential for understanding the dynamics of adsorption processes.

## Authors

- **Yesol Hyun** - School of Mathematics and Computing (Computational Science and Engineering), Yonsei University - yesol2@yonsei.ac.kr
- **Geunwoo Oh** - School of Mathematics and Computing (Computational Science and Engineering), Yonsei University - gwoh@yonsei.ac.kr
- **Jung-Il Choi** - School of Mathematics and Computing (Computational Science and Engineering), Yonsei University - jic@yonsei.ac.kr


## Installation

Clone the repository to your local machine using:

```bash
git clone https://github.com/MPMC-Lab/BPM_GasAdsorption.git
```

Alternatively, the source files can be downloaded through github menu 'Download ZIP'.

## Contributing

Contributions to this project are welcome. Please fork the repository and submit a pull request.

## Citation
Please cite the following paper if BPM_GasAdsorption is employed in your scholarly research or project implementation. This citation acknowledges the original authors and supports the ongoing academic collaboration and recognition of their work:

```bibtex
@article{hyun2024bayesian,
  title={Bayesian predictive modeling for gas purification using breakthrough curves},
  author={Hyun, Yesol and Oh, Geunwoo and Lee, Jaeheon and Jung, Heesoo and Kim, Min-Kun and Choi, Jung-Il},
  journal={Journal of Hazardous Materials},
  pages={134311},
  year={2024},
  publisher={Elsevier}
}
```
## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

Special thanks to the contributors and researchers who have provided insights and feedback on this project.

## Contact Information

Should you require assistance or have any inquiries, kindly contact Dr. Jung-Il Choi via email at [jic@yonsei.ac.kr](mailto:jic@yonsei.ac.kr).

For detailed information and further reading related to BPM_GasAdsorption, you are encouraged to consult the reference paper. Additional insights and resources are available at School of Mathematics and Computing, within the domain of Computational Science and Engineering at Yonsei University. For more details, please visit our website: [mpmc.yonsei.ac.kr](http://mpmc.yonsei.ac.kr).

