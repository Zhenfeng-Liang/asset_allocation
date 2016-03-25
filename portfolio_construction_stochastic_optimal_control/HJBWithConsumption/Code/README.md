## General description:

This is a MATLAB program for stochastic optimal control problem. This folder includes an approximation solver for HJB equation and a Monte Carlo back test program with constraints.


## Code files description:

### runPaper1.m 

#### This is the file you want to execute.

### paper1.m (This is the one you may want to edit)

#### This is the main function includes the drivers for the following six programs. 1. LogNormal model's HJB equation approximation solver test; 2. Mean Reverting model's HJB equation approximation solver test; 3. CIR model's HJB equation approximation test; 4. Monte Carlo back test for LogNormal; 5. Monte Carlo back test for Mean Reverting; 6. Monte Carlo back test for CIR. Check out the script for further parameters input instructions.

### BtEngine.m

#### This is a class implementing the back test Engine

### Constraint.m

#### This is a class implementing Leverage ratio constraint, consumption ratio constraint, maximum return constraint and maximum draw down constraint.

### ModelEvolver.m

#### This is a class implementing a data simulator according to the models specified in the program for the Monte Carlo back test.

### WKBHierarchySolver.m

#### This is a class implementing the approximation solver.

### HamiltonianSystem.m

#### This is a class implementing the Hamiltonian System to calculate all kinds of hamiltonians intermediate results.

### PortfolioCalculator.m

#### This is a class implementing the portfolio calculator to calulate all kinds of portfolio related intermediate results.

### UtilityCalculator.m

#### This is a class implementing a calculator to calculate all kinds of utility related intermediate results.

### Model.m

#### This is a class implementing different models related intermediate results.


### Questions???

#### Any questions about the code, please contact: Zhenfeng.Jose.Liang@gmail.com