# composite-optimizer

Part of the ASE student challenge.

**No extra libraries are necessary.**

## _1. Set your values_

Directly below the imports in [main.cpp](main.cpp), one can change the optimizer options.
Large forces will work; however, due to large stack sizes and more possible combinations, the variance in the output will
be larger.

```
// Step 1: Insert custom values here. There are no constraints.

// Dimensions [mm]
double a = 400.0;
double b = 150.0;
double plyThickness = 0.184;

// Force per width [kN/mm]
double Nx = -0.5;
double Ny = 0.2;
double Tau = 0.1;

// Material Properties [MPa]
double E1 = 130000;
double E2 = 10000;
double G12 = 5000;

double nu12 = 0.3;
```

## _2. Set up optimizer parameters_

These are the standard values. Several tests have shown that these parameters yield good results for large and
small stacks. On rare occasions, these values can be adapted. Should the optimizer not stop with the flag _---
Converged! ---_ either the _initial population_ or the _fillPopulation_ values should be changed. This will have a penalty
on the runtime. The _fillPopulation_ option will ensure that the population size is not reduced and only the implemented
convergence measure leads to convergence. This has a significant effect on the run time.

```
// Optimization Parameters
double RF_weight = 1;
double size_weight = 100000;
int initial_population = 10000;
int initial_max_stack_size = 200;
bool fillPopulation = false;
```

## _3. Run main_

Run the main function.

## _4. Read output_

The solution will have the following form:

```
Found following Solution:
  RF: ...
  Fitness: ...
  Stack Size: ...
  Stacking Sequence:
      .
      .
      .
  --- Mid Plane ---
```
