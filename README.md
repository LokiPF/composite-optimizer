# composite-optimizer

Part of the ASE student challenge.

**No extra libraries are necessary.**

The algorithm is based on the _Nondominated Sorting Genetic Algorithm II (NSGAII)_ [1]. This algorithm belongs to the
_multi-objective evolutionary algorithms (MOEAs)_.

## _Design of Algorithm_

1. Initialization: An equally distributed population is generated based on the maximum stack size and initial population
   size.
2. NSGAII algorithm: based on the reserve factor (strength measure) and stack size (mass) the best stacks are selected
   using a mix of tournament and elite selection. Crossover and adaptive mutation ensures a large diversity in the
   beginning and a smaller diversity in the end. The population is divided into Pareto fronts based on dominance, and crowding distance ensures that the population is spread out across the objective space. This         prevents the algorithm from optimizing only one objective at the expense of the other.

3. Convergence: Optimally, the GA finds the Pareto front in which further mutations will decrease one objective while
   increasing another. This, however, cannot always be ensured. Using the hypervolume of the Pareto front, the
   algorithm can then converge. If the hypervolume improvement falls below a certain threshold over several
   generations, the GA is said to be converged.
4. Output: In the end, all the solutions found are being output or only the best one is seen over several runs.


## _1. Set your values_

Directly below the imports in [main.cpp](main.cpp), one can change the optimizer options.
Large forces will work; however, due to large stack sizes and more possible combinations, the variance in the output
will
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

The standard values can be changed. However, these values should lead to good results.

```
// Step 2: Adapt optimization parameters if necessary.

// Initial Parameters
int initial_population = 150;
double initial_max_stack_size = 50;

// Mutation Parameters
/// Adaptive mutation is used. ///
double initial_mutation_rate = 1; // All positive values are allowed. A value of 0 deactivates ply orientation mutation.
double slope = 0.5; // Decides how fast the mutation rate is reduced. Linear Slope. 0: mutation rate remains constant.
int stack_size_variability = 2; // Changes in size of the stack. Does not depend on slope or initial mutation rate.

// Convergence
double threshold = 0.0001; // By how much the hypervolume must improve.
int previous_generations = 100; // Number of generations of the hypervolume to improve for convergence to be true.
const std::vector<double>  reference_point = {16, 1};

// Consistent Measures
constexpr int runs = 2; // To improve the consistency further, the NSGAII algorithm is being run several times.

// Output
bool show_all_solutions = false; // Only the best solution is given in the console when this is set to false.
```

## _3. Run main_

Run the main function.

## _4. Read output_

Only balanced and symmetric stacks are given as an output. The stack size is always even by design.
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

Sources:

[1]: K. Deb, A. Pratap, S. Agarwal and T. Meyarivan, "A fast and elitist multiobjective genetic algorithm: NSGA-II,"
in
IEEE Transactions on Evolutionary Computation, vol. 6, no. 2, pp. 182-197, April 2002, doi: 10.1109/4235.996017.
keywords: {Genetic algorithms;Sorting;Computational complexity;Evolutionary computation;Computational
modeling;Testing;Decision making;Associate members;Diversity reception;Constraint optimization},

[2]: OpenAI, "Consultation on small portions of composite laminate optimization and  helper function implementation," ChatGPT model, version 4, OpenAI, 2024. [Online].


