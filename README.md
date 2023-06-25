[![Build and run demo](https://github.com/abouter/gomea/actions/workflows/run_demo.yml/badge.svg?branch=main)](https://github.com/abouter/gomea/actions/workflows/run_demo.yml)

# README
The Gene-pool Optimal Mixing Evolutionary Algorithm (GOMEA) is a state-of-the-art Model-Based Evolutionary Algorithm (MBEA) with an ongoing line of research into various domains of (evolutionary) optimization.
The initial release of the GOMEA library (v1.0.0) is described in the publication:
- [Bouter, A., & Bosman, P.A.N. (2023). A Joint Python/C++ Library for Efficient yet Accessible Black-Box and Gray-Box Optimization with GOMEA. In Proceedings of the Genetic and Evolutionary Computation Conference Companion](https://arxiv.org/abs/2305.06246)

The current release of the GOMEA library supports optimization in the following domains:
- Discrete single-objective optimization
- Real-valued single-objective optimization

Support for multi-objective optimization in both of these domains is planned for a future release.

Relevant most recent publications are the following:
- [Dushatskiy, A., Virgolin, M., Bouter, A., Thierens, D., & Bosman, P.A.N. (2021). Parameterless Gene-pool Optimal Mixing Evolutionary Algorithms. arXiv preprint arXiv:2109.05259.](https://arxiv.org/abs/2109.05259)
- [Bouter, A., Alderliesten, T., & Bosman, P.A.N. (2021). Achieving highly scalable evolutionary real-valued optimization by exploiting partial evaluations. Evolutionary computation, 29(1), 129-155.](https://ieeexplore.ieee.org/abstract/document/9367090)

## INSTALL
The most straightforward installation option is through pip, as follows.
```
pip install gomea
```

To build from source, execute the following in the root GOMEA folder.
```
make install
```

## RUN
Running a simple benchmark problem, e.g., the Rosenbrock function, can be done as follows:
```
import gomea
frv = gomea.fitness.RosenbrockFunction(10,value_to_reach=1e-10)
lm = gomea.linkage.Univariate()
rvgom = gomea.RealValuedGOMEA(fitness=frv, linkage_model=lm, lower_init_range=-115, upper_init_range=-100)
result = rvgom.run()
```

Plotting a simple convergence plot can then be done using:
```
import matplotlib.pyplot as plt
plt.grid()
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Number of evaluations')
plt.ylabel('Best objective value')
plt.plot(result['evaluations'],result['best_obj_val'])
```

Rather than using a predefined fitness function, it is also possible to define your own custom fitness function for Black-Box Optimization (BBO) or Gray-Box Optimization (GBO). The definition of a custom (real-valued) BBO function can be done as follows.
```
class CustomRosenbrockFunction(gomea.fitness.BBOFitnessFunctionRealValued):
    def objective_function( self, objective_index, variables ):
        f = 0
        for i in range(len(variables)-1):
            x = variables[i]
            y = variables[i+1]
            f += 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x)
        return f

dim = 10
f_rosenbrock = CustomRosenbrockFunction(dim,value_to_reach=1e-10)
```

The definition of a GBO function requires defining the input and output of each subfunction. Thorough definitions of this are given in the [EvoSOFT workshop paper](https://arxiv.org/abs/2305.06246).
An example of the Rosenbrock function defined in a GBO setting is as follows:
```
class CustomRosenbrockFunction(gomea.fitness.GBOFitnessFunctionRealValued):
    def number_of_subfunctions( self ):
        return self.number_of_variables-1
    
    def inputs_to_subfunction( self, subfunction_index ):
        return [subfunction_index, subfunction_index+1]

    def subfunction(self, subfunction_index, variables):
        x = variables[subfunction_index]
        y = variables[subfunction_index+1]
        return 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x)
        
dim = 10
f_rosenbrock = CustomRosenbrockFunction(dim,value_to_reach=1e-10)
```

Further examples are given in the demos directory.
