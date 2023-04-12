import gomea
import matplotlib.pyplot as plt

frv = gomea.fitness.RosenbrockFunction(10,value_to_reach=1e-10)
lm = gomea.linkage.Univariate()
rvgom = gomea.RealValuedGOMEA(fitness=frv, linkage_model=lm, lower_init_range=-115, upper_init_range=-100)
result = rvgom.run()
plt.grid()
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Number of evaluations')
plt.ylabel('Best objective value')
plt.plot(result['evaluations'],result['best_obj_val'])
plt.savefig('convergence_plot.pdf')

