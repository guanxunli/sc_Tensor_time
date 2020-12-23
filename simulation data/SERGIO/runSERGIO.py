import numpy as np
from SERGIO.sergio import sergio
np.random.seed(1)
sim = sergio(number_genes=103, number_bins = 10, number_sc = 1000, 
             noise_params = 1, decays=2,
             noise_type='dpd', dynamics= True, bifurcation_matrix=[[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0,0]])
sim.build_graph('iT2.txt', 'iR2.txt', shared_coop_state=1)
np.random.seed(1)
sim.simulate_dynamics()
np.random.seed(1)
exprU, exprS = sim.getExpressions_dynamics()
np.random.seed(1)
exprU, exprS = sim.convert_to_UMIcounts_dynamics(exprS, exprU)
exprS = np.concatenate(exprS, axis = 1)
exprU = np.concatenate(exprU, axis = 1)
np.savetxt('diffOutput_S.csv',exprS, delimiter=',')
np.savetxt('diffOutput_U.csv',exprU, delimiter=',')
