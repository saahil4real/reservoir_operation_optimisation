from platypus import *
# from platypus.algorithms import NSGAII
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def get_lists():
    global reservoir_areas
    df = pd.read_excel('ravishankar.ods', engine='odf')
    inflow = np.asarray(df["Inflow \n(MCM)"])
    evap = np.asarray(df["Evaporation\n(MCM)"])
    return inflow, evap


def expressions(vars):
    global N
    global Obj_num
    global Constraints_num
    global S_max
    Inflow, Evap = get_lists()
    objs=[]
    obj1=(-1.0*np.sum(np.asarray([(vars[i]-vars[(N//2+i)]) for i in range(N//2)]))) # Obj function: Max z = ∑ S_t − Q_t
    objs.append(obj1)
    constraints = [(vars[i+1]-vars[i]-Inflow[i]+ vars[(N//2+i)] + Evap[i] + max(0,vars[i+1]-S_max)) for i in range(N//2-1) ]
    # print("DEBUG: objs.shape: ", len(objs), " constraints.len: ", len(constraints))
    return objs, constraints


N=int(324*2) # Decision Vars 324 S_t and 324 Q_t
Obj_num =1
S_min=143.6 #in mcm
S_max=910.5
Q_max = 1539.0
Constraints_num =int(N/2)-1  #S_t+1 = S_t + I_t − Q_t − E_t − O_t
reservoir_areas=[267.97]  # list of areas of reservoirs [BARGI: km^2]
constraint=[]

for i in range(N//2):
    constraint.append(Real(S_min, S_max)) # S_min <= S_t <=S_max
for i in range(N//2):
    constraint.append(Real(0.01, Q_max))

# algorithms = [NSGAII, (NSGAIII, {"divisions_outer":12})]

problem = Problem(N, Obj_num, Constraints_num)
problem.types[:] = constraint
problem.constraints[:] = "==0"
problem.function = expressions
# problem.directions[:] = Problem.MAXIMIZE

# with ProcessPoolEvaluator(4) as evaluator:
#         results = experiment(algorithms, problems, nfe=10, evaluator=evaluator)

#         hyp = Hypervolume(minimum=[0, 0, 0], maximum=[1, 1, 1])
#         hyp_result = calculate(results, hyp, evaluator=evaluator)
#         display(hyp_result, ndigits=3)


# algorithm = NSGAII(problem, population_size = 140)
algorithm = NSGAIII(problem,divisions_outer = 17, population_size = 820)
#algorithm = EpsMOEA(problem, epsilons = [0.08],population_size = 190)
algorithm.run(10)
feasible_solutions = [s for s in algorithm.result if s.feasible]
nondominated_solutions = nondominated(algorithm.result)

print("# Feasible solutions: ", len(feasible_solutions))
# print(feasible_solutions)
# print(nondominated_solutions)

result = algorithm.result

result_dict ={}
# result_dict["obj"] = []
result_dict["storage"] =[]
result_dict["outflow"]=[]

counterr=0

stor = []
outf = []

# print(feasible_solutions)

#for NSGAII
"""
for solution in result:
    print("Count: ", counterr)
    counterr+=1
    # print("Obj: ",solution.objectives)
    # print(solution.variables)
    # result_dict["obj"].append(solution.objectives[0])
    res = np.asarray(solution.variables)
    # print(len(res))
    
    if(counterr == 99):
        #print(solution.variables)
        stor.append(res[:N//2])
        outf.append(res[N//2:])
    # print(len(res))
    # result_dict["storage"].append(res[:N//2])
    # result_dict["outflow"].append(res[N//2:])
"""

#for NSGAIII

for solution in result:
    print("Count: ", counterr)
    counterr+=1
    # print("Obj: ",solution.objectives)
    # print(solution.variables)
    # result_dict["obj"].append(solution.objectives[0])
    res = np.asarray(solution.variables)
    # print(len(res))
    
    if(counterr == 3):
        print(solution.variables)
        stor.append(res[:N//2])
        outf.append(res[N//2:])
    # print(len(res))
    # result_dict["storage"].append(res[:N//2])
    # result_dict["outflow"].append(res[N//2:])


#for EMOEA
"""
for solution in result:
    print("Count: ", counterr)
    counterr+=1
    # print("Obj: ",solution.objectives)
    # print(solution.variables)
    # result_dict["obj"].append(solution.objectives[0])
    res = np.asarray(solution.variables)
    # print(len(res))
    
    # if(counterr == 3):
    # print(solution.variables)
    stor.append(res[:N//2])
    outf.append(res[N//2:])
    # print(len(res))
    # result_dict["storage"].append(res[:N//2])
    # result_dict["outflow"].append(res[N//2:])
"""

print(result_dict)

# print(stor[0])


# my_array = np.array([stor[0],outf[0]])

new_arr = []

for i in range(N//2):
    new_arr.append([outf[0][i],stor[0][i]])


df = pd.DataFrame(new_arr, columns = ['Outflow','Storage'])
#df.to_excel("RaviShankar_NSGA2_Sens.xlsx")
df.to_excel("RaviShankar_NSGA3_Sens.xlsx")
#df.to_excel("RaviShankar_EpsMOEA_Sens.xlsx")
print(df)
print(type(df))
