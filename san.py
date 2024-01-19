import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import brentq as root
from rhodium import *
# from platypus import *
"""
class CubicDPSLever(Lever):
    
    def __init__(self, name, length = 1, c_bounds = (-2, 2), r_bounds = (0, 2)):
        super(CubicDPSLever, self).__init__(name)
        self.length = length
        self.c_bounds = c_bounds
        self.r_bounds = r_bounds

    def to_variables(self):
        result = []
        
        for _ in range(self.length):
            result += [Real(self.c_bounds[0], self.c_bounds[1])] # the center
            result += [Real(self.r_bounds[0], self.r_bounds[1])] # the radius
            result += [Real(0, 1)]                               # the weight
        
        return result

    def from_variables(self, variables):
        policy = {}
        policy["length"] = self.length
        policy["rbfs"] = []

        for i in range(self.length):
            policy["rbfs"] += [{
                "center" : variables[i*3+0],
                "radius" : variables[i*3+1],
                "weight" : variables[i*3+2] }]

        weight_sum = sum([p["weight"] for p in policy["rbfs"]])
        
        for i in range(self.length):
            policy["rbfs"][i]["weight"] /= weight_sum

        return policy

def evaluateCubicDPS(policy, current_value):
    value = 0
    
    for i in range(policy["length"]):
        rbf = policy["rbfs"][i]
        value += rbf["weight"] * abs((current_value - rbf["center"]) / rbf["radius"])**3
        
    value = min(max(value, 0.01), 0.1)
    return value    
"""
def lake_problem(policy,  # the DPS policy
         b = 0.42,        # decay rate for P in lake (0.42 = irreversible)
         q = 2.0,         # recycling exponent
         mean = 0.02,     # mean of natural inflows
         stdev = 0.001,   # standard deviation of natural inflows
         alpha = 0.4,     # utility from pollution
         delta = 0.98,    # future utility discount rate
         nsamples = 100,  # monte carlo sampling of natural inflows
         steps = 100):    # the number of time steps (e.g., days)
    Pcrit = root(lambda x: x**q/(1+x**q) - b*x, 0.01, 1.5)
    X = np.zeros((steps,))
    decisions = np.zeros((steps,))
    average_daily_P = np.zeros((steps,))
    reliability = 0.0
    utility = 0.0
    inertia = 0.0

    for _ in range(nsamples):
        X[0] = 0.0
        
        natural_inflows = np.random.lognormal(
                math.log(mean**2 / math.sqrt(stdev**2 + mean**2)),
                math.sqrt(math.log(1.0 + stdev**2 / mean**2)),
                size = steps)
        
        for t in range(1,steps):
            decisions[t-1] = evaluateCubicDPS(policy, X[t-1])
            X[t] = (1-b)*X[t-1] + X[t-1]**q/(1+X[t-1]**q) + decisions[t-1] + natural_inflows[t-1]
            average_daily_P[t] += X[t]/float(nsamples)
        
        reliability += np.sum(X < Pcrit)/float(steps)
        utility += np.sum(alpha*decisions*np.power(delta,np.arange(steps)))
        inertia += np.sum(np.diff(decisions) > -0.01)/float(steps-1)
      
    max_P = np.max(average_daily_P)
    reliability /= float(nsamples)
    utility /= float(nsamples)
    inertia /= float(nsamples)
    
    return (max_P, utility, inertia, reliability)

def get_lists():
    global reservoir_areas
    df = pd.read_excel('inputbargi.ods', engine='odf')
    # print(df.head())
    # print(df.columns)
    inflow = np.asarray(df["BARGI RESERVOIR INFLOW (VIRGIN FLOW)"])
    evap = np.asarray(df["BARGI EVAP. (mm)"])
    evap = evap*(reservoir_areas[0]/1000.0) # For uniform dimensions. We want all to be in MCM. evap is in mm, area in km^2, hence divide by 1000
    # overflow = np.zeros(evap.shape)
    # print(inflow.shape)
    # print(evap.shape)
    return inflow, evap

def get_constraints(vars):
    global N
    global Obj_num
    global Constraints_num
    global S_max
    Inflow, Evap = get_lists()
    constraints = [(vars[i+1]-vars[i]-Inflow[i]+ vars[(N//2+i)] + Evap[i] + max(0,vars[i+1]-S_max)) for i in range(N//2-1) ]
    
    # print("DEBUG: objs.shape: ", len(objs), " constraints.len: ", len(constraints))
    return constraints

def Objs(vars):
    objs=[]
    obj1=(-1.0*np.sum(np.asarray([(vars[i]-vars[(N//2+i)]) for i in range(N//2)]))) # Obj function: Max z = ∑ S_t − Q_t
    objs.append(obj1)
    return obj1 

N=int(324*2) # Decision Vars 324 S_t and 324 Q_t
Obj_num =1
S_min=409.0
S_max=425.7
Q_max =302.0
Constraints_num =int(N/2)-1  #S_t+1 = S_t + I_t − Q_t − E_t − O_t
reservoir_areas=[267.97]  # list of areas of reservoirs [BARGI: km^2]
constraint=[]

    # S_min <= S_t <=S_max

for i in range(N//2):
    constraint.append(Constraint(str( "S_" + str(i) + "  >= " + str(S_min) )))
    constraint.append(Constraint(str( "S_" + str(i) + "  <= " + str(S_max) )))
        
for i in range(N//2):
    constraint.append(Constraint(str( "Q_" + str(i) + "  >= " + str(S_min) )))
    constraint.append(Constraint(str( "Q_" + str(i) + "  <= " + str(S_max) )))


model = Model(Objs)

S_t_list = [Parameter(str("S_"+str(i))) for i in range(0,N//2)]
Q_t_list = [Parameter(str("Q_"+str(i))) for i in range(0,N//2)]


model.parameters = S_t_list + Q_t_list

#model.parameters = [Parameter("policy"),Parameter("b"),Parameter("q"),Parameter("mean"),Parameter("stdev"),Parameter("delta")]

model.responses = [Response("obj1", Response.MINIMIZE)]

#model.responses = [Response("max_P", Response.MINIMIZE),Response("utility", Response.MAXIMIZE),Response("inertia", Response.MAXIMIZE),Response("reliability", Response.MAXIMIZE)]

# Use our new DPS lever
Lever_S_t = [RealLever(str("S_"+str(i)), S_min, S_max, length=10) for i in range(0,N//2)]
Levers_Q_t = [RealLever(str("Q_"+str(i)), 0.01, Q_max, length=10) for i in range(0,N//2)]
model.levers = Lever_S_t + Levers_Q_t
#model.levers = [CubicDPSLever("policy", length=3)]

model.constraints = constraint
# Some parameters are exogeneous uncertainties, and we want to better
# understand how these uncertainties impact our model and decision making
# process
model.uncertainties = [UniformUncertainty("b", 0.1, 0.45),
                       UniformUncertainty("q", 2.0, 4.5),
                       UniformUncertainty("mean", 0.01, 0.05),
                       UniformUncertainty("stdev", 0.001, 0.005),
                       UniformUncertainty("delta", 0.93, 0.99)]

setup_cache(file="example.cache")
output = cache("dps_output", lambda: optimize(model, "NSGAII", 10000))
output.save('optimization_results.csv')

print("output")



# Use Seaborn settings for pretty plots
# sns.set()

# # Plot the points in 2D space
# scatter2d(model, output)
# plt.show()
"""
# The optional interactive flag will show additional details of each point when
# hovering the mouse
scatter2d(model, output, brush="reliability >= 0.5 and utility > 0.5")
plt.show()

# Most of Rhodiums's plotting functions accept an optional expr argument for
# classifying or highlighting points meeting some condition
scatter2d(model, output, x="reliability", brush=Brush("reliability >= 0.2"))
plt.show()

# Plot the points in 3D space
scatter3d(model, output, s="reliability", show_legend=True)
plt.show()
  
# Kernel density estimation plots show density contours for samples.  By
# default, it will show the density of all sampled points
kdeplot(model, output, x="max_P", y="utility")
plt.show()

# Alternatively, we can show the density of all points meeting one or more
# conditions
kdeplot(model, output, x="max_P", y="utility",
        brush=["reliability >= 0.2", "reliability < 0.2"],
        alpha=0.8)
plt.show()

# Pairwise scatter plots shown 2D scatter plots for all outputs
pairs(model, output)
plt.show()

# We can also highlight points meeting one or more conditions
pairs(model, output,
      brush=["reliability >= 0.2", "reliability < 0.2"])
plt.show()

# Joint plots show a single pair of parameters in 2D, their distributions using
# histograms, and the Pearson correlation coefficient
joint(model, output, x="max_P", y="utility")
plt.show()

# A histogram of the distribution of points along each parameter
hist(model, output)
plt.show()
 
# A parallel coordinates plot to view interactions among responses
parallel_coordinates(model, output, colormap="rainbow", zorder="reliability", brush=Brush("reliability > 0.2"))     
plt.show()
"""