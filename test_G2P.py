from Gurobi_2_Pyomo import *
import pyomo.environ as po
from gurobipy import *

gm = Model()

x = gm.addVar(vtype=GRB.BINARY)
y = gm.addVar(vtype=GRB.INTEGER, ub=10)
z = gm.addVar(lb=0, ub=2)

gm.addConstr(x*y+2*x*z - y*z <= 5)
gm.addConstr(10*x+2*y+z <= 30)
gm.setObjective(x**2-2*x*y+z**2-z*x)

#gm.setParam('Nonconvex', 2)
#gm.optimize()
#print((x.x, y.x, z.x))

gm.update()

G2P = Gurobi_2_Pyomo(gm)
G2P.convert_model()

G2P.mdl.pprint()
