# MIQCQP_To_MILP
A class to approximate a nonconvex MIQCQP, modeled in Gurobi, with a MILP by relaxing all quadratic and bivariate terms, using any of several methods. These methods will be formally published in an upcoming paper.

The class is titled Model_qcqp.py. To obtain a MIP-relaxes model for the QCQP, construct a Model_qcqp object, passing in the Gurobi model and desired settings, then access the constructed model via the object.mdl property. For example, to construct a relaxed model with default settings:

from gurobipy import *
from Model_qcqp import Model_qcqp

# Create the model
m = Model()

x = m.addVar(lb=1, ub=3, name='x')
y = m.addVar(lb=1, ub=3, name='y')
z = m.addVar(lb=2, ub=3, name='z')

m.setObjective(-(x-.35)**2 + (y-.48)**2 - (z-.74)**2 - 5*x*y - 2*x*z - 2*y*z)
m.addConstr(x+y+z <= 5.5)
m.update()

# Convert the model
mdl_holder = Model_qcqp(m)
mdl = mdl_holder.mdl

# Optimize the model
mdl.optimize()

For more usage information, see the example code in test_model_QCQP.py. 

This repository also includes a class to convert a Gurobi model to a Pyomo model with all integer variables relaxed to continuous variables, titled Gurobi_2_Pyomo.py. This class could be useful in creating a callback function to perform local optimization within a Gurobi solve of a relaxed model, e.g. using ipopt. 

To use this class, construct a Gurobi_2_Pyomo object, passing in the original model as an argument, then access the constructed model via the object.mdl property. For example, 

from Gurobi_2_Pyomo import *
from gurobipy import *

gm = Model()

x = gm.addVar(vtype=GRB.BINARY)
y = gm.addVar(vtype=GRB.INTEGER, ub=10)
z = gm.addVar(lb=0, ub=2)

gm.addConstr(x*y+2*x*z - y*z <= 5)
gm.addConstr(10*x+2*y+z <= 30)
gm.setObjective(x**2-2*x*y+z**2-z*x)

gm.update()

G2P = Gurobi_2_Pyomo(gm)
G2P.convert_model()

G2P.mdl.pprint()
