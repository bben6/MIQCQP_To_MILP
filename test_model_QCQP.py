from glob import glob as glb
from os import path
from gurobipy import *
from Model_qcqp import Model_qcqp

fsplit = path.splitext

# Flags
UNIVAR = Model_qcqp.UNIVAR
NMDT = Model_qcqp.NMDT
DNMDT = Model_qcqp.DNMDT
TNMDT = Model_qcqp.TNMDT
TDNMDT = Model_qcqp.TDNMDT
BIN2 = Model_qcqp.BIN2
BIN3 = Model_qcqp.BIN3
MCCORMICK = Model_qcqp.MCCORMICK
UNIVAR_NOSTUB = Model_qcqp.UNIVAR_NOSTUB
UNIVAR_NOST = Model_qcqp.UNIVAR_NOST

mthd_name_strs = {
    DNMDT: 'DNMDT',
    NMDT: 'NMDT',
    TDNMDT: 'T-DNMDT',
    TNMDT: 'T-NMDT',
    UNIVAR: 'Univar',
    BIN2: 'Bin2',
    BIN3: 'Bin3',
    MCCORMICK: 'McCormick Envelopes',
    UNIVAR_NOST: 'Univar, No Satooth Relaxations',
    UNIVAR_NOSTUB: 'Univar, No Satooth Upper-Bounds'
}

# Core settings
### Model settings; user input ###
# All methods
L = 1

# DNMDT #
# ldas: list of lambda values to use for DNMT.
ldas = [0, 1]

# NMDT #
disc_left = 1
disc_right = 0

# Sawtooth #
L1 = 1

# Choose methods
mthds = [UNIVAR, NMDT, DNMDT, BIN2, BIN3, MCCORMICK]

m = Model()

x = m.addVar(lb=1, ub=3, name='x')
y = m.addVar(lb=1, ub=3, name='y')
z = m.addVar(lb=2, ub=3, name='z')

m.setObjective(-(x-.35)**2 + (y-.48)**2 - (z-.74)**2 - 5*x*y - 2*x*z - 2*y*z)
m.addConstr(x+y+z <= 5.5)

# Optimize directly
m.setParam(GRB.Param.NonConvex, 2)
m.setParam(GRB.Param.OutputFlag, 0)

print('Optimizing directly')
m.optimize()

print(f'Optimal value: {m.getObjective().getValue()}')
print(', '.join(f"{v.varName}={v.x}" for v in m.getVars() if v.varName in ['x', 'y', 'z']), end='\n\n')

m.update()


for mthd in mthds:
    print(f'Commencing method {mthd_name_strs[mthd]}')
    mdl_holder = Model_qcqp(m, L=L, L1=L1, ldas=ldas, method=mthd,
                            disc_left=disc_left, disc_right=disc_right, st_do_approx=0)

    mdl = mdl_holder.mdl

    mdl.setParam(GRB.Param.OutputFlag, 0)

    mdl.optimize()

    # Print optimized solution
    print(f'Optimal value: {mdl.getObjective().getValue()}')
    print(', '.join(f"{v.varName}={v.x}" for v in mdl.getVars() if v.varName in ['x', 'y', 'z']), end='\n\n')
