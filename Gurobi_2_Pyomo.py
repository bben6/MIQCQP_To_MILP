from gurobipy import *
import pyomo.environ as po

# The formulations in this module allow the user to provide the auxiliary variables. This give the user the option
#   To use the same auxiliary variables among different terms.

# Model assumptions:
#   ALL VARIABLES MUST HAVE BOTH UPPER AND LOWER BOUNDS, EXPLICITLY SET AS GUROBI ATTRIBUTES.
#   No general constraints.
#   Model has been updated, so that all variables and constraints have been added.


class Gurobi_2_Pyomo:
    def __init__(self, mdl_o):

        # Read in model variables and constraints
        self.mdl_o = mdl_o
        self.xs_o = mdl_o.getVars()
        self.I = range(len(self.xs_o))

        # Create the model
        self.mdl = po.ConcreteModel()

        # Var stuff
        self.xs = []
        self.xotox = {}
        self.name2x = {}

        # Constraint stuff
        self.obj = 0

    g2b = {GRB.BINARY: po.Integers,
           GRB.INTEGER: po.Integers,
           GRB.CONTINUOUS: po.Reals}

    def make_vars(self):
        mdl = self.mdl
        mdl_o = self.mdl_o
        I = self.I
        xs_o = self.xs_o

        # copy all modifiable variable attributes
        def gattr(attr, i): return mdl_o.getAttr(attr, [xs_o[i]])[0]
        def bounds(m, i): return gattr('lb', i), gattr('ub', i)
        def domain_type(m, i): return self.g2b[gattr('vtype', i)]

        mdl.I = po.Set(initialize=I)
        mdl.xs = po.Var(mdl.I, bounds=bounds, domain=domain_type)
        self.xs = mdl.xs
        xs = mdl.xs

        #raise ValueError(str([domain_type('', i) for i in I]))
        self.name2x = {gattr('varname', i): xs[i] for i in I}
        self.xotox = {xs_o[i]: xs[i] for i in I}

    def convert_model(self):
        self.make_vars()

        self.process_lconstrs()
        self.process_qconstrs()
        self.process_obj()

    def process_lconstrs(self):
        mdl = self.mdl
        mdl_o = self.mdl_o

        mdl_o_constrs = mdl_o.getConstrs()

        def make_constraint(m, j):
            lc = mdl_o_constrs[j]

            LHS = self.process_lexpr(mdl_o.getRow(lc))
            RHS = lc.rhs
            sense = lc.sense

            if sense == '=':
                return LHS == RHS
            elif sense == '>':
                return LHS >= RHS
            else:
                return LHS <= RHS

        mdl.J = po.Set(initialize=range(len(mdl_o_constrs)))
        mdl.lin_constrs = po.Constraint(mdl.J, rule=make_constraint)

    def freeze_integers(self):
        mdl_o = self.mdl_o
        I = self.I
        xs_o = self.xs_o
        xs = self.xs

        def gattr(attr, i): return mdl_o.getAttr(attr, [xs_o[i]])[0]

        for i in I:
            if gattr('vtype', i) in [GRB.INTEGER, GRB.BINARY]:
                xs[i].fixed = True

    def process_qconstrs(self):
        mdl = self.mdl
        mdl_o = self.mdl_o

        mdl_o_constrs = mdl_o.getQConstrs()

        def make_constraint(m, j):
            qc = mdl_o_constrs[j]

            LHS = self.process_qexpr(mdl_o.getQCRow(qc))
            RHS = qc.qcrhs
            sense = qc.qcsense

            if sense == '=':
                return LHS == RHS
            elif sense == '>':
                return LHS >= RHS
            else:
                return LHS <= RHS

        mdl.K = po.Set(initialize=range(len(mdl_o_constrs)))
        mdl.qp_constrs = po.Constraint(mdl.K, rule=make_constraint)

    g2p_sense = {GRB.MINIMIZE: po.minimize,
                 GRB.MAXIMIZE: po.maximize}

    def process_obj(self):
        mdl = self.mdl
        mdl_o = self.mdl_o

        obj_o = mdl_o.getObjective()
        if isinstance(obj_o, QuadExpr):
            obj = self.process_qexpr(obj_o)
        else:
            obj = self.process_lexpr(obj_o)
        sense = mdl_o.ModelSense

        self.obj = obj
        mdl.obj = po.Objective(rule=lambda m: obj, sense=self.g2p_sense[sense])

    def process_qexpr(self, qexpr):
        mdl = self.mdl
        xotox = self.xotox

        result_expr = 0

        # process quadratic part
        nqt = qexpr.size()   # number of quadratic terms

        for j in range(nqt):

            coeff = qexpr.getCoeff(j)

            x = xotox[qexpr.getVar1(j)]
            y = xotox[qexpr.getVar2(j)]

            result_expr += coeff*x*y

        # process linear expression
        lexpr = qexpr.getLinExpr()
        result_expr += self.process_lexpr(lexpr)

        return result_expr

    def process_lexpr(self, lexpr):
        xotox = self.xotox

        result_expr = 0

        # process linear part
        nl = lexpr.size()
        for j in range(nl):
            coeff = lexpr.getCoeff(j)
            var = xotox[lexpr.getVar(j)]

            result_expr += coeff*var

        # process constant part
        c = lexpr.getConstant()
        result_expr += c

        return result_expr
