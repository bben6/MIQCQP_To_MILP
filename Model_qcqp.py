from gurobipy import *

# The formulations in this module allow the user to provide the auxiliary variables. This give the user the option
#   To use the same auxiliary variables among different terms.

# Model assumptions:
#   ALL VARIABLES MUST HAVE BOTH UPPER AND LOWER BOUNDS, EXPLICITLY SET AS GUROBI ATTRIBUTES.
#   No general constraints.
#   Model has been updated, so that all variables and constraints have been added.
class Model_qcqp:
    NMDT = 0
    DNMDT = 1
    UNIVAR = 2
    BIN2 = 3
    BIN3 = 4
    MCCORMICK = 5
    TNMDT = 6
    TDNMDT = 7
    UNIVAR_NOSTUB = 8
    UNIVAR_NOST = 9

    def __init__(self, mdl_o, L=3, L1=10, method=DNMDT, ldas=(0.5,), disc_left=1, disc_right=0, st_do_approx=0):
        self.mdl_o = mdl_o
        self.L = L
        self.L1 = max(L, L1)  # Applies only to univariate formulation, or NMDT-based quadratic terms when strengthened with sawtooth LB
        self.method = method

        self.tighten_nmdt = (method in [self.TNMDT, self.TDNMDT])

        # disc_left, disc_right apply only to NMDT. Treated as binaries
        if not(disc_left or disc_right):
            self.disc_left = 1
            self.disc_right = 0
        else:
            self.disc_left = disc_left
            self.disc_right = disc_right

        # Store whether to model approximations or relaxations for quadratic terms
        self.st_do_approx = st_do_approx
        if self.method in [self.BIN2, self.BIN3]:
            self.st_do_approx = False


        # ldas apply only to D-NMDT
        if len(ldas) == 0:
            ldas = [0.5]
        self.ldas = ldas

        # Read in model variables and constraints
        self.xs_o = mdl_o.getVars()
        self.lcs_o = mdl_o.getConstrs()
        self.qcs_o = mdl_o.getQConstrs()
        self.I = range(len(self.xs_o))
        self.J = range(1, self.L+1)
        self.J1 = range(1, self.L1+1)

        #self.xtoi = {self.xs_o[i]: i for i in self.I}

        # Create the model
        self.mdl = Model()
        mdl = self.mdl
        # get data (bounds, etc) for original variables, then add to new model. For now, only bnds; names+etc later.

        # Initialize variables
        I = self.I
        xs_o = self.xs_o
        self.xs = mdl.addVars(I)
        xs = self.xs

        self.xotox = {xs_o[i]: xs[i] for i in I}

        # Copy over variables attributes
        for i in I:
            self.copy_var_attrs(xs[i], xs_o[i])

        # Update model so that variables can be used as keys
        mdl.update()

        # Get variable bounds
        gattr = lambda attr, var: mdl_o.getAttr(attr, [var])[0]
        self.xbnds = {xs[i]: (gattr("lb", xs_o[i]), gattr("ub", xs_o[i])) for i in I}

        # dictionaries for z-variables used to model z=x^2 or z=xy.
        self.xy_vdict = {}
        self.xsq_vdict = {}

        # dictionaries for auxiliary variables
        #   Binary variables, shared by all methods
        self.alphas = {}

        # sawtooth-based auxiliary definitions
        #   Stuff for single product terms only
        self.zplus = {}  # auxiliary LP-vars z for z>=(x+y)^2
        self.zplus_g = {}  # auxiliary LP-vars g for z>=(x+y)^2
        self.zminus = {}  # auxiliary LP-vars z for z>=(x-y)^2:
        self.zminus_g = {}  # auxiliary LP-vars g for z>=(x-y)^2
        #   Stuff for lots of terms
        self.gs = {}  # auxiliary vars for sawtooth-based individual quadratic terms z=xy or z<=xy

        # D-NMDT auxiliary definitions
        #   Stuff for single product terms only
        self.vs = {}  # auxiliary vars for alpha_x*(y-stuff)
        self.ws = {}  # auxiliary vars for alpha_y*(z-stuff)
        self.delzs = {}  # auxiliary vars for deltax * deltay stuff
        #   Stuff for lots of terms
        self.delxs = {}  # auxiliary vars for the deltax terms

        # Auxiliary constraints
        self.xy_constrs = {}
        self.xsq_constrs = {}

        # Process stuff to prep model for solving
        self.process_lconstrs()
        self.process_qconstrs()
        self.process_obj()

    def copy_var_attrs(self, x, xo):
        # copy all modifiable variable attributes
        mdl = self.mdl
        mdl_o = self.mdl_o
        gattr = lambda attr: mdl_o.getAttr(attr, [xo])[0]
        sattr = lambda attr, val: mdl.setAttr(attr, [x], [val])

        attr_list = 'lb,ub,vtype,varname,start'.split(',')
        for attr1 in attr_list:
            sattr(attr1, gattr(attr1))
        # omitting simplex-related stuff for now: vbasis, pstart; these may cause errors
        # Omitting 'partition' for now: setting at all seems to spark partition algorithm
        # Omitting 'start' for now. And vtag,varhintval,VarHintPri,BranchPriority

    def process_lconstrs(self):
        mdl = self.mdl
        mdl_o = self.mdl_o
        for lc in self.lcs_o:
            LHS = self.process_lexpr(mdl_o.getRow(lc))
            RHS = lc.rhs
            sense = lc.sense

            if sense == '=':
                mdl.addConstr(LHS == RHS)
            elif sense == '>':
                mdl.addConstr(LHS >= RHS)
            else:
                mdl.addConstr(LHS <= RHS)

    def process_qconstrs(self):
        mdl = self.mdl
        mdl_o = self.mdl_o
        for qc in self.qcs_o:
            LHS = self.process_qexpr(mdl_o.getQCRow(qc))
            RHS = qc.qcrhs
            sense = qc.qcsense

            if sense == '=':
                mdl.addConstr(LHS == RHS)
            elif sense == '>':
                mdl.addConstr(LHS >= RHS)
            else:
                mdl.addConstr(LHS <= RHS)

    def process_obj(self):
        mdl = self.mdl
        mdl_o = self.mdl_o

        obj_o = mdl_o.getObjective()
        if isinstance(obj_o, QuadExpr):
            obj = self.process_qexpr(obj_o)
        else:
            obj = self.process_lexpr(obj_o)
        sense = mdl_o.ModelSense

        mdl.setObjective(obj, sense)

    def process_qexpr(self, qexpr):
        xy_vdict = self.xy_vdict
        xsq_vdict = self.xsq_vdict
        mdl = self.mdl
        xotox = self.xotox
        method = self.method
        DNMDT = Model_qcqp.DNMDT
        NMDT = Model_qcqp.NMDT
        TDNMDT = Model_qcqp.TDNMDT
        TNMDT = Model_qcqp.TNMDT
        UNIVAR = Model_qcqp.UNIVAR
        UNIVAR_NOSTUB = Model_qcqp.UNIVAR_NOSTUB
        UNIVAR_NOST = Model_qcqp.UNIVAR_NOST
        BIN2 = Model_qcqp.BIN2
        BIN3 = Model_qcqp.BIN3
        MCCORMICK = self.MCCORMICK

        result_expr = 0

        # process quadratic part
        nqt = qexpr.size()   # number of quadratic terms

        isApprox = self.st_do_approx
        APPROX = Model_qcqp.ST_APPROX
        RELAX = Model_qcqp.ST_RELAX

        st_sense = (APPROX if isApprox else RELAX)

        for j in range(nqt):

            coeff = qexpr.getCoeff(j)

            x = xotox[qexpr.getVar1(j)]
            y = xotox[qexpr.getVar2(j)]

            if x is y:
                # Then this is a pure quadratic term
                if x not in xsq_vdict:
                    # TODO: model z = var1*var2 using some method. Add to dict.
                    xsq_vdict[x] = mdl.addVar(lb=-GRB.INFINITY, name=f'{x.varName}_sq')
                    if method in (DNMDT, TDNMDT):
                        self.add_DNMDTsq(x)
                    elif method in (NMDT, TNMDT):
                        self.add_NMDTsq(x)
                    elif method in [UNIVAR, BIN2, BIN3]:
                        self.sawtooth(x, sense=st_sense)
                    elif method == MCCORMICK:
                        self.add_McCormicksq(x)
                    elif method == UNIVAR_NOSTUB:
                        self.addsq_stlb(x)
                    elif method == UNIVAR_NOST:
                        self.addsq(x)

                result_expr += coeff * xsq_vdict[x]
            else:
                if (x, y) in xy_vdict:
                    result_expr += coeff * xy_vdict[x, y]
                elif (y, x) in xy_vdict:
                    result_expr += coeff * xy_vdict[y, x]
                else:
                    # TODO:
                    #     If using NMDT, we may require the user to specify a list of variables to be discretized, then
                    #       use some (simple/arbitrary) method to choose more to discretize if needed.
                    #     To discretize the fewest variables possible with NMDT for general problems, we need to
                    #       solve the maximal independent set problem, which is NP-hard (oof) to get the largest set of
                    #       variables to not discretize. (Dr. Hildebrand says it solves fast in practice.) Note that we
                    #       must always discretize variables with pure quadratic terms (i.e. diagonal of graph matrix
                    #       may have ones).
                    xy_vdict[x, y] = mdl.addVar(lb=-GRB.INFINITY, name=f'{x.varName}_times_{y.varName}')

                    if method in (NMDT, TNMDT):
                        self.add_NMDT(x, y)
                    elif method in (DNMDT, TDNMDT):
                        self.add_DNMDT(x, y)
                    elif method == UNIVAR:
                        self.add_univar(x, y)
                    elif method in (BIN2, BIN3):
                        is_plus = method == BIN2
                        self.add_bin(x, y, is_plus)
                    elif method == MCCORMICK:
                        self.add_McCormick(x, y)
                    elif method == UNIVAR_NOSTUB:
                        self.add_UNIVAR_NOSTUB(x, y)
                    elif method == UNIVAR_NOST:
                        self.add_UNIVAR_NOST(x, y)

                    result_expr += coeff * xy_vdict[x, y]

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

    # Sawtooth and univariate reformulation-related definitions

    ST_UB = 0
    ST_APPROX = 1
    ST_RELAX = 2

    def sawtooth(self, x, sense=ST_RELAX):
        # model y=x^2 with sawtooth approximation.
        # Possible types: UB, APPROX, and RELAX.
        #   UB: model y<=x^2 only
        #   APPROX: model y=x^2 using overapproximation
        #   RELAX: model sawtooth relaxation, with underapproximation accuracy L1
        mdl = self.mdl
        L = self.L
        J = self.J
        L1 = self.L1
        J1 = self.J1
        xsq_constrs = self.xsq_constrs
        xsq_vdict = self.xsq_vdict

        UB = Model_qcqp.ST_UB
        APPROX = Model_qcqp.ST_APPROX
        RELAX = Model_qcqp.ST_RELAX

        alphas = self.alphas
        gs = self.gs
        xname = x.varname

        xsq_constrs[x] = []
        xsqc = self.xsq_constrs[x]
        y = xsq_vdict[x]

        xlb = x.lb
        xub = x.ub

        lx = xub-xlb
        if lx >= 1e-10:
            if x not in alphas:
                alphas[x] = mdl.addVars(J, vtype=GRB.BINARY, name=f'{xname}_alpha')
                Jg = (J if sense in (UB, APPROX) else J1)
                gs[x] = mdl.addVars(Jg, lb=0, ub=1, name=f'{xname}_g')
                gs[x].update({0: (x-xlb)/lx})

            alpha = alphas[x]
            g = gs[x]
            yh = (y-xlb*(2*x-xlb))/lx**2
            xh = g[0]

            # Jg1: except j=0, g-indices for convex LB.
            Jg1 = (J if sense in (UB, APPROX) else J1)

            c1 = mdl.addConstrs(
                (2*(g[j-1] - alpha[j]) <= g[j]
                 for j in J)
            )
            c2 = mdl.addConstrs(
                (2*(alpha[j] - g[j-1]) <= g[j]
                 for j in J)
            )
            c3 = mdl.addConstrs(
                (g[j] <= 2*g[j-1]
                 for j in Jg1)
            )
            c4 = mdl.addConstrs(
                (g[j] <= 2*(1-g[j-1])
                 for j in Jg1)
            )

            xsqc += [c1, c2, c3, c4]
            if sense == UB:
                xsqc += [mdl.addConstr(yh <= xh - sum(2**(-2*j)*g[j] for j in J))]
            elif sense == APPROX:
                xsqc += [mdl.addConstr(yh == xh - sum(2**(-2*j)*g[j] for j in J))]
            else:
                c1 = mdl.addConstr(yh <= xh - sum(2**(-2*j)*g[j] for j in J))
                c2 = mdl.addConstrs(
                    (yh >= xh - sum(2**(-2*j)*g[j] for j in range(1, l+1)) - 2**(-2*l-2)
                     for l in range(L1+1)))
                c3 = mdl.addConstr(yh >= 0)
                c4 = mdl.addConstr(yh >= 2*xh-1)

                xsqc += [c1, c2, c3, c4]
        else:
            if sense in [APPROX, RELAX]:
                xsqc += [mdl.addConstr(y == ((xlb+xub)/2)**2)]
            elif sense == UB:
                xsqc += [mdl.addConstr(y <= ((xlb+xub)/2)**2)]

            if x not in alphas:
                alphas[x] = mdl.addVars(J, name=f'{xname}_alpha')
                Jg = (J if sense in (UB, APPROX) else J1)
                gs[x] = mdl.addVars(Jg, lb=0, ub=1, name=f'{xname}_g')
                gs[x].update({0: 0})

    def sawtooth_lb(self, x):
        # Model y>=x^2 with sawtooth relaxation.

        mdl = self.mdl
        L1 = self.L1
        J1 = self.J1
        J = self.J
        xsqc = self.xsq_constrs[x]
        xsq_vdict = self.xsq_vdict

        gs = self.gs
        xname = x.varname

        xlb = x.lb
        xub = x.ub

        lx = xub-xlb

        gdef = (x in gs)

        if lx >= 1e-10:
            # In the current implementation, gdef will always be false. But this code may support more advanced future
            #    alternatives, where the UB and LB parts can be added separately as needed, in the future.
            Jg1 = (J1 if not gdef else [j for j in J1 if j not in J])
            if not gdef:
                gs[x] = mdl.addVars(J1, lb=0, ub=1, name=f'{xname}_g')
                gs[x].update({0: (x-xlb)/lx})
            else:
                newv = mdl.addVars(Jg1, lb=0, ub=1, name=f'{xname}_g')
                gs[x].update({j: newv[j] for j in Jg1})

            y = xsq_vdict[x]
            g = gs[x]

            yh = (y-xlb*(2*x-xlb))/lx**2
            xh = g[0]

            c1 = mdl.addConstrs(
                (g[j] <= 2*g[j-1]
                 for j in Jg1)
            )
            c2 = mdl.addConstrs(
                (g[j] <= 2*(1 - g[j-1])
                 for j in Jg1)
            )
            xsqc += [c1, c2]

            c1 = mdl.addConstrs(
                (yh >= xh - sum(2**(-2*j)*g[j] for j in range(1, l+1)) - 2**(-2*l - 2)
                 for l in range(L1+1)))
            c2 = mdl.addConstr(yh >= 0)
            c3 = mdl.addConstr(yh >= 2*xh-1)

            xsqc += [c1, c2, c3]
        else:
            Jg1 = (J1 if not gdef else [j for j in J1 if j not in J])
            y = xsq_vdict[x]
            if not gdef:
                gs[x] = mdl.addVars(J1, lb=0, ub=1, name=f'{xname}_g')
                gs[x].update({0: 0})
            else:
                newv = mdl.addVars(Jg1, lb=0, ub=1, name=f'{xname}_g')
                gs[x].update({j: newv[j] for j in Jg1})
            xsqc += [mdl.addConstr(y == ((xlb+xub)/2)**2)]

    def sawtooth_p(self, key, x, xlb, xub, is_plus=True, varname=''):
        # Model x=p^2 with sawtooth formulation, where x is an expression, not a model variable.
        #    key specifies the dictionary key to store the g-variables in.
        #    Meant for p of the form p=x+y or p=x-y, where x,y are model variables, from univar-type formulations
        #    For now, used only in univar formulation
        mdl = self.mdl
        L = self.L
        L1 = self.L1
        J = self.J
        J1 = self.J1
        xyc = self.xy_constrs[key]
        alphas = self.alphas

        if is_plus:
            zs = self.zplus
            gs = self.zplus_g
        else:
            zs = self.zminus
            gs = self.zminus_g

        lx = xub - xlb

        zs[key] = mdl.addVar(lb=-GRB.INFINITY, name=f'{varname}_z')

        y = zs[key]

        if lx >= 1e-10:
            alphas[key] = mdl.addVars(J, vtype=GRB.BINARY, name=f'{varname}_alpha')
            alpha = alphas[key]
            gs[key] = mdl.addVars(J1, lb=0, ub=1, name=f'{varname}_g')
            gs[key].update({0: (x - xlb)/lx})

            g = gs[key]
            yh = (y - xlb*(2*x - xlb))/lx**2
            xh = g[0]

            c1 = mdl.addConstrs(
                (g[j] <= 2*g[j-1]
                 for j in J1)
            )
            c2 = mdl.addConstrs(
                (g[j] <= 2*(1 - g[j-1])
                 for j in J1)
            )
            c3 = mdl.addConstrs(
                (2*(g[j-1] - alpha[j]) <= g[j]
                 for j in J)
            )
            c4 = mdl.addConstrs(
                (2*(alpha[j] - g[j-1]) <= g[j]
                 for j in J)
            )
            c5 = mdl.addConstrs(
                (yh >= xh - sum(2**(-2*j)*g[j] for j in range(1, l+1)) - 2**(-2*l - 2)
                 for l in range(L1+1)))
            c6 = mdl.addConstr(yh <= xh - sum(2**(-2*j)*g[j] for j in J))
            c7 = mdl.addConstr(yh >= 0)
            c8 = mdl.addConstr(yh >= 2*xh-1)

            xyc += [c1, c2, c3, c4, c5, c6, c7, c8]
        else:
            gs[key] = mdl.addVars(range(1, L1+1), lb=0, ub=1, name=f'{varname}_g')
            gs[key].update({0: 0})
            xyc += [mdl.addConstr(y == ((xlb+xub)/2)**2)]

    def sawtooth_lb_p(self, key, x, xlb, xub, is_plus=True, varname=''):
        # Model x>=p^2 with sawtooth formulation, where x is an expression, not a model variable.
        #    key specifies the dictionary key to store the g-variables in.
        #    Meant for p of the form p=x+y or p=x-y, where x,y are model variables, from univar-type formulations
        #    For now, used only in univar formulation
        mdl = self.mdl
        L1 = self.L1
        J1 = self.J1
        xyc = self.xy_constrs[key]

        if is_plus:
            zs = self.zplus
            gs = self.zplus_g
        else:
            zs = self.zminus
            gs = self.zminus_g

        lx = xub - xlb

        zs[key] = mdl.addVar(lb=-GRB.INFINITY, name=f'{varname}_z')
        y = zs[key]

        if lx >= 1e-10:
            gs[key] = mdl.addVars(J1, lb=0, ub=1, name=f'{varname}_g')
            gs[key].update({0: (x - xlb)/lx})

            g = gs[key]
            yh = (y - xlb*(2*x - xlb))/lx**2
            xh = g[0]

            c1 = mdl.addConstrs(
                (g[j] <= 2*g[j-1]
                 for j in J1)
            )
            c2 = mdl.addConstrs(
                (g[j] <= 2*(1 - g[j-1])
                 for j in J1)
            )
            c3 = mdl.addConstrs(
                (yh >= xh - sum(2**(-2*j)*g[j] for j in range(1, l+1)) - 2**(-2*l - 2)
                 for l in range(L1+1)))
            c4 = mdl.addConstr(yh >= 0)
            c5 = mdl.addConstr(yh >= 2*xh - 1)

            xyc += [c1, c2, c3, c4, c5]
        else:
            gs[key] = mdl.addVars(range(1, L1+1), lb=0, ub=1, name=f'{varname}_g')
            gs[key].update({0: 0})
            xyc += [mdl.addConstr(y == ((xlb+xub)/2)**2)]

    def add_univar(self, x, y):
        mdl = self.mdl
        alphas = self.alphas
        xy_constrs = self.xy_constrs
        xy_vdict = self.xy_vdict
        xsq_vdict = self.xsq_vdict

        zp1s = self.zplus
        zp2s = self.zminus
        isApprox = self.st_do_approx

        xy_constrs[x, y] = []
        xyc = xy_constrs[x, y]

        APPROX = Model_qcqp.ST_APPROX
        RELAX = Model_qcqp.ST_RELAX

        sense = (APPROX if isApprox else RELAX)

        xname = x.varname
        yname = y.varname

        sawtooth = self.sawtooth
        for var in [x, y]:
            if var not in alphas:
                xsq_vdict[var] = mdl.addVar(lb=-GRB.INFINITY, name=f'{var.varName}_sq')
                sawtooth(var, sense=sense)

        xlb = x.lb
        xub = x.ub
        ylb = y.lb
        yub = y.ub

        p1 = x+y
        p2 = x-y

        p1lb = xlb+ylb
        p1ub = xub+yub
        p2lb = xlb-yub
        p2ub = xub-ylb

        sawtooth_lb = self.sawtooth_lb_p

        sawtooth_lb((x, y), p1, p1lb, p1ub, is_plus=True, varname=f'{xname}_times_{yname}_plus')
        sawtooth_lb((x, y), p2, p2lb, p2ub, is_plus=False, varname=f'{xname}_times_{yname}_minus')

        zx = xsq_vdict[x]
        zy = xsq_vdict[y]
        zp1 = zp1s[x, y]
        zp2 = zp2s[x, y]

        z = xy_vdict[x, y]

        c1 = mdl.addConstr(z >= 0.5*(zp1-zx-zy))
        c2 = mdl.addConstr(z <= 0.5*(zx+zy-zp2))

        xyc += [c1, c2]
        xyc += self.McCormick(x, y, z, [xlb, xub], [ylb, yub])

    def add_UNIVAR_NOSTUB(self, x, y):
        mdl = self.mdl
        xy_constrs = self.xy_constrs
        xy_vdict = self.xy_vdict
        xsq_vdict = self.xsq_vdict

        zp1s = self.zplus
        zp2s = self.zminus

        xy_constrs[x, y] = []
        xyc = xy_constrs[x, y]

        xname = x.varname
        yname = y.varname

        addsq_stlb = self.addsq_stlb
        for var in [x, y]:
            if var not in xsq_vdict:
                xsq_vdict[var] = mdl.addVar(lb=-GRB.INFINITY, name=f'{var.varName}_sq')
                addsq_stlb(var)

        xlb = x.lb
        xub = x.ub
        ylb = y.lb
        yub = y.ub

        p1 = x + y
        p2 = x - y

        p1lb = xlb + ylb
        p1ub = xub + yub
        p2lb = xlb - yub
        p2ub = xub - ylb

        sawtooth_lb = self.sawtooth_lb_p

        sawtooth_lb((x, y), p1, p1lb, p1ub, is_plus=True, varname=f'{xname}_times_{yname}_plus')
        sawtooth_lb((x, y), p2, p2lb, p2ub, is_plus=False, varname=f'{xname}_times_{yname}_minus')

        zx = xsq_vdict[x]
        zy = xsq_vdict[y]
        zp1 = zp1s[x, y]
        zp2 = zp2s[x, y]

        z = xy_vdict[x, y]

        c1 = mdl.addConstr(z >= 0.5 * (zp1 - zx - zy))
        c2 = mdl.addConstr(z <= 0.5 * (zx + zy - zp2))

        xyc += [c1, c2]
        xyc += self.McCormick(x, y, z, [xlb, xub], [ylb, yub])

    def add_UNIVAR_NOST(self, x, y):
        mdl = self.mdl
        xy_constrs = self.xy_constrs
        xy_vdict = self.xy_vdict
        xsq_vdict = self.xsq_vdict

        zp1s = self.zplus
        zp2s = self.zminus

        xy_constrs[x, y] = []
        xyc = xy_constrs[x, y]

        addsq = self.addsq
        for var in [x, y]:
            if var not in xsq_vdict:
                xsq_vdict[var] = mdl.addVar(lb=-GRB.INFINITY, name=f'{var.varName}_sq')
                addsq(x)

        p1 = x + y
        p2 = x - y

        zx = xsq_vdict[x]
        zy = xsq_vdict[y]

        zp1s[x, y] = mdl.addVar(lb=-GRB.INFINITY, name=f'{x.varname}_times_{y.varname}_plus_z')
        zp2s[x, y] = mdl.addVar(lb=-GRB.INFINITY, name=f'{x.varname}_times_{y.varname}_minus_z')

        zp1 = zp1s[x, y]
        zp2 = zp2s[x, y]

        c1 = mdl.addConstr(zp1 >= p1**2)
        c2 = mdl.addConstr(zp2 >= p2**2)

        z = xy_vdict[x, y]

        c3 = mdl.addConstr(z >= 0.5 * (zp1 - zx - zy))
        c4 = mdl.addConstr(z <= 0.5 * (zx + zy - zp2))

        xyc += [c1, c2, c3, c4]

    def add_bin(self, x, y, is_plus=True):
        mdl = self.mdl
        alphas = self.alphas
        xy_constrs = self.xy_constrs
        xy_vdict = self.xy_vdict
        xsq_vdict = self.xsq_vdict

        if is_plus:
            zps = self.zplus
        else:
            zps = self.zminus

        xy_constrs[x, y] = []
        xyc = xy_constrs[x, y]

        RELAX = Model_qcqp.ST_RELAX
        sense = RELAX

        xname = x.varname
        yname = y.varname

        sawtooth = self.sawtooth
        for var in [x, y]:
            if var not in alphas:
                xsq_vdict[var] = mdl.addVar(lb=-GRB.INFINITY, name=f'{var.varName}_sq')
                sawtooth(var, sense=sense)

        xlb = x.lb
        xub = x.ub
        ylb = y.lb
        yub = y.ub

        if is_plus:
            p = x+y
            plb = xlb+ylb
            pub = xub+yub
        else:
            p = x-y
            plb = xlb-yub
            pub = xub-ylb

        self.sawtooth_p((x, y), p, plb, pub, is_plus=is_plus, varname=f'{xname}_times_{yname}_zp')

        zx = xsq_vdict[x]
        zy = xsq_vdict[y]
        zp = zps[x, y]

        z = xy_vdict[x, y]

        if is_plus:
            xyc.append(mdl.addConstr(z == 0.5*(zp-zx-zy)))
        else:
            xyc.append(mdl.addConstr(z == 0.5*(zx+zy-zp)))

        xyc += self.McCormick(x, y, z, [xlb, xub], [ylb, yub])

    # Definitions related to forms of NMDT
    def McCormick(self, x, y, z, bndx, bndy):
        # McCormick envelopes for z=xy
        # get model, bounds
        mdl = self.mdl
        lx, ux = bndx
        ly, uy = bndy

        # Add McCormick envelopes
        c1 = mdl.addConstr(z >= x*ly + y*lx - lx*ly)
        c2 = mdl.addConstr(z >= x*uy + y*ux - ux*uy)
        c3 = mdl.addConstr(z <= x*ly + y*ux - ux*ly)
        c4 = mdl.addConstr(z <= x*uy + y*lx - lx*uy)

        return [c1, c2, c3, c4]

    def McCormicksq(self, x, y, bndx):
        # get model, bounds
        mdl = self.mdl
        lx, ux = bndx

        # Add McCormick envelopes
        c1 = mdl.addConstr(y >= 2*x*lx - lx**2)
        c2 = mdl.addConstr(y >= 2*x*ux - ux**2)
        c3 = mdl.addConstr(y <= x*(lx + ux) - ux*lx)

        #c1 = mdl.addConstr(y >= x**2)
        #c2 = mdl.addConstr(y <= x*(lx + ux) - ux*lx)

        return [c1, c2, c3]

    def McCormickub(self, x, y, z, bndx, bndy):
        # McCormick envelopes for z=xy
        # get model, bounds
        mdl = self.mdl
        lx, ux = bndx
        ly, uy = bndy

        # Add McCormick envelopes
        c1 = mdl.addConstr(z <= x*ly + y*ux - ux*ly)
        c2 = mdl.addConstr(z <= x*uy + y*lx - lx*uy)

        return [c1, c2]

    def McCormicksqub(self, x, y, bndx):
        # get model, bounds
        mdl = self.mdl
        lx, ux = bndx

        # Add UB McCormick envelopes
        c1 = mdl.addConstr(y <= x*(lx + ux) - ux*lx)

        return [c1]

    def add_DNMDT(self, x, y):
        # Map xh,xy to [0,1], then apply DNMDT to the mapped interval, to model z=xy
        mdl = self.mdl
        L = self.L
        J = self.J
        ldas = self.ldas
        alphas = self.alphas
        delxs = self.delxs
        xy_constrs = self.xy_constrs
        xy_vdict = self.xy_vdict

        xname = x.varname
        yname = y.varname

        # TODO: do this step in update method instead.
        #   for key in xy_constrs: for constr in xy_constrs[key]: mdl.remove(constr)
        if (x, y) in xy_constrs and 0:
            for constr in xy_constrs[x, y]:
                mdl.remove(constr)

        xy_constrs[x, y] = []
        xyc = xy_constrs[x, y]

        xadd = x not in alphas
        yadd = y not in alphas
        for var in (x, y):
            if var not in alphas:
                alphas[var] = mdl.addVars(J, vtype=GRB.BINARY, name=f'{var.varname}_alpha')
                delxs[var] = mdl.addVar(lb=0, ub=2**(-L), name=f'{var.varname}_delx')

        vs = self.vs
        ws = self.ws
        delzs = self.delzs

        delzs[x, y] = mdl.addVar(lb=-GRB.INFINITY, name=f'{xname}_{yname}_delz')

        for lda in ldas:
            # Prepare to model both lda=0 and lda=1
            vs[x, y, lda] = mdl.addVars(J, lb=-GRB.INFINITY, name=f'{xname}_{yname}_v_{lda}')
            ws[x, y, lda] = mdl.addVars(J, lb=-GRB.INFINITY, name=f'{xname}_{yname}_w_{lda}')

        delzh = delzs[x, y]

        alpha = alphas[x]
        beta = alphas[y]

        delxh = delxs[x]
        delyh = delxs[y]

        z = xy_vdict[x, y]

        xlb = x.lb
        ylb = y.lb
        xub = x.ub
        yub = y.ub

        lx = xub-xlb
        ly = yub-ylb

        tol = 1e-10
        if lx <= tol and ly <= tol:
            xmid = 0.5*(xub+xlb)
            ymid = 0.5*(yub+ylb)
            xyc += [mdl.addConstr(z == xmid*ymid)]
        elif lx <= tol:
            xmid = 0.5*(xub+xlb)
            xyc += [mdl.addConstr(z == xmid*y)]
            if yadd:
                yh = (y-ylb)/ly
                xyc += [mdl.addConstr(yh == sum(2**(-j)*beta[j] for j in J) + delyh)]
        elif ly <= tol:
            ymid = 0.5*(yub+ylb)
            xyc += [mdl.addConstr(z == ymid*x)]
            if xadd:
                xh = (x-xlb)/lx
                xyc += [mdl.addConstr(xh == sum(2**(-j)*alpha[j] for j in J) + delxh)]
        else:
            xh = (x-xlb)/lx
            yh = (y-ylb)/ly
            zh = (z-lx*xh*ylb-ly*yh*xlb-xlb*ylb)/(lx*ly)

            if xadd:
                xyc += [mdl.addConstr(xh == sum(2**(-j)*alpha[j] for j in J) + delxh)]
            if yadd:
                xyc += [mdl.addConstr(yh == sum(2**(-j)*beta[j] for j in J) + delyh)]

            McCormick = self.McCormick
            for lda in ldas:
                v = vs[x, y, lda]
                w = ws[x, y, lda]
                xyc += [mdl.addConstr(zh == sum(2**(-j)*(v[j]+w[j]) for j in J) + delzh)]

                for j in J:
                    # Apply DNMDT with the user-specified value of lambda
                    xyc += McCormick(lda*delyh+(1-lda)*yh, alpha[j], v[j], [0, lda*2**(-L) + (1-lda)], [0, 1])
                    xyc += McCormick((1-lda)*delxh+lda*xh, beta[j], w[j], [0, (1-lda)*2**(-L) + lda], [0, 1])
            xyc += McCormick(delxh, delyh, delzh, [0, 2**(-L)], [0, 2**(-L)])

    def add_DNMDTsq(self, x):
        # Apply the y=x^2 version of DNMDT. Note that there is no longer any need for lda.
        mdl = self.mdl
        L = self.L
        J = self.J
        alphas = self.alphas
        delxs = self.delxs
        xsq_constrs = self.xsq_constrs
        xsq_vdict = self.xsq_vdict

        do_tight_lb = self.tighten_nmdt
        xname = x.varname

        # TODO: do this step in update method instead.
        #   for key in xy_constrs: for constr in xy_constrs[key]: mdl.remove(constr)
        if x in xsq_constrs and 0:
            for constr in xsq_constrs[x]:
                mdl.remove(constr)

        xsq_constrs[x] = []
        xsqc = xsq_constrs[x]

        xadd = x not in alphas
        if xadd:
            alphas[x] = mdl.addVars(J, vtype=GRB.BINARY, name=f'{xname}_alpha')
            delxs[x] = mdl.addVar(lb=0, ub=2**(-L), name=f'{xname}_delx')

        vs = self.vs
        delzs = self.delzs

        vs[x] = mdl.addVars(J, lb=-GRB.INFINITY, name=f'{xname}_v')
        delzs[x] = mdl.addVar(lb=-GRB.INFINITY, name=f'{xname}_delz')

        v = vs[x]
        delz = delzs[x]
        alpha = alphas[x]
        delx = delxs[x]

        z = xsq_vdict[x]

        xlb = x.lb
        xub = x.ub

        lx = xub-xlb

        if xadd:
            xsqc += [mdl.addConstr(x == lx * sum(2**(-j)*alpha[j] for j in J) + lx*delx + xlb)]
        xsqc += [mdl.addConstr(z == lx * sum(2**(-j)*v[j] for j in J) + lx**2*delz + xlb*(x+lx*delx))]

        if not do_tight_lb:
            McCormicksq = self.McCormicksq
            McCormick = self.McCormick
        else:
            # Only use McCormick for UB; use sawtooth relaxation for LB
            McCormicksq = self.McCormicksqub
            McCormick = self.McCormickub

        for j in J:
            xsqc += McCormick(lx*delx + x, alpha[j], v[j], [xlb, xub + 2**(-L)*lx], [0, 1])
        xsqc += McCormicksq(delx, delz, [0, 2**(-L)])

        # Tighten relaxation to improve competitiveness
        if do_tight_lb:
            self.sawtooth_lb(x)

    def add_NMDT(self, x, y):
        # Apply NMDT to model z=xy. Left or right variable depends on user parameters.

        mdl = self.mdl
        L = self.L
        J = self.J
        disc_left = self.disc_left
        disc_right = self.disc_right
        alphas = self.alphas
        delxs = self.delxs
        xy_constrs = self.xy_constrs
        xy_vdict = self.xy_vdict

        xname = x.varname
        yname = y.varname

        # TODO: do this step in update method instead.
        #   for key in xy_constrs: for constr in xy_constrs[key]: mdl.remove(constr)
        if (x, y) in xy_constrs and 0:
            for constr in xy_constrs[x, y]:
                mdl.remove(constr)

        xy_constrs[x, y] = []
        xyc = xy_constrs[x, y]

        vs = self.vs
        delzs = self.delzs

        z = xy_vdict[x, y]

        McCormick = self.McCormick

        xlb = x.lb
        xub = x.ub
        ylb = y.lb
        yub = y.ub
        if disc_left:
            lx = xub-xlb

            if x not in alphas:
                alphas[x] = mdl.addVars(J, vtype=GRB.BINARY, name=f'{xname}_alpha')
                delxs[x] = mdl.addVar(lb=0, ub=lx*2**(-L), name=f'{xname}_delx')

            vs[x, y] = mdl.addVars(J, lb=-GRB.INFINITY, name=f'{xname}_{yname}_v')
            v = vs[x, y]
            alpha = alphas[x]
            delx = delxs[x]

            delzs[x, y] = mdl.addVar(lb=-GRB.INFINITY, name=f'{xname}_{yname}_delz')
            delz = delzs[x, y]

            c1 = mdl.addConstr(x == lx*sum(2**(-j)*alpha[j] for j in J) + delx + xlb)
            c2 = mdl.addConstr(z == lx*sum(2**(-j)*v[j] for j in J) + delz + xlb*y)
            xyc += [c1, c2]

            for j in J:
                # Apply McCormick envelopes
                xyc += McCormick(y, alpha[j], v[j], [ylb, yub], [0, 1])
            xyc += McCormick(delx, y, delz, [0, lx*2**(-L)], [ylb, yub])

        if disc_right:
            ly = yub-ylb

            if y not in alphas:
                alphas[y] = mdl.addVars(J, vtype=GRB.BINARY, name=f'{yname}_alpha')
                delxs[y] = mdl.addVar(lb=0, ub=ly*2**(-L), name=f'{yname}_delx')

            vs[y, x] = mdl.addVars(J, lb=-GRB.INFINITY, name=f'{yname}_{xname}_v')
            delzs[y, x] = mdl.addVar(lb=-GRB.INFINITY, name=f'{yname}_{xname}_delz')

            alpha = alphas[y]
            dely = delxs[y]
            v = vs[y, x]
            delz = delzs[y, x]

            c1 = mdl.addConstr(y == ly*sum(2**(-j)*alpha[j] for j in J) + dely + ylb)
            c2 = mdl.addConstr(z == ly*sum(2**(-j)*v[j] for j in J) + delz + ylb*x)
            xyc += [c1, c2]

            for j in J:
                # Apply McCormick envelopes
                xyc += McCormick(x, alpha[j], v[j], [xlb, xub], [0, 1])
            xyc += McCormick(dely, x, delz, [0, ly*2**(-L)], [xlb, xub])

    def add_NMDTsq(self, x):
        # Apply the y=x^2 version of NMDT.
        mdl = self.mdl
        L = self.L
        J = self.J
        alphas = self.alphas
        delxs = self.delxs
        xsq_constrs = self.xsq_constrs
        xsq_vdict = self.xsq_vdict

        xname = x.varname

        # TODO: do this step in update method instead.
        #   for key in xy_constrs: for constr in xy_constrs[key]: mdl.remove(constr)
        if x in xsq_constrs and 0:
            for constr in xsq_constrs[x]:
                mdl.remove(constr)

        xsq_constrs[x] = []
        xsqc = xsq_constrs[x]

        xlb = x.lb
        xub = x.ub

        lx = xub-xlb
        if x not in alphas:
            alphas[x] = mdl.addVars(J, vtype=GRB.BINARY, name=f'{xname}_alpha')
            delxs[x] = mdl.addVar(lb=0, ub=lx*2**(-L), name=f'{xname}_delx')

        vs = self.vs
        delzs = self.delzs

        vs[x] = mdl.addVars(J, lb=-GRB.INFINITY, name=f'{xname}_v')
        delzs[x] = mdl.addVar(lb=-GRB.INFINITY, name=f'{xname}_delz')

        v = vs[x]
        delz = delzs[x]
        alpha = alphas[x]
        delx = delxs[x]

        z = xsq_vdict[x]

        c1 = mdl.addConstr(x == lx * sum(2**(-j)*alpha[j] for j in J) + delx + xlb)
        c2 = mdl.addConstr(z == lx * sum(2**(-j)*v[j] for j in J) + delz + xlb*x)

        do_tight_lb = self.tighten_nmdt
        if not do_tight_lb:
            McCormick = self.McCormick
        else:
            McCormick = self.McCormickub
        xsqc += [c1, c2]
        for j in J:
            xsqc += McCormick(x, alpha[j], v[j], [xlb, xub], [0, 1])
        xsqc += McCormick(delx, x, delz, [0, 2**(-L)*lx], [xlb, xub])

        if do_tight_lb:
            self.sawtooth_lb(x)

    def add_McCormick(self, x, y):
        # Apply McCormick envelopes to model z=xy.
        mdl = self.mdl
        xy_constrs = self.xy_constrs
        xy_vdict = self.xy_vdict

        # TODO: do this step in update method instead.
        #   for key in xy_constrs: for constr in xy_constrs[key]: mdl.remove(constr)
        if (x, y) in xy_constrs and 0:
            for constr in xy_constrs[x, y]:
                mdl.remove(constr)

        xy_constrs[x, y] = []
        xyc = xy_constrs[x, y]

        z = xy_vdict[x, y]

        xyc += self.McCormick(x, y, z, [x.lb, x.ub], [y.lb, y.ub])

    def add_McCormicksq(self, x):
        # Apply McCormick envelopes to model y=x^2.
        mdl = self.mdl
        xsq_constrs = self.xsq_constrs
        xsq_vdict = self.xsq_vdict

        # TODO: do this step in update method instead.
        #   for key in xy_constrs: for constr in xy_constrs[key]: mdl.remove(constr)
        if x in xsq_constrs and 0:
            for constr in xsq_constrs[x]:
                mdl.remove(constr)

        xsq_constrs[x] = []
        xsqc = xsq_constrs[x]

        z = xsq_vdict[x]

        xsqc += self.McCormicksq(x, z, [x.lb, x.ub])

    def addsq_stlb(self, x):
        # Apply McCormick envelopes to model y=x^2.
        mdl = self.mdl
        xsq_constrs = self.xsq_constrs
        xsq_vdict = self.xsq_vdict

        # TODO: do this step in update method instead.
        #   for key in xy_constrs: for constr in xy_constrs[key]: mdl.remove(constr)
        if x in xsq_constrs and 0:
            for constr in xsq_constrs[x]:
                mdl.remove(constr)

        xsq_constrs[x] = []
        xsqc = xsq_constrs[x]

        z = xsq_vdict[x]

        xsqc.append(mdl.addConstr(z <= x**2))
        xsqc.append(self.sawtooth_lb(x))

    def addsq(self, x):
        # Apply McCormick envelopes to model y=x^2.
        mdl = self.mdl
        xsq_constrs = self.xsq_constrs
        xsq_vdict = self.xsq_vdict

        # TODO: do this step in update method instead.
        #   for key in xy_constrs: for constr in xy_constrs[key]: mdl.remove(constr)
        if x in xsq_constrs and 0:
            for constr in xsq_constrs[x]:
                mdl.remove(constr)

        xsq_constrs[x] = []
        xsqc = xsq_constrs[x]

        z = xsq_vdict[x]

        xsqc.append(mdl.addConstr(z == x ** 2))