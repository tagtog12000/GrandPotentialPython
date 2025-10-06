from google.colab import files
import ast
import math
import vegas
from collections import Counter
import numpy as np
from itertools import product
import numdifftools as nd

!rm -rf GrandPotentialPython
!git clone https://github.com/tagtog12000/GrandPotentialPython.git
%cd GrandPotentialPython
!g++ GP.cpp -o code
!./code

def divided_difference_num(fnl, nodes, b, tol=1e-12):
    nodes = list(map(float, nodes))
    n = len(nodes)
    if n == 1:
        return fnl(nodes[0], b)
    if all(abs(nodes[i] - nodes[0]) < tol for i in range(1, n)):
        m = n - 1
        Dm = nd.Derivative(lambda x: fnl(x, b), n=m)
        return Dm(nodes[0]) / math.factorial(m)
    top = divided_difference_num(fnl, nodes[1:], b, tol) # Recursive case
    bot = divided_difference_num(fnl, nodes[:-1], b, tol)
    return (top - bot) / (nodes[-1] - nodes[0])
fnl = lambda x, b: -0.5 + 0.5 * np.tanh(b * x / 2)
eps = 1e-12
def bit_positions(n):
    return [i+1 for i in range(n.bit_length()) if (n >> i) & 1]
def evalNum(intU: int, fun, mE, b):
    val = 1
    pos = 0
    while intU:
        if intU & 1:
            val = val * fun(mE[pos], b)
        intU >>= 1
        pos += 1
    return val
def evalDen(intU: int, mE):
    val = 0
    pos = 0
    while intU:
        if intU & 1:
            val = val + mE[pos]
        intU >>= 1
        pos += 1
    return val
def remplaceVariable(refDens: list, spTr: int):
    brnch = []  # the branches of spanning tree in integers position
    signs = []
    equs = []
    for refDen in refDens:
        posD = refDen[0] ^ refDen[1]
        sign = []
        equ = []
        if posD & spTr: #the branch in the positif side
            br = posD & spTr
            posSTr = br.bit_length()
            brnch.append(posSTr)
            posD = posD ^ br
            equN = bit_positions(posD)
            equP = bit_positions(refDen[1])
            equ = equN + equP
            sign = [-1 for i in range(len(equN))]+[1 for i in range(len(equP))]
        elif refDen[1] & spTr: #the branch in the negative side
            br = refDen[1] & spTr
            posSTr = br.bit_length()
            brnch.append(posSTr)
            refDen[1] = refDen[1] ^ br
            equN = bit_positions(refDen[1])
            equP = bit_positions(posD)
            equ = equN + equP
            sign = [-1 for i in range(len(equN))]+[1 for i in range(len(equP))]
        else:
            print("error")
        equs.append(equ)
        signs.append(sign)
    return [brnch, signs, equs]
def E(x: np.float64)-> np.float64:
    return -2*np.cos(x)
def fn(x: np.float64,b: np.float64)->np.float64:
    return -0.5 + 0.5 * np.tanh(b * x / 2)
def fp(x: np.float64,b: np.float64)->np.float64:
    return  0.5 + 0.5 * np.tanh(b * x / 2)
def kronecker_delta(i: np.int8, j: np.int8) -> np.int8:
    return np.int8(1) if i == j else np.int8(0)
def pot(s1: np.int8, s2: np.int8, s3: np.int8, s4: np.int8) -> np.int8:
    return np.int8(
        (1 - kronecker_delta(s1, s2)) * (
            kronecker_delta(s1, s3) * kronecker_delta(s2, s4)
            - kronecker_delta(s1, s4) * kronecker_delta(s2, s3)
        )
    )
def generate_valid_spins(n: int):
    pairs = [(np.int8(1), np.int8(2)), (np.int8(2), np.int8(1))]
    for combo in product(pairs, repeat=n):
        yield [s for pair in combo for s in pair]
def evaluate_sum(n: int, edges):
    total = 0
    for s in generate_valid_spins(n):
        potE = 1
        for j in range(n):
            s1 = s[2*j]
            s2 = s[2*j + 1]
            s3 = s[2*edges[2*j][1] - 2]
            s4 = s[2*edges[2*j+1][1] - 1]
            potE *= pot(s1, s2, s3, s4)
        total += potE
    return total
def evalNumT(matNSij, cfij, fn, fp, mE, b):
    resN = 0
    for j in range(len(matNSij)):
        resN += cfij[j]*evalNum(matNSij[j][1], fp, mE, b)*evalNum(matNSij[j][0] ^ matNSij[j][1], fn, mE, b)
    return resN
def evalDenT(matDSij, mE):
    resD = 1
    for j in range(len(matDSij)):
        resD *= evalDen(matDSij[j][0] ^ matDSij[j][1], mE)-evalDen(matDSij[j][1], mE)
    return resD
def evalNumTD(matNSij, fn, fp, mE, b):
    resN = 0
    for j in range(len(matNSij)):
        resN += matNSij[j][0]*evalNum(matNSij[j][1], fp, mE, b)*evalNum(matNSij[j][2] ^ matNSij[j][1], fn, mE, b)
    return resN
def evalDenTD(matDSij, fnl, mE, b):
    lisdivD = []
    for j in range(len(matDSij)):
        arg = evalDen(matDSij[j][0] ^ matDSij[j][1], mE)-evalDen(matDSij[j][1], mE)
        lisdivD.append(arg)
    return divided_difference_num(fnl, lisdivD, b)
def redDiag(lisDD0, mE, b):
    ln = len(lisDD0)//4
    matNS = []
    matDS = []
    for ij in range(ln):
        matNS.append(lisDD0[4*ij+0])
        matDS.append(lisDD0[4*ij+2])
    res = 0.0
    for ij in range(len(matNS)):
        resNN = evalNumTD(matNS[ij], fn, fp, mE, b)
        resDD = 1
        for ijk in range(len(matDS[ij])):
            resDD *= evalDenTD(matDS[ij][ijk], fnl, mE, b)
        res += resNN * resDD
    return res
with open("N.dat") as f:
    n = ast.literal_eval(f.read())
with open("matCoefF-"+str(n)+".dat") as f:
    matCoefF = ast.literal_eval(f.read())
with open("matDenNigF-"+str(n)+".dat") as f:
    matDenNigF = ast.literal_eval(f.read())
with open("matNumSignF-"+str(n)+".dat") as f:
    matNumSignF = ast.literal_eval(f.read())
with open("graphs-"+str(n)+".dat") as f:
    graphs = ast.literal_eval(f.read())
with open("spanningTrees-"+str(n)+".dat") as f:
    spanTrees = ast.literal_eval(f.read())
with open("symmetries-"+str(n)+".dat") as f:
    syms = ast.literal_eval(f.read())
with open("refDenominators-"+str(n)+".dat") as f:
    refDens = ast.literal_eval(f.read())
with open("DivDifMB-"+str(n)+".dat") as f:
    lisDD = ast.literal_eval(f.read())
dn = 2*n
ones = (1 << (dn))-1
def integrand(u):
    var = [-np.pi + 2*np.pi*ui for ui in u]
    res_total = 0.0
    iDD = 0
    for i in range(len(graphs)):
        sym = syms[i]
        spTr = spanTrees[i]
        refDen = refDens[i]
        varUInt = ones ^ spTr
        variables = bit_positions(varUInt)
        brnch, signs, equs = remplaceVariable(refDen, spTr)
        dn = 2*n
        xfull = [0.0 for _ in range(dn)]
        for k in range(n+1):
            if k < len(variables):
                xfull[variables[k]-1] = var[k]
        for k in range(n-1):
            for j in range(len(signs[k])):
                xfull[brnch[k]-1] += signs[k][j] * xfull[equs[k][j]-1]
        mE = [E(xfull[j]) for j in range(dn)]
        edges = graphs[i]
        POT = evaluate_sum(n, edges)/sym  # half-filled case
        matNS = matNumSignF[i]
        matDS = matDenNigF[i]
        cfi   = matCoefF[i]
        res = 0.0
        con = False
        for ij in range(len(matDS)):
            for j in range(len(matDS[ij])):
                bl = bit_positions(matDS[ij][j][0])
                if len(bl) == 2:
                    con = True
                    break
            if con:
                break
            resD = evalDenT(matDS[ij], mE)
            resN = evalNumT(matNS[ij], cfi[ij], fn, fp, mE, 1.0)
            if abs(resD) < eps:
                continue  # skip this term, or set resD = np.sign(resD)*eps
            res += resN / resD
        res_total += POT * res
        if con:
            res = redDiag(lisDD[iDD], mE, 1.0)
            iDD += 1
        res_total += POT * res
    return res_total  # jacobian cancels with (2π)^(n+1)
integ = vegas.Integrator([[0,1]]*(n+1)) # dimension = n+1
result = integ(integrand, nitn=5, neval=1000)  # adapt for your problem
print("Integral =", result.mean, "+/-", result.sdev)
