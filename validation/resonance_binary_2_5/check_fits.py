#!/usr/bin/env python3
"""Recompute old/new trace fits and check the exact transfer-matrix resonance."""
import csv
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent

def load(path):
    return [float(v) for v in path.read_text().split()]

def fit(y):
    n=len(y)
    mx=(n-1)/2
    my=sum(y)/n
    cov=sum((i-mx)*(v-my) for i,v in enumerate(y))
    varx=n*(n*n-1)/12
    slope=cov/varx
    residual=sum((v-(my+slope*(i-mx)))**2 for i,v in enumerate(y))
    total=sum((v-my)**2 for v in y)
    return slope, (1-residual/total if total else 1.)

def matmul(a,b):
    return ((a[0][0]*b[0][0]+a[0][1]*b[1][0], a[0][0]*b[0][1]+a[0][1]*b[1][1]),
            (a[1][0]*b[0][0]+a[1][1]*b[1][0], a[1][0]*b[0][1]+a[1][1]*b[1][1]))

def check_resonance(e):
    F=((float(e),-1.),(1.,0.))
    F3=matmul(matmul(F,F),F)
    sign=-1 if e==1 else 1
    assert all(abs(F3[i][j]-(sign if i==j else 0))<1e-12 for i in range(2) for j in range(2))
    A=((e-0.1,-1.),(1.,0.))
    B=matmul(A,F)
    tr=B[0][0]+B[1][1]
    assert abs(tr-(-1.1 if e==1 else -0.9))<1e-12
    assert abs(tr)<2
    return tr

for E in (-1,1):
    old=load(HERE/'historical_traces'/f'energy_{E}_pot_0.1.dat')
    assert len(old)==5000
    s,r2=fit(old)
    print(f'E={E:+d} archived: Nrows={len(old)}, slope={s:.12g}, xi={-2/s:.10g}, R2={r2:.8g}')
    for N in (5000,10000,100000):
        d=HERE/'new_runs'/f'N{N}_E{E:+d}'
        y=load(d/'energy_00000.dat')
        assert len(y)==N+1
        slope,r2=fit(y)
        with (d/'summary.csv').open() as f:
            row=next(csv.DictReader(f))
        assert math.isclose(slope,float(row['slope']),rel_tol=2e-5,abs_tol=1e-16)
        assert math.isclose(r2,float(row['r_squared']),rel_tol=2e-4,abs_tol=1e-10)
        print(f'   new N={N}: mean_lnT_end={y[-1]:.9g}, slope={slope:.6g}, R2={r2:.6g}')
    tr=check_resonance(E)
    print(f'   exact: tr(A*F)={tr:.2f}, |trace|<2, bounded transfer-matrix powers')
print('All trace and exact-matrix checks passed.')
