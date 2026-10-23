"""Polynomial differential reduction and exploratory high-degree systems."""

from pathlib import Path
import sys
import sympy as sp
from mpmath import mp
from generapy import Curve

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"code"))
from hyperelliptic_scalar_first.special_q.model import scalar_data, canonical_field
from hyperelliptic_scalar_first.sigma_modes.model import SigmaModes


def second_polynomials(F,x,g):
    return [sp.expand(sum((m+1-k)*F.coeff(x,m+1+k)*x**m/4
                         for m in range(k,2*g+2-k))) for k in range(1,g+1)]


def reduce_differential(B,F,x):
    """B dx/y = d(P*y)+sum c_k omega_k+sum d_k dr_k, exactly."""
    F=sp.expand(F)
    degree=sp.degree(F,x)
    if degree%2!=1 or degree<3:
        raise ValueError("Use an odd-degree curve of positive genus")
    g=(degree-1)//2
    remainder,P=sp.Poly(B,x),sp.Integer(0)
    while not remainder.is_zero and remainder.degree()>=2*g:
        k=remainder.degree()-2*g
        term=x**k
        image=sp.Poly(sp.diff(term,x)*F+term*sp.diff(F,x)/2,x)
        multiplier=remainder.LC()/image.LC()
        P+=multiplier*term
        remainder-=sp.Poly(multiplier,x)*image
    second=second_polynomials(F,x,g)
    ds=[]
    for form in second:  # descending leading degrees 2g-1,...,g
        degree=int(sp.degree(form,x))
        coefficient=remainder.nth(degree)/sp.Poly(form,x).LC()
        ds.append(coefficient)
        remainder-=sp.Poly(coefficient*form,x)
    cs=[remainder.nth(k) for k in range(g)]
    residual=sp.expand(B-(sp.diff(P,x)*F+P*sp.diff(F,x)/2)
                       -sum(cs[k]*x**k for k in range(g))
                       -sum(d*form for d,form in zip(ds,second)))
    assert residual==0
    return dict(P=sp.expand(P),c=cs,d=ds,residual=residual)


def systems():
    x=sp.Symbol("rho")
    specifications=[("g2_high_q_high_B",(-7,-5,-3,-1,1,3,5,7),2*x**4-84*x**2+x+210,5,8),
                    ("g3_fixed_curve_high_B",(-3,-2,-1,0,1,2,3),sp.Integer(1),8,4)]
    result=[]
    for name,gamma,q,K,start in specifications:
        variables=sp.symbols(f"I1:{len(gamma)+1}")
        scalar=-sum(variables)/len(gamma)
        L=variables[0]-variables[1]-(gamma[0]-gamma[1])
        Q=q.subs(x,scalar)+L*scalar**K
        data=scalar_data(Q,variables,(1,)*len(gamma),gamma,x)
        assert sp.Poly(data["F"],x).LC()==4 and data["squarefree"]
        data.update(name=name,rho_start=sp.Integer(start),rho_end=sp.Integer(start)+sp.Rational(1,500))
        result.append(data)
    return result


class PointKernels(SigmaModes):
    """Use the existing point-kernel methods without imposing B=0.

    Only initialization is local: production SigmaModes is a special-Q
    verifier. Here no correction is fitted or discarded; actual B is
    integrated separately by reduce_differential.
    """
    def __init__(self,data,dps=40,progress=print):
        self.data=data
        self.ctx=c=mp.clone(); c.dps=dps
        self.genus=g=data["genus"]; self.n=len(data["beta"])
        assert g in (2,3)
        self.scale=c.one
        self.normalization=dict(scale=sp.Integer(1),F=data["F"])
        self.sharp=tuple(int(k==1) for k in range(g))
        self.flat=(0,0) if g==2 else (0,0,1)
        self.a=[self.number(a) for a in data["a"]]
        self.b_original=self.b=[self.number(b) for b in data["ordinate"]]
        self.rho_start=self.number(data["rho_start"]); self.rho_end=self.number(data["rho_end"])
        x=data["rho"]
        self.F_original=self.F=self.compile(x,data["F"])
        self.q=self.compile(x,data["q"])
        coefficients=[self.number(data["F"].coeff(x,k)) for k in range(data["degree"]+1)]
        self.curve=Curve({(0,2):1}|{(k,0):-v for k,v in enumerate(coefficients) if v},ctx=c)
        progress(f"{data['name']}: periods and compatible fixed lifts")
        first=self.curve.periods_kind_1(); second=self.curve.periods_kind_2()
        self.periods=first
        self.sigma_options=dict(omega=first.omega,tau=first.tau,kappa=second.kappa,
                                characteristic=self.curve.riemann_constant().characteristic,
                                normalization="hyperelliptic",ctx=c)
        self.anchor_Y=c.sqrt(self.F(self.rho_start))
        anchor=(self.rho_start,self.anchor_Y)
        self.anchor_w=self.curve.abel_map_kind_1(anchor,reduce=False).value
        self.anchor_second=self.curve.abel_map_kind_2(anchor,reduce=False).value
        self.roots=[c.mpc(str(sp.re(r)),str(sp.im(r))) for r in sp.Poly(data["F"],x).nroots(n=dps+5,maxsteps=1000)]
        self.height=1+max(abs(r.imag) for r in self.roots)
        self.second_forms=[self.compile(x,p) for p in second_polynomials(data["F"],x,g)]
        self.eta=[]; self.second_at_fixed=[]
        for a,b in zip(self.a,self.b):
            dw,dR,end_y=self.fixed_path(a)
            sign=1 if abs(end_y-b)<abs(end_y+b) else -1
            assert abs(sign*end_y-b)<c.mpf(10)**(-(dps-8))*max(1,abs(b))
            self.eta.append(sign*(self.anchor_w+dw))
            self.second_at_fixed.append(sign*(self.anchor_second+dR))
        self.normals=[]; self.h=[]
        for eta in self.eta:
            jet=self.jet(eta,order=2); normal=jet[self.sharp]
            zeta=[jet[self.extend(self.sharp,k)]/normal for k in range(g)]
            h2=zeta[1]/2
            h=[zeta[0]-jet[(1,)+(0,)*(g-1)]*h2/normal,h2]
            if g==3: h.append(zeta[2])
            self.normals.append(normal); self.h.append(h)


def columns(g,d,x,y):
    if d<=g:
        return [x**k for k in range(g+1)]+[y/2]
    return [f for _,f in sorted([(2*k,x**k) for k in range(d+1)]+
            [(2*g+1+2*k,x**k*y/2) for k in range(d-g)],key=lambda item:item[0])]


def interpolation_constant(g,d,leading):
    if d<=g: return sp.Integer(1)
    m=d-g-1
    return (-1)**(1+m*(m+1)//2)*(leading/2)**m
