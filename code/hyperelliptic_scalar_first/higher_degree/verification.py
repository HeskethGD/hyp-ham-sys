"""Higher polynomial corrections and longer Onishi sigma partitions."""

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter
import sympy as sp

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
try:
    from .model import (systems,reduce_differential,second_polynomials,PointKernels,
                        columns,interpolation_constant,canonical_field)
except ImportError:
    from hyperelliptic_scalar_first.higher_degree.model import (systems,reduce_differential,second_polynomials,PointKernels,
                       columns,interpolation_constant,canonical_field)
from hyperelliptic_scalar_first.abel_lifts import integer_lattice_coordinates


def exact_checks():
    rows={}
    x,y=sp.symbols("x y")
    for g in (2,3,4):
        F=4*x**(2*g+1)-x+1
        for k in range(0,2*g+7):
            reduction=reduce_differential(x**k,F,x)
            rows[f"g={g}, B=x^{k}: differential reduction"]=reduction["residual"]==0
    for g,d in ((2,4),(2,5),(2,6),(3,5),(3,6),(4,6)):
        q=3*x**d+2*x+1
        basis=columns(g,d,x,y); n=len(basis); a=list(range(n-1))
        matrix=sp.Matrix([[f.subs({x:point,y:q.subs(x,point)},simultaneous=True) for f in basis] for point in a]+[basis])
        determinant=matrix.det(method="domain-ge")
        delta=sp.prod(a[j]-a[i] for i in range(n-1) for j in range(i+1,n-1))
        constant=interpolation_constant(g,d,sp.Poly(q,x).LC())
        rows[f"g={g}, degree q={d}: longer-column interpolation"]=sp.expand(determinant-constant*delta*(y-q)/2)==0
        opposite=sp.Matrix([[f.subs({x:point,y:-q.subs(x,point)},simultaneous=True) for f in basis] for point in a]+[basis]).det(method="domain-ge")
        opposite_constant=interpolation_constant(g,d,-sp.Poly(q,x).LC())
        rows[f"g={g}, degree q={d}: opposite interpolation"]=sp.expand(opposite-opposite_constant*delta*(y+q)/2)==0
    return rows


def verify_system(data,dps=40,progress=print):
    started=perf_counter()
    m=PointKernels(data,dps,progress)
    c,g,n,x=m.ctx,m.genus,m.n,data["rho"]
    rows=[]; exact={}
    tolerance=c.mpf(10)**(-(dps-14))
    def record(name,actual,expected):
        error=abs(actual-expected)/max(c.one,abs(actual),abs(expected))
        rows.append(dict(check=name,scaled_error=c.nstr(error,12),passed=bool(error<tolerance)))
    reductions=[reduce_differential(B,data["F"],x) for B in data["B"]]
    P=[m.compile(x,r["P"]) for r in reductions]
    B=[m.compile(x,b) for b in data["B"]]
    Bprime=[m.compile(x,sp.diff(b,x)) for b in data["B"]]
    cs=[[m.number(a) for a in r["c"]] for r in reductions]
    ds=[[m.number(a) for a in r["d"]] for r in reductions]
    exact["all reduced differentials agree"]=all(r["residual"]==0 for r in reductions)
    exact["weighted algebraic endpoint terms cancel"]=sp.expand(sum(r["P"] for r in reductions))==0
    exact["weighted second-kind coefficients cancel"]=all(sp.simplify(sum(r["d"][k] for r in reductions))==0 for k in range(g))
    args,field,H=canonical_field(data)
    field_fn=m.compile(args,field,cse=True)
    jac_fn=m.compile(args,sp.Matrix(field).jacobian(args).tolist(),cse=True)
    H_fn=m.compile(args,H,cse=True)
    samples=[]
    for sheet in (1,-1):
        progress(f"{data['name']}: full sigma modes with high-degree corrections, sheet {sheet}")
        w0=sheet*m.anchor_w
        start=m.evaluate(m.rho_start,sheet,w=w0)
        initial=m.initial(start)
        for fraction in (c.mpf("0.5"),c.one):
            rho=m.rho_start+fraction*(m.rho_end-m.rho_start)
            dw=m.abel_increment(m.rho_start,rho,sheet)
            end=m.evaluate(rho,sheet,w=w0+dw)
            y=end["y"]
            # A finite base place fixes the same local lift on either sheet;
            # independently based infinity paths can differ by a period.
            dR=m.curve.abel_map_kind_2((rho,end["Y"]),
                                     base_place=(m.rho_start,start["Y"]),reduce=False).value
            tag=f"sheet {sheet}, fraction {fraction}"
            for k in range(g):
                reference=c.quad(lambda a:m.second_forms[k](a)/(sheet*c.sqrt(m.F(a))),[m.rho_start,rho])
                record(tag+f": moving second-kind primitive {k+1}",dR[k],reference)
            state=m.solution(start,end,initial)
            lu=[]; lv=[]
            potential_u=[]; potential_v=[]
            for j in range(n):
                exponent=P[j](rho)*y-P[j](m.rho_start)*start["y"]
                exponent+=sum(cs[j][k]*dw[k]+ds[j][k]*dR[k] for k in range(g))
                direct=c.quad(lambda a:B[j](a)/(sheet*c.sqrt(m.F(a))),[m.rho_start,rho])
                record(tag+f": correction {j+1} reduced primitive versus direct quadrature",exponent,direct)
                state[j]*=c.exp(exponent); state[n+j]*=c.exp(-exponent)
                b=B[j](rho)
                ru,rv=end["rates_u"][j],end["rates_v"][j]
                lu.append(ru+b); lv.append(rv-b)
                potential_u.append(end["potentials_u"][j]+2*ru*b+b*b+y*Bprime[j](rho))
                potential_v.append(end["potentials_v"][j]-2*rv*b+b*b-y*Bprime[j](rho))
                # Independent first-order quadrature; no reduced primitives.
                logu=c.quad(lambda a:(m.b[j]-sheet*c.sqrt(m.F(a)))/(2*(m.a[j]-a)*sheet*c.sqrt(m.F(a)))+B[j](a)/(sheet*c.sqrt(m.F(a))),[m.rho_start,rho])
                logv=c.quad(lambda a:(-m.b[j]-sheet*c.sqrt(m.F(a)))/(2*(m.a[j]-a)*sheet*c.sqrt(m.F(a)))-B[j](a)/(sheet*c.sqrt(m.F(a))),[m.rho_start,rho])
                record(tag+f": u_{j+1} sigma versus direct first-order quadrature",state[j],initial[j]*c.exp(logu))
                record(tag+f": v_{j+1} sigma versus direct first-order quadrature",state[n+j],initial[n+j]*c.exp(logv))
            values=field_fn(*state)
            second=c.matrix(jac_fn(*state))*c.matrix(values)
            for j in range(n):
                u,v=state[j],state[n+j]
                record(tag+f": mode {j+1} product",u*v,m.a[j]-rho)
                record(tag+f": u_{j+1} canonical derivative",u*lu[j],values[j])
                record(tag+f": v_{j+1} canonical derivative",v*lv[j],values[n+j])
                record(tag+f": u_{j+1} sigma second derivative versus full field Jacobian",u*potential_u[j],second[j])
                record(tag+f": v_{j+1} sigma second derivative versus full field Jacobian",v*potential_v[j],second[n+j])
            record(tag+": U",c.fprod(state[:n]),(m.q(rho)+y)/2)
            record(tag+": V",c.fprod(state[n:]),(m.q(rho)-y)/2)
            record(tag+": Hamiltonian",H_fn(*state),0)
            samples.append(dict(sheet=sheet,rho=c.nstr(rho,14),z=c.nstr(dw[0],14)))
    determinant=None
    if sp.degree(data["q"],x)>g+1:
        progress(f"{data['name']}: full longer determinant and residual sigma partition")
        determinant=verify_partition(m,reductions,record)
    return dict(name=data["name"],genus=g,N=n,degree_q=int(sp.degree(data["q"],x)),
                degrees_B=[None if b == 0 else int(sp.degree(b,x)) for b in data["B"]],F=str(data["F"]),q=str(data["q"]),
                Q=str(data["Q"]),B=[str(b) for b in data["B"]],
                reductions=[dict(P=str(r["P"]),c=[str(v) for v in r["c"]],d=[str(v) for v in r["d"]]) for r in reductions],
                exact_checks=exact,determinant=determinant,samples=samples,
                all_passed=all(exact.values()) and all(r["passed"] for r in rows),
                check_count=len(rows)+len(exact),checks=rows,tolerance=c.nstr(tolerance,8),
                maximum_scaled_error=c.nstr(max(c.mpf(r["scaled_error"]) for r in rows),12),seconds=round(perf_counter()-started,3))


def verify_partition(m,reductions,record):
    c,g,N=m.ctx,m.genus,m.n
    x,y=sp.symbols("x y")
    q=m.data["q"].subs(m.data["rho"],x)
    d=int(sp.degree(q,x))
    basis=columns(g,d,x,y); size=len(basis)
    A=list(range(size-1)); residual=list(range(size-1,N))
    assert len(residual)==g
    total=sum(m.eta,c.matrix(g,1))
    integers,error=integer_lattice_coordinates(m,total)
    for k,v in enumerate(error):record(f"partition: total Abel period {k+1}",v,0)
    mm,nn=c.matrix(integers[:g]),c.matrix(integers[g:])
    period_shift=2*m.periods.omega*mm+2*m.periods.omega_prime*nn
    second=m.curve.periods_kind_2()
    h_shift=2*second.eta*mm+2*second.eta_prime*nn
    eta=[v.copy() for v in m.eta]; eta[-1]-=period_shift
    h=[v[:] for v in m.h]; h[-1]=[h[-1][k]-h_shift[k] for k in range(g)]
    for k in range(g):
        record(f"partition: balanced Abel sum {k+1}",sum(v[k] for v in eta),0)
        record(f"partition: holomorphic correction cancellation {k+1}",
               sum(m.number(r["c"][k]) for r in reductions)-sum(v[k] for v in h),0)
    # Cache the z-independent fixed--fixed sigma factors in each identity.
    def sigma(w,index=None):return m.jet(w)[index or (0,)*g]
    def setup_fs(fixed):
        count=len(fixed)+1
        sign=(-1)**((count-g-1)*(count+g*g+2*g)//2)
        fixed_factor=c.mpf(sign)
        for i in range(len(fixed)):
            for j in range(i+1,len(fixed)):fixed_factor*=sigma(fixed[i]-fixed[j],m.flat)
        fixed_factor/=c.fprod(sigma(v,m.sharp)**count for v in fixed)
        summed=sum(fixed,c.matrix(g,1))
        return lambda w:fixed_factor*sigma(w+summed)*c.fprod(sigma(v-w,m.flat) for v in fixed)/sigma(w,m.sharp)**count
    main_fs=setup_fs([eta[j] for j in A])
    star_fs=setup_fs([-eta[j] for j in A])
    aux_fs=setup_fs([-eta[j] for j in residual])
    row_functions=[m.compile((x,y),f) for f in basis]
    delta=c.fprod(m.a[j]-m.a[i] for p,i in enumerate(A) for j in A[p+1:])
    constant=m.number(interpolation_constant(g,d,sp.Poly(q,x).LC()))
    star_constant=m.number(interpolation_constant(g,d,-sp.Poly(q,x).LC()))
    E=sum((eta[j] for j in A),c.matrix(g,1))
    def raws(w,sign):
        den=sigma(w,m.sharp)
        return [sigma(w+sign*v,m.flat)/(den*sigma(v,m.sharp)) for v in eta]
    constants=[]
    for sheet in (1,-1):
        w0=sheet*m.anchor_w
        p0=m.evaluate(m.rho_start,sheet,w=w0)
        ru0,rv0=raws(w0,1),raws(w0,-1)
        initial=m.initial(p0)
        for fraction in (c.zero,c.one):
            rho=m.rho_start+fraction*(m.rho_end-m.rho_start)
            w=w0+m.abel_increment(m.rho_start,rho,sheet)
            Y=sheet*c.sqrt(m.F(rho)); qq=m.q(rho)
            tag=f"partition sheet {sheet}, fraction {fraction}"
            D=c.det(c.matrix([[f(m.a[j],m.b[j]) for f in row_functions] for j in A]+[[f(rho,Y) for f in row_functions]]))
            Dstar=c.det(c.matrix([[f(m.a[j],-m.b[j]) for f in row_functions] for j in A]+[[f(rho,Y) for f in row_functions]]))
            record(tag+": longer algebraic determinant",D,constant*delta*(Y-qq)/2)
            record(tag+": opposite longer algebraic determinant",Dstar,star_constant*delta*(Y+qq)/2)
            record(tag+": full main sigma identity",main_fs(w),D)
            record(tag+": full opposite sigma identity",star_fs(w),Dstar)
            auxiliary=c.det(c.matrix([[m.a[j]**k for k in range(g+1)] for j in residual]+[[rho**k for k in range(g+1)]]))
            record(tag+": full residual auxiliary sigma identity",aux_fs(w),auxiliary)
            ru,rv=raws(w,1),raws(w,-1)
            composite=sigma(w+E)/sigma(w,m.sharp)
            constants.append((tag,composite/c.fprod(rv[j] for j in residual),D/c.fprod(rv)))
            record(tag+": complete raw u factorization",c.fprod(ru)/c.fprod(ru0),(qq+Y)/(m.q(m.rho_start)+p0["y"]))
            record(tag+": complete raw v factorization",c.fprod(rv)/c.fprod(rv0),(qq-Y)/(m.q(m.rho_start)-p0["y"]))
            record(tag+": full U equals normalized longer determinant",c.fprod(initial[:N])*c.fprod(ru)/c.fprod(ru0),Dstar/(star_constant*delta))
            record(tag+": full V equals normalized longer determinant",c.fprod(initial[N:])*c.fprod(rv)/c.fprod(rv0),-D/(constant*delta))
    for tag,composite,factor in constants[1:]:
        record(tag+": residual splitting constant",composite,constants[0][1])
        record(tag+": raw product splitting constant",factor,constants[0][2])
    return dict(size=size,degree_q=d,extra_ordinate_powers=list(range(1,d-g)),
                fixed_rows=len(A),residual_rows=len(residual),period_coordinates=integers,
                interpolation_constant=str(interpolation_constant(g,d,sp.Poly(q,x).LC())),
                opposite_interpolation_constant=str(interpolation_constant(g,d,-sp.Poly(q,x).LC())))


def run_verification(dps=40,progress=print):
    exact=exact_checks()
    reports=[verify_system(data,dps,progress) for data in systems()]
    return dict(exact_checks=exact,examples=reports,all_passed=all(exact.values()) and all(r["all_passed"] for r in reports),
                check_count=len(exact)+sum(r["check_count"] for r in reports))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps",type=int,default=40)
    parser.add_argument("--output",type=Path,default=HERE/"results.json")
    args=parser.parse_args()
    report=run_verification(args.dps,lambda s:print(s,flush=True))
    args.output.write_text(json.dumps(report,indent=2)+"\n")
    print(report["check_count"],"checks; passed=",report["all_passed"])
    for r in report["examples"]:
        print(r["name"],r["maximum_scaled_error"])
        for row in r["checks"]:
            if not row["passed"]: print(row)
    if not report["all_passed"]: raise SystemExit(1)
