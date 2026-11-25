"""Independent mechanical invariants and mesh checks of elastic screening."""
from copy import deepcopy
import json
import numpy as np
from model import default,run,mesh,surface_loads,BASE

def verify():
    checks=[]
    def check(name,error,tolerance):
        assert np.isfinite(error)
        checks.append(dict(check=name,error=float(error),tolerance=tolerance,passed=bool(error<tolerance)))
    c=default();c['mesh']['n']=12;xyz,q,b,l=mesh(c);a,p=surface_loads(xyz,q)
    # Exact area of spherical cap is pi*(base radius²+rise²).
    exact=np.pi*(c['geometry']['radius']**2+c['geometry']['rise']**2)
    check('spherical_cap_area',abs(sum(a)/exact-1),.01)
    polygon=abs(sum(np.cross(xyz[i],xyz[j])[2] for i,j in zip(b,b[1:]+b[:1]))/2)
    check('pressure_projected_area',abs(sum(p[:,2])/polygon-1),1e-12)
    for shape in ('dome','hyperboloid'):
        for support in ('pinned','linear_springs'):
            c=default();c.update(shape=shape,support=support);c['mesh']['n']=12
            r,d=run(c)
            check(f'{shape}_{support}_force_balance',r['relative_force_residual'],1e-9)
            check(f'{shape}_{support}_moment_balance',np.linalg.norm(r['moment_residual_kNm']),1e-6)
            moon=deepcopy(c);moon['planet']='moon';rm,dm=run(moon)
            check(f'{shape}_{support}_gravity_scaling',np.max(abs(dm['displacements']-d['displacements']*1.62/9.81))/r['max_displacement_m'],1e-9)
            check(f'{shape}_{support}_mass_invariant',abs(rm['mass_tonne']-r['mass_tonne']),1e-12)
            if support=='linear_springs':
                _,_,base,lengths=mesh(c)
                weighted=sum(lengths[i]*d['displacements'][i,2] for i in base)
                check(f'{shape}_spring_settlement_equilibrium',abs(weighted*c['foundation']['k_vertical_per_m']/r['weight_kN']+1),1e-9)
                stiff=deepcopy(c);stiff['foundation']={k:1e12 for k in c['foundation']}
                rs,ds=run(stiff);pin=deepcopy(c);pin['support']='pinned';rp,dp=run(pin)
                check(f'{shape}_rigid_support_limit',np.max(abs(ds['displacements']-dp['displacements']))/rp['max_displacement_m'],1e-5)
    c=default();c['mesh']['n']=12;r,d=run(c);c['material_override']={'E':2*r['elastic_params_used']['E']}
    rr,dd=run(c)
    check('elastic_E_inverse_scaling',np.max(abs(dd['displacements']*2-d['displacements']))/r['max_displacement_m'],1e-9)
    c=default();c['mesh']['n']=12;c.update(pressure_kPa=1.,gravity_factor=0.,support='linear_springs')
    r,d=run(c);c['planet']='moon';rr,dd=run(c)
    check('pressure_independent_of_gravity',np.max(abs(d['displacements']-dd['displacements'])),1e-12)
    check('pressure_hold_down_balance',abs(r['required_hold_down_kN']/polygon-1),1e-9)
    convergence=[]
    for shape in ('dome','hyperboloid'):
        for support in ('pinned','linear_springs'):
            rows=[]
            for n in (12,24,48,72):
                c=default();c.update(shape=shape,support=support);c['mesh']['n']=n
                r,_=run(c);rows.append(dict(n=n,max_displacement_m=r['max_displacement_m'],area=r['surface_area_m2']))
            error=abs(rows[-2]['max_displacement_m']/rows[-1]['max_displacement_m']-1)
            check(f'{shape}_{support}_mesh_48_to_72',error,.02)
            convergence.append(dict(shape=shape,support=support,rows=rows))
    report=dict(checks=checks,mesh=convergence,all_passed=all(r['passed'] for r in checks),
                scope='Elastic aluminium demo. Mesh tolerance 2% is a project choice, not a design standard.')
    (BASE/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    if not report['all_passed']:raise SystemExit(2)

if __name__=='__main__':verify()
