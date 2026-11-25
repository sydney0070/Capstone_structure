"""Check the actual calibrated fiber section: exact integrals and cyclic axial response."""
from pathlib import Path
import json
import numpy as np
import openseespy.opensees as ops
import dome_in_mars_tuned as m


def main():
    section,_=m.make_section_for_target_mass(200.)
    fibers=m.shs_fibers(section)
    for target,actual in [(section['A'],sum(a for y,z,a in fibers)),(section['Iz'],sum(a*y*y for y,z,a in fibers)),(section['Iy'],sum(a*z*z for y,z,a in fibers))]:
        assert np.isclose(target,actual,rtol=1e-12)
    ops.wipe();ops.model('basic','-ndm',3,'-ndf',6)
    m.register_tuned_shs(section)
    ops.node(1,0.,0.,0.);ops.node(2,1.,0.,0.)
    ops.fix(1,1,1,1,1,1,1);ops.fix(2,0,1,1,1,1,1)
    ops.geomTransf('Linear',1,0.,0.,1.)
    ops.element('forceBeamColumn',1,1,2,1,m.STEEL_INTEGRATION_TAG,'-iter',40,1e-10)
    ops.timeSeries('Linear',1);ops.pattern('Plain',1,1);ops.load(2,1.,0.,0.,0.,0.,0.)
    ops.constraints('Transformation');ops.numberer('RCM');ops.system('UmfPack');ops.test('NormUnbalance',1e-4,50);ops.algorithm('Newton')
    ops.integrator('DisplacementControl',2,1,1e-5);ops.analysis('Static')
    assert ops.analyze(1)==0
    initial=ops.getLoadFactor(1)/(section['A']*1e-5)
    assert abs(initial/section['E']-1)<1e-6
    data=[[ops.nodeDisp(2,1),ops.getLoadFactor(1)/section['A']]]
    start=data[-1][0]
    for target in [.004,-.004,0.]:
        for u in np.linspace(start,target,51)[1:]:
            ops.integrator('DisplacementControl',2,1,float(u-ops.nodeDisp(2,1)))
            assert ops.analyze(1)==0
            data.append([ops.nodeDisp(2,1),ops.getLoadFactor(1)/section['A']])
        start=target
    a=np.array(data);fy=m.TUNED_STEEL['params']['Fy']
    assert np.max(np.abs(a[:,1]))>fy*.95
    assert abs(a[-1,1])>fy*.1
    out=Path(__file__).resolve().parent/'results_mars_dome_tuned'/'verification';out.mkdir(parents=True,exist_ok=True)
    np.savetxt(out/'section_cyclic.csv',a,delimiter=',',header='axial_strain,stress_Pa',comments='')
    result=dict(pass_=True,initial_E_Pa=initial,expected_E_Pa=section['E'],max_stress_Pa=float(np.max(abs(a[:,1]))),stress_at_final_zero_strain_Pa=float(a[-1,1]),scope='actual fiber forceBeamColumn, ideal axial member; not SHS buckling or dome validation')
    (out/'section_verification.json').write_text(json.dumps(result,indent=2));ops.wipe();print(result)

if __name__=='__main__':main()
