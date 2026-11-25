"""Reusable material factories; functions never wipe an existing domain."""
import openseespy.opensees as ops

def uniaxial(tag, c):
    p, m = c['params'], c['material']
    keys = {
        'Concrete02': ('fpc','epsc0','fpcu','epsU','lam','ft','Ets'),
        'Steel02': ('Fy','E','b','R0','cR1','cR2','a1','a2','a3','a4'),
        'PySimple1': ('soil_type','pult','y50','Cd','c'),
        'TzSimple1': ('soil_type','tult','z50','c'),
        'QzSimple1': ('qz_type','qult','z50','suction','c'),
        'ElasticPPGap': ('K','Fy','gap','eta'),
        'Viscous': ('C','alpha'),
    }
    ops.uniaxialMaterial(m, tag, *(p[k] for k in keys[m]))

def continuum(tag, c, nd=3):
    p, m = c['params'], c['material']
    if m == 'ElasticIsotropic':
        ops.nDMaterial(m,tag,p['E'],p['nu'],p['rho'])
    elif m == 'ElasticOrthotropic':
        keys=('Ex','Ey','Ez','nu_xy','nu_yz','nu_zx','Gxy','Gyz','Gzx','rho')
        ops.nDMaterial(m,tag,*(p[k] for k in keys))
    elif m == 'J2Plasticity':
        K=p['E']/(3*(1-2*p['nu'])); G=p['E']/(2*(1+p['nu']))
        ops.nDMaterial(m,tag,K,G,p['sig0'],p['sigInf'],p['delta'],p['H'])
    elif m == 'DruckerPrager':
        from math import sin,cos,radians,sqrt
        phi,psi=radians(p['phi']),radians(p['psi'])
        # Compression-meridian Mohr-Coulomb fit from OpenSees documentation.
        friction=2*sqrt(2)*sin(phi)/(sqrt(3)*(3-sin(phi)))
        dilatancy=2*sqrt(2)*sin(psi)/(sqrt(3)*(3-sin(psi)))
        sy=6*p['cohesion']*cos(phi)/(sqrt(2)*(3-sin(phi)))
        ops.nDMaterial(m,tag,p['K'],p['G'],sy,friction,dilatancy,
                       0.,0.,0.,0.,p['H'],p['theta'],p['rho'],p['pa'])
    elif m == 'PressureIndependMultiYield':
        ops.nDMaterial(m,tag,nd,p['rho'],p['G'],p['K'],p['cohesion'],
                       p['gamma_peak'],0.,100.,0.,20)
    elif m == 'PressureDependMultiYield02':
        keys = ('rho','G','K','phi','gamma_peak','p_ref','d','phi_pt',
                'contraction1','contraction3','dilation1','dilation3',
                'surfaces','contraction2','dilation2','liquefaction1',
                'liquefaction2','e','cs1','cs2','cs3','pa')
        ops.nDMaterial(m,tag,nd,*(p[k] for k in keys))
    elif m == 'PM4Sand':
        if nd != 2:
            raise ValueError('This OpenSees PM4Sand interface is PlaneStrain only')
        ops.nDMaterial(m,tag,p['Dr'],p['G0'],p['hpo'],p['rho'],p['pa'],
                       -1.,p['emax'],p['emin'],.5,.1,-1.,-1.,250.,-1.,33.,p['nu'])
    elif m == 'PM4Silt':
        if nd != 2:
            raise ValueError('This OpenSees PM4Silt interface is PlaneStrain only')
        ops.nDMaterial(m,tag,p['Su'],p['Su_ratio'],p['G0'],p['hpo'],p['rho'],
                       p['Su_factor'],p['pa'],p['nu'],.75,.5,p['e'])
    else:
        raise ValueError(m)

def pile_capacities(pult_per_length, tributary_length, shaft_stress,
                    diameter, tip_stress):
    """kN/m * m, kPa * m² => nodal spring capacities in kN."""
    from math import pi
    if min(pult_per_length,tributary_length,shaft_stress,diameter,tip_stress) <= 0:
        raise ValueError('all capacities and geometry must be positive')
    return dict(pult=pult_per_length*tributary_length,
                tult=shaft_stress*pi*diameter*tributary_length,
                qult=tip_stress*pi*diameter**2/4)

def shell_section(sec_tag,c,mat_tag=1,plate_tag=2):
    """For 3D/6DOF ShellMITC4 or a separately verified nonlinear shell.
    PlateFiber imposes sigma_zz=0. Orthotropy is along shell local material axes.
    J2Plasticity has no density argument; assign nodal mass in the full model.
    """
    p=c['params']
    continuum(mat_tag,c,3)
    ops.nDMaterial('PlateFiber',plate_tag,mat_tag)
    ops.section('PlateFiber',sec_tag,plate_tag,p['thickness'])

def rc_section(tag, cover_tag, core_tag, steel_tag, width=.4, depth=.4,
               cover=.04, bar_area=.000314, bars_per_face=3):
    """2D fiber section. cover denotes core boundary and bar-center offset here.
    Replace geometry with drawings. Equal-area fibers subtract bar area from
    concrete at each bar position to prevent area double counting.
    """
    if min(width,depth) <= 2*cover or bars_per_face < 2 or bar_area <= 0:
        raise ValueError('invalid RC geometry')
    ops.section('Fiber',tag)
    y, z = depth/2, width/2
    yc, zc = y-cover, z-cover
    ops.patch('rect',core_tag,16,16,-yc,-zc,yc,zc)
    ops.patch('rect',cover_tag,2,16,-y,-z,-yc,z)
    ops.patch('rect',cover_tag,2,16,yc,-z,y,z)
    ops.patch('rect',cover_tag,16,2,-yc,-z,yc,-zc)
    ops.patch('rect',cover_tag,16,2,-yc,zc,yc,z)
    # Bars are deliberately inside the core, rather than exactly on its boundary.
    import numpy as np
    for yy in (-yc*.95,yc*.95):
        for zz in np.linspace(-zc*.95,zc*.95,bars_per_face):
            ops.fiber(float(yy),float(zz),-bar_area,core_tag)
            ops.fiber(float(yy),float(zz),bar_area,steel_tag)

def rc_beam(element_tag, i, j, section_tag, transf_tag, integration_tag):
    """Call after creating a 2D ndm=2/ndf=3 domain and geometric transformation."""
    ops.beamIntegration('Lobatto',integration_tag,section_tag,5)
    ops.element('forceBeamColumn',element_tag,i,j,transf_tag,integration_tag)

def pile_spring(element_tag, soil_node, pile_node, material_tag, direction):
    """Nodes must be colocated; direction uses zeroLength LOCAL coordinates."""
    ops.element('zeroLength',element_tag,soil_node,pile_node,
                '-mat',material_tag,'-dir',direction)
