"""SSI contact and radiation components. No default foundation stiffness."""
import math
from .mechanics import contact_parameters,radiation_dashpots
from .materials import _tag

def register_contact(ops,element_tag,foundation_node,ground_node,parameters,*,normal_axis=3):
 """Coincident3DOF contact nodes; ground-to-foundation normal points positive axis.
 A6DOF beam/shell node needs a3DOF interface node coupled in translations.
 Kn/Kt are nodalN/m, cohesive_force isN (stress times tributary area).
 """
 for t in [element_tag,foundation_node,ground_node]:_tag(t)
 if foundation_node==ground_node:raise ValueError('distinct contact nodes required')
 if normal_axis not in [1,2,3]:raise ValueError('axis1,2,3 required; rotate geometry externally for other normals')
 for n in [foundation_node,ground_node]:
  if len(ops.nodeDisp(n))!=3:raise ValueError('release adapter requires3DOF nodes; use translational interface nodes')
 for key in ['Kn_N_m','Kt_N_m']:
  if not math.isfinite(parameters[key]) or parameters[key]<=0:raise ValueError('positive nodal penalty required')
 for key in ['mu','cohesive_force_N']:
  if not math.isfinite(parameters[key]) or parameters[key]<0:raise ValueError('nonnegative friction/cohesion')
 a=ops.nodeCoord(foundation_node);b=ops.nodeCoord(ground_node)
 if max(abs(x-y) for x,y in zip(a,b))>1e-9:raise ValueError('this adapter accepts initially coincident nodes only')
 if parameters['cohesive_force_N']!=0:raise ValueError('Selected ASDimplex interface is unbonded friction only; adhesion cannot be silently ignored')
 normal=[0.,0.,0.];normal[normal_axis-1]=1.
 ops.element('zeroLengthContactASDimplex',element_tag,ground_node,foundation_node,parameters['Kn_N_m'],parameters['Kt_N_m'],parameters['mu'],'-orient',*normal,'-intType',0)
 return dict(element_tag=element_tag,normal_axis=normal_axis,foundation_node=foundation_node,ground_node=ground_node,**parameters,contact_scope='ASDimplex implicit0 opening/closure and unbonded Coulomb slip; general solver; nonzero adhesion rejected')

def register_dashpot(ops,element_tag,material_tag,node,anchor_node,coefficient_N_s_m,dof):
 for tag in [element_tag,material_tag,node,anchor_node]:_tag(tag)
 if not math.isfinite(coefficient_N_s_m) or coefficient_N_s_m<=0:raise ValueError('positive viscous coefficient')
 if dof not in [1,2,3]:raise ValueError('translation DOF only')
 if max(abs(x-y) for x,y in zip(ops.nodeCoord(node),ops.nodeCoord(anchor_node)))>1e-9:raise ValueError('coincident zeroLength nodes required')
 ops.uniaxialMaterial('Viscous',material_tag,coefficient_N_s_m,1.)
 ops.element('zeroLength',element_tag,node,anchor_node,'-mat',material_tag,'-dir',dof)
 return dict(element_tag=element_tag,C_N_s_m=coefficient_N_s_m,dof=dof,scope='relative-velocity viscous component; anchor boundary/input-motion definition belongs to assembled model')
