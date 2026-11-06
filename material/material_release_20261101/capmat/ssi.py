"""SSI contact and radiation components. No default foundation stiffness."""
import math
from .mechanics import contact_parameters,radiation_dashpots
from .materials import _tag

def register_contact(ops,element_tag,foundation_node,ground_node,parameters,*,auxiliary_nodes,normal_material_tag,rotation_material_tag,friction_tag,normal_axis=3,restore_ndf=3):
 """Translational Coulomb/uplift interface via flatSliderBearing3D.
 The verified friction kernel acts between6DOF auxiliary nodes, with their
 rotations fixed and translations tied to physical nodes. These auxiliary
 rotational constraints do NOT restrain the physical foundation rotation.
 Auxiliary bearing P-Delta moment reactions are formulation artifacts and
 are not transmitted to the physical nodes. Use small relative slip, fixed
 normal, no rolling moment; this is a point-traction interface, not a bearing
 capacity model or finite-sliding surface contact.
 On tested OpenSees3.4, force-controlled assembly MUST use Plain constraints;
 Transformation with these auxiliary MP/SP constraints gave accumulated
 physical displacement errors. Do not transfer the prescribed-SP patch
 handler into a force-controlled model.
 """
 for t in [element_tag,foundation_node,ground_node,normal_material_tag,rotation_material_tag,friction_tag,*auxiliary_nodes]:_tag(t)
 if len(auxiliary_nodes)!=2 or len(set(auxiliary_nodes))!=2:raise ValueError('two distinct unused auxiliary node tags')
 if any(n in ops.getNodeTags() for n in auxiliary_nodes):raise ValueError('auxiliary node tags already used')
 if normal_material_tag==rotation_material_tag:raise ValueError('distinct material tags')
 if foundation_node==ground_node:raise ValueError('distinct contact nodes required')
 if normal_axis not in [1,2,3] or restore_ndf not in [3,6]:raise ValueError('normal axis1/2/3 and restore_ndf3/6 required')
 for n in [foundation_node,ground_node]:
  if len(ops.nodeDisp(n)) not in [3,6]:raise ValueError('physical3DOF/6DOF nodes required')
 for key in ['Kn_N_m','Kt_N_m']:
  if not math.isfinite(parameters[key]) or parameters[key]<=0:raise ValueError('positive nodal penalty required')
 for key in ['mu','cohesive_force_N']:
  if not math.isfinite(parameters[key]) or parameters[key]<0:raise ValueError('nonnegative friction/cohesion')
 if parameters['cohesive_force_N']!=0:raise ValueError('Selected interface is unbonded friction only; adhesion cannot be silently ignored')
 a=ops.nodeCoord(foundation_node);b=ops.nodeCoord(ground_node)
 if max(abs(x-y) for x,y in zip(a,b))>1e-9:raise ValueError('initially coincident physical nodes required')
 ground_aux,foundation_aux=auxiliary_nodes
 ops.model('basic','-ndm',3,'-ndf',6)
 for physical,aux in [(ground_node,ground_aux),(foundation_node,foundation_aux)]:
  ops.node(aux,*ops.nodeCoord(physical));ops.fix(aux,0,0,0,1,1,1);ops.equalDOF(physical,aux,1,2,3)
 ops.uniaxialMaterial('Elastic',normal_material_tag,parameters['Kn_N_m']);ops.uniaxialMaterial('Elastic',rotation_material_tag,0.)
 ops.frictionModel('Coulomb',friction_tag,parameters['mu'])
 normal=[0.,0.,0.];normal[normal_axis-1]=1.;tangent=[0.,0.,0.];tangent[1 if normal_axis==1 else 0]=1.
 ops.element('flatSliderBearing',element_tag,ground_aux,foundation_aux,friction_tag,parameters['Kt_N_m'],'-P',normal_material_tag,'-T',rotation_material_tag,'-My',rotation_material_tag,'-Mz',rotation_material_tag,'-orient',*normal,*tangent,'-iter',100,1e-10)
 ops.model('basic','-ndm',3,'-ndf',restore_ndf)
 return dict(element_tag=element_tag,normal_axis=normal_axis,foundation_node=foundation_node,ground_node=ground_node,auxiliary_nodes=list(auxiliary_nodes),**parameters,contact_scope='translational point friction/uplift; no tensile bond or cohesion; fixed normal and small relative slip; no added physical rotational stiffness; force-controlled assembly requires Plain handler with identity equalDOF and homogeneous SP; use DisplacementControl for settlement; nonzero prescribed SP needs separate verified integration')

def register_dashpot(ops,element_tag,material_tag,node,anchor_node,coefficient_N_s_m,dof):
 for tag in [element_tag,material_tag,node,anchor_node]:_tag(tag)
 if not math.isfinite(coefficient_N_s_m) or coefficient_N_s_m<=0:raise ValueError('positive viscous coefficient')
 if dof not in [1,2,3]:raise ValueError('translation DOF only')
 if max(abs(x-y) for x,y in zip(ops.nodeCoord(node),ops.nodeCoord(anchor_node)))>1e-9:raise ValueError('coincident zeroLength nodes required')
 ops.uniaxialMaterial('Viscous',material_tag,coefficient_N_s_m,1.)
 ops.element('zeroLength',element_tag,node,anchor_node,'-mat',material_tag,'-dir',dof)
 return dict(element_tag=element_tag,C_N_s_m=coefficient_N_s_m,dof=dof,scope='relative-velocity viscous component; anchor boundary/input-motion definition belongs to assembled model')
