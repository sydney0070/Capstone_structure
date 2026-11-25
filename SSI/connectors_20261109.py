"""Foundation–dome fixed, hinged and rotational-spring connectors."""
import math
import numpy as np
import openseespy.opensees as ops

def register_tension_anchor(ground_auxiliary_node,foundation_node,*,tag,stiffness_N_m,capacity_N):
    """Positive-extension tie using an existing ground-side six-DOF port.

    This preserves the prototype ElasticPPGap law; it is not an anchor pullout
    model and needs explicitly supplied stiffness and capacity.
    """
    ops.uniaxialMaterial('ElasticPPGap',tag,stiffness_N_m,capacity_N,0.,0.)
    ops.element('zeroLength',tag,ground_auxiliary_node,foundation_node,'-mat',tag,'-dir',3)

def connect_dome_port(foundation_node,dome_node,*,rotation='fixed',rotational_stiffness_Nm_rad=None,tag=900000):
    """Coincident six-DOF nodes; no rigidLink or chained MP constraints.

    fixed: identity six-DOF connection to the deformable foundation, NOT a fixed
    global support. hinged: translations only; dome rotations remain independent.
    spring: translations tied, three rotational elastic springs, supplied K only.
    """
    if rotation not in ('fixed','hinged','spring'):raise ValueError('rotation mode')
    if len(ops.nodeDisp(foundation_node))!=6 or len(ops.nodeDisp(dome_node))!=6:raise ValueError('6DOF ports required')
    if np.linalg.norm(np.array(ops.nodeCoord(foundation_node))-ops.nodeCoord(dome_node))>1e-8:raise ValueError('coincident ports')
    if foundation_node==dome_node:raise ValueError('distinct nodes')
    if rotation=='spring' and (rotational_stiffness_Nm_rad is None or not math.isfinite(rotational_stiffness_Nm_rad) or rotational_stiffness_Nm_rad<=0):raise ValueError('explicit positive rotational stiffness required')
    ops.equalDOF(foundation_node,dome_node,*(range(1,7) if rotation=='fixed' else (1,2,3)))
    if rotation=='spring':
        ops.uniaxialMaterial('Elastic',tag,rotational_stiffness_Nm_rad)
        ops.element('zeroLength',tag,foundation_node,dome_node,'-mat',tag,tag,tag,'-dir',4,5,6)
    return dict(foundation_node=foundation_node,dome_node=dome_node,rotation=rotation,rotational_stiffness_Nm_rad=rotational_stiffness_Nm_rad)
