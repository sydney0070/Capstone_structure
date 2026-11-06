"""Explicit lumped translational mass for massless material/element signatures.
Keep element self mass and added shielding mass as separate physical sources.
"""
import math

def add_translational_mass(ops,node_tags,total_mass_kg,*,source_key,ledger):
 if not node_tags or len(set(node_tags))!=len(node_tags):raise ValueError('distinct nonempty node list')
 if not math.isfinite(total_mass_kg) or total_mass_kg<=0:raise ValueError('positive finite mass')
 if not source_key or source_key in ledger:raise ValueError('mass source already added or missing key')
 sizes=[len(ops.nodeDisp(n)) for n in node_tags]
 if any(n not in [3,6] for n in sizes):raise ValueError('3D nodes3/6DOF required')
 for n,ndf in zip(node_tags,sizes):
  old=[ops.nodeMass(n,d) for d in range(1,ndf+1)]
  for d in range(3):old[d]+=total_mass_kg/len(node_tags)
  ops.mass(n,*old)
 ledger[source_key]=dict(mass_kg=total_mass_kg,nodes=list(node_tags),distribution='uniform per listed node; geometry weights must be supplied separately for nonuniform meshes')
 return ledger[source_key]

def attach_material_volume_mass(ops,node_tags,volume_m3,material_metadata,*,source_key,ledger):
 if material_metadata['mass_in_material']:raise ValueError('material already supplies density; external element mass would duplicate self mass')
 if not math.isfinite(volume_m3) or volume_m3<=0:raise ValueError('positive volume')
 return add_translational_mass(ops,node_tags,volume_m3*material_metadata['density_kg_m3'],source_key=source_key,ledger=ledger)
