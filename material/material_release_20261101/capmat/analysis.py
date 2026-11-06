"""Conservative displacement-step subdivision with committed-state checks.
No material changes and no convergence tolerance relaxation during recovery.
"""
def advance_displacement(ops,node,dof,target,*,min_increment=1e-10,max_depth=16,on_commit=None,statistics=None):
 stats=statistics if statistics is not None else {}
 stats.setdefault('accepted_steps',0);stats.setdefault('failed_attempts',0);stats.setdefault('max_subdivision_depth',0)
 def go(end,depth):
  begin=ops.nodeDisp(node,dof);before_time=ops.getTime();du=end-begin
  if abs(du)<1e-14:return True
  ops.integrator('DisplacementControl',node,dof,du)
  ok=ops.analyze(1)
  if ok==0:
   if abs(ops.nodeDisp(node,dof)-end)>1e-8:raise RuntimeError('accepted step missed displacement target')
   stats['accepted_steps']+=1;stats['max_subdivision_depth']=max(stats['max_subdivision_depth'],depth)
   if on_commit is not None:on_commit()
   return True
  stats['failed_attempts']+=1
  if abs(ops.nodeDisp(node,dof)-begin)>1e-9 or abs(ops.getTime()-before_time)>1e-7*max(abs(before_time),1):raise RuntimeError('failed analysis did not restore committed state')
  if depth>=max_depth or abs(du)/2<min_increment:return False
  return go(begin+du/2,depth+1) and go(end,depth+1)
 return go(target,0)
