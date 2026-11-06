"""Compact scientific verification figure from saved numeric evidence."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent/'results'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2})
f,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
a=np.genfromtxt(R/'steel_coupon/validation_LP2_id3.csv',delimiter=',',names=True)
ax=axes[0,0];ax.plot(a['strain']*100,a['target_Pa']/1e6,color='#555555',label='Measured coupon');ax.plot(a['strain']*100,a['model_Pa']/1e6,'--',color='#007a87',label='Steel02 replay');ax.set(xlabel='Engineering strain (%)',ylabel='Stress (MPa)',title='S355: independent specimen LP2/id3');ax.legend()
ax=axes[0,1]
for n,style in [(60,'-'),(120,'--')]:
 a=np.genfromtxt(R/f'rc_hinge/lp0.3_axial0/force_ny36_ip6_steps{n}_ne1.csv',delimiter=',',names=True,dtype=None,encoding='utf-8');p=a[a['phase']=='path'];ax.plot(p['u_m']*1000,p['force_N']/1000,style,label=f'{n} steps/segment')
ax.set(xlabel='Tip displacement (mm)',ylabel='Lateral force (kN)',title='RC: fixed 0.30 m hinge, no axial preload');ax.legend()
ax=axes[1,0]
for p in [15000,30000,60000]:
 a=np.genfromtxt(R/f'soil_material/moon_jsc1a_18_peak_{p}.csv',delimiter=',',names=True);ax.plot(a['axial_compression_increment']*100,-a['sz']/1000,label=f'Confinement {p/1000:g} kPa')
ax.set(xlabel='Axial compression increment (%)',ylabel='Axial stress (kPa)',title='JSC-1A reference: DP law verification');ax.legend()
ax=axes[1,1];a=np.genfromtxt(R/'ssi_contact/area1_adhesion0_axis3.csv',delimiter=',',names=True);ax.plot(a['ux']*1000,a['tangent_N']/1000,color='#b75d20');ax.set(xlabel='Relative slip (mm)',ylabel='Tangential force (kN)',title='SSI: slip, opening and reclosure (1 m²)')
f.suptitle('Numerical verification — source assumptions and scope remain explicit',fontsize=14)
f.savefig(R/'verification_overview.png',dpi=160);plt.close(f)
