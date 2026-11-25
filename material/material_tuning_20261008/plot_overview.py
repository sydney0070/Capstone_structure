"""Plot saved results only. Curves are demonstration seeds, not experiments."""
import os
import tempfile
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'capstone-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE=Path(__file__).resolve().parent

if __name__=='__main__':
    names=['shell_steel','shell_aluminum','shell_cfrp_elastic',
           'shell_concrete','shell_sintered_regolith','shell_geopolymer',
           'shell_sulfur_concrete','moon_regolith_dp','mars_regolith_dp']
    fig,axes=plt.subplots(3,3,figsize=(13,10),layout='constrained')
    fig.suptitle('Candidate material coupons | DEMONSTRATION, NOT CALIBRATED',fontsize=15)
    for i,(name,ax) in enumerate(zip(names,axes.ravel())):
        a=np.genfromtxt(BASE/'results'/name/'history.csv',delimiter=',',names=True)
        scale=1000. if name.startswith('shell') else 1.
        ax.plot(100*a['x'],a['y']/scale,lw=1.4,color='#0b7285' if scale==1000 else '#bd6b35')
        ax.set(title=name.replace('shell_','').replace('_',' '),xlabel='strain (%)',
               ylabel='axial stress (MPa)' if scale==1000 else 'shear stress (kPa)')
        ax.grid(alpha=.2)
    fig.savefig(BASE/'overview.png',dpi=170)
    plt.close(fig)
