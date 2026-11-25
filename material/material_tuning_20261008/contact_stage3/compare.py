"""Three gravities, two vertical restraint types, common 10 kN lateral load."""
from pathlib import Path
import json,shutil,os,tempfile
import numpy as np
from contact_model import BASE,default,save
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'capstone-mpl'))

def main():
    rows=[]
    for planet in ('earth','moon','mars'):
        for support in ('contact_only','contact_anchor'):
            c=default();c.update(planet=planet,support=support,horizontal_force_kN=10.)
            dest=BASE/'results'/f'{planet}_{support}'
            cache=BASE/'check48'/support
            if planet=='moon' and cache.exists() and json.loads((cache/'config.json').read_text())==c:
                shutil.copytree(cache,dest,dirs_exist_ok=True)
                r=json.loads((dest/'summary.json').read_text())
            else:r=save(c,dest)
            r['case']=dest.name;rows.append(r);print(dest.name,r['status'],flush=True)
    c=default();c.update(planet='moon',horizontal_force_kN=10.)
    c['anchor']['allowable_tension_per_m']=.1
    r=save(c,BASE/'results/weak_anchor_example');r['case']='weak_anchor_example';rows.append(r)
    (BASE/'comparison.json').write_text(json.dumps(rows,indent=2)+'\n')
    plot(rows)

def plot(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    colors=dict(earth='#1976a3',moon='#cd7117',mars='#35804a')
    for r in rows:
        if r['case']=='weak_anchor_example':continue
        h=json.loads((BASE/'results'/r['case']/'history.json').read_text())
        h=[x for x in h if x['phase']=='pressure' or (x['phase']=='lateral' and x['horizontal_load_factor']==1.)]
        color=colors[r['planet']];style='-' if r['support']=='contact_anchor' else '--'
        for ax,key,scale in zip(axes,('contact_length_fraction','max_gap_m','anchor_tension_kN'),(1.,1000.,1.)):
            ax.plot([x['pressure_kPa'] for x in h],[x[key]*scale for x in h],style,color=color,
                    label=r['case'].replace('contact_',''),lw=1.4)
    for ax,label in zip(axes,('Contact length / perimeter','Maximum vertical opening (mm)','Total anchor tension (kN)')):
        ax.set(xlabel='Pressure difference (kPa)',ylabel=label);ax.grid(alpha=.2)
    axes[0].legend(fontsize=7)
    fig.suptitle('DEMONSTRATION | elastic dome, unilateral vertical contact | independent lateral guides, H=10 kN')
    fig.savefig(BASE/'contact_comparison.png',dpi=160);plt.close(fig)
    a=np.load(BASE/'results/moon_contact_only/last_converged.npz')
    i=a['boundary'];theta=np.degrees(np.arctan2(a['xyz'][i,1],a['xyz'][i,0]))
    fig,axes=plt.subplots(2,1,figsize=(8,6),layout='constrained')
    axes[0].plot(theta,a['u'][i,2]*1000);axes[0].axhline(0,color='k',lw=.6)
    axes[0].set(ylabel='Vertical displacement (mm)',title='Moon, no anchors: last converged pressure step, NOT a collapse state')
    axes[1].plot(theta,a['contact']/a['boundary_lengths'])
    axes[1].set(xlabel='Base angle (degree)',ylabel='Compression per length (kN/m)')
    for ax in axes:ax.grid(alpha=.2)
    fig.savefig(BASE/'partial_contact.png',dpi=160);plt.close(fig)

if __name__=='__main__':main()
