"""72 gravity-only, equal-thickness elastic scenarios and explicit uplift example."""
from copy import deepcopy
import json, os, tempfile
from model import BASE,default,save
os.environ.setdefault('MPLCONFIGDIR',str(__import__('pathlib').Path(tempfile.gettempdir())/'capstone-mpl'))

MATERIALS=['shell_steel','shell_aluminum','shell_concrete','shell_sintered_regolith',
           'shell_geopolymer','shell_sulfur_concrete']

def main():
    results=[]
    for shape in ('dome','hyperboloid'):
        for material in MATERIALS:
            for support in ('pinned','linear_springs'):
                for planet in ('earth','moon','mars'):
                    c=default();c.update(shape=shape,material=material,support=support,planet=planet)
                    c['mesh']['n']=48
                    name=f'{shape}_{material}_{support}_{planet}'
                    r=save(c,BASE/'results'/name);results.append(r)
                    print('OK',name,flush=True)
    (BASE/'comparison.json').write_text(json.dumps(results,indent=2)+'\n')
    for planet in ('earth','moon','mars'):
        c=default();c.update(planet=planet,support='linear_springs',pressure_kPa=1.)
        c['mesh']['n']=48
        save(c,BASE/'results'/f'pressure_demo_{planet}')
    plot(results)

def plot(results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for row,shape in enumerate(('dome','hyperboloid')):
        for col,support in enumerate(('pinned','linear_springs')):
            ax=axes[row,col]
            for j,planet in enumerate(('earth','moon','mars')):
                values=[next(r['max_displacement_m']*1000 for r in results if
                        (r['shape'],r['material'],r['support'],r['planet'])==(shape,m,support,planet)) for m in MATERIALS]
                ax.bar(np.arange(6)+(j-1)*.25,values,width=.25,label=planet)
            ax.set_xticks(range(6),['steel','aluminum','concrete','sintered','geopolymer','sulfur'],rotation=25)
            ax.set(title=f'{shape} | {support}',ylabel='Max displacement (mm)')
            ax.grid(axis='y',alpha=.2);ax.legend()
    fig.suptitle('ELASTIC DEMONSTRATION | equal 20 mm thickness | gravity only | same support stiffness')
    fig.savefig(BASE/'comparison.png',dpi=160);plt.close(fig)
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig=plt.figure(figsize=(12,5),layout='constrained')
    from model import mesh
    for i,shape in enumerate(('dome','hyperboloid'),1):
        c=default();c['shape']=shape;c['mesh']['n']=12;xyz,q,b,l=mesh(c)
        ax=fig.add_subplot(1,2,i,projection='3d')
        ax.add_collection3d(Poly3DCollection(xyz[q],facecolors='#d4e8ec',edgecolors='#397a88',linewidths=.3))
        ax.auto_scale_xyz(xyz[:,0],xyz[:,1],xyz[:,2]);ax.set_box_aspect((10,10,max(xyz[:,2])))
        ax.set(xlabel='x (m)',ylabel='y (m)',zlabel='z (m)',title=f'{shape}: assumed prototype geometry')
    fig.savefig(BASE/'geometry.png',dpi=160);plt.close(fig)

if __name__=='__main__':main()
