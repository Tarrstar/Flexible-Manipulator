"""Plot terminal-link tip positions from four independent dynamics runs."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    root=Path(__file__).resolve().parents[3]
    parser.add_argument('--directory',type=Path,default=root/'data/space_7dof_all_flexible')
    args=parser.parse_args(); directory=args.directory
    keys=['matlab_flexible','exudyn_flexible','matlab_rigid','exudyn_rigid']
    raw={k:np.loadtxt(directory/(k+'_tip.csv'),delimiter=',',skiprows=1) for k in keys}
    for k,d in raw.items():
        assert d.shape[1]==4 and np.isfinite(d).all(),k
        _,inds=np.unique(d[:,0],return_index=True)
        raw[k]=d[np.sort(inds)]
    t=raw['matlab_flexible'][:,0]
    p={}
    for k,d in raw.items():
        assert d[0,0]<=t[0]+1e-10 and d[-1,0]>=t[-1]-1e-10,k
        p[k]=np.column_stack([np.interp(t,d[:,0],d[:,i]) for i in range(1,4)])
    def statistics(delta):
        norm=np.linalg.norm(delta,axis=1)
        return dict(max_position_error_m=float(norm.max()),
                    rms_position_error_m=float(np.sqrt(np.mean(norm**2))),
                    final_position_error_m=float(norm[-1]),
                    max_component_error_m=np.abs(delta).max(axis=0).tolist())
    results=dict(flexible_matlab_vs_exudyn=statistics(p[keys[0]]-p[keys[1]]),
                 rigid_matlab_vs_exudyn=statistics(p[keys[2]]-p[keys[3]]),
                 matlab_flexible_vs_rigid=statistics(p[keys[0]]-p[keys[2]]),
                 exudyn_flexible_vs_rigid=statistics(p[keys[1]]-p[keys[3]]))
    styles=[dict(color='#1565C0',ls='-',lw=2.0),dict(color='#E66B24',ls='--',lw=1.7),
            dict(color='#2F7D50',ls='-.',lw=1.6),dict(color='#353A40',ls=':',lw=1.8)]
    labels=['MATLAB: all flexible','Exudyn: all flexible','MATLAB: all rigid','Exudyn: all rigid']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.titleweight':'semibold',
                         'savefig.facecolor':'white'})
    fig,axes=plt.subplots(4,1,figsize=(10,11.2),sharex=True,
                          gridspec_kw={'height_ratios':[1,1,1,.9]})
    for j,ax in enumerate(axes[:3]):
        for k,label,style in zip(keys,labels,styles): ax.plot(t,p[k][:,j],label=label,**style)
        ax.set_ylabel(f'{"xyz"[j]} (m)'); ax.grid(alpha=.22); ax.set_xlim(t[0],t[-1])
    axes[0].legend(ncol=2,loc='best',frameon=False,fontsize=10)
    axes[0].set_title('Spatial 7-joint arm | terminal tip position',loc='left',pad=14)
    axes[3].plot(t,1e3*np.linalg.norm(p[keys[0]]-p[keys[1]],axis=1),
                 color='#1565C0',label='MATLAB vs Exudyn (flexible)',lw=1.8)
    axes[3].plot(t,1e3*np.linalg.norm(p[keys[2]]-p[keys[3]],axis=1),
                 color='#353A40',ls=':',label='MATLAB vs Exudyn (rigid)',lw=1.6)
    axes[3].set_ylabel('Position difference (mm)'); axes[3].set_xlabel('Time (s)')
    axes[3].grid(alpha=.22); axes[3].legend(frameon=False,fontsize=10)
    fig.tight_layout(pad=1.5)
    fig.savefig(directory/'endpoint_comparison.png',dpi=220)
    fig.savefig(directory/'endpoint_comparison.svg')
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(10,4.3))
    for j,color in enumerate(['#1565C0','#E66B24','#2F7D50']):
        ax.plot(t,1e3*(p[keys[0]][:,j]-p[keys[2]][:,j]),color=color,
                label=f'MATLAB: {"xyz"[j]}',lw=1.7)
        ax.plot(t,1e3*(p[keys[1]][:,j]-p[keys[3]][:,j]),color=color,ls='--',
                label=f'Exudyn: {"xyz"[j]}',lw=1.4)
    ax.set(xlabel='Time (s)',ylabel='Flexible minus rigid tip position (mm)',
           title='Effect of flexibility in all seven links',xlim=(t[0],t[-1]))
    ax.grid(alpha=.22); ax.legend(ncol=3,frameon=False,fontsize=9)
    fig.tight_layout(); fig.savefig(directory/'flexibility_effect.png',dpi=220); plt.close(fig)
    (directory/'comparison_metrics.json').write_text(json.dumps(results,indent=2))
    np.savetxt(directory/'endpoint_comparison.csv',np.column_stack([t,*[p[k] for k in keys]]),
               delimiter=',',comments='',header='time,'+','.join(k+'_'+c for k in keys for c in 'xyz'))
    print(json.dumps(results,indent=2))
    print('Figures saved in '+str(directory))


if __name__=='__main__':
    main()
