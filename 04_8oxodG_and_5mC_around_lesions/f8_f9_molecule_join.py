#!/usr/bin/env python3
"""8-oxo-dG rate in molecules grouped by their own methylation (Figure S4).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import json, matplotlib
from collections import defaultdict
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

P='/data/8oxo_project'; PREV=f'{P}/AdMSC_paper/sessions/2026-08-13/results'
OUT=f'{P}/AdMSC_paper/sessions/2026-08-14'
NT_C,PA_C='#2a78d6','#eb6834'; INK,INK2,INK3='#1a1a19','#5c5b55','#a8a69c'
plt.rcParams.update({'font.size':8.5,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
 'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':INK3,'axes.linewidth':0.8,
 'xtick.color':INK2,'ytick.color':INK2,'text.color':INK,'axes.labelcolor':INK,
 'figure.facecolor':'white','savefig.facecolor':'white'})
log=open(f'{OUT}/results/f8_f9_run.log','w')
def L(m):
    print(m,flush=True); log.write(str(m)+'\n'); log.flush()

_TJ=json.load(open(f'{P}/refs/thresholds.json'))
TH={k:v for k,v in _TJ['thresholds'].items() if k not in set(_TJ['high_fp_kmers'])}
assert len(TH)==110
THR=0.70

def cpg_p(f):
    if f[17]!='C': return None
    k=f[15].lower()
    if len(k)!=5: return None
    ok = (k[2]=='c' and k[3]=='g') if f[5]=='+' else (k[2]=='g' and k[1]=='c')
    if not ok: return None
    q=float(f[12])
    return q if f[13]=='m' else 1.0-q

res={}
for arm in ('NT','PA'):
    meth=defaultdict(lambda:[0,0])
    for line in open(f'{PREV}/panel_rerio_calls_{arm}.tsv'):
        f=line.rstrip('\n').split('\t')
        if len(f)<18: continue
        p=cpg_p(f)
        if p is None: continue
        m=meth[f[0]]; m[0]+=1; m[1]+= (p>=THR)
    oxo=defaultdict(lambda:[0,0,0])
    for b in (0,1):
        for line in open(f'{P}/panel_oxo/modcall_{arm}/batch{b}.txt'):
            f=line.rstrip('\n').split('\t')
            if len(f)<4 or f[1]=='basecalls_pos': continue
            t=TH.get(f[3].upper())
            if t is None: continue
            s=float(f[2]); o=oxo[f[0]]
            o[0]+=1; o[1]+= (s>=t); o[2]+= (s>=0.70)
    common=[r for r in meth if r in oxo]
    ncpg=np.array([meth[r][0] for r in common])
    med=int(np.median(ncpg))
    L(f'[{arm}] molecules with both marks: {len(common):,};  CpG calls per molecule '
      f'median {med}, mean {ncpg.mean():.1f}, max {ncpg.max()}')
    rows=[]
    for r in common:
        nc,nm=meth[r]; cg,lk,c70=oxo[r]
        rows.append((nm/nc, nc, cg, lk, c70))
    res[arm]=(rows,med)

L('\n=========== F9  oxidation by molecule-methylation stratum ===========')
f9=[]
for arm in ('NT','PA'):
    rows,med=res[arm]
    for label,sub in (('all molecules',rows),
                      (f'>= median CpG calls ({med})',[x for x in rows if x[1]>=med])):
        b=np.array([x[0] for x in sub])
        q1,q2=np.quantile(b,[1/3,2/3])
        L(f'  {arm}, {label}: n = {len(sub):,}   tertile bounds beta {q1:.3f} / {q2:.3f}')
        for name,sel in (('low',b<=q1),('mid',(b>q1)&(b<=q2)),('high',b>q2)):
            s=[x for x,k in zip(sub,sel) if k]
            cg=sum(x[2] for x in s); lk=sum(x[3] for x in s); c7=sum(x[4] for x in s)
            mb=float(np.mean([x[0] for x in s])) if s else float('nan')
            L(f'      {name:5} n={len(s):>6,}  mean beta {mb:.3f}  assessable G {cg:>9,}  '
              f'calls locked {lk:>4}  >0.70 {c7:>6}  rate70 {1e6*c7/cg if cg else 0:8.1f}/M')
            f9.append((arm,label,name,len(s),mb,cg,lk,c7,1e6*c7/cg if cg else float('nan')))

with open(f'{OUT}/results/f9_oxo_by_methylation_stratum.tsv','w') as o:
    o.write('arm\tsubset\tstratum\tn_molecules\tmean_beta\tassessable_G\tcalls_locked\t'
            'calls_0.70\trate_0.70_perM\n')
    for r in f9: o.write('\t'.join(str(x) for x in r)+'\n')

fig,axes=plt.subplots(1,2,figsize=(9.2,3.8),sharey=True)
for ax,subset in zip(axes,('all molecules',)*1+(f'>= median CpG calls',)*1):
    for arm,c in (('NT',NT_C),('PA',PA_C)):
        sel=[r for r in f9 if r[0]==arm and r[1].startswith(subset[:12])]
        if not sel: continue
        x=[r[4] for r in sel]; y=[r[8] for r in sel]
        ax.plot(x,y,'o-',color=c,linewidth=2.0,markersize=8,markeredgecolor='white',
                markeredgewidth=1.0,label=f'{arm}')
        for r in sel:
            ax.annotate(f'{int(r[7]):,}',xy=(r[4],r[8]),xytext=(0,9),textcoords='offset points',
                        ha='center',fontsize=6.5,color=INK2)
    ax.set_xlabel('mean methylation of the molecules in the stratum')
    ax.grid(axis='y',color=INK3,alpha=0.3,linewidth=0.6); ax.set_axisbelow(True)
    ax.set_title(subset,loc='left',color=INK,pad=6)
axes[0].set_ylabel('8-oxo-dG per million\nassessable G (cut 0.70)')
axes[0].legend(frameon=False,loc='best',fontsize=8)
fig.suptitle('8-oxo-dG counted in molecules grouped by their own methylation',y=1.0,
             fontsize=10.5,color=INK)
fig.text(0.5,0.945,'strata are tertiles of the molecules; the number above each point is the calls it '
         'rests on',ha='center',fontsize=7.5,color=INK2)
fig.tight_layout(rect=[0,0,1,0.93])
for e in ('png','pdf'): fig.savefig(f'{OUT}/figures/F9_oxo_by_methylation.{e}',dpi=300,bbox_inches='tight')
L('\n  wrote figures/F9_oxo_by_methylation.{png,pdf}')

fig,axes=plt.subplots(1,2,figsize=(9.2,3.9),sharex=True,sharey=True)
for ax,(arm,c) in zip(axes,(('NT',NT_C),('PA',PA_C))):
    rows,med=res[arm]
    sub=[x for x in rows if x[1]>=med]
    x=np.array([r[0] for r in sub]); y=np.array([1e6*r[4]/r[2] if r[2] else np.nan for r in sub])
    ax.plot(x,y,'o',color=c,markersize=2.2,alpha=0.25,markeredgecolor='none')
    bins=np.linspace(0,1,11); idx=np.digitize(x,bins)-1
    bx=[];by=[]
    for i in range(10):
        s=idx==i
        if s.sum()>20: bx.append((bins[i]+bins[i+1])/2); by.append(np.nanmedian(y[s]))
    ax.plot(bx,by,'-',color=INK,linewidth=2.0,zorder=4)
    ax.plot(bx,by,'o',color=INK,markersize=5,markeredgecolor='white',markeredgewidth=1,zorder=5)
    ax.set_title(f'{arm}  ({len(sub):,} molecules)',loc='left',color=INK,pad=6)
    ax.set_xlabel('methylation of the molecule')
    ax.grid(color=INK3,alpha=0.3,linewidth=0.6); ax.set_axisbelow(True)
axes[0].set_ylabel('8-oxo-dG per million assessable G\non that molecule (cut 0.70)')
axes[0].set_ylim(-200,12000)
fig.suptitle('Both marks read from the same molecule',y=1.0,fontsize=10.5,color=INK)
fig.text(0.5,0.945,'one point per molecule, restricted to molecules carrying at least the median '
         'number of CpG calls; the dark line is the median per methylation bin',
         ha='center',fontsize=7.5,color=INK2)
fig.tight_layout(rect=[0,0,1,0.93])
for e in ('png','pdf'): fig.savefig(f'{OUT}/figures/F8_molecule_joint.{e}',dpi=300,bbox_inches='tight')
L('  wrote figures/F8_molecule_joint.{png,pdf}')
log.close()
