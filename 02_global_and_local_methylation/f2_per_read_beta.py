#!/usr/bin/env python3
"""Distribution of reads by their own methylation fraction (Figure 9b).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import csv, matplotlib
from collections import defaultdict
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SRC='/data/8oxo_project/AdMSC_paper/sessions/2026-08-11/results/level0a_read_beta.tsv'
OUT='/data/8oxo_project/AdMSC_paper/sessions/2026-08-14'
NT_C,PA_C='#2a78d6','#eb6834'; INK,INK2,INK3='#1a1a19','#5c5b55','#a8a69c'
plt.rcParams.update({'font.size':8.5,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
 'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':INK3,'axes.linewidth':0.8,
 'xtick.color':INK2,'ytick.color':INK2,'text.color':INK,'axes.labelcolor':INK,
 'figure.facecolor':'white','savefig.facecolor':'white'})
log=open(f'{OUT}/results/f2_run.log','w')
def L(m):
    print(m,flush=True); log.write(str(m)+'\n'); log.flush()

rows=list(csv.DictReader(open(SRC),delimiter='\t'))
thr=sorted({r['threshold'] for r in rows})
L(f'  thresholds present: {thr}')
T='0.7' if '0.7' in thr else thr[0]
d=defaultdict(list)
for r in rows:
    if r['threshold']!=T: continue
    d[r['arm']].append((float(r['beta_bin_lo']),float(r['beta_bin_hi']),float(r['n_reads'])))
for a in d: d[a].sort()
tot={a:sum(x[2] for x in v) for a,v in d.items()}
L(f'  threshold used: {T}   reads  ' + '  '.join(f'{a} {int(tot[a]):,}' for a in sorted(d)))

with open(f'{OUT}/results/f2_per_read_beta.tsv','w') as o:
    o.write('arm\tthreshold\tbeta_lo\tbeta_hi\tn_reads\tfraction_of_reads\n')
    for a in sorted(d):
        for lo,hi,n in d[a]: o.write(f'{a}\t{T}\t{lo}\t{hi}\t{int(n)}\t{n/tot[a]:.6f}\n')

fig,axes=plt.subplots(1,2,figsize=(9.2,3.7))
for k,(ax,logy) in enumerate(zip(axes,(False,True))):
    for arm,c in (('NT',NT_C),('PA',PA_C)):
        v=d[arm]; x=[(lo+hi)/2 for lo,hi,_ in v]; y=[n/tot[arm] for _,_,n in v]
        ax.step(x,y,where='mid',color=c,linewidth=2.0)
        ax.plot(x,y,'o',color=c,markersize=3.5,markeredgecolor='white',markeredgewidth=0.7)
    ax.grid(axis='y',color=INK3,alpha=0.3,linewidth=0.6); ax.set_axisbelow(True)
    ax.set_xlabel('methylation of the molecule (fraction of its CpG called methylated)')
    ax.set_ylabel('fraction of reads')
    if logy:
        ax.set_yscale('log'); ax.set_title('same data, log scale',loc='left',color=INK,pad=6)
    else:
        ax.set_title('all molecules',loc='left',color=INK,pad=6)
        pass
h=[plt.Line2D([],[],color=c,linewidth=2.0,label=l) for c,l in ((NT_C,'untreated (NT)'),(PA_C,'palmitate (PA)'))]
fig.legend(handles=h,loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(0.5,-0.03))
fig.suptitle('Per-molecule CpG methylation in untreated and palmitate-treated AdMSC',
             y=1.02,fontsize=10.5,color=INK)
fig.text(0.5,0.955,'both arms are bimodal; unmethylated molecules 26.1 % -> 29.2 %, '
         'molecules above 0.6 methylation 14.8 % -> 12.2 %',ha='center',fontsize=7.5,color=INK2)
fig.tight_layout(rect=[0,0.06,1,0.93])
for e in ('png','pdf'): fig.savefig(f'{OUT}/figures/F2_per_read_beta.{e}',dpi=300,bbox_inches='tight')
L('  wrote figures/F2_per_read_beta.{png,pdf} and results/f2_per_read_beta.tsv')
log.close()
