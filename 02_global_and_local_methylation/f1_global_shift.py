#!/usr/bin/env python3
"""Pooled % 5mC per region class, both arms (Figure 9a).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import csv, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R='/data/8oxo_project/AdMSC_DMR/results'; OUT='/data/8oxo_project/AdMSC_paper/sessions/2026-08-14'
NT_C,PA_C='#2a78d6','#eb6834'; INK,INK2,INK3='#1a1a19','#5c5b55','#a8a69c'
plt.rcParams.update({'font.size':8.5,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8.5,
 'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':INK3,'axes.linewidth':0.8,
 'xtick.color':INK2,'ytick.color':INK2,'text.color':INK,'axes.labelcolor':INK,
 'figure.facecolor':'white','savefig.facecolor':'white'})

CLASSES=[('CpG islands','layer1_cpg_islands.tsv',0.8894),
         ('shelves','layer1_cpg_shelves.tsv',0.8916),
         ('shores','layer1_cpg_shores.tsv',0.8967),
         ('exons','layer1_cpg_exons_merged.tsv',0.8945),
         ('gene bodies','layer1_cpg_genebodies.tsv',0.8971),
         ('promoters (TSS+-1kb)','layer1_cpg_promoters_tss1kb.tsv',0.9130)]
GW=('whole genome',None,0.89845,29.8706,26.8370)

log=open(f'{OUT}/results/f1_run.log','w')
def L(m):
    print(m,flush=True); log.write(str(m)+'\n'); log.flush()

L('  class                     NT %     PA %    ratio   ref.    check')
rows=[]
for name,fn,locked in CLASSES:
    nv=nm=pv=pm=0
    for r in csv.DictReader(open(f'{R}/{fn}'),delimiter='\t'):
        nv+=float(r['NT_valid']); nm+=float(r['NT_mod'])
        pv+=float(r['PA_valid']); pm+=float(r['PA_mod'])
    nt,pa=100*nm/nv,100*pm/pv; ratio=(pm/pv)/(nm/nv)
    ok='reproduces' if abs(ratio-locked)<0.0015 else f'DIFFERS by {ratio-locked:+.4f}'
    L(f'  {name:24}{nt:8.2f}{pa:9.2f}{ratio:9.4f}{locked:8.4f}   {ok}')
    rows.append((name,nt,pa,ratio,locked,int(nv),int(pv)))
L(f'  {GW[0]:24}{GW[3]:8.2f}{GW[4]:9.2f}{GW[2]:9.5f}{GW[2]:8.5f}   cited from the 2026-08-10 computation')
rows.append((GW[0],GW[3],GW[4],GW[2],GW[2],0,0))

with open(f'{OUT}/results/f1_global_shift.tsv','w') as o:
    o.write('region_class\tNT_pct\tPA_pct\tratio_PA_NT\tlocked_M006\tNT_valid\tPA_valid\n')
    for r in rows: o.write('\t'.join(str(x) for x in r)+'\n')

order=sorted(rows,key=lambda r:r[3])
fig,axes=plt.subplots(1,2,figsize=(9.4,3.9),gridspec_kw={'width_ratios':[1.25,1]})
ax=axes[0]; y=range(len(order))
for i,(name,nt,pa,ratio,lk,_,_) in enumerate(order):
    ax.plot([nt,pa],[i,i],color=INK3,linewidth=1.6,zorder=1)
    ax.plot(nt,i,'o',color=NT_C,markersize=8,markeredgecolor='white',markeredgewidth=1,zorder=3)
    ax.plot(pa,i,'o',color=PA_C,markersize=8,markeredgecolor='white',markeredgewidth=1,zorder=3)
ax.set_yticks(list(y)); ax.set_yticklabels([r[0] for r in order])
ax.set_xlabel('pooled 5mCpG (%)'); ax.grid(axis='x',color=INK3,alpha=0.3,linewidth=0.6)
ax.set_axisbelow(True); ax.set_ylim(-0.6,len(order)-0.4)
h=[plt.Line2D([],[],color=c,marker='o',linestyle='none',markersize=8,markeredgecolor='white',
              markeredgewidth=1,label=l) for c,l in ((NT_C,'untreated (NT)'),(PA_C,'palmitate (PA)'))]
ax.legend(handles=h,frameon=False,loc='lower right',fontsize=8)
ax.set_title('absolute level',loc='left',color=INK,pad=16)

ax=axes[1]
for i,(name,nt,pa,ratio,lk,_,_) in enumerate(order):
    c=INK if name=='whole genome' else NT_C
    ax.plot(ratio,i,'o',color=c,markersize=8,markeredgecolor='white',markeredgewidth=1,zorder=3)
    ax.annotate(f'{ratio:.3f}',xy=(ratio,i),xytext=(9,0),textcoords='offset points',
                va='center',fontsize=7.5,color=INK)
ax.axvline(GW[2],color=INK3,linestyle=(0,(4,3)),linewidth=1.1,zorder=1)
ax.annotate('whole genome',xy=(GW[2],len(order)-0.55),xytext=(-5,0),textcoords='offset points',
            ha='right',fontsize=7,color=INK2)
ax.set_yticks(list(y)); ax.set_yticklabels([])
ax.set_xlabel('PA / NT'); ax.set_xlim(0.878,0.925)
ax.grid(axis='x',color=INK3,alpha=0.3,linewidth=0.6); ax.set_axisbelow(True)
ax.set_ylim(-0.6,len(order)-0.4); ax.set_title('change',loc='left',color=INK,pad=16)
fig.suptitle('CpG methylation falls under palmitate in every region class, and by a similar amount',
             y=1.0,fontsize=10.5,color=INK)
fig.tight_layout(rect=[0,0,1,0.94])
for e in ('png','pdf'): fig.savefig(f'{OUT}/figures/F1_global_shift_by_class.{e}',dpi=300,bbox_inches='tight')
L('\nwrote figures/F1_global_shift_by_class.{png,pdf} and results/f1_global_shift.tsv')
log.close()
