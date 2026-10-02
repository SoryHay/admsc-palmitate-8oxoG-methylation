#!/usr/bin/env python3
"""8-oxo-dG by trinucleotide context (Figure S8).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import json, matplotlib
from collections import defaultdict
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

P='/data/8oxo_project'; OUT=f'{P}/AdMSC_paper/sessions/2026-08-14'
NT_C,PA_C='#2a78d6','#eb6834'; INK,INK2,INK3='#1a1a19','#5c5b55','#a8a69c'
plt.rcParams.update({'font.size':8,'axes.labelsize':9,'xtick.labelsize':7,'ytick.labelsize':8,
 'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':INK3,'axes.linewidth':0.8,
 'xtick.color':INK2,'ytick.color':INK2,'text.color':INK,'axes.labelcolor':INK,
 'figure.facecolor':'white','savefig.facecolor':'white'})
log=open(f'{OUT}/results/f5_run.log','w')
def L(m):
    print(m,flush=True); log.write(str(m)+'\n'); log.flush()

_TJ=json.load(open(f'{P}/refs/thresholds.json'))
TH={k:v for k,v in _TJ['thresholds'].items() if k not in set(_TJ['high_fp_kmers'])}
assert len(TH)==110
L(f'[1] allowlist {len(TH)} 5-mers -> collapsed to the central 3-mer')

CUTS=['locked','0.70']
data={}
for arm in ('NT','PA'):
    opp=defaultdict(int); calls={c:defaultdict(int) for c in CUTS}
    n=0
    for b in (0,1):
        for line in open(f'{P}/panel_oxo/modcall_{arm}/batch{b}.txt'):
            f=line.rstrip('\n').split('\t')
            if len(f)<4 or f[1]=='basecalls_pos': continue
            k=f[3].upper(); t=TH.get(k)
            if t is None: continue
            n+=1; tri=k[1:4]
            opp[tri]+=1
            s=float(f[2])
            if s>=t: calls['locked'][tri]+=1
            if s>=0.70: calls['0.70'][tri]+=1
    data[arm]=(opp,calls)
    L(f'[2] {arm}: assessable G {n:,}; calls locked {sum(calls["locked"].values())}, '
      f'>0.70 {sum(calls["0.70"].values())}; contexts seen {len(opp)}')

tris=sorted(set(data['NT'][0])|set(data['PA'][0]))
with open(f'{OUT}/results/f5_oxo_context.tsv','w') as o:
    o.write('trimer\tarm\tassessable_G\tcalls_locked\tcalls_0.70\t'
            'share_of_opportunity\tshare_of_calls_0.70\tenrichment_0.70\n')
    for arm in ('NT','PA'):
        opp,calls=data[arm]; to=sum(opp.values()); tc=sum(calls['0.70'].values())
        for t in tris:
            so=opp[t]/to if to else 0; sc=calls['0.70'][t]/tc if tc else 0
            o.write(f'{t}\t{arm}\t{opp[t]}\t{calls["locked"][t]}\t{calls["0.70"][t]}\t'
                    f'{so:.6f}\t{sc:.6f}\t{(sc/so if so else float("nan")):.4f}\n')

enr={}
for arm in ('NT','PA'):
    opp,calls=data[arm]; to=sum(opp.values()); tc=sum(calls['0.70'].values())
    enr[arm]=[( (calls['0.70'][t]/tc)/(opp[t]/to) if opp[t] and tc else np.nan) for t in tris]
order=sorted(range(len(tris)),key=lambda i:-(enr['NT'][i] if enr['NT'][i]==enr['NT'][i] else -9))
tris_o=[tris[i] for i in order]

fig,axes=plt.subplots(2,1,figsize=(9.6,6.0),sharex=True,gridspec_kw={'height_ratios':[1,1.35]})
ax=axes[0]
opp_share=[100*data['NT'][0][t]/sum(data['NT'][0].values()) for t in tris_o]
ax.bar(range(len(tris_o)),opp_share,color=INK3,width=0.72)
ax.set_ylabel('share of assessable\nguanines (%)')
ax.grid(axis='y',color=INK3,alpha=0.3,linewidth=0.6); ax.set_axisbelow(True)
ax.set_title('how common each context is',loc='left',color=INK,pad=6)
ax=axes[1]
x=np.arange(len(tris_o))
ax.axhline(1.0,color=INK2,linewidth=1.0)
ax.plot(x,[enr['NT'][i] for i in order],'o-',color=NT_C,linewidth=1.8,markersize=6,
        markeredgecolor='white',markeredgewidth=0.8,label='untreated (NT)')
ax.plot(x,[enr['PA'][i] for i in order],'s-',color=PA_C,linewidth=1.8,markersize=5,
        markeredgecolor='white',markeredgewidth=0.8,label='palmitate (PA)')
ax.set_xticks(x); ax.set_xticklabels(tris_o,rotation=90)
ax.set_ylabel('enrichment of 8-oxo-dG\n(share of calls / share of opportunity)')
ax.set_xlabel('sequence context, central guanine with one base either side')
ax.grid(axis='y',color=INK3,alpha=0.3,linewidth=0.6); ax.set_axisbelow(True)
ax.legend(frameon=False,loc='upper right',fontsize=8)
fig.suptitle('8-oxo-dG by sequence context, nuclear panel',y=1.0,fontsize=10.5,color=INK)
fig.text(0.5,0.945,'calls counted at the relaxed cut 0.70 so that every context is populated; '
         'contexts ordered by untreated enrichment',ha='center',fontsize=7.5,color=INK2)
fig.tight_layout(rect=[0,0,1,0.93])
for e in ('png','pdf'): fig.savefig(f'{OUT}/figures/F5_oxo_context.{e}',dpi=300,bbox_inches='tight')
L('  wrote figures/F5_oxo_context.{png,pdf} and results/f5_oxo_context.tsv')
top=sorted(zip(tris,enr['NT']),key=lambda z:-(z[1] if z[1]==z[1] else -9))[:5]
L('  most enriched contexts in NT: '+', '.join(f'{t} {v:.2f}' for t,v in top))
log.close()
