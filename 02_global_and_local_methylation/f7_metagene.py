#!/usr/bin/env python3
"""Metagene of % 5mC along scaled genes (Figure 9d).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import os, time
import numpy as np

P='/data/8oxo_project'; OUT=f'{P}/AdMSC_paper/sessions/2026-08-14'
TMP='<SESSION_TMP>'
FLANK, NF, NB = 2000, 20, 100
MINLEN=1000
log=open(f'{OUT}/results/f7_run.log','w')
def L(m):
    print(m,flush=True); log.write(str(m)+'\n'); log.flush()

genes={}
n_all=n_short=0
for line in open(f'{P}/AdMSC_DMR/features/genebodies.bed'):
    f=line.rstrip('\n').split('\t')
    if len(f)<6: continue
    n_all+=1
    s,e=int(f[1]),int(f[2])
    if e-s<MINLEN: n_short+=1; continue
    genes.setdefault(f[0],[]).append((s,e,f[5]))
L(f'[1] genes {n_all:,}; excluded shorter than {MINLEN} bp: {n_short:,}; '
  f'used {sum(len(v) for v in genes.values()):,} on {len(genes)} contigs')

NBINS=NF+NB+NF
tot={a:(np.zeros(NBINS),np.zeros(NBINS)) for a in ('NT','PA')}
for arm in ('NT','PA'):
    t0=time.time(); nc=0
    for ent in sorted(os.scandir(f'{TMP}/{arm}'),key=lambda e:e.name):
        chrom=ent.name[:-4]
        if chrom not in genes: continue
        a=np.loadtxt(ent.path,dtype=np.int64)
        if a.size==0: continue
        if a.ndim==1: a=a.reshape(1,2)
        o=np.argsort(a[:,0]); pos=a[o,0]; mod=a[o,1].astype(np.float64)
        nc+=len(pos)
        for s,e,st in genes[chrom]:
            lo,hi=s-FLANK,e+FLANK
            i0,i1=np.searchsorted(pos,[lo,hi])
            if i1<=i0: continue
            p=pos[i0:i1]; m=mod[i0:i1]
            b=np.empty(len(p),dtype=np.int64)
            up=p<s; dn=p>=e; body=~(up|dn)
            b[up]=np.clip(((p[up]-lo)//(FLANK//NF)),0,NF-1)
            b[body]=NF+np.clip(((p[body]-s)*NB//(e-s)),0,NB-1)
            b[dn]=NF+NB+np.clip(((p[dn]-e)//(FLANK//NF)),0,NF-1)
            if st=='-': b=NBINS-1-b
            np.add.at(tot[arm][0],b,1.0)
            np.add.at(tot[arm][1],b,m)
    L(f'[2] {arm}: {nc:,} CpG calls scanned, '
      f'{int(tot[arm][0].sum()):,} placed on the metagene axis, {time.time()-t0:.0f} s')

with open(f'{OUT}/results/f7_metagene.tsv','w') as o:
    o.write('bin\tregion\trel_pos\tNT_valid\tNT_mod\tNT_pct\tPA_valid\tPA_mod\tPA_pct\tratio\n')
    for i in range(NBINS):
        if i<NF: reg,rel='upstream',(i-NF)*100
        elif i<NF+NB: reg,rel='body',(i-NF)/NB
        else: reg,rel='downstream',(i-NF-NB)*100
        nv,nm=tot['NT'][0][i],tot['NT'][1][i]; pv,pm=tot['PA'][0][i],tot['PA'][1][i]
        nt=100*nm/nv if nv else float('nan'); pa=100*pm/pv if pv else float('nan')
        o.write(f'{i}\t{reg}\t{rel}\t{int(nv)}\t{int(nm)}\t{nt:.4f}\t{int(pv)}\t{int(pm)}\t{pa:.4f}\t'
                f'{(pm/pv)/(nm/nv) if nv and pv and nm else float("nan"):.5f}\n')
L(f'  wrote results/f7_metagene.tsv')
for lbl,i in (('2 kb upstream',0),('TSS',NF),('mid-body',NF+NB//2),('TES',NF+NB-1),
              ('2 kb downstream',NBINS-1)):
    nv,nm=tot['NT'][0][i],tot['NT'][1][i]; pv,pm=tot['PA'][0][i],tot['PA'][1][i]
    L(f'    {lbl:>16}: NT {100*nm/nv:5.2f} %  PA {100*pm/pv:5.2f} %  ratio {(pm/pv)/(nm/nv):.4f}'
      f'   calls {int(nv):,}/{int(pv):,}')
log.close()
