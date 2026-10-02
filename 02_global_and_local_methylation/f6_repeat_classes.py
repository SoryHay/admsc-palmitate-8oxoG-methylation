#!/usr/bin/env python3
"""Pooled % 5mC per RepeatMasker class (Figure 9c).
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import gzip, os, subprocess, sys, time
from collections import defaultdict
import numpy as np

P='/data/8oxo_project'; OUT=f'{P}/AdMSC_paper/sessions/2026-08-14'
TMP='<SESSION_TMP>'
THR=0.70
log=open(f'{OUT}/results/f6_run.log','w')
def L(m):
    print(m,flush=True); log.write(str(m)+'\n'); log.flush()

AWK = r'''
{
  if ($18 != "C") next;
  k = tolower($16);
  if (length(k) != 5) next;
  if ($6 == "+") { if (substr(k,3,1) != "c" || substr(k,4,1) != "g") next }
  else           { if (substr(k,3,1) != "g" || substr(k,2,1) != "c") next }
  q = $13 + 0;
  p = ($14 == "m") ? q : 1 - q;
  print $3 "\t" (p >= THR ? 1 : 0) >> (DIR "/" $4 ".txt");
}
'''

for arm in ('NT','PA'):
    d=f'{TMP}/{arm}'
    if os.path.isdir(d) and any(os.scandir(d)):
        L(f'[1] {arm}: stage-1 output already present in {d}, reusing')
        continue
    os.makedirs(d,exist_ok=True)
    t0=time.time()
    L(f'[1] {arm}: stage 1 starting -> {d}')
    cmd=['awk','-F','\t','-v',f'THR={THR}','-v',f'DIR={d}',AWK,
         f'{P}/AdMSC_DMR/readlevel_calls/{arm}.calls.tsv']
    subprocess.run(cmd,check=True)
    n=sum(1 for _ in os.scandir(d))
    L(f'[1] {arm}: stage 1 done in {time.time()-t0:.0f} s, {n} contigs written')

L('[2] reading RepeatMasker')
iv=defaultdict(list)
with gzip.open(f'{P}/refs/T2T-CHM13v2.0.rmsk.out.gz','rt') as fh:
    for i,line in enumerate(fh):
        f=line.split()
        if len(f)<11 or not f[0].isdigit(): continue
        cls=f[10].split('/')[0]
        iv[f[4]].append((int(f[5]),int(f[6]),cls))
classes=sorted({c for v in iv.values() for _,_,c in v})
CID={c:i+1 for i,c in enumerate(classes)}
L(f'[2] {sum(len(v) for v in iv.values()):,} repeat intervals on {len(iv)} contigs, '
  f'{len(classes)} classes')

counts={arm:defaultdict(lambda:[0,0]) for arm in ('NT','PA')}
for arm in ('NT','PA'):
    d=f'{TMP}/{arm}'
    for ent in sorted(os.scandir(d),key=lambda e:e.name):
        chrom=ent.name[:-4]
        a=np.loadtxt(ent.path,dtype=np.int64)
        if a.size==0: continue
        if a.ndim==1: a=a.reshape(1,2)
        pos,mod=a[:,0]+1,a[:,1]
        mx=int(pos.max())+2
        cmap=np.zeros(mx,dtype=np.uint8)
        for s,e,cls in iv.get(chrom,()):
            if s<mx: cmap[s:min(e+1,mx)]=CID[cls]
        cid=cmap[pos]
        for k in np.unique(cid):
            sel=cid==k
            name='non-repeat' if k==0 else classes[k-1]
            counts[arm][name][0]+=int(sel.sum())
            counts[arm][name][1]+=int(mod[sel].sum())
        del cmap
    tot=sum(v[0] for v in counts[arm].values())
    L(f'[3] {arm}: {tot:,} CpG calls assigned to a class')

rows=[]
for name in sorted(set(counts['NT'])|set(counts['PA']),
                   key=lambda n:-counts['NT'][n][0]):
    nv,nm=counts['NT'][name]; pv,pm=counts['PA'][name]
    if nv<10000 or pv<10000: continue
    nt,pa=100*nm/nv,100*pm/pv
    rows.append((name,nv,nm,nt,pv,pm,pa,(pm/pv)/(nm/nv)))
with open(f'{OUT}/results/f6_repeat_classes.tsv','w') as o:
    o.write('repeat_class\tNT_valid\tNT_mod\tNT_pct\tPA_valid\tPA_mod\tPA_pct\tratio_PA_NT\n')
    for r in rows: o.write('\t'.join(f'{x:.6g}' if isinstance(x,float) else str(x) for x in r)+'\n')
L(f'\n  {"repeat class":<20}{"NT %":>8}{"PA %":>8}{"ratio":>9}{"NT calls":>14}')
for r in rows: L(f'  {r[0]:<20}{r[3]:>8.2f}{r[6]:>8.2f}{r[7]:>9.4f}{r[1]:>14,}')
L(f'\n  wrote results/f6_repeat_classes.tsv  (classes with >= 10,000 calls in both arms)')
log.close()
