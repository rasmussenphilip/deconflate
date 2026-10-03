# Monte Carlo reproduction of Rasmussen et al. (2024) Table 5 (20,000 draws;
# PERT central values treated as modes; negative OR draws rejected).
# Run from this directory: python3 reference_2024_montecarlo.py
import numpy as np, sys
sys.path.insert(0,'.')
from python_reference import or_to_joint
rng=np.random.default_rng(7)
D=['CK','CM','DA','DYS','LAM','MET','MF','OC','PTB','RP','SCK','SCM']; n=12
def pert(mn,mean,mx,size):
    mode=mean; a=1+4*(mode-mn)/(mx-mn); b=1+4*(mx-mode)/(mx-mn)
    return mn+(mx-mn)*rng.beta(a,b,size)
N=20000
inc={'CK':('b',13.95,441.59),'CM':('b',8.55,18.20),'DA':('b',27.50,1245.11),'DYS':('p',1.90,5.99,10.80),
 'LAM':('b',78.29,227.42),'MET':('b',17.73,167.35),'MF':('b',7.76,313.41),'OC':('p',2.70,11.46,19.07),
 'PTB':('p',1.19,10.01,21.08),'RP':('b',33.75,239.46),'SCK':('b',178.14,193.73),'SCM':('b',116.23,167.02)}
def draw_inc(spec):
    if spec[0]=='b': return rng.beta(spec[1],spec[2],N)
    return pert(spec[1]/100,spec[2]/100,spec[3]/100,N)
I={d:draw_inc(inc[d]) for d in D}
P=np.column_stack([I[d] if d=='PTB' else 1-np.exp(-I[d]) for d in D])
ORs="""CK:CM p 2.13 1.20 3.40;CK:LAM p 1.65 1.20 2.40;CK:MF n 1.60 .13;CK:OC p 1.97 1.30 4.10;CK:RP p 1.55 1.00 1.90;CK:SCK n 6.95 1.28;CK:SCM n 2.40 .41;CM:LAM f 2.10;CM:PTB n 1.89 .20;CM:RP n 2.70 .33;CM:SCK n 1.64 .20;CM:SCM p 3.05 1.30 6.50;DA:CM p 3.45 1.40 4.80;DA:MF n 2.50 .48;DA:RP p 3.50 1.60 4.60;DA:SCK n 3.87 .34;DA:SCM n 3.60 1.35;DYS:LAM n 2.09 .26;DYS:OC f 0.40;DYS:RP p 2.74 1.25 5.96;LAM:OC n 2.63 1.44;LAM:PTB n 2.70 1.22;LAM:RP n 1.50 .31;LAM:SCK n 2.01 .20;MET:CK p 2.42 1.20 10.4;MET:CM p 2.30 1.20 3.8;MET:DA p 3.40 1.60 7.60;MET:DYS p 2.95 .98 9.72;MET:LAM n 6.10 1.45;MET:MF n 1.50 .15;MET:OC p 1.94 1.20 3.00;MET:RP p 3.53 1.8 6.52;MET:SCK n 1.94 .09;MF:DYS n 9.70 1.30;MF:LAM f 3.60;MF:RP n 2.40 .20;RP:OC p 2.18 1.78 2.57;SCK:RP n 1.52 .19"""
def draw(spec):
    t=spec[0]; v=list(map(float,spec[1:]))
    if t=='f': return np.full(N,v[0])
    if t=='n': return rng.normal(v[0],v[1],N)
    return pert(v[1],v[0],v[2],N)
pairs=[]
for s in ORs.split(';'):
    p=s.split(); a,b=p[0].split(':'); pairs.append((D.index(a),D.index(b),draw(p[1:])))
Y={'CK':'p .43 .24 1.04','CM':'n 3.25 .76','DA':'p 2.84 -1.45 9.19','DYS':'n 4.92 .97','LAM':'n 4.81 .87','MET':'n 5.61 1.35','MF':'f .54','OC':'p 3.75 1.71 4.33','PTB':'n 4.30 .67','RP':'n 4.20 1.15','SCK':'n 8.40 1.19','SCM':'n 6.29 1.20'}
F={'CK':'n 1.45 .36','CM':'n 8.42 2.42','DA':'n 1.08 2.04','DYS':'n 2.40 .93','LAM':'p 3.30 1.19 10.71','MET':'n 14.67 8.54','MF':'p 2.41 2.03 3.10','OC':'p 9.69 5.04 21.43','PTB':'n 5.35 2.53','RP':'n 6.76 1.56','SCK':'n 1.12 1.82','SCM':'p 0.26 -.12 5.68'}
T5Y=[.03,1.36,1.18,3.48,2.62,2.87,.07,2.59,3.37,2.30,7.11,5.58]
T5F=[.34,6.09,.78,1.11,1.86,11.22,1.06,9.03,4.23,3.74,.39,.04]
for lab,spec,T5 in [('yield',Y,T5Y),('fert',F,T5F)]:
    M=np.column_stack([draw(spec[d].split()) for d in D])
    out_pub=np.zeros((N,n)); out_lin=np.zeros((N,n)); bad=0
    for r in range(N):
        E=np.zeros((n,n)); ok=True
        for i,k,v in pairs:
            o=v[r]
            if o<=0: ok=False; break
            pi,pk=P[r,i],P[r,k]; p11=or_to_joint(pi,pk,o)
            E[k,i]=p11/pi-(pk-p11)/(1-pi); E[i,k]=p11/pk-(pi-p11)/(1-pk)
        if not ok: bad+=1; out_pub[r]=np.nan; out_lin[r]=np.nan; continue
        m=M[r]; conf=E.T@m - 0  # E.T[i,k]=E[k,i]
        out_pub[r]=m**2/(m+conf); out_lin[r]=np.linalg.solve(np.eye(n)+E.T,m)
    print(lab,'rejected',bad)
    print(' %-4s %6s %8s %8s %8s'%('', 'T5','pubMean','pubMed','linMean'))
    for j,d in enumerate(D): print(' %-4s %6.2f %8.2f %8.2f %8.2f'%(d,T5[j],np.nanmean(out_pub[:,j]),np.nanmedian(out_pub[:,j]),np.nanmean(out_lin[:,j])))
