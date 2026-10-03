# Reference calculations for the 2022 UK dairy reproduction (published and
# simultaneous methods) and for cell-wise Shapley values (multiplicative loss).
# Run from this directory: python3 reference_2022_shapley.py
import numpy as np, sys, itertools
sys.path.insert(0,'.')
from python_reference import or_to_joint, excess_matrix, published, simultaneous, ipf
np.set_printoptions(precision=10, suppress=True, linewidth=150)
D=['CO','DA','DYS','FAS','GIN','LAM','MAS','MET','MF','NEO','PTB','RP','SCK']
P=dict(LAM=.30,MAS=.30,SCK=.22,GIN=.21,NEO=.15,MET=.10,FAS=.10,CO=.09,MF=.08,PTB=.07,RP=.05,DA=.03,DYS=.02)
ors="RP:MET 6.20;DA:SCK 4.25;RP:DYS 4.10;MET:DA 3.40;MET:DYS 3.20;LAM:PTB 2.70;MF:DA 2.50;MAS:MET 2.30;DA:RP 2.20;MAS:DA 2.10;SCK:MF 2.10;LAM:SCK 2.01;MAS:MF 1.90;MAS:PTB 1.89;MAS:CO 1.65;MAS:SCK 1.64;SCK:CO 1.60;MET:SCK 1.40;SCK:RP 1.20"
Y=dict(RP=7.38,FAS=7.33,PTB=5.90,LAM=5.54,MAS=4.57,NEO=4.20,DYS=4.05,DA=4.04,MET=3.95,GIN=3.28,SCK=3.05,MF=0.41,CO=0)
F=dict(LAM=12.47,CO=11.26,NEO=7.21,DYS=6.96,PTB=5.79,MET=4.74,RP=2.74,SCK=1.50,GIN=1.20,DA=0,MAS=0,MF=0,FAS=0)
H=dict(DA=3.83,LAM=3.40,MAS=2.78,MF=2.50,PTB=2.40,MET=2.20,SCK=2.10,DYS=1.90,NEO=1.60,CO=1,GIN=1,FAS=1,RP=1)
n=len(D); Pv=np.array([P[d] for d in D]); OR=np.ones((n,n))
for s in ors.split(';'):
    pr,v=s.split(); a,b=pr.split(':'); i,k=D.index(a),D.index(b); OR[i,k]=OR[k,i]=float(v)
c=0.27
ep=np.array([0.0 if H[d]==1 else (lambda p11: p11/P[d]-(c-p11)/(1-P[d]))(or_to_joint(c,P[d],H[d])) for d in D])
print("cull raw ep", ep)
outs={'yield':(np.array([Y[d] for d in D])/100,'dec',8737,0.3022),'fertility':(np.array([F[d] for d in D])/100,'inc',401,13*0.3022),'culling':(ep,'inc',27,13.3536)}
for meth,fn in [('published',published),('simultaneous',simultaneous)]:
  tot=0
  for o,(m,dirn,x,pr) in outs.items():
    adj=np.nan_to_num(fn(Pv,OR,m)) if meth=='simultaneous' else np.array([0 if m[i]==0 else v for i,v in enumerate(fn(Pv,OR,m))])
    L=adj@Pv; xh=x/(1-L) if dirn=='dec' else x/(1+L); gap=abs(xh-x); tot+=gap*pr
    print(meth,o,'xh',repr(xh),'value',repr(gap*pr))
    if o=='culling' and meth=='published':
      Hv=np.array([H[d] for d in D]); print('  hr_pub',np.where(ep>0,adj*Hv/np.where(ep>0,ep,1),Hv))
    if meth=='published': print('  adj',adj)
  print(meth,'TOTAL',repr(tot),'+vet',repr(tot+71.09))
# Shapley multiplicative on supplement
P3=np.array([.1,.15,.2]); OR3=np.ones((3,3)); OR3[0,1]=OR3[1,0]=2; OR3[1,2]=OR3[2,1]=3
cells,p=ipf(P3,OR3)
m=np.array([.02,.04,.06])
from math import factorial
phi=np.zeros(3); total=0
for cell,pr in zip(cells,p):
    pres=[i for i in range(3) if cell[i]]
    k=len(pres)
    if k==0: continue
    def v(S): 
        return 1-np.prod([1-m[i] for i in S])
    for i in pres:
        others=[j for j in pres if j!=i]
        for r in range(len(others)+1):
            for S in itertools.combinations(others,r):
                w=factorial(len(S))*factorial(k-len(S)-1)/factorial(k)
                phi[i]+=pr*w*(v(list(S)+[i])-v(list(S)))
    total+=pr*v(pres)
print('shapley mult',phi,'total',total, phi.sum())
