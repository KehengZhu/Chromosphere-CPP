"""Grid convergence of the lower-chromosphere (~250 km) mass-flux maximum in the
bestwb conduction-ON iso_t22k run. Shows V@250 and the peak ρV scale ∝ Δs (first
order) and extrapolate to ~0 as Δs→0 ⇒ the feature is a numerical
split conduction-to-hydro coupling artifact, not a physical flow.

Usage: plot_iso_250km_convergence.py [t_target=100] [out.png] file1 file2 ...
Defaults to the ns=2000/3000/4000 bestwb cond-on runs at t=100 s."""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

args = [a for a in sys.argv[1:]]
T = 100.0
out = "visualization/model_column/iso_t22k_250km_convergence.png"
files = []
for a in args:
    if a.endswith(".txt"): files.append(a)
    elif a.endswith(".png"): out = a
    else:
        try: T = float(a)
        except: pass
if not files:
    files = ["outputs/model_column/iso_t22k_ns2000_gamma105_bestwb_condon.txt",
             "outputs/model_column/iso_bestwb_condon_ns3000.txt",
             "outputs/model_column/iso_bestwb_condon_ns4000.txt"]

def frame_near(fn, target):
    f=open(fn); ns,neq=map(int,f.readline().split())
    h=np.array(list(map(float,f.readline().split())))
    best=None; bestdt=1e30; cur=None; t=None
    for line in f:
        if 't =' in line:
            if cur is not None and len(cur)==ns:
                U=np.array(cur)
                if U.shape==(ns,neq) and np.all(np.isfinite(U)) and abs(t-target)<bestdt:
                    best=(t,U); bestdt=abs(t-target)
            cur=[]; t=float(line.split('t =')[1].split('step')[0])
        elif cur is not None:
            v=line.split()
            if len(v)==neq: cur.append(list(map(float,v)))
    if cur is not None and len(cur)==ns:
        U=np.array(cur)
        if U.shape==(ns,neq) and np.all(np.isfinite(U)) and abs(t-target)<bestdt: best=(t,U)
    return h, best

ds_l, v250_l, rvpk_l, ns_l, tt = [], [], [], [], None
for fn in files:
    h, fr = frame_near(fn, T)
    if fr is None: continue
    t, U = fr; tt = t
    V=U[:,2]/U[:,0]; rhoV=(U[:,0]+U[:,1])*V
    m=(h>=40)&(h<=700)
    ns=U.shape[0]; ds=(h[-1]-h[0])/(ns-1)
    ds_l.append(ds); ns_l.append(ns)
    v250_l.append(float(np.interp(250,h,V)))
    rvpk_l.append(float(np.max(rhoV[m])))

ds=np.array(ds_l); v250=np.array(v250_l); rvpk=np.array(rvpk_l)
order=np.argsort(ds); ds,v250,rvpk,ns_l=ds[order],v250[order],rvpk[order],[ns_l[i] for i in order]

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(11,4.4))
# linear fit V = A + B ds (least squares)
Bv,Av=np.polyfit(ds,v250,1); Br,Ar=np.polyfit(ds,rvpk,1)
xg=np.linspace(0,ds.max()*1.05,50)
ax1.plot(ds,v250,"o-",ms=8,color="C0")
ax1.plot(xg,Av+Bv*xg,"--",color="C0",alpha=.6,label=f"fit: {Av:+.2f}+{Bv:.1f}·Δs")
ax1.axhline(0,color="k",lw=.7)
ax1.scatter([0],[Av],marker="*",s=180,color="crimson",zorder=5,label=f"Δs→0: {Av:+.2f} m/s")
for x,y,n in zip(ds,v250,ns_l): ax1.annotate(f"ns={n}",(x,y),textcoords="offset points",xytext=(6,6),fontsize=8)
ax1.set_xlabel("Δs [km]"); ax1.set_ylabel("V @ 250 km [m/s]")
ax1.set_title("V@250 km ∝ Δs → 0"); ax1.legend(fontsize=8); ax1.set_xlim(left=0)

ax2.plot(ds,rvpk*1e4,"s-",ms=8,color="C2")
ax2.plot(xg,(Ar+Br*xg)*1e4,"--",color="C2",alpha=.6,label=f"fit: {Ar*1e4:+.2f}+{Br*1e4:.1f}·Δs")
ax2.axhline(0,color="k",lw=.7)
ax2.scatter([0],[Ar*1e4],marker="*",s=180,color="crimson",zorder=5,label=f"Δs→0: {Ar*1e4:+.2f}")
for x,y,n in zip(ds,rvpk,ns_l): ax2.annotate(f"ns={n}",(x,y*1e4),textcoords="offset points",xytext=(6,6),fontsize=8)
ax2.set_xlabel("Δs [km]"); ax2.set_ylabel("peak ρV, h∈[40,700] km [10⁻⁴ kg m⁻² s⁻¹]")
ax2.set_title("lower-chromosphere ρV peak ∝ Δs → 0"); ax2.legend(fontsize=8); ax2.set_xlim(left=0)

fig.suptitle(f"iso_t22k bestwb conduction-ON: ~250 km maximum is NUMERICAL (first-order, t≈{tt:.0f} s)",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.96])
os.makedirs("visualization",exist_ok=True)
fig.savefig(out,dpi=130)
print("wrote",out,"| Δs=",list(np.round(ds,3)),"V@250=",list(np.round(v250,2)),"A_V=",round(Av,3))
