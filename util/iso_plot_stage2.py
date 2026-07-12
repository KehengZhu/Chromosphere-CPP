import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
mp=1.6726219e-27;kb=1.380649e-23;g=273.95
f=open("outputs/iso_stage2.txt");ns,neq=map(int,f.readline().split());h=np.array(list(map(float,f.readline().split())))
frames=[];cur=None;t=None
def flush():
    if cur is None:return
    rows=[r for r in cur if len(r)==neq]
    if len(rows)==ns:frames.append((t,np.array(rows)))
for line in f:
    if line.startswith('#'):flush();t=float(line.split('t =')[1].split('step')[0]);cur=[]
    else:
        try:cur.append(list(map(float,line.split())))
        except:cur.append([])
flush()
phi=g*(h-h[0])*1000
def prim(U):
    rho_i,rho_n=U[:,0],U[:,1];V=U[:,2]/rho_i
    p_i=2/3*U[:,4]-1/3*rho_i*V*V-2/3*rho_i*phi
    p_n=2/3*U[:,5]-1/3*rho_n*(U[:,3]/rho_n)**2-2/3*rho_n*phi
    rho=rho_i+rho_n;T=(p_i+p_n)/(rho*((1+1e-3)*kb/mp));return V,T,rho
times=[f[0] for f in frames]
def near(tt): return frames[int(np.argmin(np.abs(np.array(times)-tt)))]
fig,ax=plt.subplots(1,3,figsize=(15,4.5))
sel=[0.0,10,40,150,600,3000]
import matplotlib.cm as cm
cols=cm.viridis(np.linspace(0,0.92,len(sel)))
for tt,c in zip(sel,cols):
    t,U=near(tt);V,T,rho=prim(U)
    ax[0].plot(h,T,c=c,label=f't={t:.0f}s')
    ax[1].plot(h,V/1e3,c=c)
    ax[2].semilogy(h,rho,c=c)
ax[0].set_xlabel('height [km]');ax[0].set_ylabel('T [K]');ax[0].set_title('T: top heats (conductive flux from hot ghost)');ax[0].legend(fontsize=8)
ax[1].axhline(0,c='k',lw=0.5);ax[1].set_xlabel('height [km]');ax[1].set_ylabel('V [km/s]');ax[1].set_title('V: upflow + downflow front (evaporation)')
ax[2].set_xlabel('height [km]');ax[2].set_ylabel(r'$\rho$ [kg/m³]');ax[2].set_title('density')
plt.tight_layout();plt.savefig("visualization/iso_stage2_evaporation.png",dpi=130)
# report the peak transient
for tt in [2,5,10,20,40,80]:
    t,U=near(tt);V,T,rho=prim(U)
    iup=np.argmax(V);idn=np.argmin(V)
    print(f"t={t:6.1f}  Vmax(up)={V[iup]/1e3:6.2f}km/s @h={h[iup]:.0f}  Vmin(down)={V[idn]/1e3:6.2f}km/s @h={h[idn]:.0f}  Ttop={T[-5]:.0f}K")
print("saved visualization/iso_stage2_evaporation.png")
