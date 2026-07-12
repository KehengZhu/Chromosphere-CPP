import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
mp=1.6726219e-27;kb=1.380649e-23;g=273.95
def load(fn):
    f=open(fn);ns,neq=map(int,f.readline().split());h=np.array(list(map(float,f.readline().split())))
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
    flush();return h,frames
h,frames=load("outputs/iso_stage1.txt")
phi=g*(h-h[0])*1000
def prim(U):
    rho_i,rho_n=U[:,0],U[:,1];V=U[:,2]/rho_i;Un=U[:,3]/rho_n
    p_i=2/3*U[:,4]-1/3*rho_i*V*V-2/3*rho_i*phi
    p_n=2/3*U[:,5]-1/3*rho_n*Un*Un-2/3*rho_n*phi
    rho=rho_i+rho_n; T=(p_i+p_n)/(rho*((1+1e-3)*kb/mp))
    return V,T,rho,(p_i+p_n)
t0,U0=frames[0]; tf,Uf=frames[-1]
V0,T0,r0,p0=prim(U0); Vf,Tf,rf,pf=prim(Uf)
# linear fit of relaxed T(h)
A=np.polyfit(h,Tf,1); Tfit=np.polyval(A,h); R2=1-np.sum((Tf-Tfit)**2)/np.sum((Tf-Tf.mean())**2)
fig,ax=plt.subplots(1,3,figsize=(14,4.2))
ax[0].plot(h,T0,'--',c='gray',label=f'IC (t=0)')
ax[0].plot(h,Tf,'-',c='C3',label=f'relaxed (t={tf:.0f}s)')
ax[0].plot(h,Tfit,':',c='k',lw=1,label=f'linear fit R²={R2:.5f}')
ax[0].set_xlabel('height [km]');ax[0].set_ylabel('T [K]');ax[0].set_title('Temperature: linear (isentropic)');ax[0].legend(fontsize=8)
ax[1].plot(h,V0,'--',c='gray',label='IC');ax[1].plot(h,Vf,'-',c='C0',label=f'relaxed: max|V|={np.max(np.abs(Vf)):.1f} m/s')
ax[1].axhline(0,c='k',lw=0.5);ax[1].set_xlabel('height [km]');ax[1].set_ylabel('V [m/s]');ax[1].set_title('Velocity ≈ 0 (tens of m/s)');ax[1].legend(fontsize=8)
ax[2].semilogy(h,p0,'--',c='gray',label='IC');ax[2].semilogy(h,pf,'-',c='C2',label='relaxed')
ax[2].set_xlabel('height [km]');ax[2].set_ylabel('p [Pa]');ax[2].set_title('Pressure (hydrostatic, stable)');ax[2].legend(fontsize=8)
plt.tight_layout();plt.savefig("visualization/iso_stage1_relaxation.png",dpi=130)
print(f"relaxed: max|V|={np.max(np.abs(Vf)):.1f} m/s  T[{Tf.min():.0f},{Tf.max():.0f}]K  linear-fit R²={R2:.6f}  lapse={A[0]:.3f} K/km")
print("saved visualization/iso_stage1_relaxation.png")
