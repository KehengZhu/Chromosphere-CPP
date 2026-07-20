import sys
import numpy as np
sys.path.insert(0, 'util')
import animate_isentropic as anim


def test_analyze_bestwb_massflux():
    files = {
        'condon': 'outputs/model_column/iso_t22k_ns2000_gamma105_bestwb_condon.txt',
        'condoff': 'outputs/model_column/iso_t22k_ns2000_gamma105_bestwb_condoff.txt',
    }
    for label, fn in files.items():
        h, frames = anim.load(fn)
        print(f'\n[{label}] ns={h.size} frames={len(frames)} t=[{frames[0][0]:.6g},{frames[-1][0]:.6g}] dh={h[1]-h[0]:.9g} km')
        ids = sorted(set(i for i in [0,1,2,3,4,5,10,20,50,100,200,400,800,len(frames)-1] if i < len(frames)))
        for k in ids:
            t,U = frames[k]
            rho = U[:,0] + U[:,1]
            V = U[:,2] / U[:,0]
            mf = rho * V
            base = h <= 400.0
            y = mf[base]
            d2 = np.diff(y,2)
            odd_even = abs(np.mean(y[::2]) - np.mean(y[1::2]))
            imax = np.argmax(np.abs(y))
            print(f't={t:10.4f} min={y.min(): .6e} max={y.max(): .6e} maxabs_h={h[base][imax]:.3f}km V={V[base][imax]: .6e} d2rms={np.sqrt(np.mean(d2*d2)): .6e} oddeven={odd_even: .6e}')
        t,U = frames[-1]
        rho = U[:,0] + U[:,1]
        V = U[:,2] / U[:,0]
        mf = rho*V
        print('final first 40 cells: i h_km rho V mflux')
        for i in range(40):
            print(i, f'{h[i]:.6f}', f'{rho[i]:.9e}', f'{V[i]:.9e}', f'{mf[i]:.9e}')
