# Coarse `model_column` with physical-face conduction

## Status

The reduced Gamma/Saha `model_column` now uses the validated coarse physical-conduction configuration below as its canonical release model. It removes the dominant grid-dependent artificial thermal drive while using a substantially coarser mesh.

Canonical release configuration:

- domain: 1600--2153 km;
- coarse-equivalent `ISO_NS=500`;
- static outer refinement: R4, beginning at 500 km above the base with the existing 20 km transition;
- actual cells: 661;
- coarse spacing: about 1.105 km;
- finest spacing: about 0.276 km;
- upper hydro temperature decoupling retained;
- conductive boundary: fixed `T_face = 22000 K` at the physical outer face;
- physical conduction only: `ISO_NUMERICAL_DIFFUSIVITY_MULT=0` (the production default);
- production runtime CFL: 0.50.

This result supersedes the old production-shaped N=2000/R4 configuration for this reduced conduction-driven `model_column`. It does not reinterpret historical runs that used the old ghost-center thermal boundary or mesh-scaled numerical conduction.

The `model_column` scenario itself now applies this validated shape as override-preserving defaults, including `data/eos/gamma1_hydrogen_v1.dat` when neither `GAMMA_TABLE` nor an explicit diagnostic `ISO_GAMMA` is supplied. `scripts/run_chromo_realtime.sh` is runtime-only: it adds the validated CFL=0.50, OpenMP, and low-I/O defaults without carrying a second copy of the model physics. The generic `scripts/run_chromo_omp.sh` remains unchanged.

## Numerical changes

### Physical-face outer conduction boundary

The imposed 22,000 K temperature is now located at the physical domain face at 2153 km. The live top cell remains free to evolve.

For the top cell,

```text
d_boundary = 0.5 * ds_top
q_phys = kappa_face * (T_face - T_top) / d_boundary
T_face = 22000 K
```

The boundary-face conductivity is evaluated from a physical face state constructed using the fixed external hydro pressure and `T_face`, with the hydrogen ionization fraction evaluated consistently from Saha equilibrium. It is not evaluated at a fictitious reflected ghost-center temperature.

The nonlinear backward-Euler conduction residual and the outer-conduction diagnostic use the same half-cell geometry.

### Hydro ghosts remain separate from conduction

The hydro boundary was not redesigned:

- first outer hydro ghost temperature remains zero-gradient from the live top cell;
- the fixed external hydro back pressure is retained;
- the second outer hydro ghost retains the existing HSE reconstruction role required by MUSCL;
- the second hydro ghost has no physical role in the conductive boundary.

Thus the hydro reconstruction still has two ghosts, while the conductive Dirichlet datum exists only at the physical face.

### Numerical conduction removed from production physics

The active `model_column` production baseline now has

```text
q_total = q_physical
```

with no mesh-scaled artificial conductivity. `ISO_NUMERICAL_DIFFUSIVITY_MULT` defaults to zero. Setting it explicitly to one remains available only for controlled historical/diagnostic comparisons with the former `C_num=2000 m/s` behavior.

### Finite-volume initial temperature

The Gamma/Saha initial temperature is represented by a cell average of the prescribed C7/PCHIP profile instead of a cell-center point sample. A two-point Gauss-Legendre rule is used inside each finite-volume cell. The existing Saha-consistent hydrostatic construction and equilibrium-reference well balancing are otherwise unchanged.

This removes the cheap O(ds) confound in which refinement changed the upper initial atmosphere simply because the top cell center moved.

## Physical-conduction-only positivity repair

An initial N=500/R4, CFL=0.50 long run failed near t = 705.46 s with `rho_n must be positive and finite`. Investigation localized the failure to the explicit hydro predictor immediately before the single-fluid equilibrium projection.

The predictor can produce a very small finite undershoot in a trace density row even while total density, total momentum, total energy, and internal energy remain valid. The equilibrium projection then immediately reconstructs both density rows from the conserved total density and the Saha equilibrium fraction. Requiring each predictor row to be positive before that reconstruction was therefore an inconsistent precondition.

The repair is deliberately narrow:

- the equilibrium projection accepts finite individual predictor density rows;
- total density must still be positive and finite;
- row momenta and energies must remain finite;
- projected internal energy and the EOS domain remain strictly validated;
- the ordinary EOS state decoders still require each stored density row to be positive;
- the projected state still enforces the existing positive trace-fraction floor.

No global artificial diffusion, new limiter, or generalized positivity framework was added.

After this repair, the final N=500/R4, CFL=0.50 physical-conduction-only run reached 1000 s cleanly, crossing both the earlier 705 s failure and the old physical-only failure time near 921 s.

## Verification

After the numerical changes and the targeted positivity repair:

```text
OpenMP build: success
CTest:        8/8 passed
```

The focused tests cover the physical-face boundary location, live top-cell evolution, hydro/conduction temperature decoupling, second hydro ghost independence, physical-conduction-only operation, finite-volume IC averaging, and recovery of a conservative predictor state with a tiny trace-row undershoot.

For the production Gamma/Saha configuration, beam heating is off and the vacuum-floor stage is not active. The EOS debug clamp remains off. The successful long run therefore does not rely on floor/clamp rescue.

## Compact resolution study

All cases used the same 1600--2153 km domain, R4 outer refinement, fixed 22,000 K physical-face boundary, physical conduction only, and the same Gamma/Saha physics.

### Early rejection of N=250

At 100 s:

| coarse-equivalent N | actual cells | finest spacing | mean conserved face mass flux | outer physical flux | thermal-layer thickness |
|---:|---:|---:|---:|---:|---:|
| 250 | 331 | 0.552 km | 1.062e-9 kg m^-2 s^-1 | 0.758 W m^-2 | 12.15 km |
| 500 | 661 | 0.276 km | 1.341e-9 kg m^-2 s^-1 | 0.892 W m^-2 | 11.32 km |
| 1000 | 1320 | 0.138 km | 1.431e-9 kg m^-2 s^-1 | 0.941 W m^-2 | 10.90 km |

N=250 is materially under-resolved for this purpose: relative to N=500 its mean conserved face mass flux is about 21% lower and its outer physical conductive flux about 15% lower. It was rejected without further expensive runs.

### N=500 versus N=1000 at 500 s

Both runs used CFL=0.50.

| quantity | N=500 | N=1000 | practical difference |
|---|---:|---:|---:|
| actual cells | 661 | 1320 | -- |
| mean conserved face mass flux | 6.255e-10 | 6.085e-10 kg m^-2 s^-1 | ~2.8% |
| outer physical conductive flux | 0.6588 | 0.6839 W m^-2 | ~3.7% |
| integrated physical conductive input | 376.94 | 384.69 J m^-2 | ~2.0% |
| top-cell temperature | 21862 | 21929 K | ~0.3% |
| thermal-layer thickness | 13.53 | 12.14 km | ~11% |
| mean velocity, 2130--2150 km | 9.24 | 8.94 m s^-1 | ~3.4% |
| max |velocity| | 43.00 | 40.88 m s^-1 | ~5.2% |

The thermal-layer thickness remains the largest measured spatial difference. However, the N=500 layer is about 13.5 km thick and is resolved by about 50 cells, so it is not a few-cell unresolved structure. The quantities most directly tied to the thermal drive and evaporation response are within roughly 2--5% of N=1000 at 500 s.

The cell-centered `rho V` diagnostic remains visibly more oscillatory on the coarser mesh. This is not used as the convergence decision variable: the conserved face mass flux is substantially smoother and is the physically relevant continuity diagnostic.

## CFL comparison for N=500

At 500 s, CFL=0.25 and CFL=0.50 were compared on the same N=500/R4 mesh.

The CFL=0.50 result differed from CFL=0.25 by approximately:

- mean conserved face mass flux: 1.6%;
- outer physical conductive flux: 2.0%;
- integrated conductive input: 2.7%;
- thermal-layer thickness: 2.1%.

The CFL=0.25 case had a column-mass drift of about -7.8e-5 at 500 s. No qualitative change in the evaporation response was observed. CFL=0.50 is therefore accepted for this selected coarse candidate.

The solver's conservative hard-coded default remains CFL=0.25; the 0.50 recommendation is a validated runtime choice for this production-shaped model.

## Final 1000 s production-candidate run

Configuration: N=500/R4, CFL=0.50, physical conduction only, 12 OpenMP threads.

```text
termination                  end_time
final physical time          1000 s
final step                   178054
representative dt            ~5.62e-3 s
actual cells                 661
min cell width               0.2759 km
max cell width               1.1050 km
top cell width               0.2761 km
wall time                    123.68 s
q_num / q_total              0
column-mass relative drift   -1.20e-5
```

At 1000 s:

- mean conserved face mass flux, 2130--2150 km: 3.85e-10 kg m^-2 s^-1;
- outer physical conductive flux: 0.572 W m^-2;
- integrated physical conductive input: 671.6 J m^-2;
- top-cell temperature: 21881 K;
- thermal-layer thickness: 13.53 km, about 50 N=500 cells;
- top pressure: 1.02714e-2 Pa;
- mean velocity in 2130--2150 km: 5.79 m/s;
- max |velocity|: 35.84 m/s;
- conserved face-flux roughness `Q_eff`: 5.1e-4.

The run completed without a positivity exception and without numerical conduction, vacuum floors, beam-heating floors, or EOS debug clamping.

## Cost reduction

The previous production-shaped reference used N=2000/R4 and 2638 actual cells. The new candidate uses 661 cells, almost exactly 4x fewer.

Using the documented old 2000 s run (1,395,700 steps) as the reference, a 1000 s equivalent requires about 697,850 steps. The new candidate requires 178,054 steps, about 3.9x fewer. A simple cells-times-steps work proxy therefore falls by about 15.6x.

The old 2000 s run reported 1826.5 s wall time. Normalized to 1000 s of simulated time, that is about 913 s versus 123.68 s for the new candidate, an indicative wall-time improvement of about 7.4x. This is not a controlled microbenchmark: output cadence and the old numerical-conduction work differ, so the 15.6x work proxy and 7.4x observed normalized wall improvement should be treated as approximate performance evidence rather than exact kernel speedups.

## Production decision

For the intended reduced conduction-driven evaporation model, the recommended current production candidate is:

```text
N = 500 coarse-equivalent
R4 outer refinement
661 actual cells
T_face(2153 km) = 22000 K
physical conduction only
cell-averaged C7/PCHIP IC
existing hydro boundary
existing MUSCL/Rusanov hydro
CFL = 0.50
```

This meets the practical objective: the dominant mesh-dependent artificial thermal drive has been removed, the upper thermal structure remains well resolved, the conserved evaporation response is close to the next finer useful reference, and the expensive N=2000 grid is no longer supported as scientifically necessary for this reduced model.

## Remaining limitations

- No claim of formal second-order convergence is made.
- N=500 thermal-layer thickness is about 11% larger than N=1000 at 500 s, although the layer still spans about 50 cells.
- The N=1000 comparison was carried to 500 s, not 1000 s; the selected N=500 candidate alone received the 1000 s long-run stability gate.
- The cell-centered `rho V` diagnostic still contains substantial reconstruction-scale ripple; the conserved face mass flux is much smoother and is the quantity used for the physical decision.
- The 1000 s result demonstrates stable conduction-driven evolution, not a mathematically limiting or fully time-asymptotic evaporation rate.
- The measured wall-time speedup is indicative rather than an apples-to-apples benchmark.
- The narrow projection positivity repair is specific to the equilibrium single-fluid representation; strict per-row positivity remains in ordinary EOS decoding.
