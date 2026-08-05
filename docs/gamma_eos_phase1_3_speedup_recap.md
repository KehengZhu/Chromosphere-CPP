# Gamma-table solver speedup: Phases 1--3 recap

Date: 2026-07-21

## Result

Phases 1--3 are implemented without changing the CFL, spatial discretization,
EOS model, source ordering, or activation boundary. The optimized path remains
the signed-off classical pure-H Saha/CRASH-Gamma1 `model_column` configuration.

The previous 2,000-cell, 1,000 s production run required 4,417 s wall time for
85,892 steps (51.43 ms/step). Three repeated output-free Release benchmarks of
the optimized solver advance the same configuration for 10 physical seconds in
10.57, 10.61, and 10.72 s. The median is 10.61 s for 877 steps, or
12.10 ms/step: a measured **4.25x per-step speedup**.

Applying that per-step cost to the previous 85,892-step trajectory gives an
indicative 17.3-minute 1,000 s runtime. This is an extrapolation, not a newly
completed 1,000 s run.

The final profiled 877-step run reports:

```text
wall time                         10.63 s
decoded-field time                 3.616 s
gamma RHS time                     2.984 s
projection time                    2.587 s
conduction time                    3.209 s
boundary time                      0.097 s
EOS inversion calls            7,032,534
accepted exact guesses           684,630
average / maximum iterations   29.68 / 57
conduction Newton updates       4.236 average, 5 maximum
active timestep limiter          acoustic: 877 / 877
```

Decoded-field time includes the predictor decode inside the RHS, so region
timings intentionally overlap and must not be summed as exclusive percentages.

## Phase 1: profiler, benchmark controls, and inversion fast path

`CHROMO_PROFILE=1` now reports:

- wall time for CFL, decoded-field construction, gamma RHS, projection,
  nonlinear conduction, boundary updates, and output;
- inversion calls, exact-guess accepts, one-update convergence, average and
  maximum iterations, bisection fallbacks, and bracket evaluations;
- conduction Newton calls and update counts;
- the final active timestep limiter.

Counters are ordinary serial counters and are completely bypassed when
profiling is disabled. A valid temperature guess is checked against the
existing `2e-13` caloric residual tolerance before evaluating the two table
bracket endpoints. Failed or invalid guesses fall through to the unchanged
safeguarded Newton/bisection solver.

`CHROMO_OUTPUT=0` and `CHROMO_GAMMA_DIAG=0` independently disable state and
sidecar output for timing runs. No files are created when both are zero.

## Phase 2: immutable decoded fields and conservative known-T packing

Each gamma step constructs a state-bound `DecodedMixtureField` after the
boundary update. It contains all physical cells and two explicitly decoded
ghost layers. Physical cells use cell-centered potential, while both inner and
outer ghost layers use their boundary-face potentials. Center and `i+/-1,2`
stencils index this extended field instead of reinverting shifted conserved
copies.

The same immutable `U^n` field is consumed by CFL and gamma RHS. The internal
MUSCL predictor receives a different cache generation and uses the corresponding
center temperatures only as optional inversion guesses. State pointer, memory,
size, Grid identity, and all four ghost-buffer signatures are checked before a
cache can be consumed.

The nonlinear conduction solve now reuses its converged temperature while
retaining the flux-derived target total energy as authoritative. The charged
energy row is constructed analytically and the neutral row is the conservative
remainder after the charged row is rounded. The known temperature is accepted
only when it satisfies a redundant post-solve residual assertion. The
conduction solve normalizes by the old cell energy, whereas this assertion
normalizes by the flux-updated target; the assertion therefore uses twice the
`2e-11` solve tolerance so a converged cooling cell is not rejected solely
because the two denominators differ.

For the 877-step profile, the optimized path performs approximately four
temperature inversions per physical cell-step. The pre-optimization source audit
found approximately fourteen, so shifted-state caching and known-temperature
packing remove about 71% of caloric inversions.

## Phase 3: reusable scratch storage

`GammaConductionScratch` is owned by each `Grid` and reuses every nonlinear
array: density, old/target energy, temperature, composition, conductivity,
heat capacity, face coefficients, tridiagonal coefficients, residual, and
Newton update. The double-precision Thomas solver now operates in place and
does not copy or allocate its diagonal/RHS/solution inside Newton iterations.

`GammaRhsScratch` similarly reuses decoded predictor storage, mixture stencil
arrays, limiter fields, reconstructed states, face states/fluxes, predictor,
and final flux/source buffers. There is no function-static scratch, so separate
Grid instances remain isolated.

## Numerical comparison

An output-enabled 10 s run was compared with matching frames from the earlier
1,000 s PCHIP production artifact at steps 0, 200, 400, 600, and 800.
Initialization is identical. At step 800 (t approximately 9.11 s):

```text
maximum |Delta T|          0.74 K
maximum |Delta p|          0.04 Pa
maximum |Delta v|          0.2924 m/s  (peak |v| approximately 308 m/s)
largest normalized Linf    9.50e-4
largest relative L2        1.52e-3
```

The small accumulated difference is expected: the old conduction pack replaced
the flux-derived target by `e_EOS(rho,T)` after Newton convergence, whereas the
new path preserves the target energy and leaves the accepted nonlinear residual
only in the carrier-row partition. This is the conservation-preserving behavior
required by the reviewed plan.

## Regression coverage

New tests verify:

1. direct-versus-cached center and `i+/-1,2` thermodynamics;
2. both inner and outer ghost layers at their carried potentials;
3. rejection after state replacement or boundary-buffer mutation;
4. no more than one predicted physical-cell inversion plus four ghosts in the
   cached RHS;
5. exact-guess acceptance without bracket evaluation;
6. known-temperature mass, momentum, and authoritative-energy conservation;
7. rejection of a supplied temperature outside the nonlinear residual gate.

Current direct unit result: 7,206 checks, zero failures. The permanent gate is
all registered tests, not this fixed count.

## Reproduction commands

```bash
cmake -S . -B build
cmake --build build -j4
ctest --test-dir build --output-on-failure

GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat \
ISO_H_BASE=0 ISO_DH=2152.6 ISO_T_TOP=22000 \
ISO_HEAT_FLUX=1 ISO_NS=2000 CHROMO_T_END=10 \
CHROMO_OUTPUT=0 CHROMO_GAMMA_DIAG=0 CHROMO_PROFILE=1 \
./build/chromo_main /tmp/chromo_phase13_profile.txt \
full no-ionization model_column - 9.0 no-cooling
```

## Remaining measured bottlenecks

The optimized serial profile still spends most of its time in gamma RHS,
projection, conduction, and the top-level current-state decode. Temperature
inversions average many safeguarded updates in ionization-transition cells, so
further serial work should focus on the inversion algorithm only after a
separate correctness review. Parallelization and any CFL experiment remain
outside Phases 1--3.
