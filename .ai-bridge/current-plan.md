# N1000 4000s ripple resolution test

Updated: 2026-08-07T17:46:25.646Z
Workspace: /Users/zkeheng/SWMFSoftware/Chromosphere2026
Target agent: Codex (codex)

## Plan

Work in /Users/zkeheng/SWMFSoftware/Chromosphere2026. Execute the focused experiment from the user's attached task. Do not redesign the solver, reopen the full convergence study, add numerical conduction, change production defaults, or commit. Preserve all unrelated uncommitted changes.

Mission: determine whether the visible upper-chromosphere/below-TR local velocity ripple in the canonical model_column N500/R4 4000 s release is primarily a coarse-grid artifact by running exactly one N1000/R4 long case to 4000 s.

Verified current release contract: selecting model_column supplies override-preserving defaults GAMMA_TABLE=data/eos/gamma1_hydrogen_v1.dat, ISO_H_BASE=1600, ISO_DH=553, ISO_NS=500, ISO_HEAT_FLUX=1, ISO_T_TOP=22000, ISO_HYDRO_T_DECOUPLE=1, ISO_REFINE_PROFILE=outer, ISO_REFINE_FACTOR=4, ISO_REFINE_S_LO_KM=500, ISO_REFINE_TRANSITION_KM=20, ISO_NUMERICAL_DIFFUSIVITY_MULT=0, cooling/TRAC/corona/two-fluid/ionization off. scripts/run_chromo_realtime.sh adds 12-thread OpenMP, CHROMO_CFL=0.50, low-I/O policy. Therefore override only ISO_NS=1000 plus requested runtime/output diagnostics.

Run, adjusting only if the actual CLI requires it:
ISO_NS=1000 CHROMO_T_END=4000 CHROMO_OUTPUT=1 CHROMO_GAMMA_DIAG=1 CHROMO_FACE_FLUX_DIAG=1 CHROMO_OUTER_COND_DIAG=1 CHROMO_FRAME_DT=8 scripts/run_chromo_realtime.sh outputs/model_column/model_column_release_N1000_4000s.txt full no-ionization model_column - 20.0 no-cooling
Capture console to outputs/model_column/model_column_release_N1000_4000s.console.log without changing simulation semantics. Required sidecars: .gamma_diag, .faceflux, .outercond.

Acceptance: termination=end_time at 4000 s; no NaN/Inf, positivity failure, numerical conduction, EOS debug clamps, floor/vacuum rescue, or fallback. Report actual cell count, finest/coarsest spacing, final steps, representative dt, wall time, q_num/q_total, column-mass drift.

Compare existing outputs/model_column/model_column_release_4000s.txt (N500) and new N1000 at common physical times 500, 1000, 2000, 4000 s, fixed window 2130-2150 km. Reuse existing diagnostic tooling where practical; do not build a new framework. For V report mean/std/min/max, normalized high-frequency roughness Q_V=mean|D2V|/mean|V| (or exact established equivalent), and whether/local where negative-velocity reversals remain. For M=rho*V report mean, relative std, roughness. For production face flux use F_eff=F_total-F_ref, not raw F_total: report mean(F_eff), A_eff, Q_eff. Compare enough thermal structure to ensure bulk physics is unchanged: T_top, q_phys at outer face, integrated q_phys dt, thermal-layer thickness/cells, p_top, mean/peak V. Also report q_num/q_total and mass drift.

Use the existing N500 release run, not a rerun. Existing 500 s evidence for reference: N500 661 cells, dx_min~0.276 km; N1000 1320 cells, dx_min~0.138 km; at ~500 s Q_M 0.03063 vs 0.00570, A_eff 0.06984 vs 0.06120, Q_eff 7.85e-4 vs 3.65e-4, mean V 9.2408 vs 8.9384 m/s, q_phys 0.6588 vs 0.6839 W/m2, integrated heat 376.94 vs 384.69 J/m2. The decisive question is late time.

Interpretation:
A) If N1000 substantially removes late V ripple/sign reversals while conserved flux and thermal solution stay close: recommend N1000 for production/local dynamics; N500 remains okay for bulk energetics/mass transport.
B) If N1000 only modestly improves and clear cell-scale ripple/reversals remain: keep N500 for bulk physics and recommend separate hydro-flux/low-Mach/well-balanced investigation; do NOT run N2000 or alter solver.
C) If N1000 materially changes bulk evaporation/conductive/pressure/thermal-layer solution: N500 is not sufficiently converged even for long-duration bulk physics.

Optional only if cheap after successful run: generate comparable animation at visualization/model_column/model_column_release_N1000_4000s_evolution.mp4 using existing tooling and same axes/conventions as N500; if generated, update visualize_commands.md per AGENTS.md. Do not spend material effort on visualization before the scientific comparison.

Final report must include exact command, clean completion, runtime/steps/cells, metrics table at 500/1000/2000/4000, sign reversals, rhoV ripple, F_eff smoothness, conductive/thermal differences, whether N500 is acceptable for local velocity dynamics, and conclusion A/B/C. Do not commit.

## Implementation contract

- Work from this plan in small, reviewable steps.
- Keep edits scoped to the requested task and existing project conventions.
- Run focused verification before handing work back.
- Update .ai-bridge/agent-status.md with files touched, checks run, results, blockers, and review notes.
- Save the final review diff to .ai-bridge/implementation-diff.patch when practical.
- Append notable execution events to .ai-bridge/execution-log.jsonl when the implementation agent supports logging.
