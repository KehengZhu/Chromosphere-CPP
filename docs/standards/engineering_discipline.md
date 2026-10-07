# Engineering discipline

**Status: CURRENT** — moved verbatim from the root `CLAUDE.md` on 2026-10-07; `CLAUDE.md` keeps the summary.



Treat every task as a production-quality engineering task. The default priority order is:

1. Correct and reliable behavior.
2. Scientifically and numerically accurate results.
3. Reproducible evidence that supports the decision being made.
4. Simple, maintainable implementation.
5. Performance and convenience improvements.

Do not trade correctness for speed, but also do not confuse more code, more abstraction, or more validation with higher quality.

### Avoid overengineering

- Prefer the smallest change that clearly solves the actual problem.
- Use focused helper functions instead of introducing a framework, class hierarchy, policy system, or generalized subsystem without a demonstrated need.
- Do not refactor unrelated code while implementing a narrow task.
- Do not design for hypothetical future requirements unless the current task explicitly depends on them.
- Do not add layers of wrappers, configuration objects, templates, or indirection when a direct implementation is easier to review and equally safe.
- Preserve existing interfaces and behavior unless changing them is necessary for correctness or explicitly requested.

### Validate according to risk

Verification should answer material questions about correctness, reliability, numerical behavior, or performance. Do enough to establish trustworthy evidence, then stop.

- Prioritize tests that exercise the changed behavior and realistic production paths.
- Re-run broad test suites after changes that could affect them, but do not repeatedly rerun unchanged checks without a concrete reason.
- Use hashes and byte comparisons only when exact identity is itself important, such as regression preservation, deterministic output, or legacy compatibility.
- Do not generate or compare hashes for every intermediate artifact when the result cannot affect the decision.
- Do not require byte identity where numerical or physical equivalence is the correct standard.
- Use absolute errors, relative errors with defensible denominator floors, and physically meaningful observables where appropriate.
- Treat small timing variation, harmless formatting changes, and differences far below discretization or physical uncertainty as non-blocking unless they indicate a real bug.

### Focus on material issues

Investigate a discrepancy when it is scientifically meaningful, numerically suspicious, persistent, or large enough to change the production decision. Do not spend substantial effort on cosmetic details or artificial edge cases that cannot affect real use.

A newly discovered issue is a blocker only if it can materially affect one or more of:

- correctness;
- scientific interpretation;
- numerical stability or convergence;
- data integrity;
- requested output or termination behavior;
- production performance;
- reproducibility;
- the final GO/NO_GO decision.

Record unrelated minor issues briefly and leave them for a separate task.

### Use staged validation

For expensive experiments, simulations, or parameter sweeps, use a funnel:

1. Run inexpensive smoke and stability checks.
2. Reject clearly invalid candidates early.
3. Perform detailed comparisons only for plausible candidates.
4. Run long or expensive acceptance cases only for the selected candidate or the smallest set needed to resolve uncertainty.

Do not run several costly production-scale cases when one controlled run is sufficient.

### Protect existing work

- Inspect repository status and existing diffs before editing.
- Preserve unrelated uncommitted changes.
- Do not reset, revert, overwrite, or broadly reformat work outside the task scope.
- Keep patches narrow, reviewable, and attributable to the requested task.
- Keep generated outputs and build products out of source patches unless they are explicitly required deliverables.

### Stop when the decision is supported

Once the requested behavior is correct, the important tests pass, the evidence supports the conclusion, and material risks are documented, stop. Do not continue polishing, benchmarking, or redesigning merely because additional work is possible.

If a target is not achieved, report the measured gap and the most important remaining blocker. Do not silently expand into a new optimization or redesign phase.

### Reporting standard

Final recaps should be compact and decision-oriented. Emphasize:

- what changed;
- why it was necessary;
- what was validated;
- the material results;
- remaining production risks;
- the final status.

Avoid exhaustive command transcripts, repeated hashes, low-value implementation trivia, and long lists of insignificant differences.

The governing principle is:

```text
Reliable and accurate production behavior is mandatory.
Simple and sufficient evidence is preferred.
More validation is not automatically better validation.
```

