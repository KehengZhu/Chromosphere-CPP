#!/usr/bin/env python3
"""Run an ensemble of field lines in parallel — one chromo_main process per line.

WHY PROCESS-LEVEL (not threads inside the solver): each field line is a completely
independent 1-D run, and chromo_main is built without BLAS/LAPACK (CMakeLists.txt:
ARMA_DONT_USE_LAPACK/BLAS) so a single run is effectively single-threaded — its
hot loop (RHS/flux/integrator) is serial; the only OpenMP regions (src/state.cpp)
cover the cheap state conversion and give negligible intra-run speedup. The cheap,
robust, near-linear win is therefore to run MANY single-threaded processes at once,
one per core. We pin every child to ONE thread (OMP_NUM_THREADS=1 etc.) so that N
concurrent processes don't oversubscribe the N cores.

This is the embarrassingly-parallel ensemble driver (campaign Phase P4). It does not
touch the C++; it schedules chromo_main invocations across a process pool.

Two ways to specify the jobs:

  (A) Manifest (full per-line control — the production path; per-line beam params
      come from STIX + the AIA ribbon mask, see docs/event-data-to-code.md §5):

        python util/run_ensemble.py --manifest outputs/events/event_ensemble_jobs.json [--workers 16]

      jobs.json:
        {
          "binary": "build/chromo_main",
          "defaults": {
            "mode": "full", "ioniz": "ionization", "scenario": "pfss_field_line",
            "time_mult": 1.0, "cooling": "cooling",
            "env": {"PFSS_FLARE": "1", "FLARE_DELTA": "4", "FLARE_E_CUT": "20"}
          },
          "jobs": [
            {"name": "lineA", "dat": "scenarios/data/lineA.dat",
             "out": "outputs/pfss/lineA.txt",
             "env": {"FLARE_BEAM_FLUX": "3e7", "FLARE_T_ON": "2"}},
            ...
          ]
        }
      Per-job env overrides defaults.env; per-job mode/time_mult/... override defaults.

  (B) Quick glob mode (shared env for every line):

        python util/run_ensemble.py \
          --dats "scenarios/data/ensemble_demo/*.dat" --out-dir outputs/events/ensemble_demo \
          --time-mult 1.0 --env PFSS_FLARE=1 FLARE_BEAM_FLUX=3e7 FLARE_DELTA=4 FLARE_E_CUT=20

Each job writes <out>.txt and a sibling <out>.run.log; the driver prints a per-job
table and the overall wall-clock vs. the summed single-run time (the speedup).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# Force every child to a single thread so N concurrent processes map cleanly onto
# N cores (no thread oversubscription). chromo_main is single-threaded in practice,
# but be defensive about any OpenMP/BLAS the toolchain might enable.
THREAD_PIN = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",   # Apple Accelerate
    "MKL_NUM_THREADS": "1",
}


def build_jobs_from_manifest(path):
    with open(path) as f:
        spec = json.load(f)
    binary = spec.get("binary", "build/chromo_main")
    d = spec.get("defaults", {})
    d_env = d.get("env", {})
    jobs = []
    for j in spec["jobs"]:
        env = dict(d_env)
        env.update(j.get("env", {}))
        jobs.append({
            "name": j["name"],
            "binary": j.get("binary", binary),
            "dat": j["dat"],
            "out": j["out"],
            "mode": j.get("mode", d.get("mode", "full")),
            "ioniz": j.get("ioniz", d.get("ioniz", "ionization")),
            "scenario": j.get("scenario", d.get("scenario", "pfss_field_line")),
            "time_mult": j.get("time_mult", d.get("time_mult", 1.0)),
            "cooling": j.get("cooling", d.get("cooling", "cooling")),
            "env": env,
        })
    return jobs


def build_jobs_from_glob(args):
    dats = sorted(glob.glob(args.dats))
    if not dats:
        sys.exit(f"no .dat files match {args.dats!r}")
    env = {}
    for kv in args.env or []:
        if "=" not in kv:
            sys.exit(f"--env entries must be KEY=VALUE, got {kv!r}")
        k, v = kv.split("=", 1)
        env[k] = v
    os.makedirs(args.out_dir, exist_ok=True)
    jobs = []
    for dat in dats:
        name = os.path.splitext(os.path.basename(dat))[0]
        jobs.append({
            "name": name,
            "binary": args.binary,
            "dat": dat,
            "out": os.path.join(args.out_dir, f"{name}.txt"),
            "mode": args.mode,
            "ioniz": args.ioniz,
            "scenario": args.scenario,
            "time_mult": args.time_mult,
            "cooling": args.cooling,
            "env": dict(env),
        })
    return jobs


def run_one(job):
    """Run a single chromo_main process. Returns a result dict."""
    os.makedirs(os.path.dirname(job["out"]) or ".", exist_ok=True)
    log_path = os.path.splitext(job["out"])[0] + ".run.log"
    cmd = [
        job["binary"], job["out"], job["mode"], job["ioniz"],
        job["scenario"], job["dat"], str(job["time_mult"]), job["cooling"],
    ]
    env = dict(os.environ)
    env.update(THREAD_PIN)
    env.update({k: str(v) for k, v in job["env"].items()})
    t0 = time.perf_counter()
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, env=env, stdout=logf, stderr=subprocess.STDOUT)
    dt = time.perf_counter() - t0
    return {
        "name": job["name"], "rc": proc.returncode, "wall": dt,
        "out": job["out"], "log": log_path,
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", help="JSON job manifest (mode A)")
    # Glob mode (B):
    ap.add_argument("--dats", help="glob of .dat files (mode B)")
    ap.add_argument("--out-dir", default="outputs/events/event_ensemble")
    ap.add_argument("--binary", default="build/chromo_main")
    ap.add_argument("--mode", default="full")
    ap.add_argument("--ioniz", default="ionization")
    ap.add_argument("--scenario", default="pfss_field_line")
    ap.add_argument("--cooling", default="cooling")
    ap.add_argument("--time-mult", type=float, default=1.0)
    ap.add_argument("--env", nargs="*", help="KEY=VALUE shared env for every line (mode B)")
    # Pool:
    ap.add_argument("--workers", type=int, default=os.cpu_count(),
                    help="max concurrent processes (default: all cores)")
    ap.add_argument("--dry-run", action="store_true", help="list jobs, don't run")
    args = ap.parse_args()

    if args.manifest:
        jobs = build_jobs_from_manifest(args.manifest)
    elif args.dats:
        jobs = build_jobs_from_glob(args)
    else:
        ap.error("specify either --manifest or --dats")

    workers = max(1, min(args.workers, len(jobs)))
    print(f"[ensemble] {len(jobs)} field lines, {workers} concurrent workers "
          f"({os.cpu_count()} cores), 1 thread/process")
    if args.dry_run:
        for j in jobs:
            envs = " ".join(f"{k}={v}" for k, v in j["env"].items())
            print(f"  {j['name']:24s} dat={j['dat']}  out={j['out']}  [{envs}]")
        return

    t_start = time.perf_counter()
    results = []
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(run_one, j): j for j in jobs}
        for fut in as_completed(futures):
            r = fut.result()
            results.append(r)
            done += 1
            status = "ok " if r["rc"] == 0 else f"FAIL(rc={r['rc']})"
            print(f"  [{done:3d}/{len(jobs)}] {status} {r['name']:24s} {r['wall']:6.1f}s")
    wall = time.perf_counter() - t_start

    results.sort(key=lambda r: r["name"])
    n_ok = sum(1 for r in results if r["rc"] == 0)
    sum_single = sum(r["wall"] for r in results)
    print(f"\n[ensemble] done: {n_ok}/{len(results)} succeeded")
    print(f"[ensemble] wall-clock        = {wall:7.1f} s")
    print(f"[ensemble] summed single-run = {sum_single:7.1f} s  "
          f"(serial would take this long)")
    if wall > 0:
        print(f"[ensemble] parallel speedup  = {sum_single / wall:7.2f}x "
              f"on {workers} workers")
    failed = [r["name"] for r in results if r["rc"] != 0]
    if failed:
        print(f"[ensemble] FAILED lines: {', '.join(failed)} (see their .run.log)")
        sys.exit(1)


if __name__ == "__main__":
    main()
