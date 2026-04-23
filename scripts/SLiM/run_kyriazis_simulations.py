#!/usr/bin/env python3
"""
Batch runner for the Kyriazis DFE split-population SLiM simulation.

Uses kyriazis_model.slim which combines:
- Kim et al. 2017 DFE (gamma with shape=0.186) + 0.3% recessive lethals
- Kyriazis et al. 2023 inverse h-s relationship
- Fish-specific parameters (lumpfish gene size and demographic history)
- 19-population split topology (17 base + 2 recovery) for load dynamics analysis

Usage:
    python run_kyriazis_simulations.py --replicates 20 --parallel 8
    python run_kyriazis_simulations.py --replicates 1 --scaling-factor 10 --dry-run
    python run_kyriazis_simulations.py --replicates 20 --parallel 8 --burnin-load
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
SLIM_SCRIPT = SCRIPT_DIR / "kyriazis_model.slim"
RESULTS_DIR = SCRIPT_DIR / "results_kyriazis"


def stable_seed(replicate: int, base_seed: int | None = None) -> int:
    """Generate a deterministic seed across runs and platforms."""
    seed_str = f"kyriazis_split:{replicate}"
    if base_seed is not None:
        seed_str = f"{base_seed}:{seed_str}"
    digest = hashlib.blake2b(seed_str.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "little") % (2**31)


def parse_slim_csv(stdout: str) -> str | None:
    """Extract CSV data between CSV_START and CSV_END markers from SLiM stdout."""
    lines = stdout.split("\n")
    csv_lines = []
    in_csv = False
    for line in lines:
        stripped = line.strip()
        if stripped == "CSV_START":
            in_csv = True
            continue
        elif stripped == "CSV_END":
            break
        if in_csv and stripped:
            csv_lines.append(stripped)
    if csv_lines:
        return "\n".join(csv_lines) + "\n"
    return None


def run_simulation(
    replicate: int,
    output_dir: Path,
    base_seed: int | None,
    burnin_load: bool,
    extra_defines: list[str] | None = None,
    scaling_factor: int = 1,
    n_chrom: int = 25,
) -> dict:
    """Run a single Kyriazis split-population SLiM simulation."""
    seed = stable_seed(replicate, base_seed)

    # stdbuf forces line-buffered stdout so progress lines appear immediately
    cmd = ["stdbuf", "-oL", "slim"]
    cmd.extend(["-d", f"SEED={seed}"])
    cmd.extend(["-d", f"REPLICATE={replicate}"])
    cmd.extend(["-d", f'OUTPUT_DIR="{output_dir}"'])
    cmd.extend(["-d", f"Q={scaling_factor}"])
    cmd.extend(["-d", f"N_CHROM={n_chrom}"])
    cmd.extend(["-d", f"BURNIN_LOAD={1 if burnin_load else 0}"])
    cmd.extend(["-d", f"BURNIN_SAVE={0 if burnin_load else 1}"])

    if extra_defines:
        for define in extra_defines:
            cmd.extend(["-d", define])

    cmd.append(str(SLIM_SCRIPT))

    output_dir.mkdir(parents=True, exist_ok=True)

    result = {
        "replicate": replicate,
        "seed": seed,
        "output_prefix": str(output_dir / f"split_rep{replicate}"),
        "success": False,
        "error": None,
        "cmd": " ".join(cmd),
    }

    try:
        # Start SLiM in its own process group so we can kill it cleanly
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            preexec_fn=os.setsid,
        )

        # Stream stdout: print progress lines, capture CSV data
        csv_lines = []
        in_csv = False
        progress_tags = ("Burn-in", "[Gen", "Timing", "Severity", "Recovery", "=== T_SNAPSHOT", "=== Recovery")

        try:
            for line in proc.stdout:
                stripped = line.strip()
                if stripped == "CSV_START":
                    in_csv = True
                    continue
                elif stripped == "CSV_END":
                    in_csv = False
                    continue
                if any(tag in stripped for tag in progress_tags):
                    print(f"[Rep {replicate}] {stripped}", flush=True)
                elif in_csv and stripped:
                    csv_lines.append(stripped)
        except KeyboardInterrupt:
            # Kill the entire process group (stdbuf + slim)
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            proc.wait()
            raise

        proc.wait()
        stderr = proc.stderr.read()

        if proc.returncode == 0:
            result["success"] = True
            if csv_lines:
                csv_data = "\n".join(csv_lines) + "\n"
                csv_path = output_dir / f"split_rep{replicate}_slim_summary.csv"
                with open(csv_path, "w", encoding="utf-8") as f:
                    f.write(csv_data)
                result["slim_csv"] = str(csv_path)
            else:
                result["slim_csv"] = None
        else:
            result["error"] = stderr
    except KeyboardInterrupt:
        raise
    except Exception as exc:
        result["error"] = str(exc)

    return result


def _kill_slim_children() -> None:
    """
        Kill any SLiM child processes spawned by this script.
    """
    pid = os.getpid()
    # Find all descendant PIDs by walking the process tree via /proc
    descendants: list[int] = []
    to_visit = [pid]
    while to_visit:
        parent = to_visit.pop()
        try:
            children_dir = Path(f"/proc/{parent}/task/{parent}/children")
            if children_dir.exists():
                child_pids = children_dir.read_text().split()
                for cpid_str in child_pids:
                    cpid = int(cpid_str)
                    descendants.append(cpid)
                    to_visit.append(cpid)
        except (OSError, ValueError):
            pass
    # Kill descendants that are running slim (in reverse order: deepest first)
    for dpid in reversed(descendants):
        try:
            cmdline = Path(f"/proc/{dpid}/cmdline").read_bytes()
            if b"slim" in cmdline:
                os.kill(dpid, signal.SIGKILL)
        except (OSError, ValueError):
            pass


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run Kyriazis DFE split-population SLiM simulations"
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=20,
        help="Number of replicates to run (default: 20)",
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=1,
        help="Number of parallel processes (default: 1)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional base seed for deterministic runs",
    )
    parser.add_argument(
        "--burnin-load",
        action="store_true",
        help="Load per-replicate burn-in state (burnin_Q{Q}_rep{N}.slim).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: load_dynamics/results_kyriazis)",
    )
    parser.add_argument(
        "--define",
        action="append",
        default=[],
        help="Extra SLiM -d defines, e.g. --define ANCESTRAL_NE=3000",
    )
    parser.add_argument(
        "--scaling-factor",
        type=int,
        default=1,
        help="Scaling factor Q for SLiM simulation (default: 1, unscaled)",
    )
    parser.add_argument(
        "--chromosomes",
        type=int,
        default=5,
        help="Number of chromosomes to simulate (default: 5, min: 2).",
    )
    parser.add_argument(
        "--rep-ids",
        type=str,
        default=None,
        help="Comma-separated replicate IDs to run, e.g. --rep-ids 11,20",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running",
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else RESULTS_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    n_chrom = max(2, args.chromosomes)

    if args.rep_ids:
        rep_list = [int(x.strip()) for x in args.rep_ids.split(",")]
    else:
        rep_list = list(range(1, args.replicates + 1))

    jobs = [
        (rep, output_dir, args.seed, args.burnin_load, args.define, args.scaling_factor, n_chrom)
        for rep in rep_list
    ]

    genes_total = n_chrom * 500
    genome_mb = genes_total * 1840 / 1e6
    ne_scaled = 25000 // args.scaling_factor
    burnin = ne_scaled * 5
    print("=== Kyriazis Split Simulation Runs ===")
    print(f"SLiM script: {SLIM_SCRIPT}")
    print(f"Replicates: {args.replicates}")
    print(f"Parallel processes: {args.parallel}")
    print(f"Output directory: {output_dir}")
    print(f"Burn-in load: {args.burnin_load}")
    print(f"Scaling factor Q: {args.scaling_factor}")
    print(f"Chromosomes: {n_chrom} ({genes_total} genes, {genome_mb:.1f} Mb)")
    print(f"Scaled Ne: {ne_scaled}, burn-in: {burnin} gen")
    print(f"Base seed: {args.seed}")
    if args.define:
        print(f"Extra defines: {args.define}")
    print("=" * 40 + "\n")

    if args.dry_run:
        print("Dry run - would execute the following:")
        for rep, out_dir, base_seed, burnin_load, extra, sf, nc in jobs[:5]:
            seed = stable_seed(rep, base_seed)
            cmd_preview = [
                "slim",
                "-d", f"SEED={seed}",
                "-d", f"REPLICATE={rep}",
                "-d", f'OUTPUT_DIR="{out_dir}"',
                "-d", f"Q={sf}",
                "-d", f"N_CHROM={nc}",
                "-d", f"BURNIN_LOAD={1 if burnin_load else 0}",
                "-d", f"BURNIN_SAVE={0 if burnin_load else 1}",
            ]
            for define in extra:
                cmd_preview.extend(["-d", define])
            cmd_preview.append(str(SLIM_SCRIPT))
            print("  " + " ".join(cmd_preview))
        if len(jobs) > 5:
            print(f"  ... and {len(jobs) - 5} more")
        return 0

    results = []
    failed = []

    try:
        if args.parallel > 1:
            with ProcessPoolExecutor(max_workers=args.parallel) as executor:
                futures = {
                    executor.submit(run_simulation, *job): job
                    for job in jobs
                }
                try:
                    for i, future in enumerate(as_completed(futures)):
                        result = future.result()
                        results.append(result)
                        status = "OK" if result["success"] else "FAILED"
                        csv_note = ""
                        if result["success"] and result.get("slim_csv"):
                            csv_note = " [CSV saved]"
                        print(
                            f"[{i+1}/{len(jobs)}] Replicate {result['replicate']}: {status}{csv_note}"
                        )
                        if not result["success"]:
                            failed.append(result)
                except KeyboardInterrupt:
                    print("\n\nInterrupted! Cancelling pending jobs and killing SLiM processes...", flush=True)
                    for fut in futures:
                        fut.cancel()
                    executor.shutdown(wait=False, cancel_futures=True)
                    _kill_slim_children()
                    return 1
        else:
            for i, job in enumerate(jobs):
                result = run_simulation(*job)
                results.append(result)
                status = "OK" if result["success"] else "FAILED"
                csv_note = ""
                if result["success"] and result.get("slim_csv"):
                    csv_note = " [CSV saved]"
                print(
                    f"[{i+1}/{len(jobs)}] Replicate {result['replicate']}: {status}{csv_note}"
                )
                if not result["success"]:
                    failed.append(result)
    except KeyboardInterrupt:
        print("\n\nInterrupted! Killing SLiM processes...", flush=True)
        _kill_slim_children()
        return 1

    print("\n=== Summary ===")
    print(f"Completed: {len(results) - len(failed)}/{len(results)}")
    print(f"Failed: {len(failed)}")

    if failed:
        print("\nFailed jobs:")
        for f in failed:
            err = f["error"] or ""
            print(f"  Replicate {f['replicate']}: {err[:120]}...")

    manifest_path = output_dir / "split_simulation_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"\nManifest saved to: {manifest_path}")

    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
