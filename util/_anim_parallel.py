"""Shared helper: parallel frame rendering + ffmpeg stitching.

Pattern for callers:
    def _render_frame(k, tmpdir, *per_frame_args):
        # render frame k and save to tmpdir/frame_{k:06d}.png
        ...

    save_frames_parallel(
        _render_frame,
        [(k, per_frame_arg1, ...) for k in range(n_frames)],
        out_path, fps,
    )
"""
import os
import shutil
import subprocess
import tempfile
from multiprocessing import Pool


def save_frames_parallel(worker_fn, frame_args, out_path, fps, n_workers=None):
    """
    Render frames in parallel, stitch to MP4 with ffmpeg.

    worker_fn(k, tmpdir, *per_frame_args) must save tmpdir/frame_{k:06d}.png.
    frame_args: list of tuples (k, *per_frame_args), one per frame.
    """
    tmpdir = tempfile.mkdtemp(prefix="chromo_anim_")
    try:
        # Inject tmpdir as second arg after k.
        star_args = [(a[0], tmpdir) + a[1:] for a in frame_args]
        with Pool(processes=n_workers) as pool:
            pool.starmap(worker_fn, star_args)
        n = len(frame_args)
        # Encode quality is env-tunable (defaults reproduce the original output):
        #   ANIM_CRF      libx264 CRF (lower = higher quality; unset = x264 default ~23)
        #   ANIM_PIX_FMT  pixel format (default yuv420p; yuv444p avoids chroma loss on thin colored lines)
        #   ANIM_PRESET   x264 preset (default fast; slow = better compression at same CRF)
        pix_fmt = os.environ.get("ANIM_PIX_FMT", "yuv420p")
        preset = os.environ.get("ANIM_PRESET", "fast")
        crf = os.environ.get("ANIM_CRF")
        cmd = [
            "ffmpeg", "-y",
            "-framerate", str(fps),
            "-i", os.path.join(tmpdir, "frame_%06d.png"),
            "-c:v", "libx264",
            "-pix_fmt", pix_fmt,
            "-preset", preset,
        ]
        if crf is not None:
            cmd += ["-crf", crf]
        cmd += ["-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", out_path]
        subprocess.run(cmd, check=True, capture_output=True)
        print(f"wrote {out_path}  ({n} frames @ {fps} fps = {n / fps:.1f} s)")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
