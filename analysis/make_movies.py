#!/usr/bin/env python3

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# User settings.
FPS = 15.0
OUTPUT_DIR = Path("../output/B-weno3")

PLANES = ("xy", "yz", "xz")


def collect_frames(plot_dir: Path, plane: str):
    return sorted(plot_dir.glob(f"{plane}-*.png"))


def write_concat_file(frames, concat_path: Path, frame_duration: float) -> None:
    lines = []
    for frame in frames:
        lines.append(f"file '{frame.resolve()}'\n")
        lines.append(f"duration {frame_duration}\n")
    lines.append(f"file '{frames[-1].resolve()}'\n")
    concat_path.write_text("".join(lines))


def encode_movie(frames, output_path: Path, fps: float) -> None:
    frame_duration = 1.0 / fps
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as handle:
        concat_path = Path(handle.name)
    try:
        write_concat_file(frames, concat_path, frame_duration)
        cmd = [
            "ffmpeg",
            "-y",
            "-f", "concat",
            "-safe", "0",
            "-i", str(concat_path),
            "-vsync", "vfr",
            "-pix_fmt", "yuv420p",
            "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
            "-c:v", "libx264",
            "-crf", "18",
            "-preset", "medium",
            "-r", str(fps),
            str(output_path),
        ]
        subprocess.run(cmd, check=True)
    finally:
        concat_path.unlink(missing_ok=True)


def make_movie_for_plane(plot_dir: Path, plane: str, fps: float) -> bool:
    frames = collect_frames(plot_dir, plane)
    if not frames:
        print(f"no frames found for plane '{plane}' in {plot_dir}", file=sys.stderr)
        return False
    output_path = plot_dir / f"{plane}.mp4"
    encode_movie(frames, output_path, fps)
    print(f"wrote {output_path} ({len(frames)} frames)")
    return True


def main() -> None:
    if shutil.which("ffmpeg") is None:
        raise SystemExit("ffmpeg not found on PATH; install ffmpeg to encode movies")

    if not OUTPUT_DIR.is_dir():
        raise SystemExit(f"plot directory {OUTPUT_DIR} does not exist; run visualize.py first")

    any_success = False
    for plane in PLANES:
        any_success |= make_movie_for_plane(OUTPUT_DIR, plane, FPS)
    if not any_success:
        raise SystemExit("no movies were produced")


if __name__ == "__main__":
    main()
