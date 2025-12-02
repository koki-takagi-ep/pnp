#!/usr/bin/env python3
"""
Create GIF animation from transient PNP solver snapshots.
Shows the time evolution of electric double layer formation.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from pathlib import Path
import glob
import sys

# Check if imageio is available
try:
    import imageio
    HAS_IMAGEIO = True
except ImportError:
    HAS_IMAGEIO = False
    print("Warning: imageio not found. Will try using PIL/Pillow instead.")
    try:
        from PIL import Image
        HAS_PIL = True
    except ImportError:
        HAS_PIL = False


def setup_axis_style(ax):
    """Setup axis style."""
    ax.tick_params(axis='both', which='major', direction='in', length=6, width=1, labelsize=9)
    ax.tick_params(axis='both', which='minor', direction='in', length=3, width=0.5)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)


def load_snapshot(filename):
    """Load a single snapshot file."""
    data = np.loadtxt(filename, comments='#')
    return {
        'x_nm': data[:, 0],
        'phi_mV': data[:, 2],
        'c_plus_norm': data[:, 6],
        'c_minus_norm': data[:, 7],
        'phi_gc_mV': data[:, 9],
    }


def load_time_info(snapshot_dir):
    """Load time information from time_info.dat."""
    time_file = Path(snapshot_dir) / 'time_info.dat'
    if not time_file.exists():
        return None

    times = {}
    with open(time_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                idx = int(parts[0])
                time_ns = float(parts[1])
                times[idx] = time_ns
    return times


def create_frame(data, time_ns, frame_idx, output_path, x_max=50):
    """Create a single frame for the animation."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # (a) Electric potential
    ax = axes[0]
    ax.plot(data['x_nm'], data['phi_mV'], 'b-', linewidth=2, label='PNP (transient)')
    ax.plot(data['x_nm'], data['phi_gc_mV'], 'r--', linewidth=1.5, alpha=0.7, label='G-C (steady)')
    ax.set_xlabel(r'$x$ [nm]', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'$\phi$ [mV]', fontsize=11, fontweight='bold')
    ax.set_title('Electric Potential', fontsize=11, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlim([0, x_max])
    ax.set_ylim([0, data['phi_mV'].max() * 1.1])
    setup_axis_style(ax)

    # (b) Ion concentrations
    ax = axes[1]
    ax.plot(data['x_nm'], data['c_plus_norm'], 'r-', linewidth=2, label=r'$c_+/c_0$ (cation)')
    ax.plot(data['x_nm'], data['c_minus_norm'], 'b-', linewidth=2, label=r'$c_-/c_0$ (anion)')
    ax.axhline(y=1.0, color='gray', linestyle=':', linewidth=1)
    ax.set_xlabel(r'$x$ [nm]', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'$c / c_0$', fontsize=11, fontweight='bold')
    ax.set_title('Ion Concentrations', fontsize=11, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlim([0, x_max])
    ax.set_yscale('log')
    ax.set_ylim([0.01, 100])
    setup_axis_style(ax)

    # Add time annotation
    fig.suptitle(f'Time: {time_ns:.1f} ns', fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_path, dpi=100, bbox_inches='tight')
    plt.close()


def create_animation_imageio(frame_files, output_gif, fps=10):
    """Create GIF animation using imageio."""
    images = []
    for frame_file in frame_files:
        images.append(imageio.imread(frame_file))

    # Add some frames at the end for pause effect
    for _ in range(fps * 2):  # 2 second pause at end
        images.append(images[-1])

    imageio.mimsave(output_gif, images, fps=fps, loop=0)
    print(f"Animation saved to: {output_gif}")


def create_animation_pillow(frame_files, output_gif, fps=10):
    """Create GIF animation using PIL/Pillow."""
    images = []
    for frame_file in frame_files:
        img = Image.open(frame_file)
        images.append(img.copy())

    # Add pause frames at end
    for _ in range(fps * 2):
        images.append(images[-1].copy())

    duration = int(1000 / fps)  # milliseconds per frame
    images[0].save(
        output_gif,
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=0
    )
    print(f"Animation saved to: {output_gif}")


def main():
    # Parse arguments
    snapshot_dir = Path('results/snapshots')
    output_gif = Path('results/edl_evolution.gif')
    x_max = 50  # Max x to show in plots [nm]
    fps = 10

    if len(sys.argv) > 1:
        snapshot_dir = Path(sys.argv[1])
    if len(sys.argv) > 2:
        output_gif = Path(sys.argv[2])

    print(f"Loading snapshots from: {snapshot_dir}")

    # Find all snapshot files
    snapshot_files = sorted(glob.glob(str(snapshot_dir / 'snapshot_*.dat')))
    if not snapshot_files:
        print(f"Error: No snapshot files found in {snapshot_dir}")
        print("Run the solver with --animation flag first:")
        print("  ./build/pnp_solver --animation --dt 0.1 --t-final 10")
        sys.exit(1)

    print(f"Found {len(snapshot_files)} snapshots")

    # Load time information
    time_info = load_time_info(snapshot_dir)
    if time_info is None:
        print("Warning: time_info.dat not found, using frame index as time")
        time_info = {i: i * 1.0 for i in range(len(snapshot_files))}

    # Create temporary directory for frames
    frame_dir = snapshot_dir / 'frames'
    frame_dir.mkdir(exist_ok=True)

    # Determine x_max from first snapshot
    first_data = load_snapshot(snapshot_files[0])
    x_max = min(x_max, first_data['x_nm'].max())

    # Create frames
    print("Creating animation frames...")
    frame_files = []
    for idx, snapshot_file in enumerate(snapshot_files):
        data = load_snapshot(snapshot_file)
        time_ns = time_info.get(idx, idx * 1.0)

        frame_path = frame_dir / f'frame_{idx:05d}.png'
        create_frame(data, time_ns, idx, frame_path, x_max=x_max)
        frame_files.append(str(frame_path))

        if (idx + 1) % 10 == 0 or idx == len(snapshot_files) - 1:
            print(f"  Created frame {idx + 1}/{len(snapshot_files)}")

    # Create GIF animation
    print("\nCreating GIF animation...")
    output_gif.parent.mkdir(parents=True, exist_ok=True)

    if HAS_IMAGEIO:
        create_animation_imageio(frame_files, str(output_gif), fps=fps)
    elif HAS_PIL:
        create_animation_pillow(frame_files, str(output_gif), fps=fps)
    else:
        print("Error: Neither imageio nor PIL/Pillow is available.")
        print("Install one of them: pip install imageio or pip install Pillow")
        sys.exit(1)

    # Also create a MP4 if ffmpeg is available
    try:
        import subprocess
        mp4_output = output_gif.with_suffix('.mp4')
        cmd = [
            'ffmpeg', '-y', '-framerate', str(fps),
            '-i', str(frame_dir / 'frame_%05d.png'),
            '-c:v', 'libx264', '-pix_fmt', 'yuv420p',
            '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
            str(mp4_output)
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"MP4 video saved to: {mp4_output}")
    except Exception:
        pass  # ffmpeg not available

    print("\nDone!")


if __name__ == '__main__':
    main()
