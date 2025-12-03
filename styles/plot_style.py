#!/usr/bin/env python3
"""
Standard plot style for PNP project.
Reference: PhD thesis plots (plot1.ipynb style)

Usage:
    from styles.plot_style import setup_plot_style, setup_axis_style

    fig, ax = plt.subplots(figsize=(3.7, 3.7))
    setup_axis_style(ax)
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


def setup_plot_style():
    """Setup global matplotlib style for publication-quality plots."""
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 10,
        'axes.labelweight': 'bold',
        'axes.titlesize': 10,
        'axes.titleweight': 'bold',
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 8,
        'figure.figsize': (3.7, 3.7),
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
    })


def setup_axis_style(ax, minor_x=None, minor_y=5):
    """
    Setup axis style for publication-quality plots.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes to style
    minor_x : int or float, optional
        Minor tick interval for x-axis. If int, uses AutoMinorLocator(n).
        If float, uses MultipleLocator(value). Default: AutoMinorLocator(5)
    minor_y : int or float, optional
        Minor tick interval for y-axis. Same logic as minor_x.
        Default: AutoMinorLocator(5)
    """
    # Tick positions on both sides
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    # Major ticks
    ax.xaxis.set_tick_params(which='major', direction='in', length=6, width=1, labelsize=10)
    ax.yaxis.set_tick_params(which='major', direction='in', length=6, width=1, labelsize=10)

    # Minor ticks
    ax.xaxis.set_tick_params(which='minor', direction='in', length=3)
    ax.yaxis.set_tick_params(which='minor', direction='in', length=3)
    ax.tick_params(which='minor', direction='in', length=3)

    # Minor tick locators
    if minor_x is not None:
        if isinstance(minor_x, int):
            ax.xaxis.set_minor_locator(AutoMinorLocator(minor_x))
        else:
            ax.xaxis.set_minor_locator(MultipleLocator(minor_x))
    else:
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    if minor_y is not None:
        if isinstance(minor_y, int):
            ax.yaxis.set_minor_locator(AutoMinorLocator(minor_y))
        else:
            ax.yaxis.set_minor_locator(MultipleLocator(minor_y))
    else:
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))


def set_labels(ax, xlabel, ylabel, fontsize=10, fontweight='bold'):
    """Set axis labels with standard formatting."""
    ax.set_xlabel(xlabel, fontsize=fontsize, fontweight=fontweight)
    ax.set_ylabel(ylabel, fontsize=fontsize, fontweight=fontweight)


# Standard figure sizes
FIGURE_SIZES = {
    'single': (3.7, 3.7),      # Single panel
    'wide': (6, 4),            # Wide single panel
    'double': (7, 3.5),        # Two panels side by side
    '2x2': (7, 7),             # 2x2 grid
    'thesis': (3.7, 3.7),      # Thesis standard
}

# Standard colors (viridis-based)
from matplotlib import cm
def get_viridis_colors(n):
    """Get n colors from viridis colormap."""
    cmap = cm.get_cmap('viridis', n)
    return [cmap(i) for i in range(n)]


# Example usage
if __name__ == '__main__':
    import numpy as np

    setup_plot_style()

    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    x = np.linspace(0, 10, 100)
    colors = get_viridis_colors(5)
    for i, c in enumerate(colors):
        ax.plot(x, np.sin(x + i), color=c, linewidth=1)

    set_labels(ax, r'$x$ (units)', r'$y$ (units)')
    setup_axis_style(ax)

    ax.set_xlim(0, 10)
    ax.set_ylim(-1.5, 1.5)

    plt.tight_layout()
    plt.savefig('test_style.png')
    print("Saved: test_style.png")
