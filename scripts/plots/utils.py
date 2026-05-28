import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from pathlib import Path


def save_fig(fig, paths, dpi=300):
    """Save `fig` to each path in `paths` (format inferred from suffix), then close it."""
    for p in paths:
        path = Path(p)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
