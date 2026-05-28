import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from pathlib import Path


def save_fig(fig, base_path, formats=("png", "pdf"), dpi=300):
    """Save `fig` as `<base_path>.<fmt>` for each format, then close it."""
    base = Path(base_path)
    base.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(f"{base}.{fmt}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
