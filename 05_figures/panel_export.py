"""Helper: export figure panels.
Part of the analysis code of Gospodinova et al., IJMS 2026 (admsc-palmitate-8oxoG-methylation). """
import json, os
from matplotlib.transforms import Bbox

def export_panels(fig, panels, name, outdir, dpi=300, pad=0.06):
    os.makedirs(outdir, exist_ok=True)

    import re
    letters = [t for t in fig.texts if re.fullmatch(r'\([a-g]\)', t.get_text().strip())]
    for ax in fig.axes:
        letters += [t for t in ax.texts if re.fullmatch(r'\([a-g]\)', t.get_text().strip())]
    for t in letters: t.set_visible(False)
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    W, H = fig.get_size_inches(); pos = {}
    for letter, axes in panels.items():
        bb = Bbox.union([a.get_tightbbox(r) for a in axes]).transformed(fig.dpi_scale_trans.inverted())
        bb = Bbox.from_extents(bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad)
        fig.savefig(f'{outdir}/{name}_{letter}.png', dpi=dpi, bbox_inches=bb, facecolor='white')
        pos[letter] = dict(x0=bb.x0 / W, y0=1 - bb.y1 / H, w=bb.width / W, h=bb.height / H)
    for t in letters: t.set_visible(True)
    json.dump(dict(figure=name, width_in=W, height_in=H, panels=pos), open(f'{outdir}/{name}_panels.json', 'w'), indent=1)
    print(f'panels exported: {name} {sorted(pos)}')
