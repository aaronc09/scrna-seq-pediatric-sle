"""Render the corrected Table 2 panels without refitting any models.

Run from any directory: python Notebooks/render_table2_png.py
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'Results' / 'reviewer_revision'


def render():
    fig = plt.figure(figsize=(12, 12.5), facecolor='white')
    for panel, bottom, title, widths in [
        ('A', .535, 'Panel A. Interaction estimates', [.30, .20, .15, .19, .16]),
        ('B', .105, 'Panel B. Group means and contributing donor counts', [.30, .20, .125, .125, .125, .125]),
    ]:
        with (OUT / f'table2_panel_{panel}.csv').open(encoding='utf-8-sig', newline='') as f:
            records = list(csv.reader(f))
        assert len(records) == 5, 'Expected header and four significant illustrative combinations'
        headers, rows = records[0], records[1:]
        for row in rows:
            row[0] = row[0].replace('Non-classical Monocytes', 'Non-classical monocytes').replace('Classical Monocytes', 'Classical monocytes').replace(' to ', '\nto ')
            # Wrap long cell names without changing their meaning.
            row[0] = row[0].replace('\nto Non-classical monocytes', ' to\nNon-classical monocytes')
            row[0] = row[0].replace('Non-classical monocytes', 'Non-classical\nmonocytes')
            row[1] = row[1].replace('-', '-\n', 1)
            if panel == 'A':
                row[3] = row[3].replace(' to ', '\nto ')
            else:
                for col in range(2, 6):
                    row[col] = row[col].replace(' (', '\n(')
        headers = [h.replace('Ligand-receptor', 'Ligand–receptor').replace('Pediatric ', 'Pediatric\n').replace('Adult ', 'Adult\n') for h in headers]
        headers[0] = 'Source to\ntarget'
        headers[1] = 'Ligand-\nreceptor'
        ax = fig.add_axes([.03, bottom, .94, .37])
        ax.axis('off')
        ax.text(0, 1.06, title, fontsize=20, weight='bold', transform=ax.transAxes)
        table = ax.table(cellText=rows, colLabels=headers, colWidths=widths,
                         cellLoc='center', bbox=[0, 0, 1, 1])
        table.auto_set_font_size(False)
        table.set_fontsize(18)
        for (r, c), cell in table.get_celld().items():
            cell.set_edgecolor('#c1c6ca')
            cell.set_linewidth(.6)
            cell.PAD = .045
            line_count = max(value.count('\n') + 1 for value in ([headers] + rows)[r])
            cell.set_height(1.4 if line_count >= 3 else 1.0)
            if r == 0:
                cell.set_facecolor('#e7edf2')
                cell.get_text().set_weight('bold')
            elif r % 2 == 0:
                cell.set_facecolor('#f5f7f9')
            if c == 0:
                cell.get_text().set_ha('left')
    table_dir = ROOT / 'Results' / 'manuscript_tables'
    table_dir.mkdir(parents=True, exist_ok=True)
    path = table_dir / 'table2_interaction_estimates.png'
    temporary = path.with_name(path.stem + '_rendering.png')
    fig.savefig(temporary, dpi=600, facecolor='white', bbox_inches='tight', pad_inches=0.12)
    temporary.replace(path)
    plt.close(fig)
    print(path)


if __name__ == '__main__':
    render()
