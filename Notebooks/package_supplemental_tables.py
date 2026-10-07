"""Copy supplemental tables to numbered publication filenames without altering data."""
from pathlib import Path
import hashlib
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
GUIDE = ROOT / 'Results/reviewer_revision/figure_table_numbering.txt'
OUT = ROOT / 'Results/supplemental_tables'

def package():
    OUT.mkdir(parents=True, exist_ok=True)
    text = GUIDE.read_text(encoding='utf-8')
    current = None
    readme = ['SUPPLEMENTAL TABLES S1-S9', '',
              'Numbered publication copies of the analysis outputs. Data are unchanged.',
              'Letters identify panels (for example, Table_S3A.csv is Supplemental Table S3, Panel A).',
              'Files ending .csv.gz are compressed CSV files; decompress them before opening in a spreadsheet application.', '']
    updated = []
    count = 0
    for line in text.splitlines():
        heading = re.match(r'^Supplemental Table (S[1-9])\.', line)
        if heading:
            current = heading[1]
            readme.extend(['', line])
        elif line and not line.startswith(' '):
            current = None
        entry = re.fullmatch(r'  (?:(Panel ([A-F])): )?(Results/\S+\.csv(?:\.gz)?)', line)
        if current and entry:
            panel = entry[2] or ''
            source = ROOT / entry[3]
            suffix = '.csv.gz' if source.name.endswith('.csv.gz') else '.csv'
            dest = OUT / f'Table_{current}{panel}{suffix}'
            if source.resolve() != dest.resolve():
                shutil.copy2(source, dest)
            assert hashlib.sha256(source.read_bytes()).digest() == hashlib.sha256(dest.read_bytes()).digest()
            readme.append('  ' + dest.name)
            prefix = f'  Panel {panel}: ' if panel else '  '
            line = prefix + dest.relative_to(ROOT).as_posix()
            count += 1
        updated.append(line)
    assert count == 25, f'Expected 25 supplemental files, got {count}'
    (OUT/'README.txt').write_text('\n'.join(readme)+'\n', encoding='utf-8')
    GUIDE.write_text('\n'.join(updated)+'\n', encoding='utf-8')
    print(f'Packaged and verified {count} numbered files in {OUT}')

if __name__ == '__main__':
    package()
