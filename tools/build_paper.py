#!/usr/bin/env python3
"""Rebuild the authoritative v2 figures and cross-referenced PDFs from fixed raw data."""
import argparse
import csv
from fractions import Fraction as F
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import time
import zipfile

from inspect_zip import inspect
from reproduce import Replay, NAMES, ROOT, sha, write

ZIP_SHA = 'c863979b757ef5375b6a3537d8dfbe7a48608c57c51fdc4c8044f4f5e2350432'
FIGURES = ['readout_teaser', 'model_scales', 'profile_mechanism',
           'kernel_universality', 'finite_scale', 'terminal_law']


def source_intake(replay):
    archive = ROOT / 'source_inputs/critical_overflow_sicomp_v2_editable.zip'
    inventory = inspect(archive)
    replay.check('authoritative_zip', inventory['archive_size'] == 4054233 and inventory['archive_sha256'] == ZIP_SHA, inventory['archive_sha256'])
    target = replay.output / 'paper'
    target.mkdir()
    with zipfile.ZipFile(archive) as z:
        manifest = json.loads(z.read('SOURCE_FILES.sha256.json'))['files']
        names = {e['name'] for e in inventory['entries'] if not e['directory']}
        replay.check('authoritative_manifest_file_set', names == set(manifest) | {'SOURCE_FILES.sha256.json'}, {'manifest_files': len(manifest), 'all_files': len(names)})
        for name, expected in manifest.items():
            raw = z.read(name)
            if len(raw) != expected['bytes'] or hashlib.sha256(raw).hexdigest() != expected['sha256']:
                raise RuntimeError('authoritative member mismatch: ' + name)
        for entry in inventory['entries']:
            dest = target / entry['name']
            if entry['directory']:
                dest.mkdir(parents=True, exist_ok=True)
            else:
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(z.read(entry['name']))
    # The ZIP is the byte authority; editable TeX at the current entry must match it.
    for path in (ROOT / 'paper/manuscript').rglob('*'):
        if path.is_file() and 'figures' not in path.relative_to(ROOT / 'paper/manuscript').parts:
            name = 'manuscript/' + str(path.relative_to(ROOT / 'paper/manuscript'))
            replay.check('current_source/' + name, sha(path) == manifest[name]['sha256'], 'Unmodified authoritative source.')
            shutil.copyfile(path, target / name)
    # Validate and use the tracked figure sources/data as well; edits must never
    # be silently ignored in favor of the archived copies.
    for directory, prefix in [(ROOT / 'tools/paper', ''), (ROOT / 'paper/evidence', 'evidence/')]:
        for path in sorted(directory.iterdir()):
            if not path.is_file():continue
            name = prefix + path.name
            replay.check('current_source/' + name, name in manifest and sha(path) == manifest[name]['sha256'], 'Tracked source matches the authoritative baseline.')
            shutil.copyfile(path, target / name)
    write(replay.output / 'ZIP_INTAKE.json', {'inventory': inventory, 'source_manifest_entries_verified': len(manifest), 'original_PACKAGE_CHECKS_present': 'verification/PACKAGE_CHECKS.json' in names})
    # Discard only the disposable copy's retained PDFs/aux files and plot caches.
    shutil.rmtree(target / 'output')
    (target / 'output/pdf').mkdir(parents=True)
    shutil.rmtree(target / 'verification/plot_data')
    (target / 'verification/plot_data').mkdir()
    for name in FIGURES:
        for ext in ['pdf', 'svg', 'png']:
            (target / 'manuscript/figures' / (name + '.' + ext)).unlink()
    return target


def reconstruct_inputs(replay, work, paper):
    root = replay.output / 'figure-inputs'
    (root / 'analysis').mkdir(parents=True)
    (root / 'inputs').mkdir()
    reference = work / 'research' / NAMES['reference'] / 'results'
    lines = (reference / 'TABLES.hex').read_text().splitlines()
    assert lines[0] == '4097 16385 512'
    gbar = [float.fromhex(s) for s in lines[1:4098]]
    overflow = [float.fromhex(s) for s in lines[4098:]]
    g = list(csv.DictReader((paper / 'evidence/FROZEN_GBAR_NODES.csv').open()))
    s = list(csv.DictReader((paper / 'evidence/FROZEN_OVERFLOW_NODES.csv').open()))
    replay.check('authoritative_lookup_vs_frozen_reference', [float.fromhex(x['gbar_hex']) for x in g] == gbar and [float.fromhex(x['S0_hex']) for x in s] == overflow, {'gbar_nodes': len(g), 'overflow_nodes': len(s)})
    replay.check('authoritative_profile_vs_frozen_raw', sha(paper / 'evidence/raw_profiles.npz') == sha(work / 'research' / NAMES['profile'] / 'results/tail_h512_K24/raw_profiles.npz'), 'Exact NPZ bytes.')
    rows = []
    for group, batch, suffix, expected in [('finite', 'original', 'results/validation/outcomes.csv', 9216), ('prospective', 'prospective', 'results/micro/outcomes.csv', 16384)]:
        data = list(csv.DictReader((work / 'research' / NAMES[group] / suffix).open()))
        assert len(data) == expected
        for row in data:
            q = int(row['q']); d = int(row['Dcrit_integer'])
            node = F(abs(d) * 512, q)
            assert node.denominator == 1 and node < len(overflow)
            observed = int(row['Zcrit']) / q
            rows.append({'batch': batch, 'q': q, 'lambda': int(row['lambda']),
                         'kernel': 'matching' if row['kernel'] == '0' else 'product',
                         'replicate': int(row['replicate']), 'Dcrit_integer': d,
                         'observed': observed, 'imbalance_only': abs(d) / (2 * q),
                         'profile_prediction': overflow[int(node)]})
    endpoint = root / 'analysis/endpoint_comparison.csv'
    with endpoint.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    archives = {}
    for filename, group, members in [
        ('paper_v4_release_return.tar.gz', 'finite', ['results/validation/outcomes.csv']),
        ('sde_end_to_end_return.tar.gz', 'reference', ['results/TABLES.hex', 'results/reference_endpoints.csv'])]:
        path = root / 'inputs' / filename
        originals = []
        with tarfile.open(path, 'w:gz') as tar:
            for member in members:
                source = work / 'research' / NAMES[group] / member
                name = NAMES[group] + '/' + member
                data = source.read_bytes()
                info = tarfile.TarInfo(name);info.size = len(data);info.mtime = 0
                tar.addfile(info, io.BytesIO(data))
                originals.append({'workspace_path': 'research/' + name, 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()})
        archives[filename] = {'sha256': sha(path), 'members': originals,
            'scope': 'Runtime compatibility container containing selected byte-identical frozen members. It is not the complete original return archive and does not share its archive checksum.'}
    write(replay.output / 'FIGURE_INPUTS.json', {'endpoint_rows': len(rows), 'original_rows': 9216, 'prospective_rows': 16384, 'endpoint_sha256': sha(endpoint), 'archives': archives, 'new_random_words_or_trajectories': 0})
    replay.env['CRITICAL_SPILL_INPUT_ROOT'] = str(root)


def numeric_tables(replay, work, paper):
    base = lambda g: work / 'research' / NAMES[g] / 'results'
    raw = (paper / 'manuscript/sections/experiments.tex').read_text()
    means = json.loads((base('prospective') / 'MEAN_COMPARISONS.json').read_text())
    cdfs = json.loads((base('prospective') / 'CDF_COMPARISONS.json').read_text())
    matches = re.findall(r'([12])\s*&\s*(matching|product)\s*&\s*\$\[([^,]+),\\\s*([^\]]+)\]\$\s*&\s*([\d.]+)', raw)
    assert len(matches) == 4
    for lam, kernel, lo, hi, upper in matches:
        key = (int(lam), 0 if kernel == 'matching' else 1)
        a = next(x for x in means if (x['lambda'], x['kernel']) == key)
        b = next(x for x in cdfs if (x['lambda'], x['kernel']) == key)
        lower, high = map(F, a['true_mean_difference_interval'])
        replay.check('Table1/' + lam + '/' + kernel, F(lo.strip()) <= lower and F(hi.strip()) >= high and F(upper) >= F(b['true_CDF_distance_interval'][1]), {'outward_mean': [lo.strip(), hi.strip()], 'outward_CDF': upper, 'exact_mean': a['true_mean_difference_interval'], 'exact_CDF': b['true_CDF_distance_interval']})
    text = (paper / 'manuscript/sections/certificate_tables.tex').read_text()
    stages = re.findall(r'^(\d+)\s*&\s*(\d+)\s*&\s*(\d+)\s*&\s*(\d+)\s*&\s*([\d.]+)\s*&\s*([\d.]+)', text, re.M)
    assert len(stages) == 4
    for k, periods, radius, j, mass, upper in stages:
        stage = f'K{k}_B{periods}_J{j}'
        a = json.loads((base('certificate') / stage / 'CERTIFICATE.json').read_text())
        replay.check('TableSM1/' + stage, int(radius) == a['stage']['radius'] and F(mass) <= F(a['beta_lower_rational']) and F(upper) >= F(a['posterior_W1_upper_rational']), {'acceptance_met': a['acceptance_met'], 'printed_W1_upper': upper, 'exact_W1_upper': a['posterior_W1_upper_rational']})
    cert = json.loads((base('certificate') / 'CERTIFIED_READOUTS.json').read_text())
    intervals = re.findall(r'\$\[([^,]+),\s*\\?\s*([^\]]+)\]\$', text)
    assert len(intervals) == 5
    readouts = [cert['origin_alpha'], cert['feedback']['1'], cert['feedback']['2'], cert['spill_over_q']['h=0,d=0'], cert['spill_over_q']['h=0,d=1']]
    for index, ((lower, upper), item) in enumerate(zip(intervals, readouts)):
        replay.check('TableSM2/' + str(index + 1), F(lower.strip()) <= F(item['lower_rational']) and F(upper.strip()) >= F(item['upper_rational']), [lower.strip(), upper.strip()])
    primary = json.loads((base('finite') / 'PRIMARY_FAILURE_RATES.json').read_text())
    text = (paper / 'manuscript/sections/inference_details.tex').read_text()
    table = text.split('\\label{tab:primary}')[1].split('\\end{tabular}')[0]
    found = 0;metric = None
    for line in table.splitlines():
        if 'Feedback RMS' in line:metric = 'feedback_RMSE'
        if 'Overflow maximum' in line:metric = 'spill_window'
        m = re.search(r'&\s*([012]),\s*(fixed|product)\s*&\s*(\d+)\s*&\s*(\d+)\s*&\s*\$\[([^,]+),([^\]]+)\]', line)
        if not m:continue
        lam, kernel, certain, possible, lower, upper = m.groups()
        item = next(x for x in primary if x['lambda'] == int(lam) and x['kernel'] == (0 if kernel == 'fixed' else 1) and x['metric'] == metric)
        replay.check('TableSM3/' + metric + '/' + lam + '/' + kernel, int(certain) == item['certain_failures'] and int(possible) == item['possible_failures'] and F(lower) <= F(item['failure_probability_lower']) and F(upper) >= F(item['failure_probability_upper']), {'certain': certain, 'possible': possible, 'interval': [lower, upper]})
        found += 1
    assert found == 10
    replay.check('distinct_profile_error_bounds', F(cert['D_E_plus_e_rational']) < F('.000009038') and F(json.loads((base('certificate') / 'K10_B32_J128/CERTIFICATE.json').read_text())['posterior_W1_upper_rational']) < F('.000009034'), 'Uniform D_* and fourth-stage W1 quantities checked separately.')


def plot_and_compare(replay, paper):
    for script in ['prepare_figures.py', 'prepare_mechanism.py', 'prepare_teaser.py', 'prepare_kernel_figure.py']:
        replay.run('plot_' + script[:-3], [sys.executable, '-B', paper / script], paper, wall=120, cpu=120)
    checks = []
    archive = ROOT / 'source_inputs/critical_overflow_sicomp_v2_editable.zip'
    with zipfile.ZipFile(archive) as z:
        for path in sorted((paper / 'verification/plot_data').glob('*.csv')):
            name = 'verification/plot_data/' + path.name
            old = list(csv.DictReader(io.StringIO(z.read(name).decode())))
            new = list(csv.DictReader(path.open()))
            assert len(old) == len(new) and list(old[0]) == list(new[0])
            maximum = 0.0
            for a, b in zip(old, new):
                for key in a:
                    if a[key] == b[key]:continue
                    try:maximum = max(maximum, abs(float(a[key]) - float(b[key])))
                    except ValueError:raise RuntimeError('plot-data label mismatch: ' + name)
            replay.check('plot_data/' + path.name, maximum <= 5e-15, {'rows': len(new), 'largest_numeric_difference': maximum, 'sha256': sha(path)})
            checks.append({'file': path.name, 'rows': len(new), 'maximum_numeric_difference': maximum})
    for name in FIGURES:
        replay.check('figure/' + name, all((paper / 'manuscript/figures' / (name + '.' + ext)).is_file() for ext in ['pdf', 'svg', 'png']), 'Regenerated all three formats.')
    write(replay.output / 'PLOT_DATA_CHECKS.json', checks)


def compile_and_render(replay, paper):
    manuscript = paper / 'manuscript'
    docs = ['paper_sicomp', 'supplement_sicomp']
    for variable in ['BIBINPUTS', 'BSTINPUTS']:
        replay.env[variable] = str(manuscript) + os.pathsep + replay.env.get(variable, '')
    last = None;rounds = 0
    for iteration in range(1, 7):
        for doc in docs:
            replay.run(f'latex_{iteration}_{doc}', ['xelatex', '-interaction=nonstopmode', '-halt-on-error', '-no-shell-escape', '-output-directory=../output/pdf', doc + '.tex'], manuscript, wall=90, cpu=90)
            if iteration == 1:
                replay.run('bibtex_' + doc, ['bibtex', doc], paper / 'output/pdf', wall=30, cpu=30)
        current = tuple(sha(paper / 'output/pdf' / (doc + suffix)) for doc in docs for suffix in ['.aux', '.bbl'])
        rounds = iteration
        if iteration >= 3 and current == last:break
        last = current
    else:raise RuntimeError('cross-document references did not stabilize')
    reports = {}
    from PIL import Image, ImageDraw
    for doc in docs:
        pdf = paper / 'output/pdf' / (doc + '.pdf')
        log = pdf.with_suffix('.log').read_text(errors='replace')
        patterns = {'undefined_reference': r'(?:Reference|Citation).*undefined|There were undefined references', 'duplicate_label': r'multiply defined', 'overfull': r'Overfull \\[hv]box', 'missing_character': r'Missing character:'}
        counts = {key: len(re.findall(pattern, log)) for key, pattern in patterns.items()}
        replay.check('final_TeX_log/' + doc, not any(counts.values()), counts)
        info = subprocess.check_output(['pdfinfo', str(pdf)], text=True)
        pages = int(re.search(r'Pages:\s*(\d+)', info).group(1))
        fonts = subprocess.check_output(['pdffonts', str(pdf)], text=True)
        embedding = re.findall(r'\s+(yes|no)\s+(?:yes|no)\s+(?:yes|no)\s+\d+\s+\d+\s*$', fonts, re.M)
        replay.check('embedded_fonts/' + doc, bool(embedding) and all(x == 'yes' for x in embedding), {'font_count': len(embedding)})
        render = replay.output / 'render' / doc;render.mkdir(parents=True)
        replay.run('render_' + doc, ['pdftoppm', '-png', '-r', '85', pdf, render / 'page'], paper, wall=120, cpu=120)
        paths = sorted(render.glob('page-*.png'), key=lambda p: int(p.stem.split('-')[-1]))
        replay.check('rendered_pages/' + doc, len(paths) == pages, {'pages': pages})
        thumbs = []
        for index, path in enumerate(paths):
            im = Image.open(path).convert('RGB');im.thumbnail((245, 365))
            tile = Image.new('RGB', (265, 395), 'white');tile.paste(im, ((265 - im.width) // 2, 20))
            ImageDraw.Draw(tile).text((10, 378), str(index + 1), fill='black');thumbs.append(tile)
        sheet = Image.new('RGB', (265 * 4, 395 * ((len(thumbs) + 3) // 4)), '#bbbbbb')
        for i, im in enumerate(thumbs):sheet.paste(im, ((i % 4) * 265, (i // 4) * 395))
        sheet.save(replay.output / 'render' / (doc + '_overview.png'))
        replay.run('extract_' + doc, ['pdftotext', '-layout', pdf, replay.output / (doc + '.txt')], paper, wall=30, cpu=30)
        reports[doc] = {'pages': pages, 'sha256': sha(pdf), 'bytes': pdf.stat().st_size, 'fonts': len(embedding), 'final_log_counts': counts, 'rendered_pages': len(paths)}
    write(replay.output / 'DOCUMENT_CHECKS.json', {'engine': 'native XeLaTeX and BibTeX; unmodified SIAM class and TeX sources', 'engine_version': subprocess.check_output(['xelatex', '--version'], text=True).splitlines()[0], 'external_TEXINPUTS': os.environ.get('TEXINPUTS'), 'alternating_cross_document_rounds': rounds, 'documents': reports, 'page_count_difference_from_Tectonic_input': {doc: reports[doc]['pages'] - expected for doc, expected in zip(docs, [19, 12])}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True, help='New directory; existing paths are refused.')
    args = parser.parse_args()
    replay = Replay(args.out)
    replay.env.update({'SOURCE_DATE_EPOCH': '1792346031', 'MPLCONFIGDIR': str(replay.output / 'mplconfig'), 'TEXMFVAR': str(replay.output / 'texmf-var'), 'TEXMFCONFIG': str(replay.output / 'texmf-config')})
    error = None
    try:
        replay.verify()
        paper = source_intake(replay)
        work = replay.assemble('frozen-inputs')
        reconstruct_inputs(replay, work, paper)
        numeric_tables(replay, work, paper)
        plot_and_compare(replay, paper)
        compile_and_render(replay, paper)
    except Exception as exc:
        error = repr(exc)
    report = {'status': 'passed' if error is None else 'failed', 'error': error, 'wall_seconds': time.monotonic() - replay.start,
              'child_CPU_seconds': sum(x['CPU_seconds'] for x in replay.events),
              'checks_completed': len(replay.checks), 'stages_completed': len(replay.events),
              'child_max_RSS_kib_highwater': max((x['child_max_RSS_kib_highwater_so_far'] for x in replay.events), default=0),
              'all_completed_checks_passed': all(x['passed'] for x in replay.checks),
              'authoritative_zip_sha256': ZIP_SHA, 'authoritative_manuscript_read_in_this_environment': True,
              'display_figures': 6, 'tables': 4, 'new_random_words_or_trajectories': 0,
              'reuses_previous_core_stochastic_and_matrix_acceptance': True, 'proves_mathematical_theorems': False,
              'historical_dynamic_SDE_register_snapshot_available': False}
    write(replay.output / 'PACKAGE_CHECKS.json', report)
    print(json.dumps(report), flush=True)
    return 1 if error else 0


if __name__ == '__main__':
    sys.exit(main())
