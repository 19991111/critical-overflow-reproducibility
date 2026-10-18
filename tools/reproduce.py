#!/usr/bin/env python3
"""Portable replay of selected, hash-locked evidence; never replace sealed inputs."""
import argparse
import csv
from fractions import Fraction
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
NAMES = {'finite': 'paper_v4_release_20261017',
         'prospective': 'paper_micro_prospective_20261017',
         'reference': 'paper_sde_end_to_end_20261017',
         'certificate': 'paper_profile_certificate_20261017',
         'profile': 'paper_profile_numerics_20261017'}
VOLATILE = {'elapsed_analysis_budget_seconds', 'elapsed_budget_seconds',
            'elapsed_seconds', 'wall_seconds', 'CPU_seconds',
            'elapsed_matrix_budget_seconds'}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')


def scientific(value):
    if isinstance(value, dict):
        return {k: scientific(v) for k, v in value.items() if k not in VOLATILE}
    if isinstance(value, list):
        return [scientific(v) for v in value]
    return value


class Replay:
    def __init__(self, output):
        self.output = output.resolve()
        self.output.mkdir(parents=True, exist_ok=False)
        self.start = time.monotonic()
        self.core = min(os.sched_getaffinity(0))
        self.events = []
        self.checks = []
        self.files = json.loads((ROOT / 'provenance/FILES.json').read_text())['files']
        self.env = os.environ.copy()
        for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
            self.env[key] = '1'
        self.env['PYTHONDONTWRITEBYTECODE'] = '1'

    def check(self, name, passed, detail):
        record = {'name': name, 'passed': bool(passed), 'detail': detail}
        self.checks.append(record)
        write(self.output / 'checks.json', self.checks)
        if not passed:
            raise RuntimeError(name + ': ' + str(detail))

    def verify(self):
        for item in self.files:
            path = ROOT / item['stored_path']
            if sha(path) != item['stored_sha256']:
                raise RuntimeError('stored checksum mismatch: ' + item['stored_path'])
            raw = gzip.decompress(path.read_bytes()) if item['compression'] == 'gzip' else path.read_bytes()
            if len(raw) != item['bytes'] or hashlib.sha256(raw).hexdigest() != item['sha256']:
                raise RuntimeError('source checksum mismatch: ' + item['workspace_path'])
        self.check('frozen_inputs', True, {'files': len(self.files), 'all_lengths_and_SHA256_match': True})
        seeds = {}
        for group, filename, count in [('finite', 'SEEDS_V1.txt', 1024), ('prospective', 'SEEDS_V1.txt', 4096), ('reference', 'REFERENCE_SEEDS_V1.txt', 32768)]:
            rows = [line.split() for line in (ROOT / 'configs' / group / filename).read_text().splitlines()]
            values = [int(row[1], 16) for row in rows]
            self.check(group + '_seeds', len(rows) == count and [int(r[0]) for r in rows] == list(range(count)) and len(set(values)) == count, {'labels': len(rows), 'sha256': sha(ROOT / 'configs' / group / filename)})
            seeds[group] = set(values)
        self.check('seed_namespaces_disjoint', not (seeds['finite'] & seeds['prospective'] or seeds['finite'] & seeds['reference'] or seeds['prospective'] & seeds['reference']), 'No new labels generated or selected.')

    def assemble(self, name):
        work = self.output / name / 'work'
        work.mkdir(parents=True)
        for item in self.files:
            target = work / item['workspace_path']
            target.parent.mkdir(parents=True, exist_ok=True)
            source = ROOT / item['stored_path']
            if item['compression'] == 'gzip':
                with gzip.open(source, 'rb') as inp, target.open('wb') as out:
                    shutil.copyfileobj(inp, out)
            else:
                shutil.copyfile(source, target)
        for group in NAMES:
            base = work / 'research' / NAMES[group]
            (base / 'logs').mkdir(exist_ok=True)
            (base / 'results').mkdir(exist_ok=True)
        return work

    def run(self, name, command, cwd, wall=600, cpu=600):
        logs = self.output / 'logs'
        logs.mkdir(exist_ok=True)
        def limits():
            os.sched_setaffinity(0, {self.core})
            resource.setrlimit(resource.RLIMIT_CPU, (cpu, cpu + 1))
            resource.setrlimit(resource.RLIMIT_AS, (1024 ** 3, 1024 ** 3))
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        started = time.monotonic()
        code = None
        timed_out = False
        with (logs / (name + '.stdout.log')).open('x') as stdout, (logs / (name + '.stderr.log')).open('x') as stderr:
            try:
                process = subprocess.run(list(map(str, command)), cwd=cwd, env=self.env,
                                         preexec_fn=limits, timeout=wall,
                                         stdout=stdout, stderr=stderr)
                code = process.returncode
            except subprocess.TimeoutExpired:
                timed_out = True
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        event = {'stage': name, 'command': list(map(str, command)), 'cwd': str(cwd),
                 'exit_code': code, 'timed_out': timed_out,
                 'wall_seconds': time.monotonic() - started,
                 'CPU_seconds': after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
                 'child_max_RSS_kib_highwater_so_far': after.ru_maxrss,
                 'selected_cpu': self.core, 'threads': 1,
                 'wall_limit_seconds': wall, 'CPU_limit_seconds': cpu}
        self.events.append(event)
        write(self.output / 'stages.json', self.events)
        print(json.dumps({'stage': name, 'exit_code': code, 'wall_seconds': event['wall_seconds']}), flush=True)
        if code != 0:
            raise RuntimeError('stage failed: ' + name)

    def build(self, work):
        binaries = self.output / 'bin'
        binaries.mkdir()
        for group, filename in [('finite', 'finite_q.cpp'), ('prospective', 'micro.cpp'), ('reference', 'reference.cpp')]:
            flags = ['-std=c++17', '-O3', '-Wall', '-Wextra', '-Wpedantic', '-fno-fast-math']
            if group != 'finite':
                flags.append('-ffp-contract=off')
            self.run('compile_' + group, ['g++', *flags, ROOT / 'src' / group / filename, '-o', binaries / group], work, wall=90, cpu=90)
            args = [binaries / group, '--selftest']
            if group == 'finite':
                args.append(work / 'research' / NAMES[group] / 'results/WAVE_TABLE.hex')
            if group == 'reference':
                args.append(work / 'research' / NAMES[group] / 'results/TABLES.hex')
            self.run('selftest_' + group, args, work, wall=30, cpu=30)
        write(self.output / 'binary-freeze.json', {group: {'sha256': sha(binaries / group), 'source_sha256': sha(ROOT / 'src' / group / filename)} for group, filename in [('finite', 'finite_q.cpp'), ('prospective', 'micro.cpp'), ('reference', 'reference.cpp')]})

    def compare_json(self, work, group, filename):
        path = work / 'research' / NAMES[group] / 'results' / filename
        entry = next(e for e in self.files if e['workspace_path'] == str(path.relative_to(work)))
        stored = ROOT / entry['stored_path']
        expected = json.loads(gzip.decompress(stored.read_bytes()) if entry['compression'] == 'gzip' else stored.read_bytes())
        actual = json.loads(path.read_text())
        self.check(group + '/' + filename, scientific(expected) == scientific(actual), 'Scientific fields equal; explicit runtime fields excluded.')

    def recompute(self, work):
        base = lambda group: work / 'research' / NAMES[group]
        self.run('finite_raw_audit', [sys.executable, '-B', base('finite') / 'src/audit_observations.py'], work, wall=180, cpu=180)
        self.run('finite_analysis', [sys.executable, '-B', base('finite') / 'src/analyze_validation.py'], work, wall=120, cpu=120)
        for name in ['DESCRIPTIVE.json', 'PRIMARY_FAILURE_RATES.json', 'UNIVERSALITY_CDF.json', 'PAIRED_MEANS_NOMINAL.json', 'LAMBDA0_REFERENCE.json']:
            self.compare_json(work, 'finite', name)
        self.run('reference_means', [self.output / 'bin/reference', '--means', base('reference') / 'results/reference_endpoints.csv'], work, wall=180, cpu=180)
        means = json.loads((self.output / 'logs/reference_means.stdout.log').read_text())
        old = json.loads((base('reference') / 'logs/reference_means.stdout.log').read_text())
        self.check('reference_mean_grid', means == old, means)
        self.run('prospective_analysis', [sys.executable, '-B', base('prospective') / 'src/analyze.py'], work, wall=180, cpu=180)
        self.run('prospective_raw_audit', [sys.executable, '-B', base('prospective') / 'src/verify.py'], work, wall=180, cpu=180)
        for name in ['MICRO_MEAN_INTERVALS.json', 'REFERENCE_QUALITY_REUSED.json', 'MEAN_COMPARISONS.json', 'CDF_COMPARISONS.json']:
            self.compare_json(work, 'prospective', name)
        self.run('certified_readouts', [sys.executable, '-B', base('certificate') / 'src/certify_readouts.py'], work, wall=45, cpu=45)
        self.compare_json(work, 'certificate', 'CERTIFIED_READOUTS.json')
        for group, filename, script in [('finite', 'WAVE_TABLE.hex', 'export_wave.py'), ('reference', 'TABLES.hex', 'export_reference_tables.py')]:
            target = base(group) / 'results' / filename
            expected_hash = sha(target)
            target.unlink()
            self.run(group + '_table_export', [sys.executable, '-B', base(group) / 'src' / script], work, wall=120, cpu=120)
            self.check(group + '_table_bytes', sha(target) == expected_hash, {'expected': expected_hash, 'actual': sha(target)})
        self.mechanism(work)

    def mechanism(self, work):
        output = self.output / 'analysis'
        output.mkdir(exist_ok=True)
        prospective = work / 'research' / NAMES['prospective'] / 'results/micro/outcomes.csv'
        rows = [r for r in csv.DictReader(prospective.open()) if r['lambda'] == '1' and r['kernel'] == '0']
        table = (work / 'research' / NAMES['reference'] / 'results/TABLES.hex').read_text().splitlines()
        assert table[0] == '4097 16385 512'
        values = [Fraction(float.fromhex(x)) for x in table[1 + 4097:]]
        baseline = Fraction(0)
        profile = Fraction(0)
        with (output / 'observed_imbalance_readout.csv').open('w') as stream:
            writer = csv.writer(stream)
            writer.writerow(['replicate', 'D_over_q', 'observed_Z_over_q', 'mass_only_prediction', 'profile_prediction'])
            for row in rows:
                d = int(row['Dcrit_integer']); z = Fraction(int(row['Zcrit']), 512)
                naive = Fraction(abs(d), 1024)
                assert abs(d) < len(values)
                prediction = values[abs(d)]
                baseline += abs(z - naive); profile += abs(z - prediction)
                writer.writerow([row['replicate'], str(Fraction(d, 512)), str(z), str(naive), str(prediction)])
        result = {'rows': len(rows), 'q': 512, 'lambda': 1, 'kernel': 'identity_matching',
                  'both_predictions_use_observed_imbalance': True,
                  'mass_only_MAE': float(baseline / len(rows)),
                  'profile_MAE': float(profile / len(rows)),
                  'profile_lookup_uniform_error_upper': '0.000026',
                  'is_unconditional_diffusion_test': False,
                  'is_published_algorithm_comparison': False,
                  'confidence_claim_from_MAE': False}
        write(output / 'mechanism.json', result)
        desc = json.loads((work / 'research' / NAMES['finite'] / 'results/DESCRIPTIVE.json').read_text())
        with (output / 'endpoint_comparison.csv').open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(desc[0]))
            writer.writeheader();writer.writerows(desc)
        write(output / 'CSV_SCOPE.json', {'endpoint_comparison.csv': 'Reconstructed from fixed original outcomes via DESCRIPTIVE.json. This core schema is explicit; compatibility with unavailable manuscript prepare_figures.py remains unchecked.', 'mechanism': result})

    def rerun(self, work, group):
        base = work / 'research' / NAMES[group]
        output = self.output / 'rerun' / group
        output.mkdir(parents=True)
        binary = self.output / 'bin' / group
        if group == 'finite':
            args = [binary, '--run', base / 'results/WAVE_TABLE.hex', base / 'protocol/SEEDS_V1.txt', output]
            expected = base / 'results/validation'
        elif group == 'prospective':
            args = [binary, '--run', base / 'protocol/SEEDS_V1.txt', output]
            expected = base / 'results/micro'
        else:
            args = [binary, '--sample', base / 'results/TABLES.hex', base / 'protocol/REFERENCE_SEEDS_V1.txt', output / 'reference_endpoints.csv']
            expected = base / 'results'
            self.run('reference_micro_precision_gate', [sys.executable, '-B', base / 'src/analyze.py', '--micro'], work, wall=120, cpu=120)
            gate = json.loads((base / 'results/MICRO_PRECHECK.json').read_text())
            self.check('reference_micro_precision_gate', gate['all_four_precision_gates_passed'], gate)
            self.run('reference_cost_calibration', [binary, '--calibrate', base / 'results/TABLES.hex'], work, wall=60, cpu=60)
            calibration = json.loads((self.output / 'logs/reference_cost_calibration.stdout.log').read_text())
            self.check('reference_cost_gate', calibration['main_CPU_prediction_with_factor1_5'] <= 480, calibration)
        write(output / 'INPUT_FREEZE.json', {'before_trajectory_generation': True, 'source_sha256': sha(ROOT / 'src' / group / {'finite':'finite_q.cpp','prospective':'micro.cpp','reference':'reference.cpp'}[group]), 'binary_sha256': sha(binary), 'seed_sha256': sha(base / 'protocol' / ('REFERENCE_SEEDS_V1.txt' if group == 'reference' else 'SEEDS_V1.txt')), 'selected_cpu': self.core, 'threads': 1, 'new_seed_labels': 0, 'purpose': 'Replay the frozen batch, not an additional independent statistical sample.'})
        # Kernels require an empty data directory. Keep the pre-run freeze beside it.
        data = output / 'data';data.mkdir()
        args = [data if a == output else data / 'reference_endpoints.csv' if a == output / 'reference_endpoints.csv' else a for a in args]
        self.run('rerun_' + group, args, work, wall=1200 if group == 'reference' else 600, cpu=600)
        comparisons = []
        for path in sorted(data.iterdir()):
            old = expected / path.name
            equal = old.exists() and sha(path) == sha(old)
            passed = equal
            reason = 'exact raw bytes'
            if group == 'reference' and path.name.endswith('.progress.csv'):
                current = list(csv.DictReader(path.open()))
                previous = list(csv.DictReader(old.open()))
                labels = [int(r['replicate']) for r in current]
                cpus = [float(r['CPU_seconds']) for r in current]
                passed = labels == [int(r['replicate']) for r in previous] == list(range(256, 32769, 256)) and cpus == sorted(cpus)
                reason = 'all frozen progress labels match; measured CPU times are retained and expected to differ'
            comparisons.append({'file': path.name, 'sha256': sha(path), 'expected_sha256': sha(old) if old.exists() else None, 'exact_bytes_equal': equal, 'passed': passed, 'comparison': reason})
        write(output / 'BYTE_COMPARISONS.json', comparisons)
        self.check(group + '_rerun_bytes', bool(comparisons) and all(x['passed'] for x in comparisons), comparisons)

    def certificates(self):
        # Each calculation gets its own disposable copy of the sealed inputs.
        # Freeze files are not rewritten to match newly timed metadata.
        for direction, script in [('forward', 'run_frozen.py'), ('backward', 'backcheck.py')]:
            work = self.assemble('certificate_' + direction)
            base = work / 'research' / NAMES['certificate']
            if direction == 'forward':
                shutil.rmtree(base / 'results');(base / 'results').mkdir()
                shutil.rmtree(base / 'logs');(base / 'logs').mkdir()
            self.run('certificate_' + direction, [sys.executable, '-B', base / 'src' / script], work, wall=310 if direction == 'forward' else 90, cpu=310 if direction == 'forward' else 90)
            if direction == 'forward':
                for stage in ['K4_B8_J64', 'K6_B16_J64', 'K8_B24_J128', 'K10_B32_J128']:
                    self.compare_json(work, 'certificate', stage + '/CERTIFICATE.json')
                    self.compare_json(work, 'certificate', stage + '/LOWER_MATRIX.json')
                self.compare_json(work, 'certificate', 'SUMMARY.json')
            else:
                self.compare_json(work, 'certificate', 'INDEPENDENT_BACKWARD_CHECK.json')

    def finish(self, error):
        report = {'status': 'passed' if error is None else 'failed', 'error': error,
                  'wall_seconds': time.monotonic() - self.start,
                  'child_CPU_seconds': sum(x['CPU_seconds'] for x in self.events),
                  'child_max_RSS_kib_highwater': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                  'selected_cpu': self.core, 'threads': 1,
                  'platform': platform.platform(), 'python': sys.version,
                  'all_completed_checks_passed': all(x['passed'] for x in self.checks),
                  'checks_completed': len(self.checks), 'stages_completed': len(self.events),
                  'authoritative_manuscript_read_in_this_environment': False,
                  'proves_mathematical_theorems': False,
                  'historical_reference_dynamic_FP_snapshot_available': False}
        write(self.output / 'PACKAGE_CHECKS.json', report)
        print(json.dumps(report), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['quick', 'recompute', 'rerun', 'certificates', 'all'])
    parser.add_argument('--batch', choices=['finite', 'prospective', 'reference'])
    parser.add_argument('--out', type=Path, required=True, help='New output directory; existing paths are refused.')
    args = parser.parse_args()
    if args.mode == 'rerun' and not args.batch:
        parser.error('rerun requires --batch')
    replay = Replay(args.out)
    error = None
    try:
        replay.verify()
        if args.mode != 'certificates':
            work = replay.assemble('core')
            replay.build(work)
        if args.mode in ['recompute', 'all']:
            replay.recompute(work)
        if args.mode == 'rerun':
            replay.rerun(work, args.batch)
        if args.mode == 'all':
            for group in ['finite', 'prospective', 'reference']:
                replay.rerun(work, group)
        if args.mode in ['certificates', 'all']:
            replay.certificates()
    except Exception as exc:
        error = repr(exc)
    finally:
        replay.finish(error)
    return 1 if error else 0


if __name__ == '__main__':
    sys.exit(main())
