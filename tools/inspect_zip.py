#!/usr/bin/env python3
"""Bounded, read-only ZIP inventory. Never extract or execute archive content."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import stat
import sys
import zipfile

MAX_ARCHIVE = 128 * 1024 * 1024
MAX_MEMBER = 64 * 1024 * 1024
MAX_TOTAL = 256 * 1024 * 1024
MAX_ENTRIES = 10000
MAX_RATIO = 1000


def inspect(path):
    path = Path(path)
    if path.stat().st_size > MAX_ARCHIVE:
        raise ValueError('archive exceeds 128 MiB intake limit')
    archive_hash = hashlib.sha256()
    with path.open('rb') as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b''):
            archive_hash.update(chunk)
    with zipfile.ZipFile(path) as archive:
        entries = archive.infolist()
        if len(entries) > MAX_ENTRIES:
            raise ValueError('too many ZIP entries')
        seen = set()
        total = 0
        for item in entries:
            name = item.filename
            parts = PurePosixPath(name).parts
            if (not name or name.startswith('/') or '\\' in name or ':' in name
                    or '..' in parts or any(ord(c) < 32 for c in name)):
                raise ValueError('unsafe ZIP path: ' + repr(name))
            key = '/'.join(parts).casefold()
            if key in seen:
                raise ValueError('duplicate ZIP path: ' + repr(name))
            seen.add(key)
            mode = item.external_attr >> 16
            file_type = stat.S_IFMT(mode)
            if file_type not in (0, stat.S_IFREG, stat.S_IFDIR):
                raise ValueError('nonregular ZIP entry: ' + repr(name))
            if item.flag_bits & 1:
                raise ValueError('encrypted ZIP member: ' + repr(name))
            total += item.file_size
            if item.file_size > MAX_MEMBER or total > MAX_TOTAL:
                raise ValueError('uncompressed size limit exceeded')
            if item.file_size > MAX_RATIO * max(1, item.compress_size):
                raise ValueError('compression ratio limit exceeded')
        inventory = []
        for item in entries:
            digest = hashlib.sha256()
            count = 0
            with archive.open(item) as member:
                while True:
                    block = member.read(min(1024 * 1024, MAX_MEMBER - count + 1))
                    if not block:
                        break
                    count += len(block)
                    if count > MAX_MEMBER or count > item.file_size:
                        raise ValueError('actual size exceeds declared size or limit')
                    digest.update(block)
            if count != item.file_size:
                raise ValueError('actual size differs from declared size')
            inventory.append({'name': item.filename, 'size': count,
                              'compressed_size': item.compress_size,
                              'crc32': f'{item.CRC:08x}',
                              'sha256': digest.hexdigest(),
                              'directory': item.is_dir()})
    return {'archive_name': path.name, 'archive_size': path.stat().st_size,
            'archive_sha256': archive_hash.hexdigest(),
            'entry_count': len(inventory), 'uncompressed_bytes': total,
            'entries': inventory, 'status': 'inventory_passed',
            'note': 'Content has not been executed or judged safe.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive')
    args = parser.parse_args()
    try:
        result = inspect(args.archive)
    except (OSError, ValueError, zipfile.BadZipFile, RuntimeError, NotImplementedError) as exc:
        print(json.dumps({'status': 'rejected', 'reason': str(exc)}, ensure_ascii=False))
        return 2
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
