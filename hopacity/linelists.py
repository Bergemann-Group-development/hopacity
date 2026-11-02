"""Fetch the master line lists that OpenMULTI3D publishes on MPCDF Keeper.

The lists are plain files in one folder of a public Seafile share, one Turbospectrum-format
file per species (``vald_atoms.txt`` holds the atoms, the others one molecule or
isotopologue each) and ``Hlinedata`` with the hydrogen lines. The folder they are fetched
into can be given to ``HopEOS.opacity(linelist=...)`` as it is, or the files one by one.

Command line::

    python -m hopacity.linelists list
    python -m hopacity.linelists fetch DEST [--file NAME ...]

``fetch`` without ``--file`` downloads all files (about 1 GB). A file already complete in
``DEST`` is skipped, and an interrupted transfer resumes from its ``.part`` file. Only the
standard library is used.
"""
import argparse
import json
import ssl
import sys
import urllib.parse
import urllib.request
from pathlib import Path

SERVER = 'https://keeper.mpdl.mpg.de'
TOKEN = '3191475aebe340638422'
FOLDER = '/linelists/master_linelists'
PAGE = f'{SERVER}/d/{TOKEN}/?p={urllib.parse.quote(FOLDER)}'
TIMEOUT = 120
CHUNK = 1 << 20

__all__ = ['list_files', 'fetch']


def _ssl_context():
    """The python.org macOS builds ship without a CA trust store; certifi fills that gap
    where it is installed."""
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def _open(url, headers=None):
    req = urllib.request.Request(url, headers=headers or {})
    return urllib.request.urlopen(req, timeout=TIMEOUT, context=_ssl_context())


def list_files():
    """The published files as a dict ``{name: size in bytes}``, sorted by name."""
    url = (f'{SERVER}/api/v2.1/share-links/{TOKEN}/dirents/?'
           + urllib.parse.urlencode({'path': FOLDER + '/'}, quote_via=urllib.parse.quote))
    with _open(url, {'Accept': 'application/json'}) as resp:
        entries = json.load(resp)['dirent_list']
    return dict(sorted((e['file_name'], int(e['size'])) for e in entries if not e['is_dir']))


def _download(name, size, dest, progress):
    url = (f'{SERVER}/d/{TOKEN}/files/?'
           + urllib.parse.urlencode({'p': f'{FOLDER}/{name}', 'dl': '1'},
                                    quote_via=urllib.parse.quote))
    part = dest.with_name(dest.name + '.part')
    have = part.stat().st_size if part.exists() else 0
    if have >= size:
        have = 0
    with _open(url, {'Range': f'bytes={have}-'} if have else None) as resp:
        if have and resp.status != 206:     # the server ignored the range
            have = 0
        with part.open('ab' if have else 'wb') as fh:
            while buf := resp.read(CHUNK):
                fh.write(buf)
                have += len(buf)
                if progress:
                    progress(name, have, size)
    if have != size:
        raise OSError(f'{name}: got {have} of {size} bytes; run fetch again to resume')
    part.replace(dest)


def fetch(dest, names=None, progress=None):
    """Download the files ``names`` (default: all) into the folder ``dest``.

    Files already complete there are skipped. ``progress(name, done, size)`` is called as
    bytes arrive. Returns the paths of the requested files.
    """
    dest = Path(dest).expanduser()
    published = list_files()
    names = list(published) if names is None else list(names)
    unknown = [n for n in names if n not in published]
    if unknown:
        raise ValueError(f'not on the share: {", ".join(unknown)}; '
                         f'published are {", ".join(published)}')
    dest.mkdir(parents=True, exist_ok=True)
    for name in names:
        path = dest / name
        if not (path.is_file() and path.stat().st_size == published[name]):
            _download(name, published[name], path, progress)
    return [dest / n for n in names]


def _size(n):
    for unit in ('B', 'KiB', 'MiB', 'GiB'):
        if n < 1024 or unit == 'GiB':
            return f'{n:.0f} {unit}' if unit == 'B' else f'{n:.1f} {unit}'
        n /= 1024


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m hopacity.linelists',
                                 description=f'Fetch the master line lists from {PAGE}')
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('list', help='list the published files and their sizes')
    f = sub.add_parser('fetch', help='download files into a folder')
    f.add_argument('dest', type=Path, help='destination folder')
    f.add_argument('--file', action='append', dest='names', metavar='NAME',
                   help='a file to fetch (repeatable); default: all files')
    args = ap.parse_args(argv)

    if args.cmd == 'list':
        files = list_files()
        for name, size in files.items():
            print(f'{_size(size):>10}  {name}')
        print(f'{_size(sum(files.values())):>10}  total, {len(files)} files')
        return 0

    def progress(name, done, size):
        sys.stdout.write(f'\r{name:<20} {100 * done / size:5.1f}% of {_size(size)}')
        sys.stdout.flush()
        if done == size:
            sys.stdout.write('\n')
    try:
        paths = fetch(args.dest, args.names, progress)
    except ValueError as exc:
        ap.error(str(exc))
    print(f'{len(paths)} files in {args.dest}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
