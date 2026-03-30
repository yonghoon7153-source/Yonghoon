"""
Supabase Storage sync module.
Keeps uploads/, results/, archive/ synced with Supabase Storage
so data persists across Render redeploys.

Required env vars:
  SUPABASE_URL  - e.g. https://xxxxx.supabase.co
  SUPABASE_KEY  - service_role key (not anon)
  SUPABASE_BUCKET - bucket name (default: dem-data)
"""

import os
import io
import json
import gzip
import logging
import requests
from pathlib import Path

# Files larger than this (bytes) will be gzip-compressed before upload
COMPRESS_THRESHOLD = 40 * 1024 * 1024  # 40MB
COMPRESS_EXTENSIONS = {'.liggghts', '.csv', '.txt', '.md', '.json', '.stl'}
CHUNK_SIZE = 45 * 1024 * 1024  # 45MB chunks for Supabase 50MB limit

logger = logging.getLogger(__name__)

SUPABASE_URL = os.environ.get('SUPABASE_URL', '')
SUPABASE_KEY = os.environ.get('SUPABASE_KEY', '')
SUPABASE_BUCKET = os.environ.get('SUPABASE_BUCKET', 'dem-data')

_enabled = False


def is_enabled():
    return _enabled


def init():
    """Initialize storage sync. Returns True if configured."""
    global _enabled
    if SUPABASE_URL and SUPABASE_KEY:
        _enabled = True
        print(f"[Storage] Supabase enabled: {SUPABASE_URL} bucket={SUPABASE_BUCKET}")
        _ensure_bucket()
    else:
        print(f"[Storage] NOT configured - SUPABASE_URL={bool(SUPABASE_URL)} SUPABASE_KEY={bool(SUPABASE_KEY)}")
    return _enabled


def _headers():
    return {
        'apikey': SUPABASE_KEY,
        'Authorization': f'Bearer {SUPABASE_KEY}',
    }


def _ensure_bucket():
    """Create bucket if it doesn't exist."""
    url = f"{SUPABASE_URL}/storage/v1/bucket/{SUPABASE_BUCKET}"
    r = requests.get(url, headers=_headers())
    if r.status_code == 404:
        create_url = f"{SUPABASE_URL}/storage/v1/bucket"
        payload = {'id': SUPABASE_BUCKET, 'name': SUPABASE_BUCKET, 'public': False}
        r2 = requests.post(create_url, json=payload, headers=_headers())
        if r2.status_code in (200, 201):
            logger.info(f"Created bucket: {SUPABASE_BUCKET}")
        else:
            logger.error(f"Failed to create bucket: {r2.text}")


def _should_compress(local_path):
    """Check if file should be gzip-compressed (large text files)."""
    ext = os.path.splitext(local_path)[1].lower()
    size = os.path.getsize(local_path)
    return size > COMPRESS_THRESHOLD and ext in COMPRESS_EXTENSIONS


def _upload_raw(remote_path, data, ct='application/octet-stream'):
    """Upload raw bytes to Supabase Storage."""
    url = f"{SUPABASE_URL}/storage/v1/object/{SUPABASE_BUCKET}/{remote_path}"
    headers = _headers()
    headers['x-upsert'] = 'true'
    headers['Content-Type'] = ct
    r = requests.post(url, headers=headers, data=data)
    if r.status_code not in (200, 201):
        r = requests.put(url, headers=headers, data=data)
    if r.status_code not in (200, 201):
        print(f"  Upload failed {remote_path}: {r.status_code} {r.text[:200]}")
        return False
    return True


def upload_file(remote_path, local_path):
    """Upload a single file to Supabase Storage.
    Large text files are gzip-compressed. If still >45MB, split into chunks."""
    if not _enabled:
        return

    try:
        # Prepare data (compress if needed)
        if _should_compress(local_path):
            with open(local_path, 'rb') as f:
                data = gzip.compress(f.read(), compresslevel=9)
            remote_path = remote_path + '.gz'
            ct = 'application/gzip'
            orig_mb = os.path.getsize(local_path) // 1024 // 1024
            comp_mb = len(data) // 1024 // 1024
            print(f"  Compressed {os.path.basename(local_path)} ({orig_mb}MB → {comp_mb}MB)")
        else:
            with open(local_path, 'rb') as f:
                data = f.read()
            ext = os.path.splitext(local_path)[1].lower()
            ct_map = {'.json': 'application/json', '.csv': 'text/csv', '.md': 'text/markdown',
                      '.png': 'image/png', '.html': 'text/html', '.txt': 'text/plain'}
            ct = ct_map.get(ext, 'application/octet-stream')

        # If data fits in one upload
        if len(data) <= CHUNK_SIZE:
            _upload_raw(remote_path, data, ct)
        else:
            # Split into chunks
            n_chunks = (len(data) + CHUNK_SIZE - 1) // CHUNK_SIZE
            print(f"  Splitting {os.path.basename(remote_path)} into {n_chunks} chunks ({len(data)//1024//1024}MB)")
            for i in range(n_chunks):
                chunk = data[i * CHUNK_SIZE:(i + 1) * CHUNK_SIZE]
                chunk_path = f"{remote_path}.chunk{i:03d}"
                _upload_raw(chunk_path, chunk, ct)
            # Upload manifest
            manifest = json.dumps({'chunks': n_chunks, 'total_size': len(data)})
            _upload_raw(f"{remote_path}.manifest", manifest.encode(), 'application/json')

    except Exception as e:
        print(f"  Upload error {remote_path}: {e}")


def _download_raw(remote_path):
    """Download raw bytes from Supabase Storage."""
    url = f"{SUPABASE_URL}/storage/v1/object/{SUPABASE_BUCKET}/{remote_path}"
    r = requests.get(url, headers=_headers())
    if r.status_code == 200:
        return r.content
    return None


def download_file(remote_path, local_path):
    """Download a single file from Supabase Storage.
    Handles chunked files and auto-decompresses .gz files."""
    if not _enabled:
        return False
    try:
        # Check if it's a chunked file (manifest exists)
        if remote_path.endswith('.manifest'):
            return False  # Skip manifest files directly

        manifest_data = _download_raw(f"{remote_path}.manifest")
        if manifest_data:
            # Chunked download
            manifest = json.loads(manifest_data)
            n_chunks = manifest['chunks']
            print(f"  Downloading {n_chunks} chunks for {os.path.basename(remote_path)}")
            data = b''
            for i in range(n_chunks):
                chunk = _download_raw(f"{remote_path}.chunk{i:03d}")
                if chunk:
                    data += chunk
                else:
                    print(f"  Missing chunk {i} for {remote_path}")
                    return False
        else:
            # Single file download
            data = _download_raw(remote_path)
            if data is None:
                return False

        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        # Decompress .gz files
        if remote_path.endswith('.gz') and not local_path.endswith('.gz'):
            with open(local_path, 'wb') as f:
                f.write(gzip.decompress(data))
        else:
            with open(local_path, 'wb') as f:
                f.write(data)
        return True
    except Exception as e:
        print(f"  Download error {remote_path}: {e}")
        return False


def delete_path(remote_path):
    """Delete a file from Supabase Storage."""
    if not _enabled:
        return
    url = f"{SUPABASE_URL}/storage/v1/object/{SUPABASE_BUCKET}"
    try:
        r = requests.delete(url, headers=_headers(), json={'prefixes': [remote_path]})
        if r.status_code not in (200, 201):
            logger.warning(f"Delete failed {remote_path}: {r.status_code}")
    except Exception as e:
        logger.error(f"Delete error {remote_path}: {e}")


def delete_prefix(prefix):
    """Delete all files under a prefix (recursive)."""
    if not _enabled:
        return
    files = list_files(prefix)
    if files:
        paths = [f"{prefix}/{f['name']}" if prefix else f['name'] for f in files]
        # Supabase delete takes a list of paths
        url = f"{SUPABASE_URL}/storage/v1/object/{SUPABASE_BUCKET}"
        try:
            r = requests.delete(url, headers=_headers(), json={'prefixes': paths})
        except Exception as e:
            logger.error(f"Delete prefix error {prefix}: {e}")


def list_files(prefix, limit=1000):
    """List files under a prefix in Supabase Storage."""
    if not _enabled:
        return []
    url = f"{SUPABASE_URL}/storage/v1/object/list/{SUPABASE_BUCKET}"
    # Split prefix into folder path and search prefix
    parts = prefix.rstrip('/').rsplit('/', 1) if '/' in prefix else ['', prefix]
    folder = parts[0] if len(parts) > 1 else ''
    search = parts[1] if len(parts) > 1 else parts[0]

    payload = {
        'prefix': folder + '/' if folder else '',
        'limit': limit,
        'offset': 0,
        'search': search if search else '',
    }
    try:
        r = requests.post(url, headers={**_headers(), 'Content-Type': 'application/json'},
                         json=payload)
        if r.status_code == 200:
            return r.json()
        else:
            logger.warning(f"List failed {prefix}: {r.status_code} {r.text[:200]}")
            return []
    except Exception as e:
        logger.error(f"List error {prefix}: {e}")
        return []


def list_folders(prefix=''):
    """List folders (case IDs) under a prefix."""
    if not _enabled:
        return []
    url = f"{SUPABASE_URL}/storage/v1/object/list/{SUPABASE_BUCKET}"
    payload = {
        'prefix': prefix + '/' if prefix else '',
        'limit': 1000,
        'offset': 0,
    }
    try:
        r = requests.post(url, headers={**_headers(), 'Content-Type': 'application/json'},
                         json=payload)
        if r.status_code == 200:
            items = r.json()
            # Folders have id=null and name only
            return [item for item in items if item.get('id') is None]
        return []
    except Exception as e:
        logger.error(f"List folders error {prefix}: {e}")
        return []


def sync_dir_to_remote(local_dir, remote_prefix):
    """Upload all files in a local directory to remote prefix (recursive)."""
    if not _enabled or not os.path.isdir(local_dir):
        return
    count = 0
    for root, dirs, files in os.walk(local_dir):
        for fname in files:
            local_path = os.path.join(root, fname)
            rel = os.path.relpath(local_path, local_dir)
            remote_path = f"{remote_prefix}/{rel}".replace('\\', '/')
            upload_file(remote_path, local_path)
            count += 1
    print(f"[Storage] Synced {count} files to {remote_prefix}")


def sync_remote_to_dir(remote_prefix, local_dir):
    """Download all files from remote prefix to local directory (recursive)."""
    if not _enabled:
        return
    _download_recursive(remote_prefix, local_dir, remote_prefix)


def _download_recursive(prefix, local_base, top_prefix):
    """Recursively download files from a prefix."""
    url = f"{SUPABASE_URL}/storage/v1/object/list/{SUPABASE_BUCKET}"
    payload = {
        'prefix': prefix + '/' if prefix else '',
        'limit': 1000,
        'offset': 0,
    }
    try:
        r = requests.post(url, headers={**_headers(), 'Content-Type': 'application/json'},
                         json=payload)
        if r.status_code != 200:
            return
        items = r.json()
        for item in items:
            name = item.get('name', '')
            item_id = item.get('id')
            full_remote = f"{prefix}/{name}" if prefix else name

            if item_id is None:
                # It's a folder, recurse
                sub_local = os.path.join(local_base, name)
                _download_recursive(full_remote, local_base, top_prefix)
            else:
                # It's a file, download
                # Skip chunk and manifest files (handled by main file download)
                if '.chunk' in name or name.endswith('.manifest'):
                    continue
                # Calculate local path relative to top_prefix
                rel = full_remote[len(top_prefix):].lstrip('/')
                # Strip .gz suffix for local path (decompress on download)
                if rel.endswith('.gz'):
                    local_path = os.path.join(local_base, rel[:-3])
                else:
                    local_path = os.path.join(local_base, rel)
                download_file(full_remote, local_path)
    except Exception as e:
        logger.error(f"Download recursive error {prefix}: {e}")


def _remote_has_data(prefix):
    """Check if Supabase Storage has any data under prefix."""
    items = list_folders(prefix)
    return len(items) > 0


def restore_all(upload_dir, results_dir, archive_dir):
    """Bidirectional sync on startup:
    - Local has data, Supabase empty → push to Supabase (initial sync)
    - Local empty, Supabase has data → pull from Supabase (restore)
    """
    if not _enabled:
        logger.info("Storage sync not enabled, skipping restore")
        return

    has_local = os.path.isdir(upload_dir) and len(os.listdir(upload_dir)) > 0
    has_remote = _remote_has_data('uploads')

    print(f"[Storage] has_local={has_local}, has_remote={has_remote}")

    if has_local and not has_remote:
        # Initial sync: push local data to Supabase
        print("[Storage] Local data found, Supabase empty → pushing to Supabase...")
        sync_dir_to_remote(upload_dir, 'uploads')
        sync_dir_to_remote(results_dir, 'results')
        if os.path.isdir(archive_dir) and os.listdir(archive_dir):
            sync_dir_to_remote(archive_dir, 'archive')
        print("[Storage] Initial sync to Supabase complete!")

    elif not has_local and has_remote:
        # Restore: pull from Supabase
        print("[Storage] Local empty, Supabase has data → restoring...")
        sync_remote_to_dir('uploads', upload_dir)
        sync_remote_to_dir('results', results_dir)
        sync_remote_to_dir('archive', archive_dir)
        print("[Storage] Restore from Supabase complete!")

    elif has_local and has_remote:
        print("[Storage] Both local and remote have data, skipping startup sync")

    else:
        print("[Storage] No data anywhere, nothing to sync")
