import json
import os
import re
import tempfile
import time
import uuid
from pathlib import Path

import shutil
from apscheduler.schedulers.background import BackgroundScheduler
from flask import current_app, request


class FileSessionStore:
    def __init__(self, root, ttl_seconds=1800):
        self.root = Path(root)
        self.ttl_seconds = ttl_seconds

    def _path(self, sid):
        return self.root / f"{sid}.json"

    def _write(self, sid, payload):
        tmp = tempfile.NamedTemporaryFile('w', delete=False, dir=str(self.root))
        try:
            json.dump(payload, tmp)
            tmp.flush()
            os.fsync(tmp.fileno())
        finally:
            tmp.close()
        os.replace(tmp.name, self._path(sid))

    def create(self, meta=None):
        sid = uuid.uuid4().hex
        payload = {
            'created_at': time.time(),
            'updated_at': time.time(),
            'ttl': self.ttl_seconds,
            'meta': meta or {},
            'data': {}
        }
        self._write(sid, payload)
        return sid

    def create_with_id(self, sid, meta=None):
        payload = {
            'created_at': time.time(),
            'updated_at': time.time(),
            'ttl': self.ttl_seconds,
            'meta': meta or {},
            'data': {}
        }
        self._write(sid, payload)
        return sid

    def get_payload(self, sid):
        path = self._path(sid)
        if not path.exists():
            return None
        try:
            payload = json.loads(path.read_text())
        except Exception:
            return None
        ttl = int(payload.get('ttl', self.ttl_seconds))
        updated_at = float(payload.get('updated_at', 0))
        if time.time() - updated_at > ttl:
            self.delete(sid)
            return None
        return payload

    def update(self, sid, **kv):
        payload = self.get_payload(sid)
        if payload is None:
            return False
        payload['updated_at'] = time.time()
        payload.setdefault('data', {}).update(kv)
        self._write(sid, payload)
        return True

    def read(self, sid, key):
        payload = self.get_payload(sid)
        if payload is None:
            return None
        return payload.get('data', {}).get(key)

    def delete(self, sid):
        try:
            self._path(sid).unlink()
        except FileNotFoundError:
            pass

    def purge_expired(self):
        now = time.time()
        for f in self.root.glob('*.json'):
            try:
                payload = json.loads(f.read_text())
                ttl = int(payload.get('ttl', self.ttl_seconds))
                updated_at = float(payload.get('updated_at', 0))
                if now - updated_at > ttl:
                    f.unlink()
            except Exception:
                pass


def init_store(app):
    store_dir = Path(app.config['PEPMAP_STORE_DIR'])
    store_dir.mkdir(parents=True, exist_ok=True)
    store = FileSessionStore(store_dir, ttl_seconds=app.config['PEPMAP_STORE_TTL_SECONDS'])
    app.extensions['pepmap_store'] = store
    return store


def get_store():
    return current_app.extensions['pepmap_store']


def normalize_session_id(raw_session_id):
    if raw_session_id is None:
        return None
    session_id = raw_session_id.strip()
    if not session_id:
        return None
    if not re.match(r'^[A-Za-z0-9_-]{8,128}$', session_id):
        return None
    return session_id


def get_request_session_id():
    if request.method == 'GET':
        raw_session_id = request.args.get('session_id')
    else:
        raw_session_id = request.form.get('session_id')
    if raw_session_id is None:
        return None, None
    session_id = normalize_session_id(raw_session_id)
    if session_id is None:
        return None, 'Invalid session_id format.'
    return session_id, None


def clear_store_dir(store_dir):
    store_path = Path(store_dir)
    if not store_path.exists():
        return
    for filename in os.listdir(store_path):
        file_path = store_path / filename
        try:
            if file_path.is_file() or file_path.is_symlink():
                file_path.unlink()
            elif file_path.is_dir():
                shutil.rmtree(file_path)
        except Exception:
            pass


def start_scheduler(app, store):
    scheduler = BackgroundScheduler()
    scheduler.add_job(store.purge_expired, 'cron', hour=1)
    scheduler.start()
    app.extensions['pepmap_scheduler'] = scheduler
