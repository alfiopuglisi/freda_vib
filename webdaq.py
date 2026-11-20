"""Minimal client for the MCC WebDAQ-504 REST API (reverse-engineered from its web UI, firmware 2.2.12).

How live data works on this device:
  * a job is started with POST /api/v1.0/schedule (a one-job schedule, like "Start job" in the web UI);
  * GET /api/v1.0/schedule/jobs/<job>/status gives samplesAcquired;
  * GET /api/v1.0/schedule/jobs/<job>/samples/<start>/<count>/bin returns interleaved float64 (little endian),
    in the channel unit (g). This only returns data when the job's logging is enabled: samples are read back
    from the log file. So the app logs into a folder of its own (/storage/data/<job>/) and wipes it on every
    (re)start and stop, because the internal storage is small.
"""
import copy
import json

import numpy as np
import requests

H_JSON = {'Content-Type': 'application/json'}

# Job status codes (job.js)
QUEUED, STARTED, WAITING_TRIG, ACQUIRING, COMPLETED, STOPPED, CANCELED, JUMPED, ERROR, TIMEOUT = range(1, 11)
RUNNING_CODES = (QUEUED, STARTED, WAITING_TRIG, ACQUIRING)


def _channel(number, name, sensitivity, color):
    return {
        'number': number, 'type': 'acceleration',
        'www': {'productId': 316, 'color': color, 'dashboardScalar': True, 'dashboardStripChart': True,
                'dashboardFFT': True, 'showFFTPeak': True, 'fftWindowType': 'none', 'fftIntegration': 'none'},
        'name': name, 'range': '±5V', 'unit': 'g', 'couplingMode': 'ac', 'iepe': True,
        'sensorSensitivity': sensitivity,
        'customScaling': {'type': 'none', 'linear': {'multiplier': 1, 'offset': 0}},
    }


def default_job(name, rate):
    """Job descriptor used when the job does not exist on the device yet (channels are set by ensure_job)."""
    return {
        'productId': 316, 'type': 'job', 'version': 1, 'name': name,
        'channels': [],
        'acquisition': {'sample': {'rate': rate}, 'startTrigger': {'type': 'immediate', 'timeout': 0},
                        'stopTrigger': {'type': 'manual', 'timeout': 0}},
        'logging': {'enable': True, 'logFile': {'name': name, 'path': f'/storage/data/{name}/', 'appendTime': True},
                    'removeEmptyLogFile': True, 'note': ''},
        'alarms': [], 'events': [], 'www': {'status': 1, 'errorText': '', 'iterationIndex': 0},
    }


class WebDAQ:
    def __init__(self, host, job='freda', timeout=5.0):
        self.host = host
        self.job = job
        self.base = f'http://{host}/api/v1.0'
        self.php = f'http://{host}/WebDAQ.php'
        self.folder = f'/storage/data/{job}'
        self.timeout = timeout
        self.s = requests.Session()

    # --- REST helpers ---
    def _get(self, path):
        r = self.s.get(self.base + path, timeout=self.timeout, headers={'Cache-Control': 'no-cache'})
        r.raise_for_status()
        return r

    def _post(self, path, obj):
        r = self.s.post(self.base + path, data=json.dumps(obj), headers=H_JSON, timeout=self.timeout)
        if r.status_code >= 400:
            raise RuntimeError(f'POST {path}: HTTP {r.status_code} {r.text[:200]}')
        return r

    def info(self):
        return self._get('/system/info').json()

    def schedule_status(self):
        return self._get('/schedule/status').json()

    def job_descriptor(self):
        r = self.s.get(f'{self.base}/jobs/{self.job}/descriptor', timeout=self.timeout)
        return r.json() if r.status_code == 200 else None

    def free_bytes(self):
        r = self.s.get(self.php, params={'functionName': 'getDiskInfo'}, timeout=self.timeout)
        for d in r.json():
            if d.get('path', '').startswith('/storage'):
                return int(d['free'])
        return None

    # --- job management ---
    def ensure_job(self, rate, channels):
        """Create/update the job so that it has the requested rate and channels and logs into its own folder.
        channels: list of (number, name, sensitivity in V/g). Other settings of existing channels are kept.
        Returns the descriptor in use; its channels are sorted by number, which is the order of the samples."""
        self.wipe_folder()
        d = self.job_descriptor()
        if d is None:
            d = default_job(self.job, rate)
        want = copy.deepcopy(d)
        old = {c['number']: c for c in want['channels']}
        cols = ['#DD3222', '#FFC000', '#3482CB', '#8E24AA']
        chans = []
        for i, (num, name, sens) in enumerate(sorted(channels)):
            c = old.get(num) or _channel(num, name, sens, cols[i % len(cols)])
            c.update({'name': name, 'sensorSensitivity': sens, 'type': 'acceleration', 'unit': 'g', 'iepe': True,
                      'couplingMode': 'ac'})
            chans.append(c)
        want['channels'] = chans
        want['acquisition']['sample']['rate'] = rate
        want['acquisition']['stopTrigger'] = {'type': 'manual', 'timeout': 0}
        want['logging']['enable'] = True
        want['logging']['logFile']['path'] = self.folder + '/'
        want['alarms'], want['events'] = [], []
        if want != d:
            self._post(f'/jobs/{self.job}/descriptor', want)
        return want

    def wipe_folder(self):
        """Delete the job's log folder (it only ever holds this app's temporary log) and recreate it empty."""
        self.s.post(self.php, data={'functionName': 'delete', 'param1': json.dumps([self.folder])},
                    timeout=max(self.timeout, 20))
        r = self.s.post(self.php, data={'functionName': 'createFolder', 'param1': self.folder}, timeout=self.timeout)
        if 'Success' not in r.text:
            raise RuntimeError(f'cannot create {self.folder} on the WebDAQ: {r.text[:100]}')

    def start(self):
        sched = {'type': 'schedule', 'version': 1, 'jobs': [self.job], 'start': {'type': 'immediate', 'timeout': 0},
                 'repeat': {'enable': False, 'interval': 0, 'count': 0}, 'stopOnJobError': True,
                 'startOnBoot': False, 'www': {'isJobStart': True, 'status': 0, 'errorText': ''}}
        self._post('/schedule', sched)

    def stop(self):
        self._post('/schedule/status', {'run': False})

    def status(self):
        return self._get(f'/schedule/jobs/{self.job}/status').json()

    def samples(self, start, count, nch):
        """Samples [start, start+count) as an (n, nch) array in the channel unit. May return fewer rows."""
        r = self._get(f'/schedule/jobs/{self.job}/samples/{int(start)}/{int(count)}/bin')
        a = np.frombuffer(r.content, dtype='<f8')
        a = a[:len(a) - len(a) % nch]
        return a.reshape(-1, nch)
