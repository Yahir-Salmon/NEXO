"""Reproducible real-data fetcher. Must be run in an Internet-enabled environment.
Never falls back to synthetic data or pretends an unsuccessful download succeeded.
"""
import argparse, hashlib, json, os, pathlib, sys, time
from datetime import datetime, timezone
import requests

ROOT=pathlib.Path(__file__).resolve().parents[2]
SIMA_URL='https://www.observatoriodelaire.com/api/datos/sima-horario?anio=todos&estacion=todas&parametro=todos&v=2026-09'
NASA_URL='https://power.larc.nasa.gov/api/temporal/hourly/point?parameters=T2M,RH2M,WS10M,WD10M,PRECTOTCORR&community=RE&longitude=-100.316&latitude=25.6866&start=20260101&end=20260131&format=JSON&time-standard=UTC'


def fetch(url, name, timeout=90):
    target=ROOT/'data'/'raw'/name
    if target.exists():
        raise FileExistsError(f'{target} exists: immutable raw data, do not overwrite')
    response=requests.get(url,timeout=timeout,headers={'User-Agent':'NEXO-research/0.1 (educational)'})
    response.raise_for_status()
    ctype=response.headers.get('Content-Type','')
    if len(response.content)<200: raise RuntimeError('unexpectedly short payload')
    if b'<html' in response.content[:250].lower(): raise RuntimeError('server returned HTML instead of data')
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_bytes(response.content)
    digest=hashlib.sha256(response.content).hexdigest()
    meta={'retrieved_utc':datetime.now(timezone.utc).isoformat(),'url':url,'http_status':response.status_code,
          'content_type':ctype,'bytes':len(response.content),'sha256':digest,
          'authorization':'Verify upstream data terms and attribution; the server response does not grant a license.'}
    (ROOT/'data'/'metadata'/f'{name}.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding='utf8')
    print(f'OK {target}: {len(response.content)} bytes sha256={digest}')
    return target

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('dataset',choices=['sima','nasa'])
    args=parser.parse_args()
    try: fetch(*( (SIMA_URL,'observatorio_sima_horario.csv') if args.dataset=='sima' else (NASA_URL,'nasa_power_202601.json')))
    except Exception as err: print(f'FETCH FAILED: {err}',file=sys.stderr);sys.exit(1)
