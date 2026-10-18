"""Strict importer for an externally downloaded source: no guessed coordinates/units.
Keeps source raw untouched. Requires inspection if upstream schema differs.
"""
import argparse, pathlib, json, sys
import numpy as np
import pandas as pd
ROOT=pathlib.Path(__file__).resolve().parents[2]


def normalize(input_csv, output_csv):
    df=pd.read_csv(input_csv,encoding='utf-8-sig',low_memory=False)
    names={k:k.lower().strip().replace(' ','_') for k in df.columns}
    df=df.rename(columns=names)
    print('Source columns:',list(df.columns))
    aliases={'estacion':['nombre_estacion','estacion','station','station_name'],
             'fecha':['fecha','date','timestamp','fecha_hora','datetime'],
             'valor':['valor','value','concentracion','concentration'],
             'parametro':['clave_contaminante','nombre_contaminante','parametro','contaminante','parameter'],
             'latitud':['latitud','latitude','lat'],'longitud':['longitud','longitude','lon'],
             'unidades':['unidades','unit','units']}
    for target,possible in aliases.items():
        choice=next((s for s in possible if s in df),None)
        if choice and choice!=target:df=df.rename(columns={choice:target})
    required=['estacion','fecha','valor','parametro','latitud','longitud','unidades']
    missing=[s for s in required if s not in df]
    if missing:raise ValueError(f'Missing fields {missing}; inspect the real schema; no invented coordinates or units')
    if 'hora' in df and df['fecha'].astype(str).str.fullmatch(r'\d{4}-\d{2}-\d{2}').mean()>.95:
        df['fecha']=df['fecha'].astype(str)+' '+df['hora'].astype(str)
    df['timestamp_source']=pd.to_datetime(df['fecha'],errors='coerce')
    df['value']=pd.to_numeric(df['valor'],errors='coerce')
    df['latitude']=pd.to_numeric(df['latitud'],errors='coerce')
    df['longitude']=pd.to_numeric(df['longitud'],errors='coerce')
    # Input timezone must be explicitly reconciled with source documentation before using UTC
    out=df[['estacion','timestamp_source','parametro','value','unidades','latitude','longitude']].rename(columns={'estacion':'station_id','parametro':'variable','unidades':'unit'})
    counts={'rows':len(out),'missing_values_pct':round(100*out.value.isna().mean(),3),'invalid_coordinates':int((~out.latitude.between(24,27)|~out.longitude.between(-102,-98)).sum()),
            'invalid_dates':int(out.timestamp_source.isna().sum()),'timezone':'UNRESOLVED_SOURCE_TIME - NOT UTC'}
    valid=out.timestamp_source.notna()&out.latitude.between(24,27)&out.longitude.between(-102,-98)
    if (~valid).any(): print('Rows removed because of invalid dates/coordinates:',int((~valid).sum()))
    out=out.loc[valid].copy()
    output_csv=pathlib.Path(output_csv);output_csv.parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(output_csv,index=False)
    output_csv.with_suffix('.quality.json').write_text(json.dumps(counts,indent=2),encoding='utf8')
    print(json.dumps(counts,indent=2));return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=pathlib.Path);p.add_argument('target',type=pathlib.Path);a=p.parse_args();normalize(a.source,a.target)
