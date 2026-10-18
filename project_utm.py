"""Convert WGS84 lon/lat to metric UTM 14N for metropolitan Monterrey.
Coordinate transform cannot be used as a source of station observations.
"""
import argparse, pathlib
import pandas as pd
from pyproj import Transformer


def project(source,target):
    df=pd.read_csv(source)
    for s in ['latitude','longitude','station_id','timestamp_source','value','variable','unit']:
        if s not in df: raise ValueError(f'missing {s}; project only verified canonical REAL data')
    if not df.latitude.between(24,27).all() or not df.longitude.between(-102,-98).all():
        raise ValueError('coordinates outside declared local region')
    transformer=Transformer.from_crs('EPSG:4326','EPSG:32614',always_xy=True)
    df['x_m'],df['y_m']=transformer.transform(df.longitude.to_numpy(),df.latitude.to_numpy())
    pathlib.Path(target).parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(target,index=False)
    print(f'{len(df)} verified-source rows projected to EPSG:32614 (metric); time zone still unresolved')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source');parser.add_argument('target');p=parser.parse_args();project(p.source,p.target)
