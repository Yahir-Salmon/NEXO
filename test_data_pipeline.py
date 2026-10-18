import pathlib,sys
import numpy as np
import pandas as pd
import pytest
from pyproj import Transformer
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'/'preprocessing'))
from project_utm import project
from normalize_sima import normalize


def test_metric_projection_roundtrip_synthetic_fixture(tmp_path):
    # Example position in geographic bounds, NOT a monitored station.
    source=tmp_path/'input.csv';target=tmp_path/'out.csv'
    pd.DataFrame([dict(station_id='SYNTH',timestamp_source='2026-01-01 00:00:00',
                       variable='PM2.5',value=1.,unit='ug/m3',latitude=25.68,longitude=-100.31)]).to_csv(source,index=False)
    project(source,target)
    result=pd.read_csv(target)
    inverse=Transformer.from_crs('EPSG:32614','EPSG:4326',always_xy=True)
    lon,lat=inverse.transform(float(result.x_m.iloc[0]),float(result.y_m.iloc[0]))
    assert abs(lon+100.31)<1e-9 and abs(lat-25.68)<1e-9


def test_real_importer_refuses_missing_coordinates(tmp_path):
    source=tmp_path/'missing.csv'
    pd.DataFrame({'estacion':['Synthetic'],'fecha':['2026-01-01'],'valor':[23.],
                  'parametro':['PM2.5'],'unidades':['ug/m3']}).to_csv(source,index=False)
    with pytest.raises(ValueError,match='Missing fields'):
        normalize(source,tmp_path/'out.csv')
