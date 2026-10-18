# NEXO — laboratorio de reconstrucción matemática y modelación ambiental

**Versión 0.1 · 18 de octubre de 2026**

Primer prototipo de investigación científica: reconstrucción espacial, SVD, regularización de Tikhonov, validación cruzada, cálculo vectorial, conservación de masa y convergencia de un transporte-difusión periódico. Se incluyen código MATLAB, R, Python, una aplicación interactiva de navegador y un informe académico LaTeX/PDF.

> **Distinguir los datos:** todas las cifras de `results/`, las figuras, la app y la primera versión del artículo se obtuvieron con **observaciones sintéticas declaradas**. No son estimaciones observadas en Monterrey. Se identificaron fuentes públicas verificables de datos ambientales reales, pero no se pudieron descargar registros desde este contenedor; los importadores fallan de manera explícita si falta información. Ningún resultado del prototipo constituye una predicción validada del aire de Monterrey.

## Entregables y estado

| Componente | Estado de esta entrega |
|---|---|
| Exploración interactiva | `web/index.html` (JS autónomo); cómputo JS verificado con DOM de prueba, renderizado en navegador real no validado |
| App MATLAB | `matlab/app/NEXOApp.m` (uifigure, parámetros funcionales, pestañas, exportar PNG, guardar/cargar escenarios); no ejecutada por ausencia de MATLAB |
| Métodos MATLAB | IDW, RBF, SVD Tikhonov, interpolación de malla, gradientes, viento, bootstrap, EDP, convergencia; fuente disponible, no ejecutada |
| Análisis R | `R/scripts/analysis.R`, métricas y prueba de paridad; no ejecutados por ausencia de R y de datos descargados |
| Validación Python independiente | **12/12 pruebas ejecutadas satisfactoriamente** |
| Experimentos numéricos | Ejecutados en Python, semilla y tablas/figuras reproducibles |
| Pipeline de fuentes reales | Descarga HTTP con SHA-256 y esquema de datos; implementación preparada, sin descargas logradas |
| Informe PDF | `paper/main.pdf`, **14 páginas**, compilado y revisado visualmente |
| Estaciones y concentraciones reales | Pendiente (nunca sustituidas silenciosamente por datos artificiales) |

## Arranque rápido

**Solo navegador**: abre `web/index.html` desde un navegador moderno. Ajusta número de estaciones, ruido, λ, distribución, velocidad del viento, difusión y horizonte. Pulsa **Recalcular experimento** y cambia de pestaña. Exporta/importa configuraciones sintéticas en JSON. *Los datos del navegador se generan localmente; no requiere Internet.*

**Python 3.10+** (mínimos verificados: entorno actual Python 3, NumPy, SciPy, matplotlib, pandas, pytest):

```bash
python -m pip install -r requirements.txt
python -m pytest -q tests
python src/validation/run_synthetic_experiments.py
```

Los gráficos reproducibles quedarán en `results/figures/` y los datos resumidos en `results/metrics/synthetic_summary.json`.

**MATLAB R2020b o posterior** (no se requieren toolboxes comerciales en el núcleo):

```matlab
addpath(genpath('matlab'))
NEXOApp
% También:
run('matlab/run_demo.m')
run('matlab/test_core.m')
run('matlab/test_fixture_metrics.m')
```

La UI se distribuye como archivo `.m` basado en `uifigure`, **no como `.mlapp` creado con App Designer**. Algunas funciones de interfaz o exportación pueden variar con la versión. La suite en MATLAB debe ejecutarse antes de declarar compatibilidad definitiva.

**R (recomendado R >= 4.2):**

```bash
Rscript R/scripts/check_fixture.R
Rscript R/scripts/analysis.R data/processed/sima_normalized.csv results
```

El segundo comando requiere un CSV REAL normalizado, no disponible en este paquete. Si faltan variables o registros suficientes, el script no finge resultados. Se usa base R para máxima portabilidad; no hay dependencia de paquetes GIS en el script básico. La prueba compartida en `data/synthetic/` permite evaluar consistencia de RMSE/MAE entre lenguajes cuando MATLAB y R estén instalados.

**LaTeX**:

```bash
cd paper
pdflatex -interaction=nonstopmode main.tex
bibtex8 main   # o bibtex main, si está instalado
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

Necesitas una distribución TeX con `babel-spanish`, `amsmath`, `graphicx`, `microtype`, `booktabs` y dependencias comunes; ya se incluye `references.bib`. La compilación fue efectivamente ejecutada con `pdflatex` y `bibtex8`.

## Descargar datos auténticos en una máquina con Internet

```bash
python src/acquisition/fetch_sources.py sima
python src/acquisition/fetch_sources.py nasa
python src/preprocessing/normalize_sima.py data/raw/observatorio_sima_horario.csv data/processed/sima_normalized.csv
```

- Portal oficial INECC Nuevo León: [ficha datos.gob.mx](https://www.datos.gob.mx/dataset/concentraciones_validadas_calidad_aire_entidad_federativa/resource/8d674b36-46f4-4d38-a91c-6413c8b2e674), conjunto horario validado de 2022, CC BY 4.0 según el catálogo. La descarga oficial se realizará manualmente si el endpoint automatizado sigue inaccesible.
- [Observatorio del Aire](https://www.observatoriodelaire.com/datos): según su interfaz, registro de 15 estaciones SIMA, PM2.5/PM10, enero 2025 a septiembre 2026. Se documentó su enlace directo CSV en `src/acquisition/fetch_sources.py`; no se han verificado los campos internos del archivo y la licencia específica del recurso requiere revisión. **No redistribuir datos originales sin comprobar sus condiciones.**
- [NASA POWER hourly API](https://power.larc.nasa.gov/docs/services/api/temporal/hourly/): parámetros T2M/RH2M/WS10M/WD10M/PRECTOTCORR; productos meteorológicos de escala mayor, no equiparables a una medición de estación SIMA.
- [OpenAQ API v3](https://docs.openaq.org/api): documentación verificada; requiere credencial X-API-Key. No se proporciona una clave en el repositorio.

El módulo de adquisición evita sobrescribir archivos crudos existentes, registra URL, fecha, bytes y SHA-256, y marca errores de red en vez de fabricar datos. **El normalizador utiliza un contrato estricto** y puede requerir adaptación tras inspeccionar las cabeceras verdaderas de la descarga. No asigna zonas horarias falsas; `timestamp_source` permanece sin referencia UTC hasta que el proveedor documente la hora de origen. Una implementación adicional deberá almacenar CRS original WGS84/EPSG:4326 y convertir coordenadas a un CRS métrico pertinente antes de derivar campos. Tampoco se han determinado los faltantes reales.

## Primer experimento reproducible (todos los resultados son sintéticos)

Semilla 20261018; 24 ubicaciones simuladas en 1–59 km por 1–49 km; ruido N(0,1²) µg/m³; IDW con potencia 2, RBF gaussiana ancho 15 km y Tikhonov λ=0.05. Envolvente convexa para limitar mapa interpolado. Validación dejando una estación fuera (LOOCV), aparte del RMSE de malla respecto a la función sintética conocida.

| Método | RMSE LOOCV (µg/m³) | RMSE interior de hull (µg/m³) |
|---|---:|---:|
| IDW | 2.788590 | 3.103758 |
| RBF + Tikhonov SVD | 2.532853 | 2.557605 |

Número de condición de diseño RBF: aproximadamente 39435 (en esta realización). Transporte periódico u=1.5 m/s, D=10 m²/s, T=500 s; errores de convergencia con N=24/48/96: 0.083259, 0.037047 y 0.021035 µg/m³. Error relativo de masa verificado en un único paso: 0 reportado en coma flotante. Ver `results/metrics/synthetic_summary.json` para valores completos; la curva de selección λ es **exploratoria**, no una estimación anidada del error tras elegir hiperparámetros.

## Diferencias metodológicas importantes

- MATLAB y Python validan **Tikhonov de coeficientes**: `argmin ||Aβ-y||² + λ||β||²`, mediante SVD; el navegador emplea **kernel ridge** `(K+λI)α=y`. No son idénticos ni se deben esperar las mismas métricas.
- Las distribuciones espaciales sintéticas son abstractas. No se utilizaron coordenadas de estaciones reales. Nunca interpoles sobre grados lat/lon tratándolos como metros.
- La EDP implementada supone condiciones **periódicas**, constantes de transporte y difusión, y no incluye emisiones ni topografía real. No está calibrada ni validada frente a registros observados.
- RMSE de una única semilla no estima por sí solo la variación entre periodos/años o estados de la red; para esa inferencia quedan pendientes validación bloqueada, estudios de faltantes, bootstrap espacial calibrado y contraste con un modelo espacial probabilístico.

## Estructura

`src/acquisition` (descarga), `src/preprocessing` (contrato y validación), `src/integration` (algoritmos Python de referencia), `src/validation` (reproducibilidad y figuras), `matlab/` (núcleo y UI), `R/` (estadística), `web/` (UI complementaria), `tests/` (pruebas Python), `data/` (fuentes/metadata/sintéticos), `config/` (experimento), `experiments/` (logs), `results/` (métricas y figuras), `paper/` (LaTeX y PDF).

## Licencias

El código del proyecto se ofrece bajo MIT (`LICENSE`). Los datos externos conservan sus propias licencias y atribuciones. La licencia del código **no** autoriza a redistribuir datos de terceros ni imágenes externas.

## Flujo geoespacial real, preparado para cuando haya descarga

```bash
python src/preprocessing/project_utm.py data/processed/sima_normalized.csv data/processed/sima_utm.csv
```

El código establece la transformación WGS84 a UTM zona 14N (`EPSG:32614`) con `pyproj`, conservando lat/lon originales y agregando `x_m`,`y_m`. Con unidades y cobertura verificadas, en MATLAB puede llamarse `run_real_snapshot('data/processed/sima_utm.csv','PM2.5','2026-01-15')`. Es un módulo **pendiente de prueba con datos medidos** que falla explícitamente si no hay registros válidos, no una demostración de resultados reales. El informe de R incluye, cuando haya suficientes registros, resúmenes diarios, comparación de medias por estaciones, bootstrap de promedios de estación y una línea base cronológica fuera de muestra; todas estas rutinas permanecen sin ejecutar en el entorno actual.
