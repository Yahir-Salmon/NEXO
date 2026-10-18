# Estado comprobado de la entrega NEXO 0.1

Fecha: 18/10/2026. Todo experimento numérico de concentración está etiquetado como **SINTÉTICO**.

- **Implementado y probado en este entorno**: núcleo independiente Python, validación LOOCV, reconstrucción IDW y RBF regularizada, análisis de singularidad, gradientes, esquema conservativo periódico, convergencia; **12 pruebas de pytest aprobadas**. Figuras y métricas generadas. Sintaxis JS y lógica científica ejecutadas con DOM simulado en Node.
- **Implementado, no ejecutado en su plataforma de destino**: aplicación programática MATLAB, scripts MATLAB para simulación, proyección UTM en Python probado, scripts R para limpieza, estadísticas, regresiones/baselines, PCA y bootstrap (sin R disponible); aplicación web HTML/JS (navegador Chrome real bloqueado por administrador del entorno).
- **PDF verificado**: compilado con pdflatex + bibtex8, 14 páginas, sin advertencias en la compilación final; portada y páginas con figuras renderizadas y revisadas.
- **Documentado, pero no completado**: catálogo INECC, Observatorio del Aire, NASA POWER, OpenAQ; no se adquirieron observaciones reales de calidad del aire ni registros meteorológicos por limitaciones de conectividad. No se calculó porcentaje real de datos faltantes.
- **Pendiente de investigación empírica**: descargas reales, auditoría de calidad/unidades/licencias, configuración geográfica de estaciones, modelos de kriging, calibración meteorológica, validación temporal y espacial sobre datos reales, cuantificación de incertidumbre real, simulación física calibrada y comparación MATLAB/R ejecutada.

Comparación sintética LOOCV RMSE: IDW 2.788590 µg/m³; RBF 2.532853 µg/m³. Simulación periódica: error RMS 0.083259 a 0.021035 µg/m³ al refinar 24→96 celdas/eje. **No extrapolar estos valores a Monterrey.**
