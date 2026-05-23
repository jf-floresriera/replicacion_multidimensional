Aquí está el README completo para copiar y pegar directamente en GitHub:

Ve a → github.com/jf-floresriera/replicacion_multidimensional → clic en "Add a README" o crear archivo README.md → pega esto:

text
# Replicación Multidimensional – Capítulo 21.5

**Autor:** Jesus Enrique Flores Riera  
**Fecha:** 22 de mayo de 2026  
**Curso:** Programación SIG – Cubos Multidimensionales Vector-Raster

## Descripción

Réplica del análisis de vientos extremos del Capítulo 21 en tres lenguajes de programación: Python, R y Julia. Implementa la combinación probabilística de vientos huracanados y no huracanados bajo el supuesto de independencia estocástica:

$$p_c = 1 - (1 - p_h)(1 - p_{nh})$$

## Estructura del repositorio
replicacion_multidimensional/
├── analisis.qmd # Documento Quarto integrador principal
├── README.md
├── .gitignore
├── scripts/
│ ├── analisis.py # Ejercicio 1: Python (xarray + rioxarray + matplotlib)
│ ├── analisis.R # Script R de referencia (stars + sf + st_apply)
│ └── analisis.jl # Ejercicio 2: Julia (Rasters.jl + dot-broadcasting)
├── renders/
│ ├── analisis.html # Entregable HTML renderizado
│ └── analisis.pdf # Entregable PDF renderizado
└── figs/
└── mapa_rafaga_700.png # Mapa ráfaga combinada T=700 años

text

## Ejercicios implementados

### Ejercicio 1 – Python (xarray & GeoPandas)
- Ingesta de rásteres `.tif` con `xarray` y `rioxarray`
- Construcción de `xarray.Dataset` multidimensional con coordenada `mri` (Periodo de Retorno)
- Combinación probabilística con `apply_ufunc` de xarray
- Mapa exportado con `matplotlib` para T = 700 años

### Ejercicio 2 – Julia (Rasters.jl & DataFrames.jl)
- Carga de `.tif` con `Rasters.jl` y construcción de `RasterStack`
- Álgebra de mapas celda a celda con *dot-broadcasting* (`@.`)
- Visualización con `Plots.jl`

### Script R – Referencia del capítulo (stars & sf)
- Ingesta con `read_stars` y cubo multidimensional
- Conversión vectorial con `st_xy2sfc`
- Combinación probabilística con `st_apply`

## Consultas teóricas respondidas

1. **Filosofía del Apply** – Diferencia entre `for` tradicional vs `st_apply`/`apply_ufunc`/dot-broadcasting
2. **Cubos de Datos Vectoriales** – Cuándo preferir `stars` con geometría `sf` sobre ráster tradicional
3. **Independencia Estocástica** – Extensión del modelo a tercera/cuarta amenaza mediante cubos multidimensionales

## Renderizar el documento

```bash
cd replicacion_multidimensional
quarto render analisis.qmd --output-dir renders
```

## Referencia principal

Rodríguez Avellaneda, A. H. (2020). *Spatio-temporal analysis of extreme wind velocities for infrastructure design* [Master's Thesis, Erasmus Mundus MSc in Geospatial Technologies – Universitat Jaume I / University of Münster / NOVA IMS]. https://geocorp.co/thesis/mastergeotech/index.html
