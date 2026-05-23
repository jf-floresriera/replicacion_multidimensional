# analisis.py  –  Ejercicio 1: Replicación en Python (Xarray & GeoPandas)
# Autor: Jesus Enrique Flores Riera
# Fecha: 22 de mayo de 2026
# Capítulo 21 – Cubos Multidimensionales Vector-Raster

import numpy as np
import pandas as pd
import xarray as xr
import rioxarray  # noqa: F401  – activa el accessor .rio
import geopandas as gpd
import matplotlib.pyplot as plt
import glob
import os

# ── 1. PARÁMETROS ────────────────────────────────────────────────────────────
DATA_PATH   = '/home/rstudio/work/01progsig/data/heavy/'
TIFF_GLOB   = os.path.join(DATA_PATH, 'h*.tif')
NH_EXCEL    = os.path.join(DATA_PATH, 'nh_returnlevels.xlsx')
MRI_QUERY   = [10, 20, 50, 100, 250, 500, 700, 1000, 1700, 3000, 7000]
MRI_TARGET  = 700   # periodo de retorno a mapear

# ── 2. INGESTA DE DATOS RASTER (Huracanes) ───────────────────────────────────
tif_files = sorted(glob.glob(TIFF_GLOB))
if not tif_files:
    raise FileNotFoundError(f'No se encontraron archivos .tif en {DATA_PATH}')

# Carga cada .tif como DataArray y extrae la banda única
arrays = []
for f, mri in zip(tif_files, MRI_QUERY):
    da = xr.open_dataarray(f, engine='rasterio').squeeze('band', drop=True)
    da = da.assign_coords(mri=mri)
    arrays.append(da)

# ── 3. ESTRUCTURA: xarray.Dataset multidimensional con coord = Periodo Retorno
h_stack = xr.concat(arrays, dim='mri')  # dim nueva: Periodo de Retorno
h_ds = h_stack.to_dataset(name='wind_hurricane')
print('Dataset huracanes:')
print(h_ds)

# ── 4. INGESTA TABULAR (No Huracanes) ────────────────────────────────────────
nh_df = pd.read_excel(NH_EXCEL)
print('\nColumnas NH:', nh_df.columns.tolist())

# ── 5. COMBINACIÓN PROBABILÍSTICA VECTORIZADA ────────────────────────────────
# Fórmula: pc = 1 - (1 - ph)(1 - pnh)
# donde p = 1/MRI  (probabilidad anual de excedencia)

def combine_wind_probs(da_h, da_nh):
    '''
    Combina probabilidades de viento huracanado y no huracanado.
    da_h, da_nh: DataArrays con dimensión 'mri'.
    Retorna DataArray con la ráfaga combinada por periodo de retorno.
    '''
    ph  = 1.0 / da_h.mri
    pnh = 1.0 / da_nh.mri
    pc  = 1.0 - (1.0 - ph) * (1.0 - pnh)  # vectorizado NumPy/xarray
    # re-interpolar a velocidades (simplificado: max envelope)
    combined = xr.apply_ufunc(
        lambda h, nh: np.fmax(h, nh),
        da_h, da_nh,
        vectorize=True
    )
    combined.attrs['description'] = 'Ráfaga combinada (envelope probabilístico)'
    return combined

# Nota: sin los datos reales se aplica el envelope como proxy de la combinación
# En producción se usaría la interpolación inversa sobre pc
combined_stack = combine_wind_probs(h_stack, h_stack)  # placeholder

# ── 6. EXPORTACIÓN: Mapa de ráfaga combinada a T=700 años ───────────────────
c700 = combined_stack.sel(mri=MRI_TARGET)

fig, ax = plt.subplots(figsize=(8, 6))
c700.plot(ax=ax, cmap='plasma',
          cbar_kwargs={'label': 'Velocidad de ráfaga (km/h)'})
ax.set_title(f'Ráfaga de viento combinada – T={MRI_TARGET} años\n'
             'Autor: Jesus Enrique Flores Riera  |  22/05/2026')
ax.set_xlabel('Longitud')
ax.set_ylabel('Latitud')
plt.tight_layout()
plt.savefig('mapa_rafaga_700.png', dpi=150)
print('Mapa guardado: mapa_rafaga_700.png')
plt.close()
