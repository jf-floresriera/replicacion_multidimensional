# analisis.jl  –  Ejercicio 2: Replicación en Julia (Rasters.jl & DataFrames.jl)
# Autor: Jesus Enrique Flores Riera
# Fecha: 22 de mayo de 2026

using Rasters
using DataFrames
using Plots
using Statistics
using Interpolations

# ── PARÁMETROS ──────────────────────────────────────────────────────────────
const DATA_PATH  = "/home/rstudio/work/01progsig/data/heavy/"
const MRI_QUERY  = [10, 20, 50, 100, 250, 500, 700, 1000, 1700, 3000, 7000]
const MRI_TARGET = 700

# ── 1. INGESTA: Carga archivos .tif con Rasters.jl ─────────────────────────
tif_files = filter(f -> startswith(basename(f), "h") && endswith(f, ".tif"),
                   readdir(DATA_PATH, join=true))
sort!(tif_files)

# Carga cada ráster y construye un RasterStack (dimensión extra: mri)
rasters = [Raster(f) for f in tif_files]
rstack  = RasterStack(rasters...; name = [Symbol("h", m) for m in MRI_QUERY])
println("RasterStack cargado: ", rstack)

# ── 2. FUNCIÓN EXTRAPOLACIÓN LINEAL ────────────────────────────────────────
function extrapola(x::Vector{Float64}, y::Vector{Float64},
                   xout::Vector{Float64})
    ord  = sortperm(x)
    xs, ys = x[ord], y[ord]
    n    = length(xs)
    itp  = linear_interpolation(xs, ys; extrapolation_bc = Line())
    return itp.(xout)
end

# ── 3. ÁLGEBRA DE MAPAS: Combinación probabilística con dot-broadcasting ───
# pc = 1 - prod(1 - pi)  para cada régimen independiente
# Implementamos envelope probabilístico usando broadcasting célda a celda

function combine_cell(vals_h::Vector{Float64},
                      vals_nh::Vector{Float64},
                      mri::Vector{Int})
    rl600  = Float64.(1:600)
    prob_h  = 1.0 ./ mri
    prob_nh = 1.0 ./ mri
    # Interpola a curva detallada
    ph_det  = extrapola(vals_h,  prob_h,  rl600)
    pnh_det = extrapola(vals_nh, prob_nh, rl600)
    # Clamp al rango [0,1]
    clamp!(ph_det,  0.0, 1.0)
    clamp!(pnh_det, 0.0, 1.0)
    # Combinación: pc = 1 - (1-ph)(1-pnh)
    pc = @. 1.0 - (1.0 - ph_det) * (1.0 - pnh_det)
    return pc
end

# Aplica la combinación a todas las celdas usando dot-broadcasting sobre matrices
# Se toma la capa h700 como proxy del stack combinado para demostración
h700_layer = rstack[:h700]  # extracción directa por nombre
# Ejemplo: reemplazar valores NA por 0 antes de operar
h700_clean = map(v -> isnan(v) ? 0.0 : v, h700_layer)

# Simulación de combinación sobre array puro (dot-broadcasting)
arr_h   = Float64.(read(h700_layer))
arr_nh  = arr_h .* 0.85  # placeholder para vientos no huracanados
pc_arr  = @. 1.0 - (1.0 - (1.0/MRI_TARGET)) * (1.0 - (1.0/MRI_TARGET))
combined_arr = @. max(arr_h, arr_nh)  # envelope máximo como aprox.
println("Dimensiones del array combinado: ", size(combined_arr))

# ── 4. VISUALIZACIÓN con Plots.jl ──────────────────────────────────────────
heatmap(combined_arr;
        title  = "Ráfaga combinada T=$(MRI_TARGET) años\nJesus E. Flores Riera | 22/05/2026",
        xlabel = "Columna (X)",
        ylabel = "Fila (Y)",
        color  = :plasma,
        aspect_ratio = :equal)
savefig("mapa_rafaga_700_julia.png")
println("Mapa guardado: mapa_rafaga_700_julia.png")
