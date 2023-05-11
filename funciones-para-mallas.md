# Teoría, Práctica y Aplicaciones de los Elementos Finitos
## Tarea II
### Curso de Posgrado en Matemáticas-UNAM

> Profesor: Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx

> Alumno: Mario Rafael Nolasco Estrada. mnolasco@ciencias.unam.mx

# Funciones para mallas


```julia
# Cargamos paquetes 
using Plots, DataFrames, CSV, LinearAlgebra, Combinatorics

# Cargamos el ambiente del gráfico
plotlyjs()
```




    Plots.PlotlyJSBackend()



* Function:  `leer_archivo(ruta::AbstractString)`

    La función `leer_archivo` lee un archivo CSV en la ruta especificada y crea un diccionario con los datos leídos. La salida es un diccionario con los siguientes campos: 
    - `"nbNod"`: Un entero que representa el número de nodos.
    - `"POS"`: Una matriz de tipo `Float64` con las coordenadas de posición.
    - `"LINES"`: Una matriz de tipo `Int` que contiene la información de líneas.
    - `"TRIANGLES"`: Una matriz de tipo `Int` que contiene la información de triángulos.
    - `"PNT"`: Una matriz de tipo `Int` que contiene la información de puntos.
    

* Function: `read_mesh(msh::Dict)`

    La función read_mesh toma un diccionario msh que contiene información de una malla y devuelve otro diccionario mesh que contiene más información de la misma malla. La salida es un diccionario con los siguientes campos:
    Un diccionario con los siguientes campos:
    - `nb_nodes`  : número de nodos de la malla.
    - `nodes`  : matriz con coordenadas de los nodos de la malla.
    - `elems_nodes_conn`  : matriz con las conexiones nodales de los elementos de la malla.
    - `nb_elems`  : número de elementos de la malla.
    - `nb_ofaces`  : número de caras externas de la malla.
    - `sorted_triangles` : matriz temporal con nodos de cada triángulo ordenados.
    - `elems_faces_conn`  : matriz con conexiones entre los elementos y las caras.
    - `faces_elems_conn`  : matriz con conexiones entre las caras y los elementos.
    - `faces_nodes_conn`  : matriz con conexiones nodales de las caras.
    - `ofaces_bool`  : matriz booleana que indica si una cara es externa.
    - `ofaces_nodes_conn`  : matriz con conexiones nodales de las caras externas.
    - `nb_faces`  : número de caras de la malla.
    - `onodes_bool`  : matriz booleana que indica si un nodo está en el borde de la malla.
    
  

* Function: `plot_mesh(mesh::Dict)`

    La función plot_mesh toma un diccionario mesh que contiene información de una malla y la grafica en una figura. Esta función no tiene una salida.


```julia
"""
    leer_archivo(ruta::AbstractString)

La función `leer_archivo` lee un archivo CSV en la ruta especificada 
y crea un diccionario con los datos leídos.

# Argumentos
- `ruta::AbstractString`: La ruta del archivo CSV a leer.

# Salida
Un diccionario con los siguientes campos:
- `"nbNod"`: Un entero que representa el número de nodos.
- `"POS"`: Una matriz de tipo `Float64` con las coordenadas de posición.
- `"LINES"`: Una matriz de tipo `Int` que contiene la información de líneas.
- `"TRIANGLES"`: Una matriz de tipo `Int` que contiene la información de triángulos.
- `"PNT"`: Una matriz de tipo `Int` que contiene la información de puntos.
"""
function leer_archivo(ruta::AbstractString)
    # Cargar datos del archivo
    data = CSV.read(ruta, DataFrame)
    
    # Definir los indices de búsqueda
    idx_nbN = data[:,1].== "nbNod"
    idx_POS = data[:,1].== "POS"
    idx_LIN = data[:,1].== "LINES"
    idx_TRI = data[:,1].== "TRIANGLES"
    idx_PNT = data[:,1].== "PNT"
    
    # Definir las variables
    nbNod = convert.(Int, data[idx_nbN,2][1])
    POS = convert.(Float64,[data[idx_POS,2] data[idx_POS,3]])
    LINES = convert.(Int,[data[idx_LIN,2] data[idx_LIN,3] data[idx_LIN,4]])
    TRIANGLES = convert.(Int,[data[idx_TRI,2] data[idx_TRI,3] data[idx_TRI,4] data[idx_TRI,5]])
    PNT = convert.(Int,[data[idx_PNT,2] data[idx_PNT,3]])

    # Crear y retornar el diccionario
    return Dict("nbNod" => nbNod, 
                "POS" => POS, 
                "LINES" => LINES, 
                "TRIANGLES" => TRIANGLES, 
                "PNT" => PNT)
end
```




    leer_archivo




```julia
"""
    read_mesh(msh::Dict)

La función `read_mesh` toma un diccionario `msh` que contiene información de una malla
y devuelve otro diccionario `mesh` que contiene más información de la misma malla.

# Argumentos
-`msh::Dict` es una estructura de datos del tipo Dict{String, Any}.

# Salida
Un diccionario con los siguientes campos:
- `nb_nodes`  : número de nodos de la malla.
- `nodes`  : matriz con coordenadas de los nodos de la malla.
- `elems_nodes_conn`  : matriz con las conexiones nodales de los elementos de la malla.
- `nb_elems`  : número de elementos de la malla.
- `nb_ofaces`  : número de caras externas de la malla.
- `sorted_triangles` : matriz temporal con nodos de cada triángulo ordenados.
- `elems_faces_conn`  : matriz con conexiones entre los elementos y las caras.
- `faces_elems_conn`  : matriz con conexiones entre las caras y los elementos.
- `faces_nodes_conn`  : matriz con conexiones nodales de las caras.
- `ofaces_bool`  : matriz booleana que indica si una cara es externa.
- `ofaces_nodes_conn`  : matriz con conexiones nodales de las caras externas.
- `nb_faces`  : número de caras de la malla.
- `onodes_bool`  : matriz booleana que indica si un nodo está en el borde de la malla.
"""
function read_mesh(msh::Dict)
    # Crear estructura de malla
    mesh = Dict{String, Any}()

    # Obtener el número de nodos y las coordenadas de los nodos
    mesh["nb_nodes"] =  msh["nbNod"]
    mesh["nodes"] =  msh["POS"]
    
    # Obtener la conectividad de los elementos y el número de elementos
    mesh["elems_nodes_conn"] =  msh["TRIANGLES"]
    foo = size( msh["TRIANGLES"])
    mesh["nb_elems"] = foo[1]
    
    # Obtener el número de caras externas y ordenar las caras de cada elemento
    foo = size( msh["LINES"])
    mesh["nb_ofaces"] = foo[1]  # cara externa
    sorted_triangles = zeros(Int64, mesh["nb_elems"], 3)
    mesh["elems_faces_conn"] = zeros(Int64, mesh["nb_elems"], 3)
    mesh["faces_elems_conn"] = zeros(Int64, 0, 2)
    mesh["faces_nodes_conn"] = zeros(Int64, 0, 2)
    mesh["ofaces_bool"] = zeros(0, 2)   
    for elem_idx in 1:mesh["nb_elems"]
        sorted_triangles[elem_idx, :] = sort( msh["TRIANGLES"][elem_idx, 1:3])
        lfaces = zeros(Int64, 3, 2)
        lfaces[1, :] = sorted_triangles[elem_idx, 1:2]
        lfaces[2, :] = sorted_triangles[elem_idx, 2:3]
        lfaces[3, :] = sorted_triangles[elem_idx, 1:2:3]

        for lface_idx in 1:3
            cface = lfaces[lface_idx, :]
            
            # Buscar si la cara actual ya ha sido agregada anteriormente
            test1 = (mesh["faces_nodes_conn"])[:, 1] .== cface[1]
            test2 = (mesh["faces_nodes_conn"])[:, 2] .== cface[2]
            gfidx_v = (test1 ) .& (test2)
            if sum(gfidx_v)==0  # si la cara no ha sido agregada anteriormente
                mesh["faces_nodes_conn"] = vcat(mesh["faces_nodes_conn"], cface')
                gfaces_size = size(mesh["faces_nodes_conn"])
                mesh["elems_faces_conn"][elem_idx, lface_idx] = gfaces_size[1]
                # agregar nuevo elemento a la conectividad de la cara
                mesh["faces_elems_conn"] = vcat(mesh["faces_elems_conn"], [elem_idx 0])  
            else
                gfidx = findfirst(gfidx_v)
                mesh["elems_faces_conn"][elem_idx, lface_idx] = gfidx
                # agregar nuevo elemento a la conectividad de la cara
                mesh["faces_elems_conn"][gfidx, 2] = elem_idx
            end
        end
    end

    # Obtener las caras externas y su conectividad con los nodos
    mesh["ofaces_bool"] = mesh["faces_elems_conn"][:, 2] .== 0
    mesh["ofaces_nodes_conn"] = mesh["faces_nodes_conn"][mesh["ofaces_bool"], :]

    # Obtener el número de caras y crear un vector booleano para los nodos en el borde
    foo = size(mesh["faces_elems_conn"])
    mesh["nb_faces"] = foo[1]
    mesh["onodes_bool"] = falses(mesh["nb_nodes"])
    
    # Marcar nodos en la frontera
    for i in 1:mesh["nb_ofaces"]
        mesh["onodes_bool"][mesh["ofaces_nodes_conn"][i, 1]] = true
        mesh["onodes_bool"][mesh["ofaces_nodes_conn"][i, 2]] = true
    end

    return mesh
end
```




    read_mesh




```julia
"""
    plot_mesh(mesh::Dict)

La función `plot_mesh` recibe un diccionario `mesh` que representa una malla
compuesta por nodos y elementos. 

La función grafica la malla mediante la creación de una figura en la que se dibujan
los nodos como puntos y los triángulos como líneas que unen los nodos


# Argumentos
-`mesh` es una estructura de datos del tipo Dict{String, Any}.


# Salida
La función devuelve una figura que muestra la malla graficada.
"""
function plot_mesh(mesh::Dict)
    # extraer nodos y conectividad de elementos de la malla
    nb_elems = mesh["nb_elems"]                  # número de elementos en la malla
    nodes = mesh["nodes"]                 # matriz de coordenadas de nodos
    elems_nodes_conn = mesh["elems_nodes_conn"]  # matriz de conectividad de elementos
    onodes_bool = mesh["onodes_bool"]      # vector booleano que indica los nodos de contorno
    
    # graficar los nodos como puntos y los triángulos como líneas que unen los nodos
    p = plot()   # inicializar la figura de matplotlib
    for k in 1:nb_elems
        # coordenadas x de los nodos del k-ésimo elemento
        x = nodes[elems_nodes_conn[k,:][1:3],1]
        # coordenadas y de los nodos del k-ésimo elemento
        y = nodes[elems_nodes_conn[k,:][1:3],2]  
        # graficar los nodos del i-ésimo elemento
        plot!(x, y, color=:blue, legend=false )  
    end
    # graficar los nodos de contorno de la malla
    scatter!(nodes[onodes_bool,1],nodes[onodes_bool,2],
        color=:red, legend=false, aspect_ratio=:equal )  
    return p   # devolver la figura
end
```




    plot_mesh



## Leemos las mallas


```julia
msh_filenames = ["square_1.csv"
                 "square_2.csv"
                 "square_3.csv"
                 "square_4.csv"
                 "square_5.csv" 
                 "square_6.csv"];
```


```julia
msh = leer_archivo.(msh_filenames)

# Guardamos las mallas en un arreglo
MSH = read_mesh.(msh);
```

## Graficando una malla


```julia
plot_mesh(MSH[1])
```




    <div id="f12645a0-ea58-4007-ae2d-575c07fa36b9" style="width:600px;height:400px;"></div>
    <script>
        requirejs.config({
        paths: {
            plotly: 'https://cdn.plot.ly/plotly-2.6.3.min'
        }
    });
    require(['plotly'], function (Plotly) {

    Plotly.newPlot('f12645a0-ea58-4007-ae2d-575c07fa36b9', [
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            1.0,
            0.500000000002059,
            0.643750000000294
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y1",
        "legendgroup": "y1",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            1.0,
            1.0,
            0.647916666667245
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            0.0,
            0.0,
            0.374999999999159
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y2",
        "legendgroup": "y2",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            0.500000000002059,
            0.0,
            0.375
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            0.0,
            0.499999999998694,
            0.374999999999159
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y3",
        "legendgroup": "y3",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            0.0,
            0.0,
            0.375
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            1.0,
            1.0,
            0.643750000000294
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y4",
        "legendgroup": "y4",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            0.500000000002059,
            1.0,
            0.647916666667245
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            0.500000000002059,
            0.0,
            0.281250000001193
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y5",
        "legendgroup": "y5",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            1.0,
            1.0,
            0.718750000000608
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            0.0,
            0.0,
            0.281250000001193
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y6",
        "legendgroup": "y6",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            1.0,
            0.500000000002059,
            0.718750000000608
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            1.0,
            1.0,
            0.706249999999354
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y7",
        "legendgroup": "y7",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            0.0,
            0.500000000002059,
            0.293750000000806
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            0.499999999998694,
            1.0,
            0.706249999999354
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y8",
        "legendgroup": "y8",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            0.0,
            0.0,
            0.293750000000806
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            0.643750000000294,
            0.500000000002059,
            0.281250000001193
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y9",
        "legendgroup": "y9",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            0.647916666667245,
            1.0,
            0.718750000000608
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            0.499999999998694,
            0.706249999999354,
            0.374999999999159
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y10",
        "legendgroup": "y10",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            0.0,
            0.293750000000806,
            0.375
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            0.374999999999159,
            0.706249999999354,
            0.643750000000294
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y11",
        "legendgroup": "y11",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            0.375,
            0.293750000000806,
            0.647916666667245
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            0.0,
            0.374999999999159,
            0.281250000001193
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y12",
        "legendgroup": "y12",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            0.500000000002059,
            0.375,
            0.718750000000608
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            0.374999999999159,
            0.643750000000294,
            0.281250000001193
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y13",
        "legendgroup": "y13",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            0.375,
            0.647916666667245,
            0.718750000000608
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            0.706249999999354,
            1.0,
            0.643750000000294
        ],
        "showlegend": true,
        "mode": "lines",
        "name": "y14",
        "legendgroup": "y14",
        "line": {
            "color": "rgba(0, 0, 255, 1.000)",
            "shape": "linear",
            "dash": "solid",
            "width": 1
        },
        "y": [
            0.293750000000806,
            0.500000000002059,
            0.647916666667245
        ],
        "type": "scatter"
    },
    {
        "xaxis": "x",
        "colorbar": {
            "y": 0.513888888888889,
            "title": "",
            "len": 0.9525371828521435,
            "x": 0.9934383202099738
        },
        "yaxis": "y",
        "x": [
            0.0,
            1.0,
            1.0,
            0.0,
            0.499999999998694,
            1.0,
            0.500000000002059,
            0.0
        ],
        "showlegend": true,
        "mode": "markers",
        "name": "y15",
        "legendgroup": "y15",
        "marker": {
            "symbol": "circle",
            "color": "rgba(255, 0, 0, 1.000)",
            "line": {
                "color": "rgba(0, 0, 0, 1.000)",
                "width": 1
            },
            "size": 8
        },
        "y": [
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            0.500000000002059,
            1.0,
            0.500000000002059
        ],
        "type": "scatter"
    }
]
, {
    "showlegend": false,
    "xaxis": {
        "showticklabels": true,
        "gridwidth": 0.5,
        "tickvals": [
            0.0,
            0.30000000000000004,
            0.6000000000000001,
            0.9000000000000001,
            1.2000000000000002
        ],
        "range": [
            -0.27420160734787613,
            1.2742016073478761
        ],
        "domain": [
            0.0658209390492855,
            0.9934383202099738
        ],
        "mirror": false,
        "tickangle": 0,
        "showline": true,
        "ticktext": [
            "0.0",
            "0.3",
            "0.6",
            "0.9",
            "1.2"
        ],
        "zeroline": false,
        "tickfont": {
            "color": "rgba(0, 0, 0, 1.000)",
            "family": "sans-serif",
            "size": 11
        },
        "zerolinecolor": "rgba(0, 0, 0, 1.000)",
        "anchor": "y",
        "visible": true,
        "ticks": "inside",
        "tickmode": "array",
        "linecolor": "rgba(0, 0, 0, 1.000)",
        "showgrid": true,
        "title": "",
        "gridcolor": "rgba(0, 0, 0, 0.100)",
        "titlefont": {
            "color": "rgba(0, 0, 0, 1.000)",
            "family": "sans-serif",
            "size": 15
        },
        "tickcolor": "rgb(0, 0, 0)",
        "type": "-"
    },
    "paper_bgcolor": "rgba(255, 255, 255, 1.000)",
    "annotations": [],
    "height": 400,
    "margin": {
        "l": 0,
        "b": 20,
        "r": 0,
        "t": 20
    },
    "plot_bgcolor": "rgba(255, 255, 255, 1.000)",
    "yaxis": {
        "showticklabels": true,
        "gridwidth": 0.5,
        "tickvals": [
            0.0,
            0.25,
            0.5,
            0.75,
            1.0
        ],
        "range": [
            -0.030000000000000027,
            1.03
        ],
        "domain": [
            0.03762029746281716,
            0.9901574803149606
        ],
        "mirror": false,
        "tickangle": 0,
        "showline": true,
        "ticktext": [
            "0.00",
            "0.25",
            "0.50",
            "0.75",
            "1.00"
        ],
        "zeroline": false,
        "tickfont": {
            "color": "rgba(0, 0, 0, 1.000)",
            "family": "sans-serif",
            "size": 11
        },
        "zerolinecolor": "rgba(0, 0, 0, 1.000)",
        "anchor": "x",
        "visible": true,
        "ticks": "inside",
        "tickmode": "array",
        "linecolor": "rgba(0, 0, 0, 1.000)",
        "showgrid": true,
        "title": "",
        "gridcolor": "rgba(0, 0, 0, 0.100)",
        "titlefont": {
            "color": "rgba(0, 0, 0, 1.000)",
            "family": "sans-serif",
            "size": 15
        },
        "tickcolor": "rgb(0, 0, 0)",
        "type": "-"
    },
    "width": 600
}
);
    });
    </script>





