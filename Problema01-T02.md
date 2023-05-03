# Definimos funciones


```julia
using Plots, DataFrames, CSV, LinearAlgebra, Combinatorics
plotlyjs()
```


<div style="padding: 1em; background-color: #f8d6da; border: 1px solid #f5c6cb; font-weight: bold;">
<p>The WebIO Jupyter extension was not detected. See the
<a href="https://juliagizmos.github.io/WebIO.jl/latest/providers/ijulia/" target="_blank">
    WebIO Jupyter integration documentation
</a>
for more information.
</div>






    Plots.PlotlyJSBackend()




```julia
function leer_archivo(ruta::AbstractString)
    datos = CSV.read(ruta, DataFrame)

    # Definir las variables
    nbNod = convert.(Int, datos[datos[:,1].== "nbNod",2][1])
    POS = convert.(Float64,[datos[datos[:,1].== "POS",2] datos[datos[:,1].== "POS",3]])
    LINES = convert.(Int,[datos[datos[:,1].== "LINES",2] datos[datos[:,1].== "LINES",3] datos[datos[:,1].== "LINES",4]])
    TRIANGLES = convert.(Int,[datos[datos[:,1].== "TRIANGLES",2] datos[datos[:,1].== "TRIANGLES",3] datos[datos[:,1].== "TRIANGLES",4] datos[datos[:,1].== "TRIANGLES",5]])
    PNT = convert.(Int,[datos[datos[:,1].== "PNT",2] datos[datos[:,1].== "PNT",3]])

    # Crear y retornar el diccionario
    return Dict("nbNod" => nbNod, 
                "POS" => POS, 
                "LINES" => LINES, 
                "TRIANGLES" => TRIANGLES, 
                "PNT" => PNT)
end
```




    leer_archivo (generic function with 1 method)




```julia
"""
`mesh` es una estructura de datos del tipo Dict{String, Any}. 
Los datos se almacenan utilizando nombres clave de tipo String.

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
function read_mesh(msh)
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
            gfidx_v = ( (mesh["faces_nodes_conn"])[:, 1] .== cface[1]) .& ((mesh["faces_nodes_conn"])[:, 2] .== cface[2])
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
function plot_mesh(mesh::Dict)
    # extraer nodos y conectividad de elementos de la malla
    nb_elems = mesh["nb_elems"]                  # número de elementos en la malla
    nodes = mesh["nodes"]                 # matriz de coordenadas de nodos
    elems_nodes_conn = mesh["elems_nodes_conn"]  # matriz de conectividad de elementos
    onodes_bool = mesh["onodes_bool"]      # vector booleano que indica los nodos de contorno
    
    # graficar los nodos como puntos y los triángulos como líneas que unen los nodos
    p = plot()   # inicializar la figura de matplotlib
    for k in 1:nb_elems
        x = nodes[elems_nodes_conn[k,:][1:3],1]  # coordenadas x de los nodos del k-ésimo elemento
        y = nodes[elems_nodes_conn[k,:][1:3],2]  # coordenadas y de los nodos del k-ésimo elemento
        plot!(x, y, color=:blue, legend=false )  # graficar los nodos del i-ésimo elemento
    end
    # graficar los nodos de la malla
    #scatter!(nodes[:, 1], nodes[:, 2], 
     #   color=:blue, aspect_ratio=:equal, legend=false)  
    # graficar los nodos de contorno de la malla
    scatter!(nodes[onodes_bool,1],nodes[onodes_bool,2],
        color=:red, legend=false, aspect_ratio=:equal )  
    return p   # devolver la figura
end
```




    plot_mesh (generic function with 1 method)



### Funciones para cálculo de cuadratura

$$\int_{\Omega} v(x,y) dxdy = \sum_{K\in T_h} \int_{K} v(x,y) dxdy = \sum_{K\in T_h} \int_{\hat{K}}v\left( T_k(\hat{x},\hat{y}) \right) \left\vert det\left(J_K\right) \right\vert d\hat{x}d\hat{y}$$

$$\int_{\hat{K}}v\left( T_k(\hat{x},\hat{y}) \right) d\hat{x}d\hat{y} \approx  \sum_{l=1}^{l_q}\omega_l\ v\left( T_k(\hat{\xi}_l,\hat{\eta}_l) \right)  $$

\begin{equation}
\begin{array}{cc|cc|c}
\hline k_{\mathrm{q}} & l_{\mathrm{q}} & \text { Barycentric coord. } & \text { Multiplicity } & \text { Weights } \omega_l \\
\hline 1 & 1 & \left(\frac{1}{3}, \frac{1}{3}, \frac{1}{3}\right) & 1 & S \\
\hline 2 & 3 & \left(\frac{1}{6}, \frac{1}{6}, \frac{2}{3}\right) & 3 & \frac{1}{3} S \\
\hline 2 & 3 & \left(\frac{1}{2}, \frac{1}{2}, 0\right) & 3 & \frac{1}{3} S \\
\hline 3 & 4 & \left(\frac{1}{3}, \frac{1}{3}, \frac{1}{3}\right) & 1 & -\frac{9}{16} S \\
& & \left(\frac{1}{5}, \frac{1}{5}, \frac{3}{5}\right) & 3 & \frac{25}{48} S \\
\hline 3 & 7 & \left(\frac{1}{3}, \frac{1}{3}, \frac{1}{3}\right) & 1 & \frac{9}{20} S \\
& & \left(\frac{1}{2}, \frac{1}{2}, 0\right) & 3 & \frac{2}{15} S \\
& & (1,0,0) & 3 & \frac{1}{20} S \\
\hline
\end{array}
\end{equation}


```julia
function quadratures_triangle(n)
    # Definir las variables
    k_q = [1 2 2 3 3][n]
    l_q = [1 3 3 4 7][n]
    barycentric_coord = Dict( "1" => [1/3 1/3 1/3],
                              "2" =>  [1/6 1/6 2/3],
                              "3" =>  [1/2 1/2 0],
                              "4" => [1/3 1/3 1/3; 1/5 1/5 3/5],
                              "5" => [1/3 1/3 1/3; 1/2 1/2 0; 1 0 0])    
    Multiplicity = Dict( "1" => [1],
                         "2" => [3],
                         "3" => [3],
                         "4" => [1 3],
                         "5" => [1 3 3])
    Weights_q = Dict( "1" => [1],
                "2" => [1/3],
                "3" => [1/3],
                "4" => [-9/16 25/48],
                "5" => [9/20 2/15 1/20])
    # Crear y retornar el diccionario
    bar_coo_q = barycentric_coord[string(n)]    
    multi_q = Multiplicity[string(n)]
    ω_l = Weights_q[string(n)]

    return Dict("k_q" => k_q, 
                "l_q" => l_q, 
                "bar_coo_q" => bar_coo_q, 
                "multi_q" => multi_q, 
                "ω_l" => ω_l)
end
```




    quadratures_triangle (generic function with 1 method)




```julia
function gauss_points(bar_coo_q,multi_q,ω_l)
    m,n = size(bar_coo_q)
    l_q = sum(multi_q)
    points = zeros(2,l_q)
    w_points = zeros(l_q)
    
    for l in 1:m
        pares = []
        for i in 1:n
            for j in (i+1):n
                push!(pares, (bar_coo_q[l,i], bar_coo_q[l,j]))
                push!(pares, (bar_coo_q[l,j], bar_coo_q[l,i]))
            end
        end
        coo_q = unique(pares)
            
        points_m = zeros(2,multi_q[l])
        for r in 1:multi_q[l]
            points_m[:,r] = collect.(coo_q)[r]
        end
            
        if l == 1
            c = 1
        else
            c = sum(multi_q[1:l-1])+1
        end      
        points[:,c:multi_q[l]+c-1] = points_m
        w_points[c:multi_q[l]+c-1] .= ω_l[l]          
    end
    return points, w_points
end
```




    gauss_points (generic function with 1 method)



# Definición de la función y su gradiente


```julia
v = (x) -> cos(4*π*x[1]).*cos(4*π*x[2]).^2
∇v =  (x) -> [-4*π*sin(4*π*x[1])*cos(4*π*x[2])^2, -4*π*2*cos(4*π*x[1])*cos(4*π*x[2])*sin(4*π*x[2])  ]
```




    #3 (generic function with 1 method)



# Interpolador

Para $K \in \mathcal{T}_h$, la tripleta $\left\{K, P_K, \Sigma_K\right\}$ definida por
$$
\left\{\begin{array}{l}
K=T_K(\hat{K}) \\
P_K=\left\{\psi_K^{-1}(\hat{p}) ; \widehat{p} \in \hat{P}\right\} ; \\
\Sigma_K=\left\{\left\{\sigma_{K, i}\right\}_{1 \leq i \leq n_{s h}} ; \sigma_{K, i}(p)=\hat{\sigma}_i\left(\psi_K(p)\right), \forall p \in P_K\right\}
\end{array}\right.
$$
es un elemento finito. 

Las funciones de forma locales son $$\theta_{K, i}=\psi_K^{-1}\left(\hat{\theta}_i\right), 1 \leq i \leq n_{\mathrm{sh}}$$ y el operador de interpolación local $$\mathcal{I}_K: V(K) \longmapsto P_K$$ está dado por

$$
\mathcal{I}_K v =\sum_{i=1}^{n_{\mathrm{sh}}} \sigma_{K, i}(v) \theta_{K, i}
$$

Entonces el siguiente diagrama conmuta:
$$
\begin{array}{cc}
V(K) \stackrel{\psi_K}{\longrightarrow} V(\hat{K}) \quad \\
\downarrow{\mathcal{I}_K} \quad  \quad \quad \downarrow \mathcal{I}_{\hat{K}} \\
P_K\quad \stackrel{\psi_K}{\longrightarrow}\quad \hat{P} \quad
\end{array}
$$

$$ \mathcal{I}_K v (x) = \sum_{i=1}^{n_{\mathrm{sh}}} \sigma_{K,i}(v) \theta_{K,i}(x)$$
$$ \mathcal{I}_K v(x) = \sum_{i=1}^{n_{sh}} \sigma_{K,i}(v) \theta_{K,i}(x)=  \sum_{i=1}^{n_{sh}}  v\left(T_k\left(z_i\right)\right) \psi_k^{-1}\left(\hat{\theta}_i\right)\left(x\right) =\sum_{i=1}^{n_{sh}} v\left(T_k\left(z_i\right)\right)\left(\hat{\theta}_i\circ T_k^{-1}\right)(x)$$

$$T_k^{-1}\circ (x) = \hat{x}= J^{-1}_K (x-z_0) =\begin{pmatrix} \hat{\lambda}_1 \\ \hat{\lambda}_2 \end{pmatrix} $$
$$\hat{\theta}\circ T_k^{-1}(x) =\begin{pmatrix}1 -\hat{\lambda}_1-\hat{\lambda}_2 \\ \hat{\lambda}_1 \\ \hat{\lambda}_2 \end{pmatrix} $$

$$\Vert \mathcal{I}_h v-v\Vert_{L^2}^2 = \sum_{K \in \mathcal{T}_h} \int_{K} \left( \mathcal{I}_h v\vert_{K}-v\vert_{K}\right)^2(x) dK = \sum_{K \in \mathcal{T}_h} \int_{\hat{K}} \left|\operatorname{det}\left(J_K\right)\right| \left( \mathcal{I}_K v-v\vert_{K}\right)^2\circ T_K(\hat{x}) d\hat{K} $$ 
$$\Vert \mathcal{I}_h v-v\Vert_{L^2}^2 \approx \sum_{K \in \mathcal{T}_h}\left|\operatorname{det}\left(J_K\right)\right| \sum_{l=1}^{l_q} \omega_l\left[\left( \mathcal{I}_K v-v\right)^2 \circ T_K\right]\left(\hat{\xi}_l\right)$$

$$\int_{\Omega} v(x,y) dxdy = \sum_{K\in T_h} \int_{K} v(x,y) dxdy = \sum_{K\in T_h} \int_{\hat{K}}v\left( T_k(\hat{x},\hat{y}) \right) \left\vert det\left(J_K\right) \right\vert d\hat{x}d\hat{y}$$

$$\int_{\hat{K}}v\left( T_k(\hat{x},\hat{y}) \right) d\hat{x}d\hat{y} \approx  \sum_{l=1}^{l_q}\omega_l\ v\left( T_k(\hat{\xi}_l,\hat{\eta}_l) \right)  $$


```julia
function Iₖv(dofts_T, z_T,inv_Jₖ, x)
    z₀ = z_T[:,1]    
    x_hat = inv_Jₖ*(x-z₀)
    
    θ = [1-x_hat[1]-x_hat[2], x_hat[1], x_hat[2]]
    I_vT = dot(dofts_T,θ)
    return  I_vT    
end
```




    Iₖv (generic function with 1 method)



# Error en L2


```julia
function error_L2_triang(v, dofts_T, points, g_points, w_points, l_q)
    W = zeros(l_q)
    
    for i in 1:l_q
        x = points[:,i]
        #I_vT[i] = Iₖv(dofts_T, z_T, inv_Jₖ, x)
        θ = [1-g_points[1,i]-g_points[2,i],g_points[1,i], g_points[2,i]]
        res = v(x) - dot(dofts_T, θ)
        W[i] = res^2
    end
        # Cálculo del error local en L2
        err_loc = 0.5*dot(w_points,W)
    return  err_loc    
end
```




    error_L2_triang (generic function with 1 method)



# Error en H1


```julia
function error_H1_triang(g, dofts_T, points, g_points, w_points, l_q)
    W = zeros(l_q)
    
    for i in 1:l_q
        x = points[:,i]
        res = g(x) - [-dofts_T[1]+dofts_T[2], -dofts_T[1]+dofts_T[3]]
        W[i] = dot(res,res)
    end
        # Cálculo del error local en H1
        err_loc = 0.5*dot(w_points,W)
    return  err_loc    
end
```




    error_H1_triang (generic function with 1 method)




```julia
function interpolation_error(v, ∇v, mesh::Dict, n)
    # extraer nodos y conectividad de elementos de la malla
    nb_elems = mesh["nb_elems"]                  # número de elementos en la malla
    nodes = mesh["nodes"]                 # matriz de coordenadas de nodos
    elems_nodes_conn = mesh["elems_nodes_conn"]  # matriz de conectividad de elementos
    
    # Cuadratura Info
    quad = quadratures_triangle(n)
        l_q = quad["l_q"]
        ω_l = quad["ω_l"]
        bar_coo_q = quad["bar_coo_q"]
        multi_q = quad["multi_q"]
        g_points, w_points = gauss_points(bar_coo_q,multi_q,ω_l)

    err_L2 = 0.0
    err_H1 = 0.0
    # Cálculo del error de interpolación en L2
    for k in 1:nb_elems
        q₁ = nodes[elems_nodes_conn[k,:][1:3],1]  # coordenadas x de los nodos del k-ésimo elemento
        q₂ = nodes[elems_nodes_conn[k,:][1:3],2]  # coordenadas y de los nodos del k-ésimo elemento
        Q = [q₁ q₂]
        z_T = collect(Q')
        z₀ = z_T[:,1]
        z₁ = z_T[:,2]
        z₂ = z_T[:,3]
        Jₖ = [z₂-z₀ z₁-z₀]
        points = Jₖ*g_points.+z₀
        dofts_T = [v(z₀),v(z₁), v(z₂) ]
        
        # Cálculo del error local en L2
        err_loc_L2 = error_L2_triang(v, dofts_T, points, g_points, w_points, l_q)
        err_loc_H1 = error_H1_triang(∇v, dofts_T, points, g_points, w_points, l_q)
        
        err_L2 += err_loc_L2*abs(det(Jₖ))^2
        err_H1 += err_loc_H1*abs(det(Jₖ))^2
    end
    err_L2 = sqrt(err_L2)
    err_H1 = sqrt(err_H1)
    
    return err_L2, err_H1
end
```




    interpolation_error (generic function with 1 method)



## Leemos mallas


```julia
msh_filenames = ["square_1.csv" "square_2.csv" "square_3.csv" "square_4.csv" "square_5.csv" "square_6.csv"];
```


```julia
msh = leer_archivo.(msh_filenames)
MSH = read_mesh.(msh);
```

## Graficando una malla


```julia
plot_mesh(MSH[1])
```




    <div id="fb15e83a-07d9-4264-a553-cdc89c754ee3" style="width:600px;height:400px;"></div>
    <script>
        requirejs.config({
        paths: {
            plotly: 'https://cdn.plot.ly/plotly-2.6.3.min'
        }
    });
    require(['plotly'], function (Plotly) {

    Plotly.newPlot('fb15e83a-07d9-4264-a553-cdc89c754ee3', [
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




## Evaluar la función v(x) en los nodos del elemento
nb_nodes = mesh["nb_nodes"]
nodes = mesh["nodes"]
# Evaluar la función v(x) en los nodos del elemento
    dofts = zeros(nb_nodes)
    for i in 1:nb_nodes
        z = [nodes[i,1], nodes[i,2]]
        dofts[i] = v(z)
    end
dofts;
# Cálculo de la tasa de convergencia en L2


```julia
interpolation_error(v, ∇v, MSH[3], 4)
```




    (0.031203903602320456, 0.9292548743297034)




```julia
# Elección de cuadratura
qad = 5
# Cálculo de la tasa de convergencia
L2_error_vec = zeros(6)
H1_error_vec = zeros(6)


for i in 1:6
    L2_error_vec[i], H1_error_vec[i] = interpolation_error(v, ∇v, MSH[i], qad)
end

err_rate_L2 = zeros(6)
err_rate_L2[1] = 1

err_rate_H1 = zeros(6)
err_rate_H1[1] = 1

for i = 2:6
    err_rate_L2[i] = log(L2_error_vec[i]/L2_error_vec[i-1])/log(1/2)
    err_rate_H1[i] = log(H1_error_vec[i]/H1_error_vec[i-1])/log(1/2)
end
[L2_error_vec err_rate_L2 H1_error_vec err_rate_H1]
```




    6×4 Matrix{Float64}:
     0.307458     1.0      2.81811   1.0
     0.116558     1.39935  1.64786   0.774137
     0.0298725    1.96415  0.974476  0.757895
     0.0083061    1.84657  0.476321  1.03269
     0.00214957   1.95012  0.240495  0.985926
     0.000541727  1.98841  0.120564  0.996213




```julia

```


```julia

```


```julia

```


```julia

```


```julia

```


```julia
function plot_f_mesh(f::Function,mesh::Dict)
    # extraer nodos y conectividad de elementos de la malla
    nb_elems = mesh["nb_elems"]                  # número de elementos en la malla
    nodes = mesh["nodes"]                 # matriz de coordenadas de nodos
    elems_nodes_conn = mesh["elems_nodes_conn"]  # matriz de conectividad de elementos
    onodes_bool = mesh["onodes_bool"]      # vector booleano que indica los nodos de contorno
    
    # graficar los nodos como puntos y los triángulos como líneas que unen los nodos
    p = plot()   # inicializar la figura de matplotlib
    for k in 1:nb_elems
        x = nodes[elems_nodes_conn[k,:][1:3],1]  # coordenadas x de los nodos del k-ésimo elemento
        y = nodes[elems_nodes_conn[k,:][1:3],2]  # coordenadas y de los nodos del k-ésimo elemento
        plot!(x, y, color=:blue, legend=false )  # graficar los nodos del i-ésimo elemento
    end
    # graficar los nodos de la malla
    #scatter!(nodes[:, 1], nodes[:, 2], 
     #   color=:blue, aspect_ratio=:equal, legend=false)  
    # graficar los nodos de contorno de la malla
    scatter!(nodes[onodes_bool,1],nodes[onodes_bool,2],
        color=:red, legend=false, aspect_ratio=:equal )  
    return p   # devolver la figura
end
```




    plot_f_mesh (generic function with 1 method)




```julia

```


```julia

```


```julia

```

# Integral de una función dada una malla y una regla de cuadratura:


```julia
function integrate_v_mesh(v, mesh::Dict, n)
    # extraer nodos y conectividad de elementos de la malla
    nb_elems = mesh["nb_elems"]                  # número de elementos en la malla
    nodes = mesh["nodes"]                 # matriz de coordenadas de nodos
    elems_nodes_conn = mesh["elems_nodes_conn"]  # matriz de conectividad de elementos
    
    # Cuadratura Info
    quad = quadratures_triangle(n)
    l_q = quad["l_q"]
    ω_l = quad["ω_l"]
    bar_coo_q = quad["bar_coo_q"]
    multi_q = quad["multi_q"]
    g_points,w_points = gauss_points(bar_coo_q,multi_q,ω_l)

    int_glob = 0.0
    # Cálculo de la integral global
    for k in 1:nb_elems
        q₁ = nodes[elems_nodes_conn[k,:][1:3],1]  # coordenadas x de los nodos del k-ésimo elemento
        q₂ = nodes[elems_nodes_conn[k,:][1:3],2]  # coordenadas y de los nodos del k-ésimo elemento
        Q = [q₁ q₂]
        z_T = collect(Q')
        z₀ = z_T[:,1]
        z₁ = z_T[:,2]
        z₂ = z_T[:,3]
        Jₖ = [z₂-z₀ z₁-z₀]
        points = Jₖ*g_points.+z₀
    
        # Cálculo de la integral local
        int_loc = cuadratura_triang(v, points, g_points,w_points, l_q)
        int_glob += int_loc*abs(det(Jₖ))
    end
    return int_glob
end
```




    integrate_v_mesh (generic function with 1 method)




```julia
function cuadratura_triang(f, points, g_points, w_points, l_q)
    
    V = zeros(l_q)
    for i in 1:l_q
        xᵢ = [points[1,i],points[2,i]]
        V[i] = f(xᵢ) 
    end

    # Cálculo de la integral local
    int_loc = 0.5*sum(w_points.*V)
    return  int_loc    
end
```




    cuadratura_triang (generic function with 1 method)




```julia
f(x) = x[1]*x[2]^2+x[1]^3
f(x) = x[1]^4*sin(x[1])* cos(x[2])
```




    f (generic function with 1 method)




```julia
integrate_v_mesh(f,MSH[6],5)
```




    0.12340199569053975




```julia
0, 1234019955511612696102442956932126157168484339941699241312167270
```
