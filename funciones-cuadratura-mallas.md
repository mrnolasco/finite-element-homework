# Funciones de cuadratura en una malla de triángulos

La integral de una función escalar $v:\mathbb{R}^2 \rightarrow  \mathbb{R}$ sobre una región mallada $\Omega\subset \mathbb{R}^2 $, se puede calcular mediante la suma de integrales sobre cada elemento $K$ de una malla $T_h$:

$$\int_{\Omega} v(x) dx = \sum_{K\in T_h} \int_{K} v(x) dx$$

Para calcular la integral sobre un elemento $K$, se utiliza un cambio de variable para transformar la integral sobre el triángulo $\hat{K}$, que es un triángulo de referencia, al triángulo $K$. Este cambio de variable implica una transformación afín $T_K$ y una multiplicación por el determinante del Jacobiano de la transformación:
$$\int_{\Omega} v(x) dx = \sum_{K\in T_h} \int_{K} v(x) dx = \sum_{K\in T_h} \int_{\hat{K}}v\circ T_K(\hat{x}) \left\vert det\left(J_K\right) \right\vert d\hat{x}$$
donde $$x=T_K\left(\hat{x}\right) $$ 
y es tal que

$$ K=T_K\left(\hat{K}\right) $$

Ahora para aproximar la integral sobre el triángulo de referencia $\hat{K}$, se utiliza una fórmula de cuadratura que involucra un conjunto de puntos de cuadratura y pesos, que dependen del orden de la precisión deseada:
$$\int_{\hat{K}}v\circ T_K(\hat{x})  d\hat{x} \approx  \sum_{l=1}^{l_q}\omega_l\ v\circ T_K(\hat{\xi}_l)   $$

$$\therefore
\int_{\Omega} v(x) dx \approx \sum_{K\in T_h} \left(\sum_{l=1}^{l_q}\omega_l\ v\circ T_K(\hat{\xi}_l) \right)    \left\vert det\left(J_K\right) \right\vert 
$$

En la siguiente tabla se muestran diferentes conjuntos de puntos de cuadratura y pesos para distintos grados de precisión. 

Los puntos de cuadratura se expresan en coordenadas baricéntricas, que son convenientes para trabajar con triángulos. La variable $S$ es la superficie del triángulo de referencia $\hat{K}$.

$$
\begin{array}{ccc|cc|c}
\hline n & k_{\mathrm{q}} & l_{\mathrm{q}} & \text { Barycentric coord. } & \text { Multiplicity } & \text { Weights } \omega_l \\
\hline 1 & 1 & 1 & \left(\frac{1}{3}, \frac{1}{3}, \frac{1}{3}\right) & 1 & S \\
\hline 2 & 2 & 3 & \left(\frac{1}{6}, \frac{1}{6}, \frac{2}{3}\right) & 3 & \frac{1}{3} S \\
\hline 3 & 2 & 3 & \left(\frac{1}{2}, \frac{1}{2}, 0\right) & 3 & \frac{1}{3} S \\
\hline 4 & 3 & 4 & \left(\frac{1}{3}, \frac{1}{3}, \frac{1}{3}\right) & 1 & -\frac{9}{16} S \\
& & & \left(\frac{1}{5}, \frac{1}{5}, \frac{3}{5}\right) & 3 & \frac{25}{48} S \\
\hline 5 & 3 & 7 & \left(\frac{1}{3}, \frac{1}{3}, \frac{1}{3}\right) & 1 & \frac{9}{20} S \\
& & & \left(\frac{1}{2}, \frac{1}{2}, 0\right) & 3 & \frac{2}{15} S \\
& & & (1,0,0) & 3 & \frac{1}{20} S \\
\hline
\end{array}
$$

Ern, A., & Guermond, J. L. (2019). Theory and Practice of Finite Elements. Springer. Chapter 8. Quadratures, Assembling, and Storage. Table 8.2. Nodes and weights for quadratures on a triangle of area S, p 360.

Observaciones:
1. Sea $\hat{K} = \triangle OE_1E_2$ con nodos en $(0,0)$, $(1,0)$, $(0,1)$ y $K=\triangle Z_0Z_1Z_2 \in T_h$ un elemento de la malla. Entonces la transformación $T_K: \hat{K} \rightarrow K$ que lleva los vértices de $\hat{K}$ en los vértices de $K$ está dada por: 
$$x=T_K(\hat{x}) = J_K \hat{x} + z_0 \quad \hat{x}\in \hat{K} $$
Con 
$$J_K : = \left[ z_1-z_0 \  \  \vert \ \ z_2-z_0 \right]$$ y $z_i$ el vector columna con las coordenadas del punto $Z_i$

2. Si $\hat{\xi}_l$, un punto de la cuadratura, tiene coordenadas baricéntricas $(\lambda_0: \lambda_1 :\lambda_2)$, entonces 

    - Sus coordenadas en el simplejo de referencia $\hat{K}$ son: 
    
    - $$\hat{\xi}_l = \lambda_0 \begin{pmatrix} 0\\\ 0\end{pmatrix} +\lambda_1 \begin{pmatrix}1\\\ 0\end{pmatrix} + \lambda_2 \begin{pmatrix} 0\\\ 1\end{pmatrix} = \begin{pmatrix}\lambda_1\\\ \lambda_2\end{pmatrix}  $$
    - Sus coordenadas en el simplejo $K$ son: $$ T_K(\hat{\xi}_l) = J_K \hat{\xi}_l + z_0   $$

* Function: `quadratures_triangle(n::Int)`


```julia
"""
## `quadratures_triangle(n::Int)`

La función `quadratures_triangle` devuelve un diccionario que contiene los nodos
y pesos para una cuadratura en un triángulo de área uno.

## Argumentos:
- `n` : Número de la cuadratura a utilizar. Debe ser un entero en el rango `1` a `5`.

## Salida:
- Un diccionario que contiene los siguientes elementos:
    - `k`: Orden de la función polinómica exacta que se integra.
    - `l`: Número de puntos en la cuadratura.
    - `bary_coord`: Coordenadas baricéntricas de los puntos de cuadratura.
    - `multi`: Multiplicidad de los puntos de cuadratura.
    - `ω`: Pesos de los puntos de cuadratura.
"""
function quadratures_triangle(n::Int)
    # Definir las variables
    
    # grado máximo de los polinomios que se pueden integrar exactamente
    K = [1 2 2 3 3][n]  

    # número de puntos de cuadratura en la regla n
    L = [1 3 3 4 7][n]  

    # Coordenadas baricéntricas de los puntos de cuadratura. 
    Barycentric_Coord = Dict( "1" => [1//3 1//3 1//3],
                              "2" =>  [1//6 1//6 2//3],
                              "3" =>  [1//2 1//2 0//1],
                              "4" => [1//3 1//3 1//3; 1//5 1//5 3//5],
                              "5" => [1//3 1//3 1//3; 1//2 1//2 0//1; 1//1 0//1 0//1])    
    
    # Multiplicidades de los puntos de cuadratura
    Multiplicity = Dict( "1" => [1],
                         "2" => [3],
                         "3" => [3],
                         "4" => [1 3],
                         "5" => [1 3 3])
    
    # Pesos de los puntos de cuadratura
    Weights = Dict( "1" => [1],
                    "2" => [1//3],
                    "3" => [1//3],
                    "4" => [-9//16 25//48],
                    "5" => [9//20 2//15 1//20])
    
    # Crear y retornar el diccionario
        # coordenadas baricéntricas de los puntos de cuadratura
        bary_coord = Barycentric_Coord[string(n)]

        # multiplicidades de los puntos de cuadratura
        multi = Multiplicity[string(n)] 
        
        # pesos de los puntos de cuadratura
        ω = Weights[string(n)]  

    return Dict("k" => K, 
                "l" => L, 
                "bary_coord" => bary_coord, 
                "multi" => multi, 
                "ω" => ω)
end
```




    quadratures_triangle



* Function: `gauss_points(bary_coord::Array, multi::Array, weights::Array)`


```julia
"""
## `gauss_points(bary_coord::Array, multi::Array, weights::Array)`

La función `gauss_points(bary_coord, multi, weights)` calcula los puntos 
y pesos de Gauss para la cuadratura en un triángulo.

## Argumentos:
- `coord`: matriz que contiene las triadas baricéntricas de la cuadratura. 
- `multi`: vector que especifica la multiplicidad de cada triada baricéntrica.
- `weights`: vector que especifica los pesos de cada triada baricéntrica.

## Salida:
- `gauss_coord`: matriz que contiene las coordenadas baricéntricas de los nodos de Gauss.
- `gauss_weights`: vector que contiene los pesos correspondientes a los nodos de Gauss.
"""
function gauss_points(bary_coord::Array, multi::Array, weights::Array)
    # Obtener las dimensiones de la matriz bary_coord.
    n_points, n_dim = size(bary_coord)
    n_q = sum(multi)

    # Crear matrices vacías para almacenar los puntos y pesos de Gauss.
    gauss_coord = zeros(Rational{Int64}, n_dim, n_q)
    gauss_weights = zeros(Rational{Int64}, n_q)

    idx_start = 1
    for i in 1:n_points
        # Encontrar todas las permutaciones únicas de los nodos de la triada actual.
        perms = unique(permutations(bary_coord[i,:]))
            
        # Crear una matriz para almacenar los puntos de Gauss para esta triada.
        sub_coord = zeros(Rational{Int64}, n_dim, multi[i])
        
        # Generar puntos de Gauss para esta barra.
        for j in 1:multi[i]
            sub_coord[:,j] = perms[j]
        end
        
        idx_end = idx_start + multi[i] - 1
        
        # Agregar los puntos de Gauss para esta triada a la matriz gauss_coord.
        gauss_coord[:,idx_start:idx_end] = sub_coord
        # Agregar los pesos de Gauss para esta triada a la matriz gauss_weights.        
        gauss_weights[idx_start:idx_end] .= weights[i]
        
        idx_start = idx_end + 1
    end
    
    return gauss_coord, gauss_weights
end
```




    gauss_points



## Integral de una función dada una malla y una regla de cuadratura:

* Function: `quadrature_triangle_loc(f::Function, points::Matrix, g_points::Matrix, w_points::Vector)`


```julia
"""
## `quadrature_triangle_loc(f::Function, points::Matrix, g_points::Matrix, w_points::Vector)`

La función `cuadratura_triang_loc` calcula la cuadratura de una función en un triángulo
definido por los nodos de entrada usando el método de Gauss-Legendre.

## Argumentos
- `f`: Función a integrar
- `points`: Matriz de coordenadas globales de nodos de cuadratura.
            Cada columna corresponde a un nodo de cuadratura en el triángulo de la malla.
- `g_points`: Matriz de coordenadas locales de nodos de cuadratura.
            Cada columna corresponde a un nodo de cuadratura en el triángulo de referencia.
- `w_points`: Vector que contiene los pesos correspondientes a los nodos de la cuadratura.

## Salida
- `int_loc`: Valor de la integral de `f` en el triángulo.

"""
function quadrature_triangle_loc(f, points, g_points, w_points)
# Evaluación de la función en los puntos de cuadratura
    fval = [f(points[:,i]) for i in 1:size(points,2)]
    
    # Cálculo de la integral local
    int_loc = 0.5 * dot(w_points, fval)
    return int_loc  
end
```




    quadrature_triangle_loc



* Function: `integrate_f_mesh(f::Function, mesh::Dict, n::Int)`


```julia
"""
## `integrate_f_mesh(f::Function, mesh::Dict, n::Int)`

Aproxima la integral de la función `f` en una malla triangular usando cuadratura Gaussiana.

## Argumentos
- `f`: función a integrar
- `mesh`: diccionario que contiene la información de la malla con las siguientes llaves:
  - `"nb_elems"`: número de elementos en la malla
  - `"nodes"`: matriz de coordenadas de los nodos
  - `"elems_nodes_conn"`: matriz de conectividad de los nodos de los elementos
- `n` : Número de la cuadratura a utilizar. Debe ser un entero en el rango `1` a `5`.

## Salida
- `int_glob`: aproximación al valor de la integral de `f`

"""
function integrate_f_mesh(f::Function, mesh::Dict, n::Int)
    # extraer nodos y conectividad de elementos de la malla
    nb_elems = mesh["nb_elems"]           # número de elementos en la malla
    nodes = mesh["nodes"]                 # matriz de coordenadas de nodos
    elems_nodes_conn = mesh["elems_nodes_conn"]  # matriz de conectividad de elementos
    
    # Generación de los puntos de la cuadratura
    quad = quadratures_triangle(n)
        ωl = quad["ω"]
        bary_coord = quad["bary_coord"]
        multi = quad["multi"]
    g_points, w_points  = gauss_points(bary_coord, multi, ωl)
    
    # Obtener las coordenadas locales de los puntos de la cuadratura
    g_points = g_points[2:3,:]

    # Cálculo de la integral global
    int_glob = 0.0
    for k in 1:nb_elems
        # Extracción de nodos y conectividad de elementos
        elem_nodes = elems_nodes_conn[k, 1:3] # coordenadas de los nodos del k-ésimo elemento
        q = collect(nodes[elem_nodes, :]')
        Jk = q[:,2:end] .- q[:,1]
        points = Jk * g_points .+ q[:,1]
    
        # Cálculo de la integral local
        int_loc = quadrature_triangle_loc(f, points, g_points,w_points)*abs(det(Jk))
        int_glob += int_loc
    end
    return int_glob
end
```




    integrate_f_mesh



### Ejemplo

Si $f(x,y) = x^4 sen(x) cos(y) $, entonces:

$$\int_0^1 \int_0^1 x^4 \sin(x) \cos(y) dx dy = -20 \sin^2(1) + 24 \sin(1) - \frac{13}{2}\sin(2) \approx 0.12340199555116005 $$




```julia
# Definimos la función
f(x) = x[1]^4*sin(x[1])* cos(x[2]) 

# Definimos la malla a utilizar
meshf = MSH[6]

# Definimos la cuadratura a utilizar
n_quad = 5;

# Aproximación de la integral
cuadf = integrate_f_mesh(f, meshf, n_quad)
```




    0.12340199569053976




```julia
# Valor exacto de la integral
intf = -20*sin(1)^2+24*sin(1)-13*0.5*sin(2)
```




    0.12340199555116005




```julia
# Error absoluto de la aproximación
abs(cuadf-intf)
```




    1.3937971610200606e-10

