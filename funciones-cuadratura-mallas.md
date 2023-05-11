# Funciones para cálculo de cuadratura en una malla de triángulos

La integral de una función escalar $v:\mathbb{R}^2 \rightarrow  \mathbb{R}$ sobre una región mallada $\Omega\subset \mathbb{R}^2 $, se puede calcular mediante la suma de integrales sobre cada elemento $K$ de una malla $T_h$.

$$\int_{\Omega} v(x,y) dxdy = \sum_{K\in T_h} \int_{K} v(x,y) dxdy = \sum_{K\in T_h} \int_{\hat{K}}v\left( T_k(\hat{x},\hat{y}) \right) \left\vert det\left(J_K\right) \right\vert d\hat{x}d\hat{y}$$

Para calcular la integral sobre un elemento $K$, se utiliza un cambio de variable para transformar la integral sobre el triángulo $\hat{K}$, que es un triángulo de referencia, al triángulo $K$. Este cambio de variable implica una transformación afín $T_K$ y una multiplicación por el determinante del Jacobiano de la transformación.

$$\int_{\hat{K}}v\left( T_k(\hat{x},\hat{y}) \right) d\hat{x}d\hat{y} \approx  \sum_{l=1}^{l_q}\omega_l\ v\left( T_k(\hat{\xi}_l,\hat{\eta}_l) \right)  $$

Para aproximar la integral sobre el triángulo de referencia $\hat{K}$, se utiliza una fórmula de cuadratura que involucra un conjunto de puntos de cuadratura y pesos, que dependen del orden de la precisión deseada.

En la siguiente tabla se muestran diferentes conjuntos de puntos de cuadratura y pesos para distintos grados de precisión. Los puntos de cuadratura se expresan en coordenadas baricéntricas, que son convenientes para trabajar con triángulos. La variable $S$ es la superficie del triángulo de referencia $\hat{K}$.

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


```julia
"""
    quadratures_triangle(n::Int)

La función `quadratures_triangle(n)` devuelve un diccionario que contiene los nodos
y pesos para una cuadratura en un triángulo de área `S`.

# Argumentos:
- `n` : Número de la cuadratura a utilizar. Debe ser un entero en el rango `1` a `5`.

# Salida:
    Un diccionario que contiene los siguientes elementos:

- `k_q`: Orden de la función polinómica exacta que se integra.
- `l_q`: Número de puntos en la cuadratura.
- `bar_coo_q`: Coordenadas baricéntricas de los puntos de cuadratura.
- `multi_q`: Multiplicidad de los puntos de cuadratura.
- `ω_l`: Pesos de los puntos de cuadratura.
"""
function quadratures_triangle(n::Int)
    # Definir las variables
    k_q = [1 2 2 3 3][n]  # grado máximo de los polinomios que se pueden integrar exactamente
    l_q = [1 3 3 4 7][n]  # número de puntos de cuadratura en la regla n
    
    # Coordenadas baricéntricas de los puntos de cuadratura. 
    barycentric_coord = Dict( "1" => [1//3 1//3 1//3],
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
    Weights_q = Dict( "1" => [1],
                "2" => [1//3],
                "3" => [1//3],
                "4" => [-9//16 25//48],
                "5" => [9//20 2//15 1//20])
    # Crear y retornar el diccionario
    bar_coo_q = barycentric_coord[string(n)]  # coordenadas baricéntricas de los puntos de cuadratura
    multi_q = Multiplicity[string(n)]  # multiplicidades de los puntos de cuadratura
    ω_l = Weights_q[string(n)]  # pesos de los puntos de cuadratura

    return Dict("k_q" => k_q, 
                "l_q" => l_q, 
                "bar_coo_q" => bar_coo_q, 
                "multi_q" => multi_q, 
                "ω_l" => ω_l)
end
```




    quadratures_triangle




```julia
function gauss_points(coord, multiplicity, weights)
    # Obtener las dimensiones de la matriz coord.
    n_points, n_dim = size(coord)
    n_q = sum(multiplicity)

    # Crear matrices vacías para almacenar los puntos y pesos de Gauss.
    gauss_coord = zeros(Float64, n_dim, n_q)
    gauss_weights = zeros(Float64, n_q)

    idx_start = 1
    for i in 1:n_points
        # Encontrar todas las permutaciones únicas de los nodos de la triada actual.
        perms = unique(permutations(coord[i,:]))
            
        # Crear una matriz para almacenar los puntos de Gauss para esta triada.
        sub_coord = zeros(Rational{Int64}, n_dim, multiplicity[i])
        
        # Generar puntos de Gauss para esta barra.
        for j in 1:multiplicity[i]
            sub_coord[:,j] = collect(perms[j])
        end
        
        idx_end = idx_start + multiplicity[i] - 1
        
        # Agregar los puntos de Gauss para esta triada a la matriz gauss_coord.
        gauss_coord[:,idx_start:idx_end] = sub_coord
        # Agregar los pesos de Gauss para esta triada a la matriz gauss_weights.        
        gauss_weights[idx_start:idx_end] .= weights[i]
        
        idx_start = idx_end + 1
    end
    
    return gauss_coord, gauss_weights
end
```




    gauss_points (generic function with 1 method)



### Integral de una función dada una malla y una regla de cuadratura:

$$
\int_{\Omega} v(x,y) dxdy 
= \sum_{K\in T_h} \int_{K} v(x,y) dxdy
= \sum_{K\in T_h} \int_{\hat{K}}v\left( T_k(\hat{x},\hat{y}) \right) \left\vert det\left(J_K\right) \right\vert d\hat{x}d\hat{y}
= \sum_{K\in T_h} \left(\sum_{l=1}^{l_q}\omega_l\ v\left( T_K(\hat{\xi}_l,\hat{\eta}_l) \right)   \right) \left\vert det\left(J_K\right) \right\vert 
$$

Sea $\hat{K} = \triangle OE_1E_2$ con nodos en el origen y los puntos $(1,0)$, $(0,1)$, entonces las coordenadas $(\hat{\xi}, \hat{\eta})$ en términos de las coordenadas baricéntricas $(\lambda_1: \lambda_2 :\lambda_3)$ son:

$$(\hat{\xi},\hat{\eta}) = \lambda_1 (0,0) +\lambda_2 (1,0) + \lambda_2 (0,1) = (\lambda_1,\lambda_2)  $$

Además si $K=\triangle Z_0Z_1Z_2 \in T_h$ es un elemento de la malla $K=\triangle Z_0Z_1Z_2$ la transformación $T_K: \hat{K} \rightarrow K$ que lleva los vértices de $\hat{K}$ en los vértices de $K$ está dada por:

$$
T_K(\hat{\xi},\hat{\eta}) = \left[ z_1-z_0 \  \ | \begin{matrix} \ \\ \ \end{matrix} \ z_2-z_0 \right] \left(\begin{matrix} \hat{\xi} \\\ \hat{\eta} \end{matrix}\right) + z_0
$$

Con $z_i$ el vector columna con las coordenadas del punto $Z_i$


```julia
"""
integrate_f_mesh(f::Function, mesh::Dict, n::Int)

Aproxima la integral de la función f en una malla triangular usando cuadratura Gaussiana.

# Argumentos
- `f`: función a integrar
- `mesh`: diccionario que contiene la información de la malla. Se espera que tenga las siguientes llaves:
  - `"nb_elems"`: número de elementos en la malla
  - `"nodes"`: matriz de coordenadas de los nodos
  - `"elems_nodes_conn"`: matriz de conectividad de los nodos de los elementos
- `n` : Número de la cuadratura a utilizar. Debe ser un entero en el rango `1` a `5`.

# Salida
- `int_glob`: valor de la integral calculada

"""
function integrate_f_mesh(f::Function, mesh::Dict, n::Int)
    # extraer nodos y conectividad de elementos de la malla
    nb_elems = mesh["nb_elems"]           # número de elementos en la malla
    nodes = mesh["nodes"]                 # matriz de coordenadas de nodos
    elems_nodes_conn = mesh["elems_nodes_conn"]  # matriz de conectividad de elementos
    
    # Generación de puntos de cuadratura
    quad = quadratures_triangle(n)
        ω_l = quad["ω_l"]
        bar_coo_q = quad["bar_coo_q"]
        multi_q = quad["multi_q"]
    g_points, w_points = gauss_points(bar_coo_q, multi_q, ω_l)
        # Seleccionar los dos primeros puntos de la matriz g_points como gpoints.
        g_points = 1.0*g_points[2:3,:]
        # Devolver los puntos y pesos de Gauss.
        w_points = 1.0*w_points

    int_glob = 0.0
    # Cálculo de la integral global
    for k in 1:nb_elems
        # Extracción de nodos y conectividad de elementos
        elem_nodes = elems_nodes_conn[k, 1:3] # coordenadas de los nodos del k-ésimo elemento
        q = collect(nodes[elem_nodes, :]')
        Jk = q[:,2:end] .- q[:,1]
        points = Jk * g_points .+ q[:,1]
    
        # Cálculo de la integral local
        int_loc = cuadratura_triang(f, points, g_points,w_points)*abs(det(Jk))
        int_glob += int_loc
    end
    return int_glob
end
```




    integrate_f_mesh




```julia
"""
cuadratura_triang(f, points, g_points, w_points)

Calcula la cuadratura de una función en un triángulo definido por los puntos 
de entrada usando el método de Gauss-Legendre.

# Argumentos
- `f::Function`: Función a integrar
- `points::Matrix`: Matriz de coordenadas de los puntos del triángulo. 
                    La matriz debe ser de tamaño (2,3) donde 
                    la primera fila corresponde a las coordenadas `x` de los puntos 
                    y la segunda fila a las coordenadas `y`.
- `g_points::Matrix`: Matriz de puntos de cuadratura. Cada columna corresponde a 
                    un punto de cuadratura en el triángulo de referencia.
- `w_points::Vector`: Vector de pesos de cuadratura.

# Salida
- `int_loc::Float64`: Valor de la integral de `f` en el triángulo definido por `points`.

"""
function cuadratura_triang(f, points, g_points, w_points)
# Evaluación de la función en los puntos de cuadratura
    fval = [f(points[:,i]) for i in 1:size(points,2)]
    
    # Cálculo de la integral local
    int_loc = 0.5 * dot(w_points, fval)
    return int_loc  
end
```




    cuadratura_triang



### Ejemplo

Si $f(x,y) = x^4 sen(x) cos(y) $, entonces:

$$\int_0^1 \int_0^1 x^4 \sin(x) \cos(y) dx dy = -20 \sin^2(1) + 24 \sin(1) - \frac{13}{2}\sin(2) \approx 0.12340199555116005 $$




```julia
# Definimos la función
f(x) = x[1]^4*sin(x[1])* cos(x[2]) # 0.12340199555116126961024429

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


