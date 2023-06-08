# Teoría, Práctica y Aplicaciones de los Elementos Finitos
## Tarea II
### Curso de Posgrado en Matemáticas-UNAM

> Profesor: Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx

> Alumno: Mario Rafael Nolasco Estrada. mnolasco@ciencias.unam.mx

# Solución numéricamente de una Ecuación Diferencial Parcial
Aproximar numéricamente utilizando los elementos finitos de Lagrange de primer orden la solución débil del siguiente problema con valores en la frontera:

$$
\begin{aligned}
-\nabla \cdot(k(x) \nabla u) & =f(x) & & \text { en } \Omega:=(0,1) \times(0,1), \\
u & =0 & & \text { en } \partial \Omega,
\end{aligned}
$$


donde

$$
f(x)=8 \pi \sin \left(4 \pi x_2\right)\left(4 \pi x_1^2 \sin \left(4 \pi x_1\right)+4 \pi \sin \left(4 \pi x_1\right)-x_1 \cos \left(4 \pi x_1\right)\right), \quad k(x)=1+x_1^2
$$

* Verificar que $u(x)=\sin \left(4 \pi x_1\right) \sin \left(4 \pi x_2\right)$ es la solución del problema. 

* Obtener entonces la tasa de convergencia para el error $e:=u_h-u$ en las normas $L^2(\Omega)$ y la seminorma $H^1(\Omega)$ para la familia de mallas proporcionada en la página del curso.

## La formulación débil del problema:

Multiplicando la ecuación diferencial por una función de prueba $v(x)$

$$
-\nabla \cdot(\kappa(x) \nabla u) v(x)=f(x) v(x)
$$

Integrando por partes en el dominio $\Omega$:

$$
\begin{array}{c}
\int_{\Omega} \nabla \cdot(\kappa(x) \nabla u) v(x) \mathrm{d} x=\int_{\Omega} f(x) v(x) \mathrm{d} x \\
-\int_{\partial \Omega} \kappa(x) \nabla u \cdot \mathbf{n} v(x) \mathrm{d} s+\int_{\Omega} \kappa(x) \nabla u \cdot \nabla v(x) \mathrm{d} x=\int_{\Omega} f(x) v(x) \mathrm{d} x
\end{array}
$$

Aplicando la condición de frontera homogénea $u=0$ en $\partial \Omega$:

$$
\int_{\Omega} \kappa(x) \nabla u \cdot \nabla v(x) \mathrm{d} x=\int_{\Omega} f(x) v(x) \mathrm{d} x
$$


## Aproximación numérica de la solución


```julia
κ(x) = 1+x[1]^2
f(x) = 8*π*sin(4π*x[2])*( 4*π*x[1]^2*sin(4π*x[1]) + 4*π*sin(4π*x[1]) - x[1]*cos(4π*x[1]) )
```




    f (generic function with 1 method)



## Solución de la EDP


```julia
function quad_tri_loc_2(f, X, Θ, ξ, ω)
# Evaluación de la función en los puntos de cuadratura
    fval = [f( X[:,i]) for i in 1:size(X,2)]
    Θval = [Θ(ξ[:,i]) for i in 1:size(ξ,2)]
    
    # Cálculo de la integral local
    int_loc = 0.5 * dot(ω, fval.*Θval)
    return int_loc  
end
```




    quad_tri_loc_2 (generic function with 1 method)




```julia
function solve_edp(mesh, f, κ, n)
    nb_nodes = mesh["nb_nodes"]                  # número de nodos/vértices en la malla
    mode_on = zeros(Int64, nb_nodes,2)           # creación de matriz de conectividad
    nb_on = .!(mesh["onodes_bool"])              # indices booleanos de nodos en el interior
    mode_on[:,1] = nb_on                         # llenar primera columna de matriz mode_on
    nb_elems = mesh["nb_elems"]                  # número de elementos/simplejos en la malla
    nodes = mesh["nodes"]                        # matriz con coordenadas de nodos/vértices
    elems_nodes_conn = mesh["elems_nodes_conn"]  # matriz de conectividad de elementos

    # vector de solución
    u_sol = zeros(nb_nodes)

    # llenar segunda columna de matriz mode_on
    for k in 1:nb_nodes
        mode_on[k,2]= sum(mode_on[1:k])
    end

    # número de elementos en el interior del dominio
    nb_nodes_int = mode_on[nb_nodes,2]

    # creación de matriz global de nodos interiores
    A_int = zeros(nb_nodes_int,nb_nodes_int)
    b_int = zeros(nb_nodes_int);

    # matriz booleana de nodos en el interior del dominio
    M = nb_on[elems_nodes_conn[:,1:3]]
    
    # Generación de los puntos de la cuadratura
        quad = quadratures_triangle(n)
            l = quad["l"]
            ω_l = quad["ω"]
            bary_coord = quad["bary_coord"]
            multi = quad["multi"]
            g_points, ω = gauss_points(bary_coord, multi, ω_l)
            # Obtener las coordenadas locales de los puntos de la cuadratura
            ξ = 1.0*g_points[2:3,:]
        # Funciones de forma
            θ₁(x) = (1-x[1]-x[2])
            θ₂(x) = x[1]
            θ₃(x) = x[2]
            Θ = [θ₁ θ₂ θ₃] 
    
        # Gradientes de funciones de forma
            ∇θ = [-1 -1; 
                   1  0;
                   0  1]
    
    for k in 1:nb_elems # pasando sobre todos los elemementos/simplex de la malla
        T = M[k,:]*M[k,:]' # matriz booleana de conexión en el simplex k-ésimo
        
        # Extracción de nodos y conectividad de elementos
            elem_nodes = elems_nodes_conn[k, 1:3] # indice de los nodos del k-ésimo elemento
            K = collect(nodes[elem_nodes, :]')

            # Definiendo la transformación Tₖ
            Jk = K[:,2:3] .- K[:,1]
            X = Jk * ξ .+ K[:,1]
            detJₖ = det(Jk)
            invJₖ = inv(Jk)   

        # Ensable Local
        for i in 1:3 
                if (M[k,i] == 1)
                    I = elem_nodes[i] # indice de nodos globales
                    b_int[mode_on[I,2]] += quad_tri_loc_2(f, X, Θ[i], ξ, ω)*abs(detJₖ)  
                  for j in 1:3
                    if (T[i,j] == 1) 
                        ∇θᵢ = ∇θ[i,:]'*invJₖ
                        ∇θⱼ = ∇θ[j,:]'*invJₖ
                        ∇θᵢ∇θⱼ = dot(∇θᵢ, ∇θⱼ)
                        J = elem_nodes[j] # indice de nodos globales
                        # Ensable Global
                        A_int[mode_on[I,2],mode_on[J,2]] += quad_tri_loc(κ, X , ω)*∇θᵢ∇θⱼ*abs(detJₖ)
                    end
                end
            end
        end
    end
    u_int = A_int\b_int
    u_sol[nb_on] = u_int;
    return u_sol
end
```




    solve_edp (generic function with 1 method)



## Solución Exacta y Gradiente


```julia
u_exact(x) = sin(4π*x[1])*sin(4π*x[2])
Du_exact(x) = [4π*cos(4π*x[1])*sin(4π*x[2]) , 4π*sin(4π*x[1])*cos(4π*x[2])]
```




    Du_exact (generic function with 1 method)



## Cálculo del Error Global de la Solución


```julia
function error_global_epd(v::Function, ∇v::Function, u_sol::Array, mesh::Dict, n::Int)
    # Extraer nodos y conectividad de elementos de la malla
    nb_elems = mesh["nb_elems"]           # número de elementos en la malla
    nodes = mesh["nodes"]                 # matriz de coordenadas de nodos
    elems_nodes_conn = mesh["elems_nodes_conn"]  # matriz de conectividad de elementos
    
    # Generación de los puntos de la cuadratura
    quad = quadratures_triangle(n)
        ω_l = quad["ω"]
        bary_coord = quad["bary_coord"]
        multi = quad["multi"]
        g_points, ω = gauss_points(bary_coord, multi, ω_l)
        # Obtener las coordenadas locales de los puntos de la cuadratura
        ξ = g_points[2:3,:]

    err_L2 = 0.0
    err_H1 = 0.0
    # Cálculo del error de interpolación global en L2 y H1
    for k in 1:nb_elems
        # Extracción de nodos y conectividad de elementos
        elem_nodes = elems_nodes_conn[k, 1:3] # coordenadas de los nodos del k-ésimo elemento
        K = collect(nodes[elem_nodes, :]')
        
        # Definiendo la transformación Tₖ
        Jk = K[:,2:end] .- K[:,1]
        X = Jk * ξ .+ K[:,1]
        invJₖ = inv(Jk)
        detJₖ = det(Jk)
        
        # Obteniendo los grados de libertad
        σ_T = u_sol[elem_nodes]
        
        # Cálculo del error de interpolación local en L2 y H1
        err_loc_L2 = error_local_L2(v, σ_T, X, ξ, ω, 1)  
        err_loc_H1 = error_local_H1(∇v, σ_T, invJₖ, X, ξ, ω, 1)
        
        err_L2 += err_loc_L2*abs(detJₖ)
        err_H1 += err_loc_H1*abs(detJₖ)
    end
    err_L2 = sqrt(err_L2)
    err_H1 = sqrt(err_H1)
    
    return err_L2, err_H1
end
```




    error_global_epd (generic function with 1 method)



## Refinamiento


```julia
function refina_edp(v::Function, ∇v::Function, MSH::Array, n::Int)
    # Elemento finito
        println("Elemento finito: Lagrange Orden 1")
    
    # Selección de Cuadratura
        if n == 1
            println("Cuadratura: Gaussiana-01 de 1 punto")
        elseif n == 2
            println("Cuadratura: Gaussiana-02 de 3 puntos")
        elseif n == 3
            println("Cuadratura: Gaussiana-03 de 3 puntos")
        elseif n == 4
            println("Cuadratura: Gaussiana-04 de 4 puntos")
        elseif n == 5
            println("Cuadratura: Gaussiana-05 de 7 puntos")
        else
            throw(ArgumentError("Elección de cuadratura no válida"))
        end    
    
    # Inicializar vector con error en L² para cada refinamiento
    L2_error_vec = zeros(6)
    # Inicializar vector con error en H¹ para cada refinamiento
    H1_error_vec = zeros(6)
    # Inicializar vector con parámetro global h
    h_global = zeros(6)
            
    # Ciclos de refinamiento
    for i in 1:6
        # Cálculo de solución numérica
        u_sol  = solve_edp(MSH[i],f,κ, n);
        # Calculo del error en la norma L² y en la seminorma H¹
        L2_error_vec[i], H1_error_vec[i] = error_global_epd(v, ∇v, u_sol, MSH[i], n)
        h_global[i] = get_h_global(MSH[i])
    end
    
    # Inicializar vector con tasa de convergencia en L²
    err_rate_L2 = zeros(6)
    err_rate_L2[1] = 1
            
    # Inicializar vector con tasa de convergencia en H¹
    err_rate_H1 = zeros(6)
    err_rate_H1[1] = 1

    for i = 2:6
        # Cálculo de la tasa de convergencia en L² y en H¹
        err_rate_L2[i] = log(L2_error_vec[i]/L2_error_vec[i-1])/log(1/2)
        err_rate_H1[i] = log(H1_error_vec[i]/H1_error_vec[i-1])/log(1/2)
    end
            
    # Creación de tabla con resultados 
    data = [ h_global L2_error_vec err_rate_L2 H1_error_vec err_rate_H1]    
    table= pretty_table(data; formatters = ft_printf("%1.7f") ,
        header = (["Parámetro h", "Error L²_norm", "Rate L²_norm", 
                    "Error H¹_norm", "Rate H¹_norm"]))
    
    return table
end
```




    refina_edp (generic function with 1 method)



## Presentación de Resultados


```julia
@time refina_edp(u_exact, Du_exact, MSH, 5)
```
    Elemento finito: Lagrange Orden 1
    Cuadratura: Gaussiana-05 de 7 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.7470386 │    1.0000000 │     8.5358844 │    1.0000000 │
    │   0.3423854 │     0.4516211 │    0.7260698 │     7.6397173 │    0.1600214 │
    │   0.1520212 │     0.1643357 │    1.4584667 │     4.5371990 │    0.7517173 │
    │   0.0792050 │     0.0444655 │    1.8858881 │     2.4067575 │    0.9147112 │
    │   0.0411692 │     0.0113692 │    1.9675567 │     1.2267077 │    0.9722993 │
    │   0.0201700 │     0.0028422 │    2.0000276 │     0.6148821 │    0.9964098 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
      1.831627 seconds (6.12 M allocations: 1005.700 MiB)
    


```julia
@time refina_edp(u_exact, Du_exact, MSH, 4)
```
    Elemento finito: Lagrange Orden 1
    Cuadratura: Gaussiana-04 de 4 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.7650589 │    1.0000000 │     9.1633256 │    1.0000000 │
    │   0.3423854 │     0.3828708 │    0.9987131 │     7.4438317 │    0.2998258 │
    │   0.1520212 │     0.1369886 │    1.4828024 │     4.5883176 │    0.6980802 │
    │   0.0792050 │     0.0379160 │    1.8531751 │     2.4170930 │    0.9246923 │
    │   0.0411692 │     0.0097760 │    1.9554961 │     1.2280833 │    0.9768646 │
    │   0.0201700 │     0.0024476 │    1.9978714 │     0.6150588 │    0.9976121 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
      2.067812 seconds (5.15 M allocations: 922.463 MiB, 24.10% gc time)

    


```julia
@time refina_edp(u_exact, Du_exact, MSH, 3)
```

    Elemento finito: Lagrange Orden 1
    Cuadratura: Gaussiana-03 de 3 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     1.2051372 │    1.0000000 │     9.3915216 │    1.0000000 │
    │   0.3423854 │     0.4468801 │    1.4312379 │     8.0776386 │    0.2174253 │
    │   0.1520212 │     0.1720198 │    1.3773126 │     4.5578195 │    0.8255898 │
    │   0.0792050 │     0.0461160 │    1.8992370 │     2.4097274 │    0.9194739 │
    │   0.0411692 │     0.0117522 │    1.9723367 │     1.2270990 │    0.9736183 │
    │   0.0201700 │     0.0029390 │    1.9995557 │     0.6149296 │    0.9967584 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
      1.360753 seconds (4.83 M allocations: 891.954 MiB, 5.28% gc time)


```julia
@time refina_edp(u_exact, Du_exact, MSH, 2)
```

    Elemento finito: Lagrange Orden 1
    Cuadratura: Gaussiana-02 de 3 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.7162942 │    1.0000000 │    10.1441457 │    1.0000000 │
    │   0.3423854 │     0.3859097 │    0.8922889 │     7.3620846 │    0.4624612 │
    │   0.1520212 │     0.1381731 │    1.4817866 │     4.5823493 │    0.6840269 │
    │   0.0792050 │     0.0382249 │    1.8538920 │     2.4159954 │    0.9234697 │
    │   0.0411692 │     0.0098531 │    1.9558596 │     1.2279379 │    0.9763802 │
    │   0.0201700 │     0.0024664 │    1.9981523 │     0.6150405 │    0.9974843 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
      1.382713 seconds (4.83 M allocations: 891.957 MiB, 5.19% gc time)    


```julia
@time refina_edp(u_exact, Du_exact, MSH, 1)
```
  
      Elemento finito: Lagrange Orden 1
    Cuadratura: Gaussiana-01 de 1 punto
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     1.9254634 │    1.0000000 │    15.5624702 │    1.0000000 │
    │   0.3423854 │     0.5675281 │    1.7624419 │     8.2209190 │    0.9206995 │
    │   0.1520212 │     0.2218581 │    1.3550545 │     3.4292011 │    1.2614272 │
    │   0.0792050 │     0.0618449 │    1.8429099 │     1.4969976 │    1.1958006 │
    │   0.0411692 │     0.0159408 │    1.9559304 │     0.7239684 │    1.0480733 │
    │   0.0201700 │     0.0039895 │    1.9984467 │     0.3563881 │    1.0224776 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
      1.331759 seconds (4.18 M allocations: 835.108 MiB, 5.34% gc time)
    
    
