# Solución numéricamente de una Ecuación Diferencial Parcial  


Aproximar numéricamente utilizando los elementos finitos de:
* Lagrange de primer orden
* Crouzeix–Raviart de primer orden
* Lagrange de segundo orden 

la solución débil del siguiente problema con valores en la frontera:
$$
\begin{aligned}
-\nabla \cdot(k(x) \nabla u) & =f(x) & & \text { en } \Omega:=(0,1) \times(0,1), \\
u & =0 & & \text { en } \partial \Omega,
\end{aligned}
$$


donde
$$
f(x)=8 \pi \sin \left(4 \pi x_2\right)\left(4 \pi x_1^2 \sin \left(4 \pi x_1\right)+4 \pi \sin \left(4 \pi x_1\right)-x_1 \cos \left(4 \pi x_1\right)\right), \quad k(x)=1+x_1^2$$


* Obtener entonces la tasa de convergencia para el error $e:=u_h-u$ en las normas $L^2(\Omega)$ y la seminorma $H^1(\Omega)$ para la familia de mallas proporcionada en la página del curso.



```julia
function update_mesh(mesh::Dict)
 # Obtener datos relevantes
    nb_elems = mesh["nb_elems"]           # número de elementos en la malla
    nodes = mesh["nodes"]                 # matriz de coordenadas de nodos
    elems_nodes_conn = mesh["elems_nodes_conn"]  # matriz de conectividad de elementos
    nb_nodes = mesh["nb_nodes"] # Número de nodos de la malla
    nb_faces = mesh["nb_faces"] # Número total de caras / puntos medios en las caras
    faces_nodes_conn = mesh["faces_nodes_conn"] # Matriz con conexiones nodales de las caras  
    onodes_bool = mesh["onodes_bool"]      # vector booleano que indica los nodos de contorno
    #ofaces_nodes_conn = mesh["ofaces_nodes_conn"] # Matriz con conexiones nodales de las caras externas.
    #ofaces_bool = mesh["ofaces_bool"] # Matriz booleana que indica si una cara es externa.
    
    # Calcular numero total de nuevos nodos
    nb_total = nb_nodes + nb_faces

    # Crear una matriz de conexión elems_nodes_conn_med para los nuevos nodos
    elems_nodes_conn_med = [mesh["elems_nodes_conn"][:,1:3] mesh["elems_faces_conn"].+nb_nodes]
    
    for i in 1:nb_elems
        elems_nodes_conn_med[i,1:3] = sort(elems_nodes_conn_med[i,1:3])
    end

    # Crear una matriz nodes_med para las nuevas coordenadas
    nodes_med = zeros(nb_total, 2)
    nodes_med[1:nb_nodes,:] = nodes

    index1 = mesh["faces_nodes_conn"][:,1]
    index2 = mesh["faces_nodes_conn"][:,2]

    nodes_med[nb_nodes+1:end,:] = [0.5*(nodes[index1,1]+nodes[index2,1]) 0.5*(nodes[index1,2]+nodes[index2,2])]
    nodes_med
    onodes_bool_med = [onodes_bool; mesh["ofaces_bool"]]

    
    # Crear un diccionario para el nuevo conjunto de datos de la malla
    updated_mesh = Dict{String, Any}()

    # Copiar los campos existentes en el diccionario original
    for field in keys(mesh)
            updated_mesh[field] = mesh[field]
    end
    
    # h global
    h_global = get_h_global(mesh)
    
    # Actualizar el diccionario con los nuevos datos
    updated_mesh["elems_nodes_conn_med"] = elems_nodes_conn_med
    updated_mesh["nodes_med"] = nodes_med
    updated_mesh["onodes_bool_med"] = onodes_bool_med
    updated_mesh["nb_total"] = nb_total
    updated_mesh["h_global"] = h_global
    
    return updated_mesh
end
```




    update_mesh (generic function with 1 method)




```julia
@time MSHU = update_mesh.(MSH);
```

      0.140577 seconds (554.86 k allocations: 28.968 MiB, 79.65% compilation time: 100% of which was recompilation)
    


```julia
function get_θ_functions_edp(ef::Int)
    # Funciones de forma
    if ef == 1 
        θ₁(x) = 1-x[1]-x[2]
        θ₂(x) = x[1]
        θ₃(x) = x[2]
        Θ = [θ₁ θ₂ θ₃] 
    elseif ef == 2
        θc₁(x) = 1.0-2*(1-x[1]-x[2])
        θc₂(x) = 1.0-2*x[1]
        θc₃(x) = 1.0-2*x[2]
        Θ = [θc₁ θc₂ θc₃] 
    elseif ef == 3
        θ2₁(x) = (1-x[1]-x[2])*(2*(1-x[1]-x[2])-1)
        θ2₂(x) = x[1]*(2*x[1]-1)
        θ2₃(x) = x[2]*(2*x[2]-1)
        θ2₄(x) = 4*(1-x[1]-x[2])*x[1]
        θ2₅(x) = 4*x[1]*x[2]
        θ2₆(x) = 4*(1-x[1]-x[2])*x[2]
        Θ = [θ2₁ θ2₂ θ2₃ θ2₄ θ2₅ θ2₆] 
    else
        throw(ArgumentError("Tipo de elemento finito no válido"))
    end
    return Θ
end
```




    get_θ_functions_edp (generic function with 1 method)




```julia
function get_Dθ_functions_edp(ef::Int)
    # Calcular las funciones de forma según el tipo de elemento finito
    if ef == 1
        Dθ₁(x) = [-1 -1]
        Dθ₂(x) = [1 0]
        Dθ₃(x) = [0 1]
        Dθ = [Dθ₁; Dθ₂; Dθ₃]
    elseif ef == 2
        Dcθ₁(x) = [2 2]
        Dcθ₂(x) = [-2 0]
        Dcθ₃(x) = [0 -2]
        Dθ = [Dcθ₁; Dcθ₂; Dcθ₃]
    elseif ef == 3
        Dθ2₁(x) = [1-4*(1-x[1]-x[2]) 1-4*(1-x[1]-x[2])]
        Dθ2₂(x) = [4*x[1]-1 0]
        Dθ2₃(x) = [0 4*x[2]-1] 
        Dθ2₄(x) = [4*(1-x[1]-x[2])-4*x[1] -4*x[1]]
        Dθ2₅(x) = [4*x[2] 4*x[1]]
        Dθ2₆(x) = [-4*x[2] 4*(1-x[1]-x[2])-4x[2]] 
        Dθ = [Dθ2₁; Dθ2₂; Dθ2₃; Dθ2₄; Dθ2₅; Dθ2₆]              
    else
        throw(ArgumentError("Tipo de elemento finito no válido"))
    end
    return Dθ
end
```




    get_Dθ_functions_edp (generic function with 1 method)




```julia
function solve_edp2(mesh, f, κ, n, ef)
    if ef==1
        nb = mesh["nb_nodes"]                  # número de nodos/vértices en la malla
        nb_on = .!(mesh["onodes_bool"])          # indices booleanos de nodos en el interior
        nodes = mesh["nodes"]                    # matriz con coordenadas de nodos/vértices
        elems_nodes_conn = mesh["elems_nodes_conn"]  # matriz de conectividad de elementos
        elems_conn = elems_nodes_conn[:,1:3]        
        num = 3
        itera = collect(1:num)
    elseif ef==2
        nb = mesh["nb_faces"]                  # número de caras en la malla
        nb_on = .!(mesh["ofaces_bool"])  
        nodes = mesh["nodes"]                    # matriz con coordenadas de nodos/vértices   
        elems_nodes_conn = mesh["elems_nodes_conn_med"]  # matriz de conectividad de elementos
        elems_nodes_conn = elems_nodes_conn[:,1:3]
        elems_faces_conn = mesh["elems_faces_conn"]
        elems_conn = elems_faces_conn[:,1:3]
        num = 3
        itera = [3,1,2]
    elseif ef==3
        nb = mesh["nb_total"]                  # número de nodos/vértices en la malla
        nb_on = .!(mesh["onodes_bool_med"])          # indices booleanos de nodos en el interior   
        nodes = mesh["nodes_med"]                    # matriz con coordenadas de nodos/vértices
        elems_nodes_conn = mesh["elems_nodes_conn_med"]  # matriz de conectividad de elementos
        elems_conn = elems_nodes_conn[:,1:6]
        num = 6
        itera = collect(1:num)
    end 
    
    # Creación de índices para nodos/caras en el interior del dominio
        mode_on = zeros(Int64, nb,2)           # creación de matriz de conectividad
        mode_on[:,1] = nb_on                   # llenar primera columna de matriz mode_on
    
    for k in 1:nb
        mode_on[k,2]= sum(mode_on[1:k])
    end
    # número de elementos en el interior del dominio
    nb_int = mode_on[nb,2]

    # creación de matriz global de nodos interiores
    A_int = zeros(nb_int,nb_int)
    b_int = zeros(nb_int)   
    # vector de solución
    u_sol = zeros(nb)

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
    Θ = get_θ_functions_edp(ef)
    
    # Gradientes de funciones de forma
    ∇θ = get_Dθ_functions_edp(ef)
    
    # LLenado de la matriz global
        nb_elems = mesh["nb_elems"]        # número de elementos/simplejos en la malla
     # matriz booleana de nodos en el interior del dominio
        M = nb_on[elems_conn] 
    
    # Inicializando variables
        X = similar(1.0.*ξ)
        Jk = similar(zeros(2,2))
    
    for k in 1:nb_elems # pasando sobre todos los elemementos/simplex de la malla
        T = M[k,:]*M[k,:]' # matriz booleana de conexión en el simplex k-ésimo

        # Extracción de nodos y conectividad de elementos
            elem_nodes = elems_conn[k, :]
            K = collect(nodes[elems_nodes_conn[k, 1:3], :]')
            # Definiendo la transformación Tₖ
            Jk = K[:,2:3] .- K[:,1]
            X = Jk * ξ .+ K[:,1]
            detJₖ = det(Jk)
            invJₖ = inv(Jk)   
        # Ensable 
        for i in 1:num 
                if (M[k,i] == 1)
                   I = elem_nodes[i] # indice de nodos globales
                   b_int[mode_on[I,2]] += quad_tri_loc_2(f, X, Θ[itera[i]], ξ, ω)*abs(detJₖ)
                  for j in 1:num 
                    if (T[i,j] == 1) 
                        ∇θᵢ(x) = (∇θ[itera[i],:])[1](x)*invJₖ 
                        ∇θⱼ(x) = (∇θ[itera[j],:])[1](x)*invJₖ
                        ∇θᵢ∇θⱼ(x) = dot(∇θᵢ(x), ∇θⱼ(x))
                       J = elem_nodes[j] # indice de nodos globales
                       A_int[mode_on[I,2],mode_on[J,2]] += quad_tri_loc_2(κ, X, ∇θᵢ∇θⱼ, ξ, ω)*abs(detJₖ)
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




    solve_edp2 (generic function with 1 method)




```julia
function error_global_epd2(v::Function, ∇v::Function, u_sol::Array, mesh::Dict, n::Int, ef::Int)
    # Extraer nodos y conectividad de elementos de la malla
    nb_elems = mesh["nb_elems"]           # número de elementos en la malla
    nodes = mesh["nodes"]                 # matriz de coordenadas de nodos    
    elems_nodes_conn = mesh["elems_nodes_conn_med"]
    elems_faces_conn = mesh["elems_faces_conn"]   
    
    # Generación de los puntos de la cuadratura
    quad = quadratures_triangle(n)
        ω_l = quad["ω"]
        bary_coord = quad["bary_coord"]
        multi = quad["multi"]
        g_points, ω = gauss_points(bary_coord, multi, ω_l)
        # Obtener las coordenadas locales de los puntos de la cuadratura
        ξ = g_points[2:3,:]
    
    # Extracción de nodos y conectividad de elementos
    if ef==3
        num = 6
        elems_conn = elems_nodes_conn[:, 1:6] # indice de los nodos del k-ésimo elemento
    elseif ef== 2
        num = 3
        elems_conn = elems_faces_conn[:, 1:3] # indice de los nodos del k-ésimo elemento
    elseif ef== 1
        num = 3
        elems_conn = elems_nodes_conn[:, 1:3] # indice de los nodos del k-ésimo elemento
    end
    σ_T = similar(u_sol[elems_conn])
    X = similar(1.0.*ξ)
    Jk = similar(zeros(2,2))
    
    err_L2 = 0.0
    err_H1 = 0.0
    # Cálculo del error de interpolación global en L2 y H1
    for k in 1:nb_elems
        K = collect(nodes[elems_nodes_conn[k, 1:3], :]')
        elem_nodes = elems_conn[k,:]
        # Definiendo la transformación Tₖ
        Jk = K[:,2:end] .- K[:,1]
        X = Jk * ξ .+ K[:,1]
        invJₖ = inv(Jk)
        detJₖ = det(Jk)

        # Obteniendo los grados de libertad   
        if ef ==2
            σ_T = u_sol[elem_nodes]
            σ_T = [σ_T[2], σ_T[3], σ_T[1]]
        else
            σ_T = u_sol[elem_nodes]
        end
        
        # Cálculo del error de interpolación local en L2 y H1
        err_loc_L2 = error_local_L2(v, σ_T, X, ξ, ω, ef)  
        err_loc_H1 = error_local_H1(∇v, σ_T, invJₖ, X, ξ, ω, ef)
        
        err_L2 += err_loc_L2*abs(detJₖ)
        err_H1 += err_loc_H1*abs(detJₖ)
    end
    err_L2 = sqrt(err_L2)
    err_H1 = sqrt(err_H1)
    
    return err_L2, err_H1
end
```




    error_global_epd2 (generic function with 1 method)




```julia
function refina_edp2(v::Function, ∇v::Function, f::Function, κ::Function, MSH::Array, n::Int, ef::Int )
    println("Solución numérica de una EDP\n")
    # Selección de Elemento finito
    if ef == 1
        println("Elemento finito: Lagrange Orden 1")
    elseif ef == 2
        println("Elemento finito: Crouzeix–Raviart Orden 1")
    elseif ef == 3
        println("Elemento finito: Lagrange Orden 2")
    else
        throw(ArgumentError("Tipo de elemento finito no válido"))
    end       
    
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
        u_sol  = solve_edp2(MSH[i], f, κ, n, ef);
        # Calculo del error en la norma L² y en la seminorma H¹
        L2_error_vec[i], H1_error_vec[i] = error_global_epd2(v, ∇v, u_sol, MSH[i], n, ef)
        h_global[i] = MSH[i]["h_global"]
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




    refina_edp2 (generic function with 1 method)



## Crouzeix–Raviart Orden 1


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 5, 2)
```

    Solución numérica de una EDP
    
    Elemento finito: Crouzeix–Raviart Orden 1
    Cuadratura: Gaussiana-05 de 7 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.7853244 │    1.0000000 │     7.8176018 │    1.0000000 │
    │   0.3423854 │     0.3490673 │    1.1697834 │     6.8140141 │    0.1982212 │
    │   0.1520212 │     0.1239018 │    1.4943078 │     4.4677392 │    0.6089599 │
    │   0.0792050 │     0.0344652 │    1.8459853 │     2.3953354 │    0.8993173 │
    │   0.0411692 │     0.0088404 │    1.9629611 │     1.2217941 │    0.9712265 │
    │   0.0201700 │     0.0022298 │    1.9872124 │     0.6148680 │    0.9906526 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
     35.929381 seconds (23.74 M allocations: 35.901 GiB, 17.96% gc time, 0.16% compilation time)
    


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 4, 2)
```

    Solución numérica de una EDP
    
    Elemento finito: Crouzeix–Raviart Orden 1
    Cuadratura: Gaussiana-04 de 4 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.9942653 │    1.0000000 │    16.1918205 │    1.0000000 │
    │   0.3423854 │     0.1696763 │    2.5508455 │     7.1941821 │    1.1703626 │
    │   0.1520212 │     0.0953432 │    0.8315839 │     4.5383665 │    0.6646575 │
    │   0.0792050 │     0.0264458 │    1.8500896 │     2.4060339 │    0.9155161 │
    │   0.0411692 │     0.0067340 │    1.9735054 │     1.2231848 │    0.9760145 │
    │   0.0201700 │     0.0017024 │    1.9838942 │     0.6150451 │    0.9918784 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
     32.329607 seconds (15.18 M allocations: 22.386 GiB, 13.60% gc time)
    


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 3, 2)
```

    Solución numérica de una EDP
    
    Elemento finito: Crouzeix–Raviart Orden 1
    Cuadratura: Gaussiana-03 de 3 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     1.3635098 │    1.0000000 │    18.8993423 │    1.0000000 │
    │   0.3423854 │     0.2528120 │    2.4311883 │     8.8960259 │    1.0871031 │
    │   0.1520212 │     0.0681935 │    1.8903595 │     4.8523598 │    0.8744744 │
    │   0.0792050 │     0.0161805 │    2.0753744 │     2.4471585 │    0.9875790 │
    │   0.0411692 │     0.0039293 │    2.0419281 │     1.2283624 │    0.9943712 │
    │   0.0201700 │     0.0009815 │    2.0012005 │     0.6156970 │    0.9964439 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
     26.877344 seconds (12.34 M allocations: 17.877 GiB, 13.36% gc time)
    


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 2, 2)
```

    Solución numérica de una EDP
    
    Elemento finito: Crouzeix–Raviart Orden 1
    Cuadratura: Gaussiana-02 de 3 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.8616370 │    1.0000000 │    12.7880487 │    1.0000000 │
    │   0.3423854 │     0.2062548 │    2.0626526 │     6.8412413 │    0.9024661 │
    │   0.1520212 │     0.1034131 │    0.9960085 │     4.4903015 │    0.6074458 │
    │   0.0792050 │     0.0289408 │    1.8372427 │     2.3994867 │    0.9040865 │
    │   0.0411692 │     0.0073981 │    1.9678851 │     1.2223522 │    0.9730658 │
    │   0.0201700 │     0.0018700 │    1.9840996 │     0.6149398 │    0.9911429 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
     26.409029 seconds (12.34 M allocations: 17.877 GiB, 12.73% gc time)
    


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 1, 2)
```

    Solución numérica de una EDP
    
    Elemento finito: Crouzeix–Raviart Orden 1
    Cuadratura: Gaussiana-01 de 1 punto
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     2.4454667 │    1.0000000 │    16.8702102 │    1.0000000 │
    │   0.3423854 │     0.4350511 │    2.4908530 │     6.1413189 │    1.4578575 │
    │   0.1520212 │     0.0506653 │    3.1021155 │     3.0546696 │    1.0075322 │
    │   0.0792050 │     0.0128148 │    1.9831830 │     1.4463762 │    1.0785735 │
    │   0.0411692 │     0.0032762 │    1.9677237 │     0.7115519 │    1.0234019 │
    │   0.0201700 │     0.0008002 │    2.0335835 │     0.3558518 │    0.9996925 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
     20.309239 seconds (6.66 M allocations: 8.865 GiB, 3.08% gc time)
    

## Lagrange Orden 1


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 5, 1)
```

    Solución numérica de una EDP
    
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
     13.943192 seconds (22.52 M allocations: 31.188 GiB, 48.21% gc time)
    


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 4, 1)
```

    Solución numérica de una EDP
    
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
      4.251549 seconds (14.39 M allocations: 18.157 GiB, 30.43% gc time)
    


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 3, 1)
```

    Solución numérica de una EDP
    
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
      3.597463 seconds (11.68 M allocations: 13.808 GiB, 30.01% gc time)
    


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 2, 1)
```

    Solución numérica de una EDP
    
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
      3.521531 seconds (11.68 M allocations: 13.808 GiB, 29.47% gc time)
    


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 1, 1)
```

    Solución numérica de una EDP
    
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
      2.309475 seconds (6.27 M allocations: 5.117 GiB, 26.47% gc time)
    

##  Lagrange Orden 2


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 5, 3)
```

    Solución numérica de una EDP
    
    Elemento finito: Lagrange Orden 2
    Cuadratura: Gaussiana-05 de 7 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.7322130 │    1.0000000 │    10.7551141 │    1.0000000 │
    │   0.3423854 │     0.1746292 │    2.0679684 │     5.8688186 │    0.8738808 │
    │   0.1520212 │     0.0118039 │    3.8869589 │     1.6354488 │    1.8433835 │
    │   0.0792050 │     0.0014338 │    3.0413985 │     0.4219979 │    1.9543787 │
    │   0.0411692 │     0.0001753 │    3.0317297 │     0.1071211 │    1.9779928 │
    │   0.0201700 │     0.0000214 │    3.0356819 │     0.0266886 │    2.0049497 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
     50.753511 seconds (79.29 M allocations: 10.647 GiB, 2.28% gc time)
    


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 4, 3)
```

    Solución numérica de una EDP
    
    Elemento finito: Lagrange Orden 2
    Cuadratura: Gaussiana-04 de 4 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.6743218 │    1.0000000 │     8.3992770 │    1.0000000 │
    │   0.3423854 │     0.1764792 │    1.9339388 │     3.4770391 │    1.2724058 │
    │   0.1520212 │     0.0201499 │    3.1306535 │     0.5983805 │    2.5387242 │
    │   0.0792050 │     0.0030025 │    2.7465186 │     0.1311798 │    2.1895172 │
    │   0.0411692 │     0.0003992 │    2.9109322 │     0.0310936 │    2.0768552 │
    │   0.0201700 │     0.0000499 │    3.0007831 │     0.0074728 │    2.0568961 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
     49.478052 seconds (49.72 M allocations: 9.377 GiB, 3.01% gc time)
    


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 3, 3)
```

    Solución numérica de una EDP
    
    Elemento finito: Lagrange Orden 2
    Cuadratura: Gaussiana-03 de 3 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     1.4943068 │    1.0000000 │    19.4658289 │    1.0000000 │
    │   0.3423854 │     0.1778278 │    3.0709233 │     7.1750572 │    1.4398815 │
    │   0.1520212 │     0.0083136 │    4.4188659 │     1.9126476 │    1.9074192 │
    │   0.0792050 │     0.0005418 │    3.9397440 │     0.4511283 │    2.0839615 │
    │   0.0411692 │     0.0000353 │    3.9382581 │     0.1103770 │    2.0310983 │
    │   0.0201700 │     0.0000026 │    3.7531389 │     0.0271131 │    2.0253772 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
     48.045605 seconds (39.86 M allocations: 8.938 GiB, 2.85% gc time)
    


```julia
@time refina_edp2(u_exact, Du_exact, f, κ, MSHU, 2, 3)
```

    Solución numérica de una EDP
    
    Elemento finito: Lagrange Orden 2
    Cuadratura: Gaussiana-02 de 3 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.7133601 │    1.0000000 │     7.6411378 │    1.0000000 │
    │   0.3423854 │     0.1684147 │    2.0826121 │     3.4273511 │    1.1566935 │
    │   0.1520212 │     0.0195473 │    3.1069757 │     0.7404250 │    2.2106685 │
    │   0.0792050 │     0.0029062 │    2.7497634 │     0.1696743 │    2.1255858 │
    │   0.0411692 │     0.0003845 │    2.9181991 │     0.0414110 │    2.0346831 │
    │   0.0201700 │     0.0000480 │    3.0017382 │     0.0101380 │    2.0302453 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
     49.523763 seconds (39.86 M allocations: 8.938 GiB, 1.84% gc time)
    
