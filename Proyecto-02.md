# Ejercicio Allen-Chan

Aproximar numéricamente utilizando los elementos finitos de Lagrange de primer
orden la solucion débil del problema estacionario de Allen-Chan con valores en la frontera:
$$
\begin{cases} - \Delta u -u+u^3=f(x) & \text{ en } \Omega:= (0,1)\times (0,1)\\\ u=0 & \text{ en } \partial \Omega  \end{cases} 
$$

Usando $u(x)= \sin\left(4\pi x_1\right) \sin\left(4\pi x_2\right) $ como solución, entonces calcular $f(x)$.

$$
\begin{array}{l}\frac{\partial u}{\partial x_1}=4 \pi \cos \left(4 \pi x_1\right) \sin \left(4 \pi x_2\right), \\ \frac{\partial u}{\partial x_2}=4 \pi \sin \left(4 \pi x_1\right) \cos \left(4 \pi x_2\right) .\end{array}
$$

$$
\begin{array}{l}\frac{\partial^2 u}{\partial x_1^2}=-(4 \pi)^2 \sin \left(4 \pi x_1\right) \sin \left(4 \pi x_2\right) \\ \frac{\partial^2 u}{\partial x_2^2}=-(4 \pi)^2 \sin \left(4 \pi x_1\right) \sin \left(4 \pi x_2\right)\end{array}
$$

$$
f(x) = \sin\left(4\pi x_1\right) \sin\left(4\pi x_2\right) \left( 2(4 \pi)^2-1 + \sin^2\left(4\pi x_1\right) \sin^2\left(4\pi x_2\right) \right)
$$


Observar que el problema es una EDP no lineal, por lo tanto se tendrá que resolver un
sistema no-lineal discreto, para ello entonces utilizar el algoritmo de Newton (ver libro de
Burden-Faires seccion 10.2). Obtener entonces la tasa de convergencia para el error en las normas $L^2\left(\Omega\right)$ y $H^1\left(\Omega\right)$ para la familia de mallas proporcionada en la pagina del curso.

## Forma débil del problema

Encontrar $u \in H_0^1(\Omega)$ tal que para todo $v \in H_0^1(\Omega)$ se cumpla:

$$\int_{\Omega} \nabla u \cdot \nabla v d x-\int_{\Omega} u v d x+\int_{\Omega} u^3 v d x=\int_{\Omega} f(x) v d x$$

Aquí, $H_0^1(\Omega)$ denota el espacio de funciones de prueba $v$ que se anulan en la frontera $\partial \Omega$ y tienen una derivada en $L^2(\Omega)$.


```julia
f(x) = sin(4π*x[1])*sin(4π*x[2])*(2*(4π)^2-1+(sin(4π*x[1])*sin(4π*x[2]))^2)
```




    f (generic function with 1 method)




```julia
function solve_edp3(mesh, f, n)
    nb_nodes = mesh["nb_nodes"]                  # número de nodos/vértices en la malla
    mode_on = zeros(Int64, nb_nodes,2)           # creación de matriz de conectividad
    nb_on = .!(mesh["onodes_bool"])              # indices booleanos de nodos en el interior
    mode_on[:,1] = nb_on                         # llenar primera columna de matriz mode_on
    nb_elems = mesh["nb_elems"]                  # número de elementos/simplejos en la malla
    nodes = mesh["nodes"]                        # matriz con coordenadas de nodos/vértices
    elems_nodes_conn = mesh["elems_nodes_conn"]  # matriz de conectividad de elementos

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
    
    # inicialización de vector solución
    u_k = ones(nb_nodes)
    u_k[.!nb_on] .= 0.0
    
    # inicialización del paso de Newton
    s = u_k
    while norm(s)> 1.0e-3
        # inicialización de función no lineal y Jacobiano
        F_k = zeros(nb_nodes_int)
        J_k = zeros(nb_nodes_int,nb_nodes_int)
        for k in 1:nb_elems # pasando sobre todos los elemementos de la malla
            T = M[k,:]*M[k,:]' 
            # Extracción de nodos y conectividad de elementos
                elem_nodes = elems_nodes_conn[k, 1:3] # indice de los nodos del k-ésimo elemento
                K = collect(nodes[elem_nodes, :]')

                # Definiendo la transformación Tₖ
                Jk = K[:,2:3] .- K[:,1]
                X = Jk * ξ .+ K[:,1]
                detJₖ = det(Jk)
                invJₖ = inv(Jk)   

        # Cálculo de función no lineal y Jacobiano en el paso k-ésimo
        for i in 1:3 
                if (M[k,i] == 1)
                    u = u_k[elem_nodes]
                    I = elem_nodes[i] # indice de nodos globales
                    bij = quad_tri_loc_2(f, X, Θ[i], ξ, ω)*abs(detJₖ)
                    U³ = x->Θ[i](x)*(u[1]*Θ[1](x)+u[2]*Θ[2](x)+u[3]*Θ[3](x))^3
                    Ui = quad_tri_loc_2(x->1.0, X, U³ , ξ, ω)*abs(detJₖ)
                    F_k[mode_on[I,2]] += Ui - bij  
                  for j in 1:3
                    if (T[i,j] == 1) 
                        ∇θᵢ = ∇θ[i,:]'*invJₖ
                        ∇θⱼ = ∇θ[j,:]'*invJₖ
                        ∇θᵢ∇θⱼ = dot(∇θᵢ, ∇θⱼ)
                        J = elem_nodes[j] # indice de nodos globales
                        # Ensable Global
                        Rij = 0.5*∇θᵢ∇θⱼ*abs(detJₖ)-quad_tri_loc_2(x->1.0, X, x -> Θ[i](x)*Θ[j](x), ξ, ω)*abs(detJₖ)
                        N = x->3*Θ[i](x)*Θ[j](x)*(u[1]*Θ[1](x)+u[2]*Θ[2](x)+u[3]*Θ[3](x))^2
                        Nij = quad_tri_loc_2(x->1.0, X, N, ξ, ω)*abs(detJₖ)
                        J_k[mode_on[I,2],mode_on[J,2]] += Rij + Nij
                        F_k[mode_on[I,2]] +=  Rij*u[j]
                    end
                end
            end
        end  
        end   
    # Método de Newton
    s =  -J_k\F_k
    u_k[nb_on] += s 
    end
    
    return u_k
end
```




    solve_edp3 (generic function with 1 method)




```julia
function refina_edp3(v::Function, ∇v::Function, MSH::Array, n::Int)
    # Elemento finito
        println("Ejercicio Allen-Chan")
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
        u_sol  = solve_edp3(MSH[i], f, n);
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




    refina_edp3 (generic function with 1 method)




```julia
@time refina_edp3(u_exact, Du_exact, MSH, 1)
```

    Ejercicio Allen-Chan
    Elemento finito: Lagrange Orden 1
    Cuadratura: Gaussiana-01 de 1 punto
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     1.5924317 │    1.0000000 │    14.2729140 │    1.0000000 │
    │   0.3423854 │     0.5762711 │    1.4664118 │     8.2489209 │    0.7910026 │
    │   0.1520212 │     0.2223024 │    1.3742241 │     3.4342264 │    1.2642203 │
    │   0.0792050 │     0.0618313 │    1.8461135 │     1.4967273 │    1.1981737 │
    │   0.0411692 │     0.0159328 │    1.9563391 │     0.7239227 │    1.0479039 │
    │   0.0201700 │     0.0039859 │    1.9990319 │     0.3563796 │    1.0224208 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
      9.788198 seconds (42.09 M allocations: 3.948 GiB, 9.43% gc time, 14.83% compilation time)
    


```julia
@time refina_edp3(u_exact, Du_exact, MSH, 2)
```

    Ejercicio Allen-Chan
    Elemento finito: Lagrange Orden 1
    Cuadratura: Gaussiana-02 de 3 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.7303728 │    1.0000000 │    10.1627344 │    1.0000000 │
    │   0.3423854 │     0.3882873 │    0.9115087 │     7.3649906 │    0.4645330 │
    │   0.1520212 │     0.1379786 │    1.4926803 │     4.5819330 │    0.6847273 │
    │   0.0792050 │     0.0381214 │    1.8557725 │     2.4159292 │    0.9233782 │
    │   0.0411692 │     0.0098254 │    1.9560081 │     1.2279287 │    0.9763514 │
    │   0.0201700 │     0.0024589 │    1.9985076 │     0.6150393 │    0.9974762 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
      8.772389 seconds (62.40 M allocations: 4.488 GiB, 7.62% gc time)
    


```julia
@time refina_edp3(u_exact, Du_exact, MSH, 3)
```

    Ejercicio Allen-Chan
    Elemento finito: Lagrange Orden 1
    Cuadratura: Gaussiana-03 de 3 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     1.1689687 │    1.0000000 │     9.2637414 │    1.0000000 │
    │   0.3423854 │     0.4453615 │    1.3921877 │     8.0548104 │    0.2017443 │
    │   0.1520212 │     0.1720038 │    1.3725364 │     4.5571784 │    0.8217098 │
    │   0.0792050 │     0.0460372 │    1.9015674 │     2.4096507 │    0.9193168 │
    │   0.0411692 │     0.0117294 │    1.9726773 │     1.2270892 │    0.9735840 │
    │   0.0201700 │     0.0029326 │    1.9998484 │     0.6149284 │    0.9967496 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
      8.300437 seconds (62.40 M allocations: 4.488 GiB, 7.60% gc time)
    


```julia
@time refina_edp3(u_exact, Du_exact, MSH, 4)
```

    Ejercicio Allen-Chan
    Elemento finito: Lagrange Orden 1
    Cuadratura: Gaussiana-04 de 4 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.7717938 │    1.0000000 │     9.1703430 │    1.0000000 │
    │   0.3423854 │     0.3871316 │    0.9953915 │     7.4499327 │    0.2997483 │
    │   0.1520212 │     0.1367855 │    1.5009089 │     4.5879374 │    0.6993817 │
    │   0.0792050 │     0.0378107 │    1.8550497 │     2.4170288 │    0.9246111 │
    │   0.0411692 │     0.0097480 │    1.9556188 │     1.2280741 │    0.9768370 │
    │   0.0201700 │     0.0024400 │    1.9982221 │     0.6150577 │    0.9976039 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
      9.330733 seconds (74.12 M allocations: 4.865 GiB, 12.03% gc time)
    


```julia
@time refina_edp3(u_exact, Du_exact, MSH, 5)
```

    Ejercicio Allen-Chan
    Elemento finito: Lagrange Orden 1
    Cuadratura: Gaussiana-05 de 7 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.7589085 │    1.0000000 │     8.5510135 │    1.0000000 │
    │   0.3423854 │     0.4513920 │    0.7495450 │     7.6275988 │    0.1648665 │
    │   0.1520212 │     0.1642889 │    1.4581456 │     4.5365049 │    0.7496478 │
    │   0.0792050 │     0.0443851 │    1.8880883 │     2.4066756 │    0.9145395 │
    │   0.0411692 │     0.0113458 │    1.9679221 │     1.2266980 │    0.9722616 │
    │   0.0201700 │     0.0028357 │    2.0003578 │     0.6148809 │    0.9964011 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
     11.798472 seconds (109.29 M allocations: 5.890 GiB, 12.34% gc time)
    
