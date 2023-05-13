# Problema 1

Obtener el interpolador $v_I(x)$ de la funcion 
$$v(x,y)= \cos(4\pi x)\cdot cos^2(4\pi y)$$
para el dominio $\Omega ≔ (0, 1) \times (0, 1)$ utilizando

* Los elementos finitos de **Lagrange de primer orden**
* Los elementos finitos de **Crouzeix–Raviart de primer orden**
* Los elementos finitos de **Lagrange de segundo orden**

Obtener entonces la tasa de convergencia para el error $e:= v_I - v$ en la norma $L^2(\Omega)$ y la seminorma $H^1(\Omega)$ para la familia de mallas proporcionada.

## Definición de la función y su gradiente


```julia
# Definición de la función a interpolar
v = (x) -> cos(4*π*x[1])*(cos(4*π*x[2]))^2

# Definición del gradiente
∇v =  (x) -> [-4*π*sin(4*π*x[1])*(cos(4*π*x[2]))^2
              -8*π*cos(4*π*x[1])*cos(4*π*x[2])*sin(4*π*x[2])]
```




    #5 (generic function with 1 method)



# Interpolador
## El Interpolador local

Para $K \in \mathcal{T}_h$, la tripleta $\left\lbrace K, P_K, \Sigma_K \right\rbrace$ definida por

$$
\begin{cases}
K=T_K(\hat{K}) \\
P_K=\left\lbrace \psi_K^{-1}(\hat{p}) ; \ \widehat{p} \in \hat{P}\right\rbrace ; \\
\Sigma_K=\left\lbrace \left\lbrace \sigma_{K, i}\right\rbrace_{1 \leq i \leq n_{s h}} ;\  \sigma_{K, i}(p)=\hat{\sigma}_i\left(\psi_K(p)\right), \ \forall p \in P_K\right\rbrace
\end{cases}
$$


es un elemento finito. 

Las funciones de forma locales son 

$$ \theta_{K,i} = \psi_K^{-1} \left( \hat{\theta}_i \right), 1 \leq i \leq n_{sh} $$

$$\theta_{K,i} = \psi_K ^{-1} \left( \hat{\theta}i \right), 1 \leq i \leq n{sh}  $$

y el operador de interpolación local $$\mathcal{I}_K: V(K) \longmapsto P_K$$ está dado por

$$
\mathcal{I}_K v =\sum_{i=1}^{n_{\mathrm{sh}}} \sigma_{K, i}(v) \theta_{K, i}
$$

- Para el elemento finito de Lagrange:
    * $K$ es un simplex en $\mathbb{R}^d$
    * $P_K=\mathbb{P}_k$, $n_{s h}=dim\left( \mathbb{P}_k\right), k\geq 1 $
    * $\sigma_{K, i}(p)=p(x_i)$, para $\left\lbrace x_i \right\rbrace_{i=1}^{n_{sh}}$ los nodos.
        
        Usando las coordenadas baricéntricas: $(\lambda_0, \lambda_1 ,\ldots , \lambda_d)$
        - Si $k=1$ las funciones de forma son: $$\theta_i = \lambda_i, \quad 0\leq i \leq d$$
        - Si $k=2$ las funciones de forma son: $$\theta = \begin{cases}  \lambda_i \left(2\lambda_i-1\right) & \quad 0\leq i \leq d \\   4\lambda_i\lambda_j & \quad 0\leq i < j \leq   d \\  \end{cases}$$
- Para el elemento finito de Crouzeix-Raviart:
    * $K$ es un simplex en $\mathbb{R}^d$
    * $P_K=\mathbb{P}_1$
    * $\sigma_{K, i}(v),\quad  0\leq i \leq d$ dado por: $$\sigma_{K, i}(v)=\frac{1}{\vert F_i\vert }\int_{F_i}v $$
    
    Usando las coordenadas baricéntricas: $(\lambda_0, \lambda_1 ,\ldots , \lambda_d)$
    $$\theta_i = d\left(\frac{1}{d}- \lambda_i \right), \quad 0\leq i \leq d$$


Como el siguiente diagrama conmuta:
$$
\begin{array}{cc}
V(K) \stackrel{\psi_K}{\longrightarrow} V(\hat{K}) \quad \\
\downarrow{\mathcal{I}_K} \quad  \quad \quad \downarrow \mathcal{I}_{\hat{K}} \\
P_K\quad \stackrel{\psi_K}{\longrightarrow}\quad \hat{P} \quad
\end{array}
$$

con 
$$\mathcal{I}_{\hat{K}}: V(\hat{K}) \longmapsto \hat{P}$$ 
$$
\mathcal{I}_{\hat{K}} \hat{v} =\sum_{i=1}^{n_{\mathrm{sh}}} \hat{\sigma}_{i}(\hat{v}) \hat{\theta}_i
$$ el operador de interpolación de referencia.


Entonces podemos expresar el interpolador local en términos del elemento finito de referencia:

$$ \mathcal{I}_K v(x) = \sum_{i=1}^{n_{sh}} \sigma_{K,i}(v)\ \theta_{K,i}(x)=  \sum_{i=1}^{n_{sh}}  \sigma_{K,i}(v)\ \psi_K^{-1}\left(\hat{\theta}_i\right)\left(x\right) =\sum_{i=1}^{n_{sh}} \sigma_{K,i}(v)\ \left(\hat{\theta}_i\circ T_K^{-1}\right)(x)$$


Para nuestro problema $d=2$, tenemos

- Elementos finitos de Lagrange de primer orden
$$\theta = \begin{pmatrix} \theta_0 \\ \theta_1 \\ \theta_2 \end{pmatrix}= \begin{pmatrix} \lambda_0 \\ \lambda_1 \\ \lambda_2 \end{pmatrix} $$
- Elementos finitos de Crouzeix–Raviart de primer orden
$$\theta = \begin{pmatrix} \theta_0 \\ \theta_1 \\ \theta_2 \end{pmatrix}= \begin{pmatrix} 1-2\lambda_0 \\ 1-2\lambda_1 \\ 1-2\lambda_2 \end{pmatrix} $$
- Elementos finitos de Lagrange de segundo orden
$$\theta = \begin{pmatrix} \theta_0 \\ \theta_1 \\ \theta_2 \\ \theta_3 \\ \theta_4 \\ \theta_5 \end{pmatrix}= \begin{pmatrix} \lambda_0(2\lambda_0-1) \\ \lambda_1(2\lambda_1-1)\\ \lambda_2(2\lambda_2-1) \\ 4\lambda_0\lambda_1 \\ 4\lambda_1\lambda_2 \\ 4\lambda_0\lambda_2 \end{pmatrix} $$

* Function: `get_θ_functions(ξ::Vector, ef::Int)`


```julia
"""
## `get_θ_functions(ξ::Vector, ef::Int)`

Calcula las funciones de forma locales en el punto `ξ` para el elemento finito `ef`, 

## Argumentos
- `ξ` : Coordenadas locales de un punto en el simplejo de referencia. 
- `ef` : Índice del Elemento finito
    - `ef` = 1 : Elementos finitos de Lagrange de primer orden
    - `ef` = 2 : Elementos finitos de Crouzeix–Raviart de primer orden
    - `ef` = 3 : Elementos finitos de Lagrange de segundo orden

## Salida
- `θ` : Funciones de forma locales en el punto `ξ`

"""
function get_θ_functions(ξ::Vector, ef::Int) 
    # Inicializar variable para almacenar las funciones de forma
    λ = [1 - ξ[1] - ξ[2], ξ[1], ξ[2]]
    
    # Calcular las funciones de forma según el tipo de elemento finito
    if ef == 1 # Elementos de Lagrange de primer orden
        θ = λ
    elseif ef == 2 # Elementos de Crouzeix-Raviart de primer orden
        θ = -2*λ .+ 1            
    elseif ef == 3 # Elementos de Lagrange de segundo orden
        θ = zeros(6)
        θ[1] = λ[1]*(2λ[1]-1)
        θ[2] = λ[2]*(2λ[2]-1)
        θ[3] = λ[3]*(2λ[3]-1)
        θ[4] = 4λ[1]*λ[2]
        θ[5] = 4λ[2]*λ[3]
        θ[6] = 4λ[1]*λ[3]
    else
        throw(ArgumentError("Tipo de elemento finito no válido"))
    end
    
    return θ
end
```




    get_θ_functions



## El Interpolador global

Usando la familia de elementos finitos $\left\lbrace K, P_K, \Sigma_K \right\rbrace_{K \in \mathcal{T}_h}$, se puede construir un operador de interpolación global $\mathcal{I}_h$ de la siguiente manera: primero, se elige su dominio como
$$
D\left(\mathcal{I}_h\right)=\left\{v \in\left[L^1\left(\Omega_h\right)\right]^m ; \forall K \in \mathcal{T}_h, v_{\mid K} \in V(K)\right\}
$$

donde $\Omega_h$ es la interpolación geométrica de $\Omega$. Para una función $v \in D\left(\mathcal{I}h\right)$, las cantidades $\sigma_{K, i}\left(v_{\mid K}\right)$ tienen sentido en todos los elementos de la malla y para todo $1 \leq i \leq n_{\text {sh }}$. Luego, el interpolante global $\mathcal{I}_h v$ se puede especificar elemento por elemento utilizando los operadores de interpolación local, es decir,

$$
\forall K \in \mathcal{T}_h, \quad\left(\mathcal{I}_h v\right)_{\mid K}=\mathcal{I}_K\left(v_{\mid K}\right)=\sum_{i=1}^{n_{\text {sh }}} \sigma_{K, i}\left(v_{\mid K}\right) \theta_{K, i}
$$

El operador de interpolación global se define de la siguiente manera:
$$
\mathcal{I}_h: D\left(\mathcal{I}_h\right)  \longmapsto W_h,
$$

$$
\mathcal{I}_h (v) := \sum_{K \in \mathcal{T}_h} \sum_{i=1}^{n_{\mathrm{ah}}} \sigma_{K, i}\left(v_{\mid K}\right) \theta_{K, i} 
$$
donde $W_h$, el codominio de $\mathcal{I}_h$, es
$$
W_h=\left\{v_h \in\left[L^1\left(\Omega_h\right)\right]^m ; \forall K \in \mathcal{T}_h, v_{\mid K} \in P_K\right\}
$$

El espacio $W_h$ se llama espacio de aproximación.

## Error en $L^2$ del Interpolador global

Para el error del interpolador global en $L^2$ tenemos:

$$\Vert \mathcal{I}_h v-v\Vert_{L^2}^2 = \sum_{K \in \mathcal{T}_h} \int_{K} \left( \mathcal{I}_h v\vert_{K}-v\vert_{K}\right)^2(x) dK = \sum_{K \in \mathcal{T}_h} \int_{\hat{K}} \left|\operatorname{det}\left(J_K\right)\right| \left( \mathcal{I}_K v-v\vert_{K}\right)^2\circ T_K(\hat{x}) d\hat{K} $$ 

$$\therefore \Vert \mathcal{I}_h v-v\Vert_{L^2}^2 \approx \sum_{K \in \mathcal{T}_h}\left|\operatorname{det}\left(J_K\right)\right| \sum_{l=1}^{l_q} \omega_l\left[\left( \mathcal{I}_K v-v\right)^2 \circ T_K\right]\left(\hat{\xi}_l\right)$$

* Function: `error_local_L2(v::Function, σ_T::Array, X::Array, ξ::Array, ω::Array, ef::Int)`


```julia
"""
## `error_local_L2(v::Function, σ_T::Vector, X::Array, ξ::Array, ω::Array, ef::Int)`

Calcula el error local en L² del interpolador de la función v sobre el triángulo T,
en el elemento finito `ef`:

## Argumentos
- `v` : Solución exacta
- `σ_T` : Vector con grados de libertad.
- `X` : Matriz con Coordenadas globales de los gauss points.
- `ξ` : Matriz con Coordenadas locales de los gauss points. 
- `ω` : Vector con pesos de los gauss points. 
- `ef`: Índice del Elemento finito
    - `ef` = 1 : Elementos finitos de Lagrange de primer orden
    - `ef` = 2 : Elementos finitos de Crouzeix–Raviart de primer orden
    - `ef` = 3 : Elementos finitos de Lagrange de segundo orden

## Salida
- Error local en L² del interpolador de la función v
"""
function error_local_L2(v::Function, σ_T::Array, X::Array, ξ::Array, ω::Array, ef::Int)
    # Número de puntos de la cuadratura
    l_q = length(ω)
    
    # Inicialización del vector de residuos al cuadrado
    squared_residuals = zeros(l_q)
    
    # Bucle sobre los puntos de cuadratura
    @inbounds for i in 1:l_q
        # Punto de cuadratura y funciones de forma del elemento finito
        x = X[:,i]
        θ = get_θ_functions(ξ[:,i], ef)
        
        # Interpolación y cálculo del residuo al cuadrado
        Iv_T = dot(σ_T, θ)
        res = v(x) - Iv_T
        squared_residuals[i] = res^2
    end
    
    # Cálculo del error local en L2
    err_loc = 0.5 * dot(1.0*ω, squared_residuals)
    return err_loc    
end
```




    error_local_L2



## Error en $H^1$ del Interpolador local

Para el error del interpolador local en $H^1$ tenemos:

$$\Vert  \nabla \mathcal{I}_h v-\nabla v\Vert_{L^2}^2 = \sum_{K \in \mathcal{T}_h} \int_{K} \left( \nabla \mathcal{I}_h v\vert_{K}-\nabla v\vert_{K}\right)^2(x) dK = \sum_{K \in \mathcal{T}_h} \int_{\hat{K}} \left|\operatorname{det}\left(J_K\right)\right| \left( \nabla  \mathcal{I}_K v-\nabla v\vert_{K}\right)^2\circ T_K(\hat{x}) d\hat{K} $$ 

$$\therefore \Vert \mathcal{I}_h v-v\Vert_{L^2}^2 \approx \sum_{K \in \mathcal{T}_h}\left|\operatorname{det}\left(J_K\right)\right| \sum_{l=1}^{l_q} \omega_l\left[\left( \nabla \mathcal{I}_K v-\nabla v\right)^2 \circ T_K\right]\left(\hat{\xi}_l\right)$$



Ya que $$ \mathcal{I}_K v(x) = \sum_{i=1}^{n_{sh}} \sigma_{K,i}(v)\ \theta_{K,i}(x)=  \sum_{i=1}^{n_{sh}}  \sigma_{K,i}(v)\ \psi_K^{-1}\left(\hat{\theta}_i\right)\left(x\right) =\sum_{i=1}^{n_{sh}} \sigma_{K,i}(v)\ \left(\hat{\theta}_i\circ T_K^{-1}\right)(x)$$

Entonces podemos obtener el gradiente como sigue:
$$\nabla_x \mathcal{I}_K v(x)  =\sum_{i=1}^{n_{sh}} \sigma_{K,i}(v)\nabla_x  \left(\hat{\theta}_i\circ T_k^{-1}\right)(x) =\sum_{i=1}^{n_{sh}} \sigma_{K,i}(v) \nabla_\hat{x} \hat{\theta}_i  (\hat{x}) \frac{\partial T_k^{-1}(x)}{\partial x}  =\sum_{i=1}^{n_{sh}} \sigma_{K,i}(v)  \nabla_\hat{x} \hat{\theta}_i (\hat{x}) \frac{\partial \hat{x} }{\partial x}$$
$$ \therefore \nabla_x \mathcal{I}_K v(x) =  \sum_{i=1}^{n_{sh}} \sigma_{K,i}(v) \nabla_\hat{x} \hat{\theta}_i  (\hat{x})J_K^{-1}  $$

Que en notación vectorial podemos escribir como
$$\nabla_x \mathcal{I}_K v(x) = \sigma_K^t(v) \cdot D_\hat{x} \theta (\hat{x})\cdot J_K^{-1}$$

Para nuestro problema $d=2$, tenemos

- Elementos finitos de Lagrange de primer orden
$$\theta = \begin{pmatrix} \theta_0 \\ \theta_1 \\ \theta_2 \end{pmatrix}= \begin{pmatrix} \lambda_0 \\ \lambda_1 \\ \lambda_2 \end{pmatrix} = \begin{pmatrix} 1-\hat{x}_1-\hat{x}_2 \\ \hat{x}_1\\ \hat{x}_2 \end{pmatrix} \implies D_\hat{x}\theta = \begin{pmatrix}-1 & -1 \\ \ 1&\ 0 \\ \ 0&\ 1\end{pmatrix} $$
- Elementos finitos de Crouzeix–Raviart de primer orden
$$\theta = \begin{pmatrix} \theta_0 \\ \theta_1 \\ \theta_2 \end{pmatrix}= \begin{pmatrix} 1-2\lambda_0 \\ 1-2\lambda_1 \\ 1-2\lambda_2 \end{pmatrix}\implies D_\hat{x}\theta = \begin{pmatrix}\ \ 2 & \ 2 \\ -2 &\ 0 \\ \ \ 0&-2 \end{pmatrix}   $$
- Elementos finitos de Lagrange de segundo orden
$$\theta = \begin{pmatrix} \theta_0 \\ \theta_1 \\ \theta_2 \\ \theta_3 \\ \theta_4 \\ \theta_5 \end{pmatrix}= \begin{pmatrix} \lambda_0(2\lambda_0-1) \\ \lambda_1(2\lambda_1-1)\\ \lambda_2(2\lambda_2-1) \\ 4\lambda_0\lambda_1 \\ 4\lambda_1\lambda_2 \\ 4\lambda_0\lambda_2 \end{pmatrix}\implies D_\hat{x}\theta =D_{\lambda}  \theta\cdot D_{\hat{x}} \lambda = \begin{pmatrix}
1-4\lambda_0&1-4\lambda_0\\
4\lambda_1-1&0\\
0 &4\lambda_2-1\\
4(\lambda_0-\lambda_1)&-4\lambda_1\\
4\lambda_2 &4\lambda_1\\
-4\lambda_2 &4(\lambda_0-\lambda_2)   \end{pmatrix}  $$

* Function: `get_Dθ_functions(ξ::Vector, ef::Int)`


```julia
"""
## `get_Dθ_functions(ξ::Vector, ef::Int)`

Calcula la derivada de funciones de forma locales en el punto `ξ`
para el elemento finito `ef`.

## Argumentos
- `ξ` : Coordenadas locales de un punto en el simplejo de referencia. 
- `ef` : Índice del Elemento finito
    - `ef` = 1 : Elementos finitos de Lagrange de primer orden
    - `ef` = 2 : Elementos finitos de Crouzeix–Raviart de primer orden
    - `ef` = 3 : Elementos finitos de Lagrange de segundo orden

## Salida
- `Dθ` : Derivada de funciones de forma locales en el punto `ξ`
"""
function get_Dθ_functions(ξ::Vector, ef::Int)
    # Calcular las funciones de forma según el tipo de elemento finito
        if ef == 1
            Dθ = [-1 -1; 1 0; 0 1]
        elseif ef == 2
            Dθ = -2*[-1 -1; 1 0; 0 1]
        elseif ef == 3
            λ = [1-ξ[1]-ξ[2], ξ[1], ξ[2]]
            Dθ = [1-4λ[1] 4λ[2]-1 0 4λ[1]-4λ[2] 4λ[3] -4λ[3];
                  1-4λ[1] 0 4λ[3]-1 -4λ[2] 4λ[2] 4λ[1]-4λ[3]]'    
    else
            throw(ArgumentError("Tipo de elemento finito no válido"))
        end
    return Dθ
end
```




    get_Dθ_functions



* Function: `error_local_H1(∇v::Function, σ_T::Vector, invJₖ::Matrix, X::Array, ξ::Array, ω::Array, ef::Int)`


```julia
"""
## `error_local_H1(∇v::Function, σ_T::Vector, invJₖ::Matrix, X::Array, ξ::Array, ω::Array, ef::Int)`

Calcula el error local en H¹ del interpolador de la función v sobre el triángulo T,
en el elemento finito `ef`:

## Argumentos
- `∇v` : Gradiente de la Solución exacta
- `σ_T` : Vector con grados de libertad.
- `invJₖ` : Matriz inversa de Jₖ
- `X` : Matriz con Coordenadas globales de los gauss points.
- `ξ` : Matriz con Coordenadas locales de los gauss points. 
- `ω` : Vector con pesos de los gauss points. 
- `ef`: Índice del Elemento finito
    - `ef` = 1 : Elementos finitos de Lagrange de primer orden
    - `ef` = 2 : Elementos finitos de Crouzeix–Raviart de primer orden
    - `ef` = 3 : Elementos finitos de Lagrange de segundo orden

## Salida
- Error local en H¹ del interpolador de la función v
"""
function error_local_H1(∇v::Function, σ_T::Vector, invJₖ::Matrix, X::Array, ξ::Array, ω::Array, ef::Int)
    # Número de puntos de la cuadratura
    l_q = length(ω)
    
    # Inicialización del vector de residuos al cuadrado
    squared_residuals = zeros(l_q)
    
    # Bucle sobre los puntos de cuadratura
    @inbounds for i in 1:l_q
        x = X[:, i]
        
        # Derivada del Interpolador  
        Dθ = get_Dθ_functions(ξ[:,i], ef)
        ∇Iv_T = σ_T'*Dθ*invJₖ
        
        #  Cálculo del residuo al cuadrado 
        res = ∇v(x) - vec(∇Iv_T)
        squared_residuals[i] = dot(res, res)
    end
    
    # Cálculo del error local en H1
    err_loc = 0.5 * dot(1.0*ω, squared_residuals)
    return err_loc    
end
```




    error_local_H1



# Grados de Libertad

* Function: `cuadratura_gauss(f::Function, a::Number, b::Number)`


```julia
"""
### `cuadratura_gauss(f::Function, a::Number, b::Number)`

Aproxima el valor de la integral de la función `f` en el intervalo `[a,b]` utilizando
la fórmula de cuadratura de Gauss-Legendre con tres nodos.

### Argumentos

- `f`: Función a integrar.
- `a`: Límite inferior del intervalo de integración.
- `b`: Límite superior del intervalo de integración.

### Salida

- `val`: Aproximación del valor de la integral de `f` en el intervalo `[a,b]`.

"""
function cuadratura_gauss(f::Function, a, b)
        # Nodos de la cuadratura
        x₁ = -0.5*(b-a)*sqrt(3/5) + 0.5*(a+b)   # Nodo x₁
        x₂ = 0.5*(a + b)                        # Nodo x₂
        x₃ = 0.5*(b-a)*sqrt(3/5) + 0.5*(a+b)    # Nodo x₃
    
        # Cálculo de integral utilizando la fórmula de cuadratura de Gauss-Legendre
        val = 5*f(x₁) + 8*f(x₂) + 5*f(x₃)       # Suma ponderada de f en los nodos
        val = 0.5*(b-a)*val/9                   # Multiplicación por el factor de escala
    
    return val  # Retorna el valor aproximado de la integral de f en [a,b]
end
```




    cuadratura_gauss



* Function: `get_σ_interpolator(v::Function, K::Matrix, ef::Int)`


```julia
"""
## `get_σ_interpolator(v::Function, K::Matrix, ef::Int)`

Devuelve el vector de grados de libertar en el simplejo `K` para el elemento finito `ef`.

## Argumentos
- `v` : Función a Interpolar
- `K` : Coordenadas los vértices del simplejo K. 
- `ef` : Índice del Elemento finito
    - `ef` = 1 : Elementos finitos de Lagrange de primer orden
    - `ef` = 2 : Elementos finitos de Crouzeix–Raviart de primer orden
    - `ef` = 3 : Elementos finitos de Lagrange de segundo orden

## Salida
- `σ_T` : Vector de grados de libertar en el simplejo `K`.
"""
function get_σ_interpolator(v::Function, K::Matrix, ef::Int)
        z₀ = K[:,1]
        z₁ = K[:,2]
        z₂ = K[:,3]
        
        if ef == 1 
            σ_T = zeros(3)
            σ_T = [v(z₀), v(z₁), v(z₂)]
        elseif ef == 2   
            σ_T = zeros(3)
            f₁(t) = v((1-t)*z₀+t*z₁)
            f₂(t) = v((1-t)*z₁+t*z₂)
            f₃(t) = v((1-t)*z₂+t*z₀)
            σ_T[1] = cuadratura_gauss(f₂, 0, 1)
            σ_T[2] = cuadratura_gauss(f₃, 0, 1)
            σ_T[3] = cuadratura_gauss(f₁, 0, 1)
        elseif ef == 3 
            σ_T = zeros(6)
            σ_T[1] = v(z₀)
            σ_T[2] = v(z₁)
            σ_T[3] = v(z₂)
            σ_T[4] = v(0.5*z₀+0.5*z₁)            
            σ_T[5] = v(0.5*z₂+0.5*z₁)                   
            σ_T[6] = v(0.5*z₀+0.5*z₂)
        else
            throw(ArgumentError("Tipo de elemento finito no válido"))
        end
    return σ_T
end
```




    get_σ_interpolator



# Tasa de convergencia del error en $L^2$ y $H^1$

* Function: `error_global(v::Function, ∇v::Function, mesh::Dict, n::Int, ef::Int)`


```julia
"""
## `error_global(v::Function, ∇v::Function, mesh::Dict, n::Int, ef::Int)`

Calcula el error global de interpolación en la norma L² y la seminorma H¹
en una malla de elementos finitos para una función de interpolación `v` 
y su gradiente `∇v`. Utiliza un esquema de cuadratura de Gauss para evaluar
la función a través de su interpolación y compara los valores obtenidos con
la solución exacta. 

## Argumentos
- `v` : Solución Exacta.
- `∇v` : Gradiente de la Solución Exacta.
- `mesh` : Diccionario que contiene la información de la malla de elementos finitos. 
- `n` : Índice de la Cuadratura para la integración
    - `n` = 1 : Gaussiana-01 de 1 punto
    - `n` = 2 : Gaussiana-02 de 3 puntos
    - `n` = 3 : Gaussiana-03 de 3 puntos
    - `n` = 4 : Gaussiana-04 de 4 puntos
    - `n` = 5 : Gaussiana-05 de 7 puntos
- `ef` : Índice del Elemento finito
    - `ef` = 1 : Elementos finitos de Lagrange de primer orden
    - `ef` = 2 : Elementos finitos de Crouzeix–Raviart de primer orden
    - `ef` = 3 : Elementos finitos de Lagrange de segundo orden

## Salida
- `err_L2`: error global de interpolación en la norma L².
- `err_H1`: error global de interpolación en la seminorma H¹.

"""
function error_global(v::Function, ∇v::Function, mesh::Dict, n::Int, ef::Int)
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
        σ_T = get_σ_interpolator(v, K, ef)
        
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




    error_global



# Cálculo de la tasa de convergencia en L2 y H1

* Function: `refina(v::Function, ∇v::Function, MSH::Array, ef::Int, n::Int)`


```julia
"""
## `refina(v::Function, ∇v::Function, MSH::Array, ef::Int, n::Int)`

Función que calcula los errores de aproximación en la norma L² y la seminorma H¹,
realiza un análisis de convergencia de error para el interpolador utilizando
diferentes refinamientos en una malla y diferentes elementos finitos. 

## Argumentos
- `v` : Solución Exacta.
- `∇v` : Gradiente de la Solución Exacta.
- `MESH` : Arreglo de diccionarios que contiene la información de las mallas del refinamiento.
- `ef` : Índice del Elemento finito
    - `ef` = 1 : Elementos finitos de Lagrange de primer orden
    - `ef` = 2 : Elementos finitos de Crouzeix–Raviart de primer orden
    - `ef` = 3 : Elementos finitos de Lagrange de segundo orden
- `n` : Índice de la Cuadratura para la integración
    - `n` = 1 : Gaussiana-01 de 1 punto
    - `n` = 2 : Gaussiana-02 de 3 puntos
    - `n` = 3 : Gaussiana-03 de 3 puntos
    - `n` = 4 : Gaussiana-04 de 4 puntos
    - `n` = 5 : Gaussiana-05 de 7 puntos
## Salida
- La función retorna la tabla de resultados.
"""
function refina(v::Function, ∇v::Function, MSH::Array, ef::Int, n::Int)
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
        # Calculo del error en la norma L² y en la seminorma H¹
        L2_error_vec[i], H1_error_vec[i] = error_global(v, ∇v, MSH[i], n, ef)
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




    refina



# Presentación de resultados


```julia
refina(v, ∇v, MSH, 1, 5)
```

    Elemento finito: Lagrange Orden 1
    Cuadratura: Gaussiana-05 de 7 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.6596904 │    1.0000000 │     8.2485125 │    1.0000000 │
    │   0.3423854 │     0.4496413 │    0.5530146 │     6.9950336 │    0.2378030 │
    │   0.1520212 │     0.1473130 │    1.6098901 │     5.0794873 │    0.4616481 │
    │   0.0792050 │     0.0558648 │    1.3988724 │     3.0743890 │    0.7243832 │
    │   0.0411692 │     0.0145517 │    1.9407507 │     1.5905570 │    0.9507677 │
    │   0.0201700 │     0.0036800 │    1.9833937 │     0.7998152 │    0.9917935 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
    


```julia
refina(v, ∇v, MSH, 2, 5)
```

    Elemento finito: Crouzeix–Raviart Orden 1
    Cuadratura: Gaussiana-05 de 7 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.4314069 │    1.0000000 │     6.7387760 │    1.0000000 │
    │   0.3423854 │     0.2985520 │    0.5310671 │     5.8420809 │    0.2060042 │
    │   0.1520212 │     0.0919788 │    1.6986096 │     4.3226416 │    0.4345691 │
    │   0.0792050 │     0.0326415 │    1.4945924 │     2.5023374 │    0.7886369 │
    │   0.0411692 │     0.0084430 │    1.9508759 │     1.2922393 │    0.9534030 │
    │   0.0201700 │     0.0021287 │    1.9878123 │     0.6517931 │    0.9873874 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
    


```julia
refina(v, ∇v, MSH, 3, 5)
```

    Elemento finito: Lagrange Orden 2
    Cuadratura: Gaussiana-05 de 7 puntos
    ┌─────────────┬───────────────┬──────────────┬───────────────┬──────────────┐
    │ Parámetro h │ Error L²_norm │ Rate L²_norm │ Error H¹_norm │ Rate H¹_norm │
    ├─────────────┼───────────────┼──────────────┼───────────────┼──────────────┤
    │   0.5303301 │     0.2841278 │    1.0000000 │     8.0323610 │    1.0000000 │
    │   0.3423854 │     0.0928317 │    1.6138510 │     6.9346450 │    0.2120021 │
    │   0.1520212 │     0.0217665 │    2.0925059 │     3.1538145 │    1.1367242 │
    │   0.0792050 │     0.0017820 │    3.6105501 │     0.8026173 │    1.9743137 │
    │   0.0411692 │     0.0001447 │    3.6226096 │     0.2081855 │    1.9468424 │
    │   0.0201700 │     0.0000118 │    3.6175815 │     0.0523020 │    1.9929310 │
    └─────────────┴───────────────┴──────────────┴───────────────┴──────────────┘
    

