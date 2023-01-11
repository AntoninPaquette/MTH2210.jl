"""
Résolution d'EDOs de type `` \\frac{dY}{dt}(t) = F(t,Y(t))`` avec les
conditions initiales ``Y(t_0) = Y_0`` à l'aide de la méthode d'Euler modifié
jusqu'au temps ``t_f`` avec pas constant ``h``:

``\\tilde{Y} = Y_{n} + h F(t_n,Y_n)
\\\\ Y_{n+1} = Y_n + \\frac{h}{2} \\left(F(t_n,Y_n) + F(t_n +h, \\tilde{Y}) \\right)``

# Syntaxe
```julia
(t,Y) = eulermod(fct , tspan , Y0 , nbpas)F(t_n,Y_n) +
```

# Entrée
    1.  fct     -   Fonction décrivant le système de N EDOs
    2.  tspan   -   (Array{Float,1}) Vecteur contenant le temps initial et final (tspan=[t0,tf])
    3.  Y0      -   (Array{Float,1}) Vecteur contenant les N conditions initiales
    4.  nbpas   -   (Integer) Nombre de pas de temps

# Sortie
    1.  temps   -   (Array{Float,1}) Vecteur contenant les pas de temps
    2.  Y       -   (Array{Float,2}) Matrice de dimension N x (nbpas+1) contenant les approximations

# Exemples d'appel
```julia
function my_edo(t,z)
    f = zeros(length(z))
    f[1] = z[2]
    f[2] = -z[1]
    return f
end
(t,y)   =   eulermod(my_edo , [0;10] , [1;0] , 1000)
```
```julia
(t,y)   =   eulermod((t,y) -> cos(t) , [0;2] , [1] , 1000)
```
```julia
(t,y)   =   eulermod((t,y) -> [y[2];-y[1]] , [0;10] , [1;0] , 1000)
```
"""
function eulermod(fct::Function, tspan::AbstractVector{T}, Y0::AbstractVector{T} , nbpas::Integer) where {T<:AbstractFloat}

    # Vérification des arguments d'entrée
    args_check(fct, tspan, Y0, nbpas)

    (N, Y, temps, h) = edo_init(tspan, Y0, nbpas)
    
    y_tilde =   zeros(T,N)

    for t=1:nbpas
        y_tilde .= view(Y,:,t) .+ h .* fct(temps[t], view(Y,:,t))
        Y[:,t+1] .= view(Y,:,t) .+ (h ./ 2) .* (fct(temps[t], view(Y,:,t)) .+ fct(temps[t] + h, y_tilde))
    end

    return  temps , Y

end

@inline eulermod(fct::Function , tspan::AbstractVector{<:Real} , Y0::AbstractVector{<:Real} ,nbpas::Integer) = eulermod(fct, Float64.(tspan), Float64.(Y0), nbpas)

@inline eulermod(fct::Function, tspan::AbstractVector{<:Real}, Y0::Real, nbpas::Integer) = eulermod(fct, Float64.(tspan), [Float64(Y0)], nbpas)
