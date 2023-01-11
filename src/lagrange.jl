"""
Interpolation de Lagrange

# Syntaxe
```julia
Lx = lagrange(xi , yi , x)
```

# Entrée
    1.  xi          -   (Vector{Float}) Abscisses des points d'interpolation
    2.  yi          -   (Vector{Float}) Ordonnées des points d'interpolation
    3.  x           -   (Vector{Float}) Points où le polynôme d'interpolation est évalué

# Sortie
    1.  Lx          -   (Vector{Float}) Valeur du polynôme aux points x

# Exemples d'appel
```julia
Lx = lagrange([-1,0,1] , [1,0,1] , LinRange(-1, 1, 200))
```
"""
function lagrange(xi::AbstractVector{T}, yi::AbstractVector{T}, x::AbstractVector{T}) where {T<:AbstractFloat}

    N = length(xi)

    if length(yi) != N
            error("Les vecteurs xi et yi doivent avoir la meme taille")
    elseif xi == x
            println("Le polynôme d''interpolation est évalué exactement au points d''interpolation")
    end

    # Calcul des poids barycentriques
    w   =   Vector{T}(undef, N)
    ind =   1:1:N

    for t=1:N
        w[t] = 1 ./ prod(xi[ind .!= t] .- xi[t])
    end

    # Calcul de l'interpolant aux points x
    test        =   Array{Bool,1}(undef,length(xi))
    Lx          =   Vector{T}(undef, length(x))
    vec_diff    =   Vector{T}(undef, N)
    ratio       =   Vector{T}(undef, N)

    for t=1:length(x)
        test .= (x[t] .== xi)
        if any(test)
            Lx[t] = yi[test][1]
        else
            vec_diff .= x[t] .- xi
            ratio    .= w ./ vec_diff
            Lx[t]    = sum(ratio .* yi) / sum(ratio)
        end
    end

    return Lx
end

@inline lagrange(xi::AbstractVector{<:Real}, yi::AbstractVector{<:Real}, x::AbstractVector{<:Real}) = lagrange(Float64.(xi), Float64.(yi), Float64.(x))
