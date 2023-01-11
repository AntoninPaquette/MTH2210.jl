function args_check(fct::Function, tspan::AbstractVector{T}, Y0::AbstractVector{T}, nbpas::Integer) where {T<:AbstractFloat}

    length(tspan) == 2 || throw(ArgumentError("Le vecteur tspan doit contenir 2 composantes, [t0 , tf]"))
    nbpas>=0 || throw(DomainError(nbpas, "L'argument nbpas doit être un entier > 0."))

    try
        fct(tspan[1],Y0)
    catch y
        if isa(y,BoundsError)
            error("Le nombre de composantes de Y0 et f ne concorde pas")
        else
            error(y)
        end
    end

    isa(fct(tspan[1],Y0), Union{Real, Vector{<:Real}}) || throw(DimensionMismatch("La fonction f ne retourne pas un scalaire ou un vecteur"))

    length(Y0) == length(fct(tspan[1],Y0)) || throw(DimensionMismatch("Le nombre de composantes de Y0 et f ne concorde pas"))
    
end

function edo_init(tspan::AbstractVector{T}, Y0::AbstractVector{T} , nbpas::Integer) where {T<:AbstractFloat}
    N       =   length(Y0)
    Y       =   zeros(T,N,nbpas+1)
    Y[:,1]  .=  Y0
    temps   =   LinRange{T}(tspan[1], tspan[2] , nbpas+1)
    h       =   temps[2] - temps[1]

    return (N, Y, temps, h)
end