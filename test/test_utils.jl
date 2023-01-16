struct algo
    name::Function
    order::Real
end

struct ode_prob
    fct::Function
    y0::Union{Real, Vector{<:Real}}
    tspan::Vector{<:Real}
    ex_sol::Function
    degre::Real
end

function order_computation(err::Vector{<:Real}, ratio_h::Real, tol::Real=0.2)

    N = length(err)
    N>=4 || throw(ArgumentError("Not enough data (minimum of 3)"))

    order_app = log.(err[1:end-1]./err[2:end]) ./ log(ratio_h)

    stable_region = (order_app.>0) .& (abs.(grad_vec(order_app)).<tol) #.& abs.(del2_vec(order_app)).<2*tol
    ind_stable_region = findall(stable_region)

    # Sanity check
	if isempty(ind_stable_region)
		println("No asymptotic zone")
	elseif length(ind_stable_region) < 2
		println("Asymptotic zone not large")
    elseif any(diff(ind_stable_region) .!= 1)
        print(ind_stable_region)
        println("Broken asymptotic zone")
	end
	
	order = mean(order_app[ind_stable_region])

    return (order, order_app)

end


function order_computation_nl(err::Vector{<:Real}, tol::Real=0.2)

    N = length(err)
    N>=4 || throw(ArgumentError("Not enough data (minimum of 3)"))

    order_app = log.(err[2:end-1] ./ err[3:end]) ./ log.(err[1:end-2] ./ err[2:end-1])

    stable_region = (order_app.>0) .& (abs.(grad_vec(order_app)).<tol) #.& abs.(del2_vec(order_app)).<2*tol
    ind_stable_region = findall(stable_region)

    # Sanity check
	if isempty(ind_stable_region)
		println("No asymptotic zone")
	elseif length(ind_stable_region) < 2
		println("Asymptotic zone not large")
    elseif any(diff(ind_stable_region) .!= 1)
        print(ind_stable_region)
        println("Broken asymptotic zone")
	end
	
	order = mean(order_app[ind_stable_region])

    return (order, order_app)

end

function rate_computation_nl(err::Vector{<:Real}, order::Real, tol::Real=0.2)

    N = length(err)
    N>=3 || throw(ArgumentError("Not enough data (minimum of 3)"))

    rate_app = err[2:end] ./ (err[1:end-1] .^ order)

    stable_region = (rate_app.>0) .& (abs.(grad_vec(rate_app)).<tol) #.& abs.(del2_vec(order_app)).<2*tol
    ind_stable_region = findall(stable_region)

    # Sanity check
	if isempty(ind_stable_region)
		println("No asymptotic zone")
	elseif length(ind_stable_region) < 2
		println("Asymptotic zone not large")
    elseif any(diff(ind_stable_region) .!= 1)
        print(ind_stable_region)
        println("Broken asymptotic zone")
	end
	
	rate = mean(rate_app[ind_stable_region])

    return (rate, rate_app)

end


function order_computation_bissec(err)
    # Least-square fit approximation of the order for the bissec method
    
    nb_iter	=	length(err)
    A		=	hcat(ones(nb_iter-2), log.(err[1:end-2]))
    b		=	log.(err[2:end-1])
    coeff	=	A\b;
    
    order	=	coeff[2];

    return order
    
end


function grad_vec(vec::Vector{<:Real})
    
    N = length(vec)
    N>=3 || throw(ArgumentError("Not enough data (minimum of 3)"))

    diag_center = vcat(-3/2, zeros(N-2), 3/2)
    diag_below = vcat(-1/2*ones(N-2), -2)
    diag_below2 = vcat(zeros(N-3), 1/2)
    diag_above = vcat(2, 1/2*ones(N-2))
    diag_above2 = vcat(-1/2, zeros(N-3))

    M = spdiagm(-2 => diag_below2, -1 => diag_below, 0 => diag_center, 1 => diag_above, 2 => diag_above2)

    grad = M*vec

    return grad
end

function del2_vec(vec::Vector{<:Real})
    
    N = length(vec)
    N>=3 || throw(ArgumentError("Not enough data (minimum of 3)"))

    diag_center = vcat(1, -2*ones(N-2), 1)
    diag_below = vcat(ones(N-2), -2)
    diag_below2 = vcat(zeros(N-3), 1)
    diag_above = vcat(-2, ones(N-2))
    diag_above2 = vcat(1, zeros(N-3))

    M = spdiagm(-2 => diag_below2, -1 => diag_below, 0 => diag_center, 1 => diag_above, 2 => diag_above2)

    del2 = M*vec

    return del2
end