@testset "edo" begin
    
    euler_algo = algo(euler,1)
    eulermod_algo = algo(eulermod,2)
    ptmilieu_algo = algo(ptmilieu,2)
    rk4_algo = algo(rk4,4)

    ode_integrator = (euler_algo, eulermod_algo, ptmilieu_algo, rk4_algo)

    tspan = [0,5]

    fct_degre1 = ode_prob((t,y) -> 3, 1, tspan, t -> 3*t .+ 1, 1)
    fct_degre2 = ode_prob((t,y) -> -4*t, 1, tspan, t -> -2*t.^2 .+ 1, 2)
    fct_degre4 = ode_prob((t,y) -> 5*t.^3, 1, tspan, t -> 5/4 * t.^4 .+ 1, 4)
    fct_scalar = ode_prob((t,y) -> 2*y .- t .+ 4, 1, tspan, t -> -7/4 + 1/2*t + 11/4*exp(2*t), Inf)
    fct_system = ode_prob((t,y) -> [-2 1;1 -2]*y .+ [2*exp(-t),3*t], [2,3], tspan, t -> -7/6*[1,-1] .* exp.(-3*t) .+ 4*[1,1] .* exp.(-t) .+ 1/2*[1,-1] .* exp.(-t) .+ [1,1] .* t .* exp(-t) .+ [1,2] .* t .- 1/3*[4,5], Inf)

    ode_problems_ex = (fct_degre1, fct_degre2, fct_degre4)


    for ode_int in ode_integrator

        @testset "Verification of $(String(Symbol(ode_int.name)))" begin
            
            @testset "Exact integration" begin
                
                for prob in ode_problems_ex
                    if prob.degre <= ode_int.order
                        
                        (temps , y)	=	ode_int.name(prob.fct, prob.tspan, prob.y0, 100);
                        sol_ex      =   prob.ex_sol.(temps)
                        erreur_rel	=	norm(y[:] .- sol_ex, Inf)/norm(sol_ex, Inf);
                        
                        @test erreur_rel < 1e-14

                    end
                end

            end

            @testset "Order of convergence (scalar)" begin
                
                nb_eval = 6
                nb_pas_init = 100
                nb_pas = 2 .^ (0:nb_eval-1) * nb_pas_init
                erreur = Inf*ones(nb_eval)

                tol = 0.2

                for t=1:nb_eval
                    (temps, y) = ode_int.name(fct_scalar.fct, fct_scalar.tspan, fct_scalar.y0, nb_pas[t])
                    erreur[t] = norm(y[:] - fct_scalar.ex_sol.(temps),Inf)
                end

                (ordre_app, _) = order_computation(erreur,2,tol);
                
                @test abs(ordre_app - ode_int.order) < tol

            end

            @testset "Order of convergence (system)" begin
                
                nb_eval = 6
                nb_pas_init = 100
                nb_pas = 2 .^ (0:nb_eval-1) * nb_pas_init
                erreur = Inf*ones(nb_eval)

                tol = 0.2

                for t=1:nb_eval
                    (temps, y) = ode_int.name(fct_system.fct, fct_system.tspan, fct_system.y0, nb_pas[t])
                    erreur[t] = norm(y - reduce(hcat,fct_system.ex_sol.(temps)),Inf)
                end

                (ordre_app, _) = order_computation(erreur,2,tol);
                
                @test abs(ordre_app - ode_int.order) < tol

            end

            @testset "Float and Vector(Float)" begin
                
                nb_pas = 17

                (temps, y) = ode_int.name(fct_scalar.fct, fct_scalar.tspan, fct_scalar.y0, nb_pas)
    
                @test length(temps) == (nb_pas+1)
                @test size(y) == (1,nb_pas+1)
    
                (temps, y) = ode_int.name(fct_system.fct, fct_system.tspan, fct_system.y0, nb_pas)
    
                @test size(y) == (2,nb_pas+1)
        
            end


        end
    end

end