@testset "Non linear" begin

    fct = [x -> x^2 - 10, x -> exp(x) - x^3 - (exp(pi) - pi^3)]
    dfct = [x -> 2*x, x -> exp(x) - 3*x^2]
    x0 = [2, 2.65]
    xb = [4, 3.25]
    root = [sqrt(10), pi]

    g = [x -> -x^2/10 + x + 1, x -> -x^2/6 + x + 9/6]
    pt_fix = [sqrt(10), 3]
    x0_g = [1,1]
    order_g = [1,2]
    rate_g = [-2*sqrt(10)/10 + 1, 1/6]

    fct_system(x) = [5*sin(0.1*x[1]*x[2]) - x[3] - (5*sin(-0.2)-5),
                      x[1]^2 + x[2]^2 + x[3]^2 - 30,
                      x[1] - x[2] - x[3] + 2]
    jac_fct_system(x) = [5*0.1*x[2]*cos(0.1*x[1]*x[2])  5*0.1*x[1]*cos(0.1*x[1]*x[2])  -1 ;
                  2*x[1]  2*x[2]  2*x[3] ;
                  1  -1  -1]
    root_system = [1, -2, 5]		  
    x0_system	= [1.5, -3, 4]
    
    @testset "Bissection" begin

        for t=1:length(fct)

            @testset "Function $t" begin

                (app, err) = bissec(fct[t], x0[t], xb[t], 200 ,1e-12)

                @test abs(app[end] - root[t])/(abs(root[t])) < 1e-12

                order = order_computation_bissec(err)

                @test abs(order - 1) < 0.1

            end
            
        end
        
    end 


    @testset "Secante" begin

        for t=1:length(fct)

            @testset "Function $t" begin

                (app, err) = secante(fct[t], x0[t], xb[t], 20 ,1e-14)

                @test abs(app[end] - root[t])/(abs(root[t])) < 1e-12

                (order, order_app) = order_computation_nl(err)
                
                @test abs(order - (1 + sqrt(5))/2) < 0.2

            end
            
        end
        
    end 


    @testset "newton1D" begin

        for t=1:length(fct)

            @testset "Function $t" begin

                (app, err) = newton1D(fct[t], dfct[t], x0[t], 20 ,1e-13)

                @test abs(app[end] - root[t])/(abs(root[t])) < 1e-12

                (order, order_app) = order_computation_nl(err)

                @test any(abs.(order_app .- 2) .< 0.1)

            end
            
        end

        @testset "Multiple roots" begin
            
            fct_mult(x) = x*sin(x)^2
            dfct_mult(x) = sin(x)^2 + 2*x*sin(x)*cos(x)

            (approx_r1, err_r1) = newton1D(fct_mult, dfct_mult, 1, 200, 1e-12)
            (approx_r2, err_r2) = newton1D(fct_mult, dfct_mult, 3, 200, 1e-12)

            @test abs(approx_r1[end] - 0) < 1e-10
            @test abs(approx_r2[end] - pi) < 1e-10

            (order_r1, _) = order_computation_nl(err_r1)
            (order_r2, _) = order_computation_nl(err_r2)

            @test abs(order_r1 - 1) < 0.1
            @test abs(order_r2 - 1) < 0.1

            (rate_r1, _) = rate_computation_nl(err_r1,1)
            (rate_r2, _) = rate_computation_nl(err_r2,1)

            @test abs(rate_r1 - (3-1)/3) < 0.1
            @test abs(rate_r2 - (2-1)/2) < 0.1
        end
        
    end 


    @testset "ptfixes" begin

        for t=1:length(g)

            @testset "Function $t" begin

                (app, err) = ptfixes(g[t], x0_g[t], 100, 1e-13)

                @test abs(app[end] - pt_fix[t])/(abs(pt_fix[t])) < 1e-12

                (order, order_app) = order_computation_nl(err)

                @test abs(order - order_g[t]) .< 0.1

                (rate, rate_app) = rate_computation_nl(err, order_g[t])

                @test abs(rate - rate_g[t]) .< 0.1

            end
            
        end
        
    end 

    @testset "newtonND" begin

        (app, err) = newtonND(fct_system, x0_system, 20, 1e-12)

        @test norm(app[end,:] .- root_system) < 1e-10

        (order, order_app) = order_computation_nl(err)

        @test any(abs.(order_app .- 2) .< 0.1)

    end

    @testset "newtonNDder" begin

        (app, err) = newtonNDder(fct_system, jac_fct_system, x0_system, 20, 1e-12)

        @test norm(app[end,:] .- root_system) < 1e-10

        (order, order_app) = order_computation_nl(err)

        @test any(abs.(order_app .- 2) .< 0.1)    

    end

end