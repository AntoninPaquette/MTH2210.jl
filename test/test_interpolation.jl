@testset "Interpolation" begin

    poly = [x -> (x-1) * (x+2.5) * (x-pi) * (x+11),
            x -> (x-2)^2 * (x+5)^3 * (x-exp(1))^4]
    degre = [4,9]

    function spline_exact(x)
        if x>=0 && x<1
            f = x^2
        elseif x>=1 && x<3
            f = -x^3 + 4*x^2 - 3*x + 1
        elseif x>=3 && x<=4
            f =  -5*x^2 + 24*x - 26
        else
            f = NaN
        end

        return f
    end

    spline_type_bc1 = [2,2,2,3,3,3,4,4,4]
    spline_type_bc2 = [2,3,4,2,3,4,2,3,4]
    spline_bc1 = [2,2,2,NaN,NaN,NaN,0,0,0]
    spline_bc2 = [-10,NaN,-16,-10,NaN,-16,-10,NaN,-16]



    @testset "Lagrange exact interpolation" begin

        for t=1:length(poly)
            xi = range(-5, 5, length=degre[t]+1)
            yi = poly[t].(xi)
            x_discr = range(minimum(xi), maximum(xi), length=1000)
            y_inter = lagrange(xi,yi,x_discr)
            y_ex = poly[t].(x_discr)
            rel_err = norm(y_ex .- y_inter)/norm(y_ex)

            @test rel_err < 1e-14 
        end


    end

    @testset "Spline exact interpolation" begin

        for t=1:length(spline_type_bc1)
            xi = [0,1,3,4]
            yi = spline_exact.(xi)
            x_discr = range(minimum(xi), maximum(xi), length=1000)
            y_inter = splinec(xi, yi, x_discr, [spline_type_bc1[t], spline_type_bc2[t]], [spline_bc1[t], spline_bc2[t]])
            y_ex = spline_exact.(x_discr)
            rel_err = norm(y_ex .- y_inter)/norm(y_ex)

            @test rel_err < 1e-14 
        end


    end

    @testset "Order Lagrange" begin
        
        fct(x) = exp(2*x)
        
        degre_for_order = [1,2,3,4,5]
        nb_loop = 10
        x_interest = (1/3)^nb_loop
        y_ex = fct.(x_interest)
        a = 0

        tol = 0.2

        for d=1:length(degre_for_order)

            nb_pts = degre_for_order[d] + 1
            err = NaN*ones(nb_loop)


            for t=1:nb_loop
                b = (1/2)^(t-1)
                k = 1:nb_pts
                xi = 1/2*(a+b) .+ 1/2*(b-a)*cos.((2*k .- 1)*pi/(2*nb_pts))
                yi = fct.(xi)
                y_inter = lagrange(xi, yi, [x_interest])
                err[t] = abs(y_ex - y_inter[1])
            end
           
            (order,_) = order_computation(err, 2, tol)

            @test abs(order - (degre_for_order[d] + 1)) < tol
        end


    end

    @testset "Natural curvature spline" begin

        fct(x) = -4*x^3 + pi*x^2 + 11*x - exp(1)
        d2_fct(x) = -24*x + 2*pi
        a = pi/12
        b = 4
        nb_pts = 10
        xi = range(a, b, length=nb_pts)
        yi = fct.(xi)
        x = range(a, b, length=1000)
        y_ex = fct.(x)
        
        Sx = splinec(xi, yi, x, [1,2], [NaN,d2_fct(b)])
        
        rel_err	=	norm(y_ex .- Sx)/norm(y_ex)
        @test rel_err < 1e-14

    end

    @testset "Natural spline with derivative" begin
        
        fct(x) = -4*x^3 + pi*x^2 + 11*x - exp(1)
        d_fct(x) = -12*x^2 + 2*pi*x + 11
        a = -10
        b = pi/12
        nb_pts = 10
        xi = range(a, b, length=nb_pts)
        yi = fct.(xi)
        x = range(a, b, length=1000)
        y_ex = fct.(x)
        
        Sx = splinec(xi, yi, x, [4,1], [d_fct(a),NaN])
        
        rel_err	=	norm(y_ex .- Sx)/norm(y_ex)
        @test rel_err < 1e-14
        
    end


end
