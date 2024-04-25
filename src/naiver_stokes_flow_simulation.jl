module naiver_stokes_flow_simulation

using CairoMakie
using SymbolicUtils
using LaTeXStrings
using Symbolics

    greet() = print("Hello World!")


# Newton-Raphson Method
    function newton_raphson(;
        f ::Function,          # The function f 
        df ::Function,         # The derivative of f
        x0 ::Float64,          # The initial guess
        eps ::Float64,         # The error tolerance 
        delta:: Float64,       # The derivative tolerance
        max_iter ::Int64 )     # The maximum number of iterations to do

        fp                                 # For holding the value of f'(x) in each iteration
        fx = f(x0)                         # For holding the value of f(x) in each iteration
        
        x = x0                                 # The current value of x
        for i in 1:max_iter    
            
            fp = fd(x)                         # Get the value of f'(x) in the current of x
            if abs(df(x0)) < delta             # Check if the derivative is small enough
                println("Derivada pequeña")    
                return (i, x, fx)
            end
            d = f(x) / df(x)                   # Calculate the next increment for x
            x = x - d                          # Update the value of x with the new difference   
            fx = f(x)                          # Evaluate f(x) in the current of x

            if abs(fx) < eps                   # Check if the aproximation is close enough
                println("Convergencia")
                return (i, x, fx)              # Return the number of iterations, the value of x and the value of f(x)
            end
        end
    end    



    function guessing_initial_values_equation(x; a = (-1/2401), b = (2/2401), c = (2400/2401))
       
        return a*x^2 + b*x + c
    end   


    """
        bezier_curve(n, p)

        Generate a Bezier curve that decays from 1 to 0 in n steps with a curvature parameter p.
        The curvature parameter p should be between 0 and 1. A value of 0 gives a linear decay,
        while a value of 1 gives the maximum curvature.
    """
    function bezier_curve(n, p)
        t = range(0, 1, length=n)
        x = t
        y = (1 .- t) .^ 3 .* 1.0 + 3 * (1 .- t) .^ 2 .* t .* p + 3 * (1 .- t) .* t .^ 2 .* (1 .- p) + t .^ 3 .* 0.0
        return x, y
    end

    function bezier_curve_function(t, p)
     
        y = (1 - t) ^ 3 * 1 + 3 * (1 - t) ^ 2 * t * p + 3 * (1 - t) * t ^ 2 * (1 - p) + t ^ 3 * 0
        return y
    end
     

    function generate_initial_matrix()              # Generate the initial matrix XO

        p = 1
        matrix = Matrix{Float64}(undef, 5, 50)      # Create a 5x50 matrix with undefined values
        for i in 1:5                                # Populate the matrix
            for j in 1:50
                if i == 1 && i == 5    # Set the first row and last row values to 0
                    matrix[i, j] = 0.1
                else
                    matrix[i, j] = guessing_initial_values_equation(j)  # Set the rest of the values into the matrix to quadratic functions   
                end
            end
        end
        println(matrix)

        scale = ReversibleScale(x -> asinh(x / 2) / log(10), x -> 2sinh(log(10) * x))
        fig, ax, hm = heatmap(transpose(matrix); fontsize = 28, colorscale = scale,  axis = (; xlabel = L"v_x~~donde~~v_{x_0} = 1", ylabel = L"v_y~~donde~~v_{x_i} = 0", title = L"\frac{-1}{2401}x^2 + \frac{2}{2401} + \frac{2400}{2401}") )
        Colorbar(fig[1, 2], hm)

        fig
    end    

    function generate_initial_matrix_with_bezier()              # Generate the initial matrix XO

        p = 0.9
        x, y = bezier_curve(50, p) 
        matrix = Matrix{Float64}(undef, 5, 50)      # Create a 5x50 matrix with undefined values
        for i in 1:5                                # Populate the matrix         
            matrix[i, :] = y # Set the rest of the values into the matrix to quadratic functions   
        end
        println(matrix)

        scale = ReversibleScale(x -> asinh(x / 2) / log(10), x -> 2sinh(log(10) * x))
        fig, ax, hm = heatmap(transpose(matrix); fontsize = 28, colorscale = scale,  axis = (; xlabel = L"v_x~~donde~~v_{x_0} = 1", ylabel = L"v_y~~donde~~v_{x_i} = 0", title = L"\frac{-1}{2401}x^2 + \frac{2}{2401} + \frac{2400}{2401}") )
        Colorbar(fig[1, 2], hm)

        fig
    end 


# Newton-Raphson Method
    function secant(;
        f ::Function,          # The function f 
        a ::Float64,           # The initial point of the interval
        b ::Float64,           # The final point of the interval
        eps ::Float64,         # The error tolerance
        max_iter ::Int64 )     # The maximum number of iterations to do

        fa = f(a)                          # For holding the value of f(a) in each iteration
        fb = f(b)                          # For holding the value of f(b) in each iteration                               
        fx = f(x0)                         # For holding the value of f(x) in each iteration
        delta

        function swapping()
            if abs(fa) > abs(fb)                   # Swap the values of a and b if the absolute value of f(a) is greater than f(b) 
                a, b = b, a
                fa, fb = fb, fa
            end    
        end
        
        swapping()

        for i in 2:max_iter    
            swapping()
            delta = (b - a) / (fb - fa)
            b = a
            fb = fa
            d = d * fa

            if abs(delta) < eps                   # Check if the aproximation is close enough
                println("Convergencia")
                return (i, a, fa)              # Return the number of iterations, the value of x and the value of f(x)
            end

            a = a - delta
            fa = f(a)
        end
    end   


    function plot_distribution_heatmap()
        heatmap(matrix, c=:viridis, aspect_ratio=:equal, colorbar_title="Values", xlabel="X-axis", ylabel="Y-axis", title="Heatmap Plot")
    end


    system_equation = L"u_{i,j} = \frac{u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - \frac{u_{i+1,j}u_{i,j}}{2} + \frac{u_{i-1,j}u_{i,j}}{2} + \frac{u_{i,j-1}w_{i,j}}{2} - \frac{u_{i,j+1}w_{i,j}}{2}}{4}"


    function create_system_equation()

        matrix = Matrix{LaTeXString}(undef, 5, 50) 
        for i in 1:5
            for j in 1:50
                matrix[i, j] = latexstring("u_{$i,$j} = \\frac{u_{$(i+1),$j)} + u_{$(i-1),$j} + u_{$i,$(j+1)} + u_{$i,$(j-1)} - \\frac{u_{$(i+1),$j}u_{$i,$j}}{2} + \\frac{u_{(i-1),$j}u_{$i,$j}}{2} + \\frac{u_{$i,$(j-1)}w_{$i,$j}}{2} - \\frac{u_{$i,$(j+1)}w_{$i,$j}}{2}}{4}")
            end
        end
        println(matrix)        
    end
    


    function do_the_jacobian_matrix()
        # Create a 50x50 matrix of symbolic variables
        n = 50
        jacobian = Symbolics.variables(Symbol("u"), 1:n, 1:n)
        matrix = Matrix{LaTeXString}(undef, 5, 50) 
        for i in 1:5
            for j in 1:50
                matrix[i, j] = latexstring(Differential(jacobian[3, 4])(system_equation))  # Compute ∂h/∂v_3_4
            end
        end
        println(jacobian)   
    end


end # module naiver_stokes_flow_simulation