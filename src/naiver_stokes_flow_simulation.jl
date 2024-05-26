module naiver_stokes_flow_simulation

# Backend dfor graphics generation.
using CairoMakie

# For symbolic computation
using SymbolicUtils
using Symbolics

# For latex notation handles in code
using LaTeXStrings


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
                if i == 1 || i == 5    # Set the first row and last row values to 0
                    matrix[i, j] = 0.5
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

        p = 0.5
        x, y = bezier_curve(50, p) 
        matrix = Matrix{Float64}(undef, 5, 50)      # Create a 5x50 matrix with undefined values
        for i in 1:5                                # Populate the matrix         
            matrix[i, :] = y # Set the rest of the values into the matrix to quadratic functions   
        end
        println(matrix)

        scale = ReversibleScale(x -> asinh(x / 2) / log(10), x -> 2sinh(log(10) * x))
        fig, ax, hm = heatmap(transpose(matrix); fontsize = 28, colorscale = scale,  axis = (; xlabel = L"v_x~~donde~~v_{x_0} = 1", ylabel = L"v_y~~donde~~v_{x_i} = 0") )
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


    system_equation = L"  "


    function create_system_equation()

        matrix = Matrix{LaTeXString}(undef, 5, 50) 
        for i in 1:5
            for j in 1:50
                matrix[i, j] = latexstring("u_{$i,$j} = \\frac{u_{$(i+1),$j} + u_{$(i-1),$j} + u_{$i,$(j+1)} + u_{$i,$(j-1)} - \\frac{u_{$(i+1),$j}u_{$i,$j}}{2} + \\frac{u_{$(i-1),$j}u_{$i,$j}}{2} + \\frac{u_{$i,$(j-1)}w_{$i,$j}}{2} - \\frac{u_{$i,$(j+1)}w_{$i,$j}}{2}}{4}")
            end
        end
        #println(matrix)
        matrix        
    end
    


    function do_the_jacobian_matrix()
        # Create a 50x50 matrix of symbolic variables
        n = 50
        jacobian = Symbolics.variables(Symbol("u"), 1:5, 1:50)
        matrix = Matrix{LaTeXString}(undef, 5, 50) 
        system_equation_matrix = create_system_equation() 
        for i in 1:5
            for j in 1:50
                matrix[i, j] = latexstring(expand_derivatives(Differential(jacobian[i, j])(system_equation_matrix[ i, j])) )
                println("La derivada parcial en in index $i, $j = $(matrix[i, j])") 
                f1Exp = expand_derivatives(Differential(jacobian[i, j] )(system_equation_matrix[ i, j]))
                println( "Expandido $f1Exp")
            end
        end
        #println(jacobian)   
    end



    @variables i j u[i, j] w[i, j]
    function derivate_partial()
        # Define symbolic variables
    

        # Define the equation as a LaTeX string
        eq_latex = "u_{i,j} = \\frac{u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - \\frac{u_{i+1,j}u_{i,j}}{2} + \\frac{u_{i-1,j}u_{i,j}}{2} + \\frac{u_{i,j-1}w_{i,j}}{2} - \\frac{u_{i,j+1}w_{i,j}}{2}}{4}"

        # Convert the LaTeXString to a symbolic expression
       # eq_symbolic = latexstr(eq_latex)

        # Parse the string into a symbolic expression
        eq = SymbolicUtils.parse(SymbolicUtils.Sym, eq_latex)

        # Compute partial derivatives with respect to i and j
        du_di = diff(eq, i)
        du_dj = diff(eq, j)

        # Print the partial derivatives
        println("Partial derivative with respect to i:")
        display(du_di)
        println("\nPartial derivative with respect to j:")
        display(du_dj)
    end


    function other()

        # Define the expression
        @variables i j u[i, j] w[i, j]
        
        # Given expression
        expr = (u[i + 1, j] + u[i - 1, j] + u[i, j + 1] + u[i, j - 1] -
                u[i + 1, j] * u[i, j] / 2 + u[i - 1, j] * u[i, j] / 2 +
                u[i, j - 1] * w[i, j] / 2 - u[i, j + 1] * w[i, j] / 2) / 4
        
        # Compute partial derivatives
        ∂u_∂i = diff(expr, u[i, j], i)
        ∂u_∂j = diff(expr, u[i, j], j)
        
        # Display the results
        println("Partial derivative with respect to i:")
        display(∂u_∂i)
        
        println("\nPartial derivative with respect to j:")
        display(∂u_∂j)
        
    end


    function other_f()
        # Define symbolic variables
        @variables i j u[i, j] w[i, j]

        # Define the equation as a LaTeX string
        eq_latex = L"u_{i,j} = \frac{u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - \frac{u_{i+1,j}u_{i,j}}{2} + \frac{u_{i-1,j}u_{i,j}}{2} + \frac{u_{i,j-1}w_{i,j}}{2} - \frac{u_{i,j+1}w_{i,j}}{2}}{4}"

        expr = string(eq_latex)
        # Convert the LaTeXString to a symbolic expression
        #eq_symbolic = parse(expr)

        # Parse the string into a symbolic expression
        eq = Symbolics.operation(eq_latex)

        # Compute the Jacobian matrix
        jacobian = Symbolics.jacobian([eq], [u[i, j] for i in 1:5, j in 1:50])

        # Print the Jacobian matrix
        display(jacobian)
    
    end     

    @variables i j u[i, j] w[i, j] ip1 im1 jp1 jm1

    @variables ip1(i) = i + 1
    @variables im1(i) = i - 1
    @variables jp1(j) = j + 1
    @variables jm1(j) = j - 1

    function f_symbolic()
        return (u[ip1, j] + u[im1, j] + u[i, jp1] + u[i, jm1] -
        (u[ip1, j] * u[i, j]) / 2 + (u[im1, j] * u[i, j]) / 2 +
        (u[i, jm1] * w[i, j]) / 2 - (u[i, jp1] * w[i, j]) / 2) / 4
    end

    function ot()
        eq = f_symbolic()
        jacobian = Symbolics.jacobian([eq], [u[i, j] for i in 1:5, j in 1:50])
        display(jacobian)
    end    

    using ForwardDiff

    # Define the discretized function
    function f(u, w, i, j)
        ui = u[i, j]
        uip1 = u[i+1, j]
        uim1 = u[i-1, j]
        ujp1 = u[i, j+1]
        ujm1 = u[i, j-1]
        wij = w[i, j]
    
        return (uip1 + uim1 + ujp1 + ujm1 - (uip1 * ui) / 2 + (uim1 * ui) / 2 + (ujm1 * wij) / 2 - (ujp1 * wij) / 2) / 4
    end
   function d() 
    # Define the grid size
        grid_size = (5, 50)
        
        # Initialize arrays for u and w
        u = rand(grid_size)
        w = rand(grid_size)
        
        # Compute the Jacobian matrix numerically
        jacobian = zeros(grid_size)
        for i in 1:grid_size[1], j in 1:grid_size[2]
            jacobian[i, j] = ForwardDiff.derivative(x -> f(u, w, i, j), u[i, j])
        end
        
        # Print the Jacobian matrix
        display(jacobian)

    end


    function finite_difference_jacobian(u_function, u; h=1e-6,)
        m, n = size(u)
        J = zeros(m*n, m*n) 
    
        for i in 1:m
            for j in 1:n
                u_plus_h = copy(u)
                u_plus_h[i, j] += h
    
                u_minus_h = copy(u)
                u_minus_h[i, j] -= h
    
                F_plus_h = u_function(u_plus_h)
                F_minus_h = u_function(u_minus_h)
    
                J[:, (i-1)*n+j] = (F_plus_h[:] - F_minus_h[:]) / (2 * h)
            end
        end
    
        return J
    end


    function forward_difference_jacobian(u_function, u; h=1e-6)
        m, n = size(u)
        J = zeros(m*n, m*n)
    
        for i in 1:m
            for j in 1:n
                u_plus_h = copy(u)
                u_plus_h[i, j] += h
    
                F_plus_h = u_function(u_plus_h)
    
                J[:, (i-1)*n+j] = (F_plus_h[:] - u_function(u)[:]) / h
            end
        end
    
        return J
    end

    function backward_difference_jacobian(u_function, u; h=1e-6)
        m, n = size(u)
     
        p = 0.5
        x, y = bezier_curve(50, p) 
        for i in 1:5                                # Populate the matrix         
            u[i, :] = y # Set the rest of the values into the matrix to quadratic functions   
        end

        println(u)
        J = zeros(m*n, m*n)
    
        for i in 2:m-1
            for j in 2:n-1
                u_minus_h = copy(u)
                u_minus_h[i, j] -= h
    
                F_minus_h = u_function(u_minus_h, 1, i, j)
    
                J[:, (i-1)*n+j] = (u_function(u)[:] - F_minus_h[:]) / h
            end
        end
        return J
    end

    using FiniteDifferences

# Define the discretized function
function f(u, w, i, j)
    ui = u[i, j]
    uip1 = u[i+1, j]
    uim1 = u[i-1, j]
    ujp1 = u[i, j+1]
    ujm1 = u[i, j-1]
    wij = 1 #w[i, j]

    return (uip1 + uim1 + ujp1 + ujm1 - (uip1 * ui) / 2 + (uim1 * ui) / 2 + (ujm1 * wij) / 2 - (ujp1 * wij) / 2) / 4
end


function jacobian_finite_differenecs()
# Define the grid size
grid_size = (5, 50)
n = grid_size[1] * grid_size[2]

# Initialize arrays for u and w
u = rand(grid_size)
w = rand(grid_size)

# Compute the Jacobian matrix using finite differences
jacobian = zeros(n, n)
for i in 1:grid_size[1], j in 1:grid_size[2]
    idx = (i - 1) * grid_size[2] + j
    jacobian[:, idx] = FiniteDifferences.finite_difference_jacobian(u -> f(u, w, i, j), u)
end

# Print the Jacobian matrix
display(jacobian)

# Get the machine epsilon
epsilon = eps()
end


# Define a function to calculate the step size based on machine epsilon
function calculate_step(epsilon)
    return sqrt(epsilon)
end
function cal_eps()
# Calculate the step size
    step_size = calculate_step(epsilon)
end
# Use the step size in your computations


end # module naiver_stokes_flow_simulation