using Plots

module naiver_stokes_flow_simulation

greet() = print("Hello World!")


# Newton-Raphson Method
function newton_raphson(;
    f ::Function,          # The function f 
    df ::Function,         # The derivative of f
    x0 ::Float64,          # The initial guess
    eps ::Float64,         # The error tolerance 
    delta:: Float64,       # The derivative tolerance
    max_iter ::Int64 )     # The maximum number of iterations to do

    var fp                                 # For holding the value of f'(x) in each iteration
    var fx = f(x0)                         # For holding the value of f(x) in each iteration
    
    x = x0                                 # The current value of x
    for i in 1:max_iter    
        
        fp = fd(x)                         # Get the value of f'(x) in the current of x
        if abs(df(x0)) < delta             # Check if the derivative is small enough
            println("Derivada pequeña")    
            return
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
    
    f(x) = ax^2 + bx + c
end   

     

function generate_initial_matrix()              # Generate the initial matrix XO

    matrix = Matrix{Float64}(undef, 5, 50)      # Create a 5x50 matrix with undefined values
    for i in 1:5                                # Populate the matrix
        for j in 1:50
            if i != 1 && (j == 1 || j == 5)     # Set the first row and last row values to 0
                matrix[i, j] = 0.0
            else
                matrix[i, j] = guessing_initial_values_equation(j)  # Set the rest of the values into the matrix to quadratic functions   
            end
        end
    end
end    

# Newton-Raphson Method
function secant(;
    f ::Function,          # The function f 
    a ::Float64,           # The initial point of the interval
    b ::Float64,           # The final point of the interval
    eps ::Float64,         # The error tolerance
    max_iter ::Int64 )     # The maximum number of iterations to do

    var fa = f(a)                          # For holding the value of f(a) in each iteration
    var fb = f(b)                          # For holding the value of f(b) in each iteration                               
    var fx = f(x0)                         # For holding the value of f(x) in each iteration
    var delta

    function swapping ()
        if abs(fa) > abs(fb)                   # Swap the values of a and b if the absolute value of f(a) is greater than f(b) 
            a, b = b, a
            fa, fb = fb, fa
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


end # module naiver_stokes_flow_simulation