using Plots

# Function definitions (same as before)
dx(params, x, v) = params.r*x*(1-x) - params.b*x*v
dy(params, y, x, v) = params.b*x*v - params.a*y
dv(params, v, x, y) = params.c*y - params.d*x*v - params.k*v

# Simulation function (same as before)
function sim(p, x_0, y_0, v_0; delta_t = 0.1, stop_t=800)
    xs = []
    ys = []
    vs = []
    x = x_0
    y = y_0
    v =v_0

    for step in 0:delta_t:stop_t
        push!(xs, x)
        push!(ys, y)
        push!(vs, v)

        x = x + dx(p, x, v)*delta_t
        y = y + dy(p, y, x, v)*delta_t
        v = v + dv(p, v, x, y)*delta_t
    end

    return xs, ys, vs
end

# Define parameter ranges for b
b_values = 0:0.05:3

# Arrays to store final x, y, v values for each b
final_xs = []
final_ys = []
final_vs = []

# Loop over different values of b
for b in b_values
    # Update the parameter set with the current value of b
    params = (r = 0.5, b = b, a = 0.13, c = 0.1, d = 0.01, k = 0.2)
    
    # Run the simulation
    xs, ys, vs = sim(params, 0.4, 0.3, 0.1)
    
    # Store the final values of x, y, v
    push!(final_xs, xs[end])
    push!(final_ys, ys[end])
    push!(final_vs, vs[end])
end

r_values = 0:0.05:30

for r in r_values
    # Update the parameter set with the current value of b
    params = (r = 1.5, b = r, a = 1.13, c = 0.2, d = 0.01, k = 0.1)
    
    # Run the simulation
    xs, ys, vs = sim(params, 0.4, 0.3, 0.1)
    
    # Store the final values of x, y, v
    push!(final_xs, xs[end])
    push!(final_ys, ys[end])
    push!(final_vs, vs[end])
end

# Plot b vs final x, y, v values
display(Plots.plot(r_values, [final_xs final_ys final_vs], xlabel="b", ylabel="Final values", 
     label=["Final healthy brain cells (x)" "Final infected brain cells (y)" "Final harmful effector cells (v)"],
     title="Effect of varying parameter b on final states", linewidth=3))

Plots.savefig("rrms_project_plots/b_new_params_og_b=0-30.png")