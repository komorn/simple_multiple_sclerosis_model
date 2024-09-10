using Plots

# ODEs definitions (same as before)
dx(params, x, v) = params.r*x*(1-x) - params.b*x*v
dy(params, y, x, v) = params.b*x*v - params.a*y
dv(params, v, x, y) = params.c*y - params.d*x*v - params.k*v

# simulation function (same as before)
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

# define parameter ranges for b, r, c, a
b_values = 0:0.05:30
r_values = 0:0.05:30
c_values = 0:0.05:30
a_values = 0:0.05:30

# Arrays to store final x, y, v values for each b
final_xs = []
final_ys = []
final_vs = []

for b in b_values
    # update the parameter set with the current value of b
    params = (r = 1.5, b = b, a = 1.13, c = 0.2, d = 0.01, k = 0.1) #r = 1.5, b = 0.25, c = 0.2, a = 1.13
    params2 = (r = 0.5, b = b, a = 0.13, c = 0.1, d = 0.01, k = 0.1) #r = 0.5, b = 0.28, c= 0.1, a = 0.13

    
    # simulate
    xs, ys, vs = sim(params, 0.4, 0.3, 0.1)
    
    # store the final values of x, y, v
    push!(final_xs, xs[end])
    push!(final_ys, ys[end])
    push!(final_vs, vs[end])
end



# plottinf
display(Plots.plot(b_values, [final_xs final_ys final_vs], xlabel="b", ylabel="Final values", 
     label=["Final healthy brain cells (x)" "Final infected brain cells (y)" "Final harmful effector cells (v)"],
     title="Effect of varying parameter b on final states", linewidth=3))

Plots.savefig("rrms_project_plots/b_params_og_c=0-30.png")