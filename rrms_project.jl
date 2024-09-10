using Plots
using PyPlot
using Colors

# r = growth rate of healthy cells
# b = attack rate on x(healthy) by v(virus) secreted by y(infected)
# c = secrete rate how y secrete v
# a = die rate of y
# d = rate at which x affect v
# k = die rate of v

# ODEs that define the behaviour of the system

# healthy brain cells - x
dx(params, x, v) = params.r*x*(1-x)-params.b*x*v 
# infected brain cells - y
dy(params, y, x, v) = params.b*x*v - params.a*y
# effector cells produced by y - v
dv(params, v, x, y) = params.c*y - params.d*x*v - params.k*v

# 1st parameters used in the model
params1= (
    r = 1.5, b = 0.25, a = 1.13, c = 0.2, d = 0.01, k =0.1
)

# params extreme
params_ex = (
    r = 1, b = 1, a = 1, c = 1, d = 0.1, k =0.1
)

# original initial values used in the model
inits1 = (x = 0.4, y = 0.3, v = 0.1)

# definition of some extreme initial values to verify the the statement about the initial values
inits2 = (x = 0.1, y = 0.4, v = 0.2)
inits3 = (x = 1, y = 1, v = 1)
inits4 = (x = 0.05, y = 0.8, v = 0.5)

# function simulate which models the dynamics of the biological system in time based on the previously defined ODEs
function simulate(p, i, stop_t; delta_t = 0.1)
    xs = []
    ys = []
    vs = []
    x = i.x
    y = i.y
    v = i.v


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

# here we need to change the parameters of the function simulate to use either different initial values or different parameters defined further in the code
xs_demo, ys_demo, vs_demo = simulate(params1, inits1, 200)
display(Plots.plot([xs_demo, ys_demo, vs_demo], xlabel="time", title="Number of different cell types in time", label=["healthy brain cells - x(t)" "infected brain cells - y(t)" "harmful effector cells - v(t)"], linewidth=3))

Plots.savefig("rrms_project_plots/basic_inits1.png")

# run simulation with extreme parameters
xs_ex, ys_ex, vs_ex = simulate(params_ex, inits1, 200)
display(Plots.plot([xs_ex, ys_ex, vs_ex], xlabel="time", title="Number of different cell types in time", label=["healthy brain cells - x(t)" "infected brain cells - y(t)" "harmful effector cells - v(t)"], linewidth=3))
Plots.savefig("ms_model/rrms_project_plots/extreme_params_inits1.png")

params2= (
    r = 0.5, b = 0.28, a = 0.13, c = 0.1, d = 0.01, k =0.1
)

xs_1, ys_1, vs_1 = simulate(params2, inits4, 400)
display(Plots.plot([xs_1, ys_1, vs_1], xlabel="time", title="Number of different cell types in time ", label=["healthy brain cells - x(t)" "infected brain cells - y(t)" "harmful effector cells - v(t)"], linewidth=3))
Plots.savefig("rrms_project_plots/obratene_inits4.png")


params3= (
    r = 0.5, b = 1.251, a = 0.13, c = 0.1, d = 0.01, k =0.1
)

xs_2, ys_2, vs_2 = simulate(params3, inits3, 800)
display(Plots.plot([xs_2, ys_2, vs_2], xlabel="time", title="RRMS ", label=["healthy brain cells" "infected brain cells" "harmful effector cells"], linewidth=3))


Plots.savefig("rrms_project_plots/oscillations_inits3.png")
# Plots.plot!(ts_demo, es_demo)

b_values

# xs_3, ys_3, vs_3, b_values = sim_b( 0.4, 0.3, 0.1)
# display(Plots.plot(b_values, [xs_3 ys_3 vs_3], xlabel="b", ylabel="Final values", 
#      label=["Final healthy brain cells (x)" "Final infected brain cells (y)" "Final harmful effector cells (v)"],
#      title="Effect of varying parameter b on final states", linewidth=3))
     
     
# Plots.savefig("rrms_project_plots/b_varying.png")