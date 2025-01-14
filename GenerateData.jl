## Generating Data from Lorenz System
using OrdinaryDiffEq
using JLD2, FileIO
import ReservoirComputing: rand_sparse 

# Define the system
function lorenz!(du, u, p, t)
    du[1] = p[1]*(u[2] - u[1])
    du[2] = u[1]*(p[2] - u[3]) - u[2]
    du[3] = u[1]*u[2] - p[3]*u[3]
end

# Define the initial condition
u0 = [6.0,8.0,10.0]
p = [10.0, 28.0, 8.0/3.0]

# Define the time span
tfinal = 40.0
tspan = (0.0, tfinal)

# Solve the ODE
prob = ODEProblem(lorenz!, u0, tspan, p)
sol = solve(prob)

# Extract the solution
dt = 0.01
tdata = 0:dt:tfinal
data = reduce(hcat,sol(tdata).u)

@save "data.jld2" data tdata dt

