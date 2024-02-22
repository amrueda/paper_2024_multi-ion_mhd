
using OrdinaryDiffEq
using Trixi

###############################################################################
# semidiscretization of the ideal multi-ion MHD equations
equations = IdealMhdMultiIonEquations2D(gammas = (5/3, 1.4),
                                        charge_to_mass = (1, 1/2))

"""
Initial (and exact) solution for the the manufactured solution test. Runs with 
* gammas = (2.0, 4.0),
* charge_to_mass = (2.0, 1.0)
* Domain size: [-1,1]²
"""
function initial_condition_khi(x, t, equations::IdealMhdMultiIonEquations2D)
    (; gammas) = equations
    ca = 0.1 # Alfvén speed
    theta = pi / 3 # Angle of initial magnetic field
    M = 1 # Mach number
    y0 = 1 / 20 #steepness of the shear
    v20 = 0.01 # Parameter of perturbation
    sigma = 0.1 # Parameter of perturbation

    rho1 = rho2 = 0.5
    p1 = 0.5 / gammas[1]
    p2 = 0.5 / gammas[2]
    B1 = ca * cos(theta)
    B2 = 0
    B3 = ca * sin(theta)
    v1 = 0.5 * M * tanh(x[2] / y0)
    v2 = v20 * sin(2 * pi * x[1]) * exp(-x[2]^2 / sigma^2)
    v3 = 0


    prim = SVector{nvariables(equations), real(equations)}([B1, B2, B3, rho1, v1, v2, v3, p1, rho2, v1, v2, v3, p2, 0])
    return prim2cons(prim, equations)
end

@inline function poloidal_field(u, equations::IdealMhdMultiIonEquations2D)
    return (u[1]^2 + u[2]^2)
end

initial_condition = initial_condition_khi
source_terms = source_terms_standard

volume_flux = (flux_central, flux_nonconservative_central)
surface_flux = (flux_lax_friedrichs, flux_nonconservative_central)

solver = DGSEM(polydeg=3, surface_flux=surface_flux,
               volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = (-1.0, -1.0)
coordinates_max = (1.0, 1.0)

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=7,
                n_cells_max=100_000_000,
                periodicity = (true, false))

# Simple implementation of a perfectly conducting slip wall
@inline function boundary_condition_conducting_slip_wall(u_inner, orientation,
                                              direction, x, t,
                                              surface_flux_function,
                                              equations::IdealMhdMultiIonEquations2D)
    # get the appropriate normal vector from the orientation
    if orientation == 1
        u_outer = SVector(-u_inner[1], u_inner[2], u_inner[3], u_inner[4], -u_inner[5], u_inner[6], u_inner[7], u_inner[8], u_inner[9], -u_inner[10], u_inner[11], u_inner[12], u_inner[13], u_inner[14])
    else # orientation == 2
        u_outer = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4], u_inner[5], -u_inner[6], u_inner[7], u_inner[8], u_inner[9], u_inner[10], -u_inner[11], u_inner[12], u_inner[13], u_inner[14])
    end

    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_outer, orientation, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_outer, u_inner, orientation, equations)
    end
    
    return flux
end

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_conducting_slip_wall,
                       y_pos = boundary_condition_conducting_slip_wall)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 20.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, 
                                     interval=analysis_interval, 
                                     extra_analysis_integrals=(poloidal_field, ),
                                     save_analysis=true,
                                     output_directory=joinpath(@__DIR__, "out"))
alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(dt = 0.5,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim,
                                     output_directory = joinpath(@__DIR__, "out"))

cfl=0.5

stepsize_callback = StepsizeCallback(cfl=cfl)

glm_speed_callback = GlmSpeedCallback(glm_scale = 0.5, cfl = cfl)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback,
                        glm_speed_callback)


###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary
