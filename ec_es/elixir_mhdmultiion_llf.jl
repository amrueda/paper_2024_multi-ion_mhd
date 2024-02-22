
using OrdinaryDiffEq
using Trixi

###############################################################################
"""
  electron_pressure_alpha(u, equations::IdealMhdMultiIonEquations2D)
Returns a fraction (alpha) of the total ion pressure for the electron pressure.
"""
function electron_pressure_alpha(u, equations::IdealMhdMultiIonEquations2D)
    alpha = 0.2
    prim = cons2prim(u, equations)
    p_e = zero(u[1])
    for k in eachcomponent(equations)
        _, _, _, _, p_k = get_component(k, prim, equations)
        p_e += p_k
    end
    return alpha * p_e
end
"""
Convert conservative variables to entropy
"""
@inline function entropy_density(u, equations::IdealMhdMultiIonEquations2D)
    @unpack gammas = equations
    B1, B2, B3 = magnetic_field(u, equations)
    psi = divergence_cleaning_field(u, equations)

    prim = cons2prim(u, equations)
    entropy_density = zero(eltype(u))
    rho_p_plus = zero(real(equations))
    for k in eachcomponent(equations)
        rho, v1, v2, v3, p = get_component(k, prim, equations)
        s = log(p) - gammas[k] * log(rho)
        
        entropy_density -= rho * s / (gammas[k] - 1) 
    end

    return entropy_density
end

# semidiscretization of the ideal MHD equations
equations = IdealMhdMultiIonEquations2D(gammas = (2.0, 4.0),
                                        charge_to_mass = (2.0, 1.0),
                                        electron_pressure = electron_pressure_alpha)

initial_condition = initial_condition_weak_blast_wave

volume_flux = (flux_ruedaramirez_etal, flux_nonconservative_ruedaramirez_etal)
surface_flux = (flux_lax_friedrichs, flux_nonconservative_central)

solver = DGSEM(polydeg = 3, surface_flux = surface_flux,
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = (-2.0, -2.0)
coordinates_max = (2.0, 2.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 1_000_000)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_standard)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.4)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1

for cfl in [0.400, 0.200, 0.100, 0.050, 0.025]
    analysis_callback = AnalysisCallback(semi, 
                                        interval = analysis_interval,
                                        extra_analysis_integrals=(entropy_density, ),
                                        save_analysis=true,
                                        output_directory=joinpath(@__DIR__, "out_llf_cfl" * string(cfl)))
    alive_callback = AliveCallback(analysis_interval = analysis_interval)

    save_solution = SaveSolutionCallback(dt = 0.1, # interval=100,
                                        save_initial_solution = true,
                                        save_final_solution = true,
                                        solution_variables = cons2prim,
                                        output_directory=joinpath(@__DIR__, "out_llf_cfl" * string(cfl)))

    stepsize_callback = StepsizeCallback(cfl = cfl)

    glm_speed_callback = GlmSpeedCallback(glm_scale = 0.5, cfl = cfl)

    callbacks = CallbackSet(summary_callback,
                            analysis_callback, alive_callback,
                            save_solution,
                            stepsize_callback,
                            glm_speed_callback)

    ###############################################################################
    # run the simulation

    sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
                dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
                save_everystep = false, callback = callbacks);
    summary_callback() # print the timer summary
end
