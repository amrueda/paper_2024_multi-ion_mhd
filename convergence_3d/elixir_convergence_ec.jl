
using OrdinaryDiffEq
using Trixi

###############################################################################
"""
  electron_pressure_alpha(u, equations::IdealGlmMhdMultiIonEquations3D)
Returns a fraction (alpha) of the total ion pressure for the electron pressure.
"""
function electron_pressure_alpha(u, equations::IdealGlmMhdMultiIonEquations3D)
    alpha = 0.2
    prim = cons2prim(u, equations)
    p_e = zero(u[1])
    for k in eachcomponent(equations)
        _, _, _, _, p_k = get_component(k, prim, equations)
        p_e += p_k
    end
    return alpha * p_e
end
# semidiscretization of the ideal multi-ion MHD equations
equations = IdealGlmMhdMultiIonEquations3D(gammas = (2.0, 4.0),
                                            charge_to_mass = (2.0, 1.0),
                                            electron_pressure = electron_pressure_alpha)


"""
Initial (and exact) solution for the the manufactured solution test. Runs with 
* gammas = (2.0, 4.0),
* charge_to_mass = (2.0, 1.0)
* Domain size: [-1,1]²
"""
function initial_condition_manufactured_solution(x, t, equations::IdealGlmMhdMultiIonEquations3D)
    am = 0.1
    om = π
    h = am * sin(om * (x[1] + x[2] + x[3] - t)) + 2
    hh1 = am * 0.4 * sin(om*(x[1] + x[2] + x[3] - t)) + 1
    hh2 = h-hh1

    rho_1 = hh1
    rhou_1 = hh1
    rhov_1 = hh1
    rhow_1 = 0.1*hh1
    rhoe_1 = 2*hh1^2+hh1
    rho_2 = hh2
    rhou_2 = hh2
    rhov_2 = hh2
    rhow_2 = 0.1*hh2
    rhoe_2 = 2*hh2^2+hh2
    B1 = 0.5*h
    B2 = -0.25*h
    B3 = -0.25*h

    return SVector{nvariables(equations), real(equations)}([B1, B2, B3, rho_1, rhou_1, rhov_1, rhow_1, rhoe_1, rho_2, rhou_2, rhov_2, rhow_2, rhoe_2, 0])
end

"""
Source term that corresponds to the manufactured solution test. Runs with 
* gammas = (2.0, 4.0),
* charge_to_mass = (2.0, 1.0)
* Domain size: [-1,1]²
"""
function source_terms_manufactured_solution_pe(u, x, t, equations::IdealGlmMhdMultiIonEquations3D)
    am = 0.1
    om = pi
    h1 = am * sin(om * (x[1] + x[2] + x[3] - t))
    hx = am * om * cos(om * (x[1] + x[2] + x[3] - t))

    s1 = (11*hx)/25;
    s2 = (30615*hx * h1^2+156461*hx * h1+191990*hx)/(35000*h1+75000);
    s3 = (30615*hx * h1^2+156461*hx * h1+191990*hx)/(35000*h1+75000);
    s4 = (30615*hx * h1^2+142601*hx * h1+162290*hx)/(35000*h1+75000);
    s5 = (4735167957644739545*hx * h1^2+22683915240114795103*hx * h1+26562869799135852870*hx)/(1863579030125050000*h1+3993383635982250000);
    s6 = (33*hx)/50;
    s7 = (63915*hx * h1^2+245476*hx * h1+233885*hx)/(17500*h1+37500);
    s8 = (239655*hx * h1^2+903512*hx * h1+838345*hx)/(35000*h1+75000);
    s9 = (16005*hx * h1^2+57602*hx * h1+52645*hx)/(35000*h1+75000);
    s10 = (1619415*hx * h1^2+6083946*hx * h1+5629335*hx)/(175000*h1+375000);
    s11 = (11*hx)/20;
    s12 = -((11*hx)/40);
    s13 = -((11*hx)/40);

    s = SVector{nvariables(equations), real(equations)}([s11, s12, s13, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, 0])
    S_std = source_terms_lorentz(u, x, t, equations::IdealGlmMhdMultiIonEquations3D)

    return SVector{nvariables(equations), real(equations)}(S_std .+ s)
end

initial_condition = initial_condition_manufactured_solution
source_terms = source_terms_manufactured_solution_pe

volume_flux = (flux_ruedaramirez_etal, flux_nonconservative_ruedaramirez_etal)
surface_flux = (flux_ruedaramirez_etal, flux_nonconservative_ruedaramirez_etal)

solver = DGSEM(polydeg=3, surface_flux=surface_flux,
               volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = (-1.0, -1.0, -1.0)
coordinates_max = (1.0, 1.0, 1.0)

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=3,
                n_cells_max=1_000_000)


semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)
alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=100,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim
                                     output_directory=joinpath(@__DIR__, "out"))

cfl = 0.5
stepsize_callback = StepsizeCallback(cfl = cfl)

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
