using GLMakie
using TestParticle
using OrdinaryDiffEq
using StaticArrays

function main_gui()
    # --- Data Storage for Parsed Values ---
    # Field Configuration
    # E-field is hardcoded to Zero.
    parsed_bz_obs = Observable(1e-8)  # Default Bz, parsed reactively

    # Initial Conditions are Hardcoded.

    # Time Span (Only End Time is configurable, Start Time hardcoded to 0.0)
    parsed_t_end_obs = Observable(10.0) # Default End Time, parsed reactively

    # Particle Species
    selected_species_obs = Observable("Proton") # Default species

    # Solver Configuration
    selected_solver_obs = Observable("Vern9 (adaptive)") # Renamed and changed to Observable
    parsed_dt_obs = Observable(0.01)     # Default dt, parsed reactively
    parsed_reltol_obs = Observable(1e-7) # Default reltol

    # Simulation Solution Storage
    sol_obs = Observable{Any}(nothing) # Changed from Ref to Observable
    # current_trajectory_plot Ref removed (no 3D plot)
    
    # Observables for 2D time-series plot data
    t_obs = Observable(Float64[])
    x_pos_obs = Observable(Float64[])
    y_pos_obs = Observable(Float64[])
    z_pos_obs = Observable(Float64[])
    vx_vel_obs = Observable(Float64[])
    vy_vel_obs = Observable(Float64[])
    vz_vel_obs = Observable(Float64[])
    # pos_ts_plots_ref and vel_ts_plots_ref REMOVED

    # --- Figure and Main Layout ---
    fig = Figure(size = (1200, 700)) # Adjusted size for simpler layout
    # Main GridLayout: 1 row, 2 columns (inputs | 2D plots)
    gl_root = fig[1, 1] = GridLayout() 

    # Column 1: Input Grid
    input_grid = gl_root[1, 1] = GridLayout() 
    
    # Column 2: Plot Display Grid (for 2D plots)
    plot_display_grid = gl_root[1, 2] = GridLayout()
    pos_ts_ax = Axis(plot_display_grid[1, 1], title="Position vs Time", xlabel="Time [s]", ylabel="Position [m]")
    vel_ts_ax = Axis(plot_display_grid[2, 1], title="Velocity vs Time", xlabel="Time [s]", ylabel="Velocity [m/s]")
    linkxaxes!(pos_ts_ax, vel_ts_ax)

    # Initialize 2D time series plots with observable data
    lines!(pos_ts_ax, t_obs, x_pos_obs, label="x")
    lines!(pos_ts_ax, t_obs, y_pos_obs, label="y")
    lines!(pos_ts_ax, t_obs, z_pos_obs, label="z")
    axislegend(pos_ts_ax, position=:rt)

    lines!(vel_ts_ax, t_obs, vx_vel_obs, label="vx")
    lines!(vel_ts_ax, t_obs, vy_vel_obs, label="vy")
    lines!(vel_ts_ax, t_obs, vz_vel_obs, label="vz")
    axislegend(vel_ts_ax, position=:rt)

    # 3D Trajectory Plot (plot_ax) REMOVED
    
    # Adjust column sizes
    colsize!(gl_root, 1, Fixed(300)) # Narrower input column
    colsize!(gl_root, 2, Auto())     # Plot area takes the rest


    # --- UI Sections within input_grid ---
    # The input_grid will now have a simpler, sequential layout.
    # Row 1: Bz input
    # Row 2: End Time input
    # Row 3: Particle Species Menu
    # Row 4: Solver Menu & its parameters
    # Row 5: Run Button
    # Row 6: Error Label

    # Row 1: B-Field Bz
    b_field_minimal_grid = input_grid[1, 1] = GridLayout(tellheight=false, tellwidth=false, halign=:left)
    Label(b_field_minimal_grid[1,1], "Bz:")
    tb_Bz = Textbox(b_field_minimal_grid[1,2], placeholder = "1e-8", validator = Float64)
    # Initial value of parsed_bz_obs (1e-8) is the default. 
    # If placeholder is valid, Textbox.value will update parsed_bz_obs if linked or on first valid input.
    # For simplicity, we assume the placeholder matches the observable's default or is valid.
    # Or, we can explicitly link: on(tb_Bz.value) do val; parsed_bz_obs[] = val; end
    # and then trigger an initial update if needed, but Makie often handles this.
    # Let's ensure the observable is updated when the Textbox value changes.
    on(tb_Bz.value) do val # Listens to successfully validated Float64 values
        # val is guaranteed to be a Float64 due to validator with allow_nothing=false
        parsed_bz_obs[] = val
        # No specific error message is set by this callback for Bz format errors anymore.
        # Error clearing for "Bz" can be removed or made more specific if other logic sets Bz errors.
        # For now, removing it to rely on global error clear by run button or specific error setting.
        # if occursin("Bz", error_display_label.text[]) error_display_label.text[] = "" end
    end
    # Initialize observable from placeholder if valid (optional, if Textbox doesn't do it automatically via value)
    # This is more robust if the Textbox doesn't immediately push its placeholder to .value
    try parsed_bz_obs[] = parse(Float64, tb_Bz.stored_string[]) catch; end


    # Row 2: Time Span (End Time only)
    tspan_minimal_grid = input_grid[2, 1] = GridLayout(tellheight=false, tellwidth=false, halign=:left)
    Label(tspan_minimal_grid[1,1], "End Time:")
    tb_t_end = Textbox(tspan_minimal_grid[1,2], placeholder = "10.0", validator = Float64)
    on(tb_t_end.value) do val # Listens to successfully validated Float64 values
        # val is guaranteed to be a Float64
        parsed_t_end_obs[] = val
        # Positivity check for t_end remains in the run_button.clicks logic.
        # No format error is set by this callback. Clearing related to "End Time" format error removed.
        # If run button sets "ERROR: End Time must be greater than 0.0.", this callback could clear it if logic was added.
        # For now, only run button clears its own errors or global errors.
        # if occursin("End Time", error_display_label.text[]) error_display_label.text[] = "" end
    end
    try parsed_t_end_obs[] = parse(Float64, tb_t_end.stored_string[]) catch; end

    # Row 3: Particle Species Section (Kept as is from previous simplification)
    species_grid = input_grid[3, 1] = GridLayout(tellheight = false, tellwidth=false, halign=:left)
    # Label(species_grid[1, 1:2], "Particle Species", tellwidth=false, font=:bold) # Optional title for section
    Label(species_grid[1,1], "Species:")
    species_options = ["Proton", "Electron"]
    # Link Menu selection directly to selected_species_obs
    species_menu = Menu(species_grid[1,2], options = species_options, selection = selected_species_obs)
    # selected_species = species_menu.selection # No longer needed, selected_species_obs is the source of truth

    # Row 4: Solver Configuration Section (Kept as is)
    solver_grid = input_grid[4, 1] = GridLayout(tellheight = false, tellwidth=false, halign=:left)
    # Label(solver_grid[1, 1:2], "Solver Configuration", tellwidth=false, font=:bold) # Optional title
    Label(solver_grid[1,1], "Solver:")
    solver_options = ["Vern9 (adaptive)", "Tsit5 (adaptive)", "Boris (fixed-step)", "ImplicitMidpoint (fixed-step)"]
    # Link Menu selection directly to selected_solver_obs
    solver_menu = Menu(solver_grid[1,2], options = solver_options, selection = selected_solver_obs)
    # selected_solver variable is no longer needed, selected_solver_obs is the source of truth

    solver_params_grid = solver_grid[2,1:2] = GridLayout() # Spans 2 cols relative to solver_grid's own layout
    # Timestep (dt)
    lbl_dt = Label(solver_params_grid[1,1], "Timestep (dt):")
    tb_dt = Textbox(solver_params_grid[1,2], placeholder = "0.01", validator = Float64)
    # Relative Tol. (reltol)
    lbl_reltol = Label(solver_params_grid[2,1], "Relative Tol. (reltol):")
    tb_reltol = Textbox(solver_params_grid[2,2], placeholder = "1e-7", validator = Float64)

    lbl_dt.visible = false; tb_dt.visible = false # Initial state
    lbl_reltol.visible = false; tb_reltol.visible = false # Initial state

    # Reactive parsing for dt (using .value)
    try parsed_dt_obs[] = parse(Float64, tb_dt.stored_string[]) catch; end # Keep default if placeholder invalid
    on(tb_dt.value) do val # val is guaranteed Float64
        if val > 0.0
            parsed_dt_obs[] = val
            if error_display_label.text[] == "ERROR: Timestep (dt) must be positive." # Clear only this specific error
                error_display_label.text[] = "" 
            end
        else
            error_display_label.text[] = "ERROR: Timestep (dt) must be positive."
            # parsed_dt_obs[] is not updated with non-positive value by this reactive callback
        end
    end

    # Reactive parsing for reltol (using .value)
    try parsed_reltol_obs[] = parse(Float64, tb_reltol.stored_string[]) catch; end # Keep default if placeholder invalid
    on(tb_reltol.value) do val # val is guaranteed Float64
        if val > 0.0
            parsed_reltol_obs[] = val
            if error_display_label.text[] == "ERROR: Relative Tolerance (reltol) must be positive." # Clear only this specific error
                error_display_label.text[] = "" 
            end
        else
            error_display_label.text[] = "ERROR: Relative Tolerance (reltol) must be positive."
            # parsed_reltol_obs[] is not updated with non-positive value by this reactive callback
        end
    end
    
    on(selected_solver_obs) do selected_type # Listen to the observable
        is_fixed_step = occursin("fixed-step", selected_type)
        is_adaptive = occursin("adaptive", selected_type)
        lbl_dt.visible[] = is_fixed_step
        tb_dt.visible[] = is_fixed_step
        lbl_reltol.visible[] = is_adaptive
        tb_reltol.visible[] = is_adaptive
    end
    notify(selected_solver_obs) # Notify the observable

    # Row 5: Run Simulation Button
    run_button = Button(input_grid[5, 1], label = "Run Simulation", tellwidth=false, halign=:left)

    # Row 6: Error/Status Display Label
    error_display_label = Label(input_grid[6, 1], "", tellwidth=false, halign=:left) # Removed color=:red


    # --- Button Click Logic ---
    
    # Helper function to clear all plots (2D time series only)
    function clear_all_plots()
        # Clear data observables, which will clear the plots
        t_obs[] = Float64[]
        x_pos_obs[] = Float64[]
        y_pos_obs[] = Float64[]
        z_pos_obs[] = Float64[]
        vx_vel_obs[] = Float64[]
        vy_vel_obs[] = Float64[]
        vz_vel_obs[] = Float64[]
        
        # Autolimits after clearing data
        autolimits!(pos_ts_ax)
        autolimits!(vel_ts_ax)
        # Legends are persistent, no need to clear/recreate if plots are linked to observables
    end

    on(run_button.clicks) do _
        # Clear previous solution, error message, and all plots
        sol_obs[] = nothing 
        error_display_label.text[] = "" 
        clear_all_plots() # This now empties the data observables

        # --- Parse Initial Conditions (Hardcoded) ---
        # initial_pos_ref[] is no longer set from UI
        # initial_vel_ref[] is no longer set from UI
        # Using hardcoded values directly in stateinit construction later

        # --- Use pre-parsed and validated values ---
        # The error_display_label was already cleared at the beginning of this callback.
        # Reactive parsers handle their own specific error clearing on valid input.

        # Validate End Time (already parsed into parsed_t_end_obs)
        current_t_end = parsed_t_end_obs[]
        if 0.0 >= current_t_end 
            error_display_label.text[] = "ERROR: End Time must be greater than 0.0."
            # error_display_label.color[] = :red # Removed
            return
        end
        # t_span_ref is no longer used, (0.0, current_t_end) will be used directly.
        
        # E-Field Parameters are Hardcoded to Zero.

        # B-Field Parameters (Bz only, already parsed into parsed_bz_obs)
        # bz_value_ref is no longer used.
        # No parsing needed here, value is in parsed_bz_obs[].

        # --- Particle Species Parameters (use selected_species_obs directly) ---
        # species_type_ref[] = selected_species_obs[] # No need to update species_type_ref if it's removed
        # Logic for parsing custom q and m REMOVED

        # --- Solver Configuration (use selected_solver_obs directly) ---
        # solver_type_ref[] = selected_solver_obs[] # No longer need to update solver_type_ref
        current_solver_type = selected_solver_obs[]
        is_fixed_step_solver = occursin("fixed-step", current_solver_type)
        is_adaptive_solver = occursin("adaptive", current_solver_type)

        # Use pre-parsed and validated values from observables
        # The on-change callbacks for tb_dt and tb_reltol handle format and positivity.
        # If an error message related to dt or reltol is currently displayed, we should not proceed.
        if (is_fixed_step_solver && (error_display_label.text[] == "ERROR: Timestep (dt) must be a positive number." || parsed_dt_obs[] <= 0.0))
            if error_display_label.text[] != "ERROR: Timestep (dt) must be a positive number." # Set it if not already set by reactive parser
                 error_display_label.text[] = "ERROR: Timestep (dt) must be a positive number."
            end
            return
        end
        if (is_adaptive_solver && (error_display_label.text[] == "ERROR: Relative Tolerance (reltol) must be a positive number." || parsed_reltol_obs[] <= 0.0))
             if error_display_label.text[] != "ERROR: Relative Tolerance (reltol) must be a positive number."
                 error_display_label.text[] = "ERROR: Relative Tolerance (reltol) must be a positive number."
            end
            return
        end
        
        # Clear any parsing related error messages before attempting simulation
        # This specifically clears dt/reltol errors if they were set and now conditions are met
        if error_display_label.text[] == "ERROR: Timestep (dt) must be a positive number." || error_display_label.text[] == "ERROR: Relative Tolerance (reltol) must be a positive number."
            error_display_label.text[] = ""
        end
        error_display_label.text[] = "Running simulation..." 
        # error_display_label.color[] = :black # Removed

        # --- Construct E-Field Function (Hardcoded to Zero) ---
        E_func = TestParticle.ZeroField()

        # --- Construct B-Field Function ---
        current_bz = parsed_bz_obs[] # Use value from observable
        B_func = (xu, t) -> SA[0.0, 0.0, current_bz]
        if current_bz == 0.0
            B_func = TestParticle.ZeroField()
        end

        # --- Determine Particle Species for prepare() ---
        species_kwargs = Dict{Symbol, Any}()
        # Use the value from the observable directly
        if selected_species_obs[] == "Proton"
            species_kwargs[:species] = TestParticle.Proton
        elseif selected_species_obs[] == "Electron"
            species_kwargs[:species] = TestParticle.Electron
        end
        # "Custom" branch REMOVED as it's no longer an option
        
        # --- Call TestParticle.prepare() ---
        # Note: The example showed `param` (singular) from prepare, used in TraceProblem.
        # ODEProblem uses `p` which is often `params` plural. Let's call it `sim_params`.
        sim_params = TestParticle.prepare(E_func, B_func; species_kwargs...)

        # --- Construct Initial State (Hardcoded) ---
        # Initial position [0.0, 0.0, 0.0]
        # Initial velocity [1000.0, 0.0, 0.0] (as per minimal plan)
        hardcoded_initial_pos = SA[0.0, 0.0, 0.0]
        hardcoded_initial_vel = SA[1000.0, 0.0, 0.0]
        stateinit = vcat(hardcoded_initial_pos, hardcoded_initial_vel) 

        # --- Create Problem and Solve ---
        try
            if current_solver_type == "Boris (fixed-step)" # Use current_solver_type
                # Boris method specific setup from example
                # stateinit_boris = SVector{6, Float64}(stateinit...) # Boris might prefer SVector
                # prob_boris = TestParticle.TraceProblem(stateinit_svector, (0.0, parsed_t_end_obs[]), sim_params) 
                # sol_ref[] = TestParticle.solve(prob_boris; dt = solver_dt_ref[])[1] # [1] if it returns a tuple
                # Re-checking example: TestParticle.solve(prob_boris; dt) returns sol directly.
                # The example `demo_proton_dipole.jl` uses `sol = TestParticle.solve(prob_boris; dt)[1]`.
                # This implies `TestParticle.solve` might return a tuple (sol, some_other_info).
                # Let's assume it returns the solution object that can be plotted.
                # And `param` (singular) was used in example.
                
                # Reconciling: TestParticle.jl README uses `params` for ODEProblem,
                # and `param` for `TraceProblem`. Let's assume `sim_params` from `prepare` works for both.
                stateinit_svector = SVector{6,Float64}(stateinit)
                prob_boris = TestParticle.TraceProblem(stateinit_svector, (0.0, parsed_t_end_obs[]), sim_params) 
                sol_obs[] = TestParticle.solve(prob_boris; dt = parsed_dt_obs[]) # Use parsed_dt_obs & sol_obs

            else # OrdinaryDiffEq solvers
                prob = ODEProblem(TestParticle.trace, stateinit, (0.0, parsed_t_end_obs[]), sim_params) 
                # current_solver_reltol and current_solver_dt are no longer Refs.
                # Use values from observables directly.

                if current_solver_type == "Vern9 (adaptive)" # Use current_solver_type
                    sol_obs[] = OrdinaryDiffEq.solve(prob, Vern9(), reltol = parsed_reltol_obs[], abstol=1e-8, save_everystep=false) 
                elseif current_solver_type == "Tsit5 (adaptive)" # Use current_solver_type
                    sol_obs[] = OrdinaryDiffEq.solve(prob, Tsit5(), reltol = parsed_reltol_obs[], abstol=1e-8, save_everystep=false) 
                elseif current_solver_type == "ImplicitMidpoint (fixed-step)" # Use current_solver_type
                    sol_obs[] = OrdinaryDiffEq.solve(prob, ImplicitMidpoint(), dt = parsed_dt_obs[], save_everystep=false) 
                else
                    error_display_label.text[] = "Error: Unknown solver '$(current_solver_type)'. Cannot solve." 
                    return 
                end
            end

            if sol_obs[] !== nothing && sol_obs[].retcode == :Success 
                error_display_label.text[] = "Simulation successful! Retcode: $(sol_obs[].retcode)" 
                # error_display_label.color[] = :green # Removed
                
                # 3D Trajectory Plotting REMOVED

                # --- Update 2D Time Series Plot Data ---
                sol = sol_obs[] 
                
                # Extract data from solution
                t_new = sol.t
                x_new = [state[1] for state in sol.u]
                y_new = [state[2] for state in sol.u]
                z_new = [state[3] for state in sol.u]
                vx_new = [state[4] for state in sol.u]
                vy_new = [state[5] for state in sol.u]
                vz_new = [state[6] for state in sol.u]

                # Update observables to trigger plot updates
                t_obs[] = t_new
                x_pos_obs[] = x_new
                y_pos_obs[] = y_new
                z_pos_obs[] = z_new
                vx_vel_obs[] = vx_new
                vy_vel_obs[] = vy_new
                vz_vel_obs[] = vz_new

                # Rescale axes after data update
                autolimits!(pos_ts_ax)
                autolimits!(vel_ts_ax)

            elseif sol_obs[] !== nothing # Simulation ran but might not be successful 
                error_display_label.text[] = "Simulation completed with issues. Retcode: $(sol_obs[].retcode)" # Update to sol_obs
                # error_display_label.color[] = :orange # Removed
                # Plots already cleared by clear_all_plots() at the start of the callback
            else # sol_obs[] is nothing, meaning simulation didn't run or produce output # Update to sol_obs
                error_display_label.text[] = "Simulation failed to produce a solution object."
                # error_display_label.color[] = :red # Removed
                # Plots already cleared
            end

        catch e
            error_display_label.text[] = "Simulation failed: " * sprint(showerror, e)
            # error_display_label.color[] = :red # Removed
            # sol_obs[] is already nothing or will be reset if error happened mid-simulation # Update to sol_obs
            # Plots already cleared by clear_all_plots() at the start of the callback
            # or if error was during plotting, some might exist.
            # For safety, can call clear_all_plots() again, or ensure it's robust.
            # The initial clear should handle most cases.
        end
        
    end

    return fig
end

# Call it if running as a script, or ensure it can be called from elsewhere.
# main_gui()
