using GLMakie
using TestParticle
using OrdinaryDiffEq
using StaticArrays

function main_gui()
    # --- Data Storage for Parsed Values ---
    # Field Configuration
    e_field_type_ref = Ref("Zero") # String type
    e_uniform_params_ref = Ref(SA[0.0, 0.0, 0.0]) # Ex, Ey, Ez

    b_field_type_ref = Ref("Zero") # String type
    b_uniform_params_ref = Ref(SA[0.0, 0.0, 0.0]) # Bx, By, Bz
    b_dipole_params_ref = Ref(SA[0.0, 0.0, 1.0])  # Mx, My, Mz (default Mz=1.0)
    b_cs_params_ref = Ref(SA[1.0, 0.1, 0.0])      # B₀, L, Bn (defaults)

    # Initial Conditions
    initial_pos_ref = Ref(SA[0.0, 0.0, 0.0]) # x₀, y₀, z₀
    initial_vel_ref = Ref(SA[0.0, 0.0, 0.0]) # vx₀, vy₀, vz₀

    # Time Span
    t_span_ref = Ref((0.0, 10.0)) # t_start, t_end

    # Particle Species
    species_type_ref = Ref("Proton") # String type
    custom_species_q_ref = Ref(1.0)
    custom_species_m_ref = Ref(1.0)

    # Solver Configuration
    solver_type_ref = Ref("Vern9 (adaptive)") # String type
    solver_dt_ref = Ref(0.01)
    solver_reltol_ref = Ref(1e-6)

    # Simulation Solution Storage
    sol_ref = Ref{Any}(nothing)
    # Current Trajectory Plot Object Storage
    current_trajectory_plot = Ref{Any}(nothing)
    # Refs for 2D time-series plot lines (stores plot objects like LineSegments)
    pos_ts_plots_ref = Ref{Vector{Any}}(Any[]) # For x(t), y(t), z(t) lines
    vel_ts_plots_ref = Ref{Vector{Any}}(Any[]) # For vx(t), vy(t), vz(t) lines

    # --- Figure and Main Layout ---
    fig = Figure(size = (1600, 900)) # Increased resolution for more plots
    # Main GridLayout: 1 row, 3 columns (inputs | 2D plots | 3D plot)
    gl_root = fig[1, 1] = GridLayout() 

    # Column 1: Input Grid
    input_grid = gl_root[1, 1] = GridLayout() 
    
    # Column 2: Time Series Plots (stacked vertically)
    middle_plot_grid = gl_root[1, 2] = GridLayout()
    pos_ts_ax = Axis(middle_plot_grid[1, 1], title="Position vs Time", xlabel="Time [s]", ylabel="Position [m]")
    vel_ts_ax = Axis(middle_plot_grid[2, 1], title="Velocity vs Time", xlabel="Time [s]", ylabel="Velocity [m/s]")
    # Link x-axes of time series plots for synchronized zooming/panning if desired
    linkxaxes!(pos_ts_ax, vel_ts_ax)


    # Column 3: 3D Trajectory Plot (existing plot_ax, now moved)
    plot_ax = Axis3(gl_root[1, 3], title="Particle Trajectory", aspect=:data)
    plot_ax.xlabel = "X [m]" # Keep existing labels
    plot_ax.ylabel = "Y [m]"
    plot_ax.zlabel = "Z [m]"
    
    # Adjust column sizes
    colsize!(gl_root, 1, Fixed(400))      # Input column
    colsize!(gl_root, 2, Relative(0.35))  # Middle column for 2D plots (takes 35% of remaining space after fixed)
    colsize!(gl_root, 3, Relative(0.65))  # Right column for 3D plot (takes 65% of remaining space)
    # Alternative for Auto sizing:
    # colsize!(gl_root, 1, Fixed(400))
    # colsize!(gl_root, 2, Auto())
    # colsize!(gl_root, 3, Auto())


    # --- UI Sections within input_grid ---

    # Section 0: Error/Status Display Label (added at the top of input_grid for visibility)
    # Spanning 2 columns to align with parameter grids if they use 2 columns.
    # The subsequent grids (field_cfg_grid etc.) are added from row 1 downwards.
    # Let's make this input_grid[1,1] and shift others down, or add it at the end.
    # For prominence, let's try to place it conceptually "above" the first config block.
    # This means other sections will start from input_grid[2,1] etc.
    # Let's add it at the very end of the input_grid (e.g., row 7) for simplicity first.
    # The run_button is at input_grid[6,1]. So error label at input_grid[7,1].

    # Field Configuration section (now added to input_grid)
    # Row 1 of input_grid
    field_cfg_grid = input_grid[1, 1] = GridLayout(tellheight = false)
    Label(field_cfg_grid[1, 1:2], "Field Configuration", tellwidth=false, font=:bold)

    # E-Field Section
    Label(field_cfg_grid[2, 1], "E-Field Type:")
    e_field_options = ["Zero", "Uniform"]
    e_menu = Menu(field_cfg_grid[2, 2], options = e_field_options, default = "Zero")
    selected_e_field = e_menu.selection

    e_params_grid = field_cfg_grid[3, 1:2] = GridLayout()

    # Uniform E-Field Parameters
    # Define and parent e_uniform_params_grid directly
    e_uniform_params_grid = e_params_grid[1, 1] = GridLayout()
    Label(e_uniform_params_grid[1,1], "Ex:")
    tb_Ex = Textbox(e_uniform_params_grid[1,2], placeholder = "0.0")
    Label(e_uniform_params_grid[2,1], "Ey:")
    tb_Ey = Textbox(e_uniform_params_grid[2,2], placeholder = "0.0")
    Label(e_uniform_params_grid[3,1], "Ez:")
    tb_Ez = Textbox(e_uniform_params_grid[3,2], placeholder = "0.0")
    e_uniform_params_grid.layoutobservables.visible[] = false # Initially hidden
    # The line `e_params_grid[1,1] = e_uniform_params_grid` is now redundant due to direct parenting.

    # B-Field Section
    Label(field_cfg_grid[4, 1], "B-Field Type:")
    b_field_options = ["Zero", "Uniform", "Dipole", "Current Sheet"]
    b_menu = Menu(field_cfg_grid[4, 2], options = b_field_options, default = "Zero")
    selected_b_field = b_menu.selection

    b_params_grid = field_cfg_grid[5, 1:2] = GridLayout()

    # Uniform B-Field Parameters
    b_uniform_params_grid = GridLayout()
    Label(b_uniform_params_grid[1,1], "Bx:")
    tb_Bx = Textbox(b_uniform_params_grid[1,2], placeholder = "0.0")
    Label(b_uniform_params_grid[2,1], "By:")
    tb_By = Textbox(b_uniform_params_grid[2,2], placeholder = "0.0")
    Label(b_uniform_params_grid[3,1], "Bz:")
    tb_Bz = Textbox(b_uniform_params_grid[3,2], placeholder = "0.0")
    b_uniform_params_grid.layoutobservables.visible[] = false # Initially hidden
    
    # Dipole B-Field Parameters
    b_dipole_params_grid = GridLayout()
    Label(b_dipole_params_grid[1,1], "Dipole Mx:")
    tb_Mx = Textbox(b_dipole_params_grid[1,2], placeholder = "0.0")
    Label(b_dipole_params_grid[2,1], "Dipole My:")
    tb_My = Textbox(b_dipole_params_grid[2,2], placeholder = "0.0")
    Label(b_dipole_params_grid[3,1], "Dipole Mz:")
    tb_Mz = Textbox(b_dipole_params_grid[3,2], placeholder = "1.0") # Default Mz
    b_dipole_params_grid.layoutobservables.visible[] = false # Initially hidden
    # Add to b_params_grid, but ensure only one is visible at a time.
    # We can place it in the same cell as b_uniform_params_grid, visibility handles the rest.

    # Current Sheet B-Field Parameters
    b_cs_params_grid = GridLayout()
    Label(b_cs_params_grid[1,1], "CS B₀:")
    tb_CS_B0 = Textbox(b_cs_params_grid[1,2], placeholder = "1.0")
    Label(b_cs_params_grid[2,1], "CS L:")
    tb_CS_L = Textbox(b_cs_params_grid[2,2], placeholder = "0.1")
    Label(b_cs_params_grid[3,1], "CS Bn:")
    tb_CS_Bn = Textbox(b_cs_params_grid[3,2], placeholder = "0.0")
    b_cs_params_grid.layoutobservables.visible[] = false # Initially hidden
    # Add to b_params_grid, same cell.

    # Dynamic E-Field Parameter Display Logic
    on(selected_e_field) do selected_type
        e_uniform_params_grid.layoutobservables.visible[] = (selected_type == "Uniform")
    end
    notify(selected_e_field) # Trigger initial update

    # Dynamic B-Field Parameter Display Logic
    on(selected_b_field) do selected_type
        b_uniform_params_grid.layoutobservables.visible[] = (selected_type == "Uniform")
        b_dipole_params_grid.layoutobservables.visible[] = (selected_type == "Dipole")
        b_cs_params_grid.layoutobservables.visible[] = (selected_type == "Current Sheet")

        # Ensure correct grid is in the parent if they share the same cell
        if selected_type == "Uniform"
            b_params_grid[1,1] = b_uniform_params_grid
        elseif selected_type == "Dipole"
            b_params_grid[1,1] = b_dipole_params_grid
        elseif selected_type == "Current Sheet"
            b_params_grid[1,1] = b_cs_params_grid
        else # Zero or other cases
            # Hide all, or clear the cell if necessary
            # For now, setting visible to false is handled above,
            # but explicitly clearing the cell might be safer if grids overlap
            delete!(b_params_grid, 1, 1)
        end
    end
    notify(selected_b_field) # Trigger initial update

    # Initial Conditions Section
    init_cond_grid = input_grid[2, 1] = GridLayout(tellheight = false)
    Label(init_cond_grid[1, 1:2], "Initial Conditions", tellwidth=false, font=:bold)

    Label(init_cond_grid[2,1], "x₀:")
    tb_x0 = Textbox(init_cond_grid[2,2], placeholder = "0.0")
    Label(init_cond_grid[3,1], "y₀:")
    tb_y0 = Textbox(init_cond_grid[3,2], placeholder = "0.0")
    Label(init_cond_grid[4,1], "z₀:")
    tb_z0 = Textbox(init_cond_grid[4,2], placeholder = "0.0")

    Label(init_cond_grid[5,1], "vx₀:")
    tb_vx0 = Textbox(init_cond_grid[5,2], placeholder = "0.0")
    Label(init_cond_grid[6,1], "vy₀:")
    tb_vy0 = Textbox(init_cond_grid[6,2], placeholder = "0.0")
    Label(init_cond_grid[7,1], "vz₀:")
    tb_vz0 = Textbox(init_cond_grid[7,2], placeholder = "0.0")

    # Time Span Section
    tspan_grid = input_grid[3, 1] = GridLayout(tellheight = false)
    Label(tspan_grid[1, 1:2], "Time Span", tellwidth=false, font=:bold)

    Label(tspan_grid[2,1], "Start Time:")
    tb_t_start = Textbox(tspan_grid[2,2], placeholder = "0.0")
    Label(tspan_grid[3,1], "End Time:")
    tb_t_end = Textbox(tspan_grid[3,2], placeholder = "10.0")

    # Particle Species Section
    species_grid = input_grid[4, 1] = GridLayout(tellheight = false)
    Label(species_grid[1, 1:2], "Particle Species", tellwidth=false, font=:bold)

    Label(species_grid[2,1], "Species:")
    species_options = ["Proton", "Electron", "Custom"]
    species_menu = Menu(species_grid[2,2], options = species_options, default = "Proton")
    selected_species = species_menu.selection

    custom_species_params_grid = species_grid[3,1:2] = GridLayout()
    Label(custom_species_params_grid[1,1], "Charge (q):")
    tb_q = Textbox(custom_species_params_grid[1,2], placeholder = "1.0")
    Label(custom_species_params_grid[2,1], "Mass (m):")
    tb_m = Textbox(custom_species_params_grid[2,2], placeholder = "1.0")
    custom_species_params_grid.layoutobservables.visible[] = false # Initially hidden

    on(selected_species) do selected_type
        custom_species_params_grid.layoutobservables.visible[] = (selected_type == "Custom")
    end
    notify(selected_species)

    # Solver Configuration Section
    solver_grid = input_grid[5, 1] = GridLayout(tellheight = false)
    Label(solver_grid[1, 1:2], "Solver Configuration", tellwidth=false, font=:bold)

    Label(solver_grid[2,1], "Solver:")
    solver_options = ["Vern9 (adaptive)", "Tsit5 (adaptive)", "Boris (fixed-step)", "ImplicitMidpoint (fixed-step)"]
    solver_menu = Menu(solver_grid[2,2], options = solver_options, default = "Vern9 (adaptive)")
    selected_solver = solver_menu.selection

    solver_params_grid = solver_grid[3,1:2] = GridLayout()
    
    # Timestep (dt) - for fixed-step solvers
    row_dt = 1 # Current row in solver_params_grid
    lbl_dt = Label(solver_params_grid[row_dt,1], "Timestep (dt):")
    tb_dt = Textbox(solver_params_grid[row_dt,2], placeholder = "0.01")
    
    # Relative Tol. (reltol) - for adaptive solvers
    row_reltol = row_dt + 1 # Place it in the next row or manage layout better
    lbl_reltol = Label(solver_params_grid[row_reltol,1], "Relative Tol. (reltol):")
    tb_reltol = Textbox(solver_params_grid[row_reltol,2], placeholder = "1e-6")

    # Initially hide all solver parameters, then show based on default selection
    lbl_dt.visible = false
    tb_dt.visible = false
    lbl_reltol.visible = false
    tb_reltol.visible = false

    on(selected_solver) do selected_type
        is_fixed_step = occursin("fixed-step", selected_type)
        is_adaptive = occursin("adaptive", selected_type)

        lbl_dt.visible[] = is_fixed_step
        tb_dt.visible[] = is_fixed_step
        
        lbl_reltol.visible[] = is_adaptive
        tb_reltol.visible[] = is_adaptive
    end
    notify(selected_solver) # Trigger initial update

    # Run Simulation Button (Row 6 of input_grid)
    run_button = Button(input_grid[6, 1], label = "Run Simulation", tellwidth=false)

    # Error/Status Display Label (Row 7 of input_grid)
    error_display_label = Label(input_grid[7, 1], "", color=:red, tellwidth=false, justification=:left)
    # Make it span across if input_grid's main content is effectively 2 columns wide
    # The individual sections like field_cfg_grid span 1:2 internally, but are placed in input_grid[X,1]
    # So, error_display_label at input_grid[7,1] is fine.

    # --- Button Click Logic ---
    
    # Helper function to clear all plots (3D trajectory and 2D time series)
    function clear_all_plots()
        # Clear 3D trajectory plot
        if !isnothing(current_trajectory_plot[])
            delete!(plot_ax, current_trajectory_plot[])
            current_trajectory_plot[] = nothing
        end
        autolimits!(plot_ax) # Reset limits for 3D plot

        for plt in pos_ts_plots_ref[] delete!(pos_ts_ax, plt) end
        empty!(pos_ts_plots_ref[])
        autolimits!(pos_ts_ax)

        for plt in vel_ts_plots_ref[] delete!(vel_ts_ax, plt) end
        empty!(vel_ts_plots_ref[])
        autolimits!(vel_ts_ax)
        
        # Clear legends if they exist (Makie might handle this with plot deletion, but to be safe)
        # This might require storing legend objects if they need explicit deletion.
        # For now, assume deleting lines clears relevant legend entries or legend is auto-updated/rebuilt.
        # If explicit legend object deletion were needed:
        #   legend_pos = findfirst(x -> x isa Legend, pos_ts_ax.scene.plots)
        #   if legend_pos !== nothing delete!(pos_ts_ax.scene, pos_ts_ax.scene.plots[legend_pos]) end
        # Makie's `axislegend` is generally robust to being called again.
    end

    on(run_button.clicks) do _
        # Clear previous solution, error message, and all plots
        sol_ref[] = nothing
        error_display_label.text[] = "" 
        clear_all_plots()

        # --- Parse Initial Conditions ---
        val_x0 = tryparse(Float64, tb_x0.stored_string[])
        val_y0 = tryparse(Float64, tb_y0.stored_string[])
        val_z0 = tryparse(Float64, tb_z0.stored_string[])

        if isnothing(val_x0) error_display_label.text[] = "Error: Invalid x₀. Must be a number."; return; end
        if isnothing(val_y0) error_display_label.text[] = "Error: Invalid y₀. Must be a number."; return; end
        if isnothing(val_z0) error_display_label.text[] = "Error: Invalid z₀. Must be a number."; return; end
        initial_pos_ref[] = SA[val_x0, val_y0, val_z0]

        val_vx0 = tryparse(Float64, tb_vx0.stored_string[])
        val_vy0 = tryparse(Float64, tb_vy0.stored_string[])
        val_vz0 = tryparse(Float64, tb_vz0.stored_string[])

        if isnothing(val_vx0) error_display_label.text[] = "Error: Invalid vx₀. Must be a number."; return; end
        if isnothing(val_vy0) error_display_label.text[] = "Error: Invalid vy₀. Must be a number."; return; end
        if isnothing(val_vz0) error_display_label.text[] = "Error: Invalid vz₀. Must be a number."; return; end
        initial_vel_ref[] = SA[val_vx0, val_vy0, val_vz0]

        # --- Parse Time Span ---
        val_t_start = tryparse(Float64, tb_t_start.stored_string[])
        val_t_end = tryparse(Float64, tb_t_end.stored_string[])

        if isnothing(val_t_start) error_display_label.text[] = "Error: Invalid Start Time. Must be a number."; return; end
        if isnothing(val_t_end) error_display_label.text[] = "Error: Invalid End Time. Must be a number."; return; end
        if val_t_start >= val_t_end error_display_label.text[] = "Error: Start Time must be less than End Time."; return; end
        t_span_ref[] = (val_t_start, val_t_end)
        
        # --- Parse E-Field Parameters ---
        e_field_type_ref[] = selected_e_field[]
        if e_field_type_ref[] == "Uniform"
            val_Ex = tryparse(Float64, tb_Ex.stored_string[]); if isnothing(val_Ex) error_display_label.text[] = "Error: Invalid Ex."; return; end
            val_Ey = tryparse(Float64, tb_Ey.stored_string[]); if isnothing(val_Ey) error_display_label.text[] = "Error: Invalid Ey."; return; end
            val_Ez = tryparse(Float64, tb_Ez.stored_string[]); if isnothing(val_Ez) error_display_label.text[] = "Error: Invalid Ez."; return; end
            e_uniform_params_ref[] = SA[val_Ex, val_Ey, val_Ez]
        else
            e_uniform_params_ref[] = SA[0.0, 0.0, 0.0]
        end

        # --- Parse B-Field Parameters ---
        b_field_type_ref[] = selected_b_field[]
        if b_field_type_ref[] == "Uniform"
            val_Bx = tryparse(Float64, tb_Bx.stored_string[]); if isnothing(val_Bx) error_display_label.text[] = "Error: Invalid Bx."; return; end
            val_By = tryparse(Float64, tb_By.stored_string[]); if isnothing(val_By) error_display_label.text[] = "Error: Invalid By."; return; end
            val_Bz = tryparse(Float64, tb_Bz.stored_string[]); if isnothing(val_Bz) error_display_label.text[] = "Error: Invalid Bz."; return; end
            b_uniform_params_ref[] = SA[val_Bx, val_By, val_Bz]
        elseif b_field_type_ref[] == "Dipole"
            val_Mx = tryparse(Float64, tb_Mx.stored_string[]); if isnothing(val_Mx) error_display_label.text[] = "Error: Invalid Dipole Mx."; return; end
            val_My = tryparse(Float64, tb_My.stored_string[]); if isnothing(val_My) error_display_label.text[] = "Error: Invalid Dipole My."; return; end
            val_Mz = tryparse(Float64, tb_Mz.stored_string[]); if isnothing(val_Mz) error_display_label.text[] = "Error: Invalid Dipole Mz."; return; end
            b_dipole_params_ref[] = SA[val_Mx, val_My, val_Mz]
        elseif b_field_type_ref[] == "Current Sheet"
            val_CS_B0 = tryparse(Float64, tb_CS_B0.stored_string[]); if isnothing(val_CS_B0) error_display_label.text[] = "Error: Invalid CS B₀."; return; end
            val_CS_L  = tryparse(Float64, tb_CS_L.stored_string[]);  if isnothing(val_CS_L)  error_display_label.text[] = "Error: Invalid CS L."; return; end
            val_CS_Bn = tryparse(Float64, tb_CS_Bn.stored_string[]); if isnothing(val_CS_Bn) error_display_label.text[] = "Error: Invalid CS Bn."; return; end
            b_cs_params_ref[] = SA[val_CS_B0, val_CS_L, val_CS_Bn]
        else
            b_uniform_params_ref[] = SA[0.0, 0.0, 0.0]
        end

        # --- Parse Particle Species Parameters ---
        species_type_ref[] = selected_species[]
        if species_type_ref[] == "Custom"
            val_q = tryparse(Float64, tb_q.stored_string[]); if isnothing(val_q) error_display_label.text[] = "Error: Invalid Charge (q)."; return; end
            val_m = tryparse(Float64, tb_m.stored_string[]); if isnothing(val_m) error_display_label.text[] = "Error: Invalid Mass (m)."; return; end
            if val_m <= 0.0 error_display_label.text[] = "Error: Mass (m) must be positive."; return; end
            custom_species_q_ref[] = val_q
            custom_species_m_ref[] = val_m
        else
            custom_species_q_ref[] = 1.0
            custom_species_m_ref[] = 1.0
        end

        # --- Parse Solver Configuration ---
        solver_type_ref[] = selected_solver[]
        is_fixed_step_solver = occursin("fixed-step", solver_type_ref[])
        is_adaptive_solver = occursin("adaptive", solver_type_ref[])

        if is_fixed_step_solver
            val_dt = tryparse(Float64, tb_dt.stored_string[]); if isnothing(val_dt) error_display_label.text[] = "Error: Invalid Timestep (dt)."; return; end
            if val_dt <= 0.0 error_display_label.text[] = "Error: Timestep (dt) must be positive."; return; end
            solver_dt_ref[] = val_dt
            solver_reltol_ref[] = NaN
        elseif is_adaptive_solver
            val_reltol = tryparse(Float64, tb_reltol.stored_string[]); if isnothing(val_reltol) error_display_label.text[] = "Error: Invalid Relative Tol."; return; end
            if val_reltol <= 0.0 error_display_label.text[] = "Error: Relative Tol. must be positive."; return; end
            solver_reltol_ref[] = val_reltol
            solver_dt_ref[] = NaN
        else
            error_display_label.text[] = "Warning: Unknown solver type. Using default dt/reltol." # Changed from println
            solver_dt_ref[] = NaN 
            solver_reltol_ref[] = NaN
        end
        
        # Clear any parsing related error messages before attempting simulation
        error_display_label.text[] = "Running simulation..." 
        error_display_label.color[] = :black # Neutral color for status

        # --- Construct E-Field Function ---
        E_func = TestParticle.ZeroField() # Default to zero field
        if e_field_type_ref[] == "Uniform"
            # Capture current parameter values for the closure
            current_e_params = e_uniform_params_ref[]
            E_func = (xu, t) -> SVector{3, Float64}(current_e_params...)
            # Or if TestParticle.jl expects E_func(xu)
            # E_func = xu -> SVector{3, Float64}(current_e_params...)
        end

        # --- Construct B-Field Function ---
        B_func = TestParticle.ZeroField() # Default to zero field
        if b_field_type_ref[] == "Uniform"
            current_b_uniform_params = b_uniform_params_ref[]
            B_func = (xu, t) -> SVector{3, Float64}(current_b_uniform_params...)
        elseif b_field_type_ref[] == "Dipole"
            # Assuming TestParticle.dipole(r, M) is the correct function as per worker note
            current_b_dipole_params = b_dipole_params_ref[]
            # The dipole function in TestParticle.jl might not expect time `t`.
            # If it's (xu, M), then: B_func = (xu, t) -> TestParticle.dipole(xu[1:3], current_b_dipole_params)
            # Or if it's just (xu) and M is configured elsewhere (less likely given GUI)
            # For now, assuming it can take (xu, t) and ignore t, or (xu) and dipole handles it.
            # Let's use (xu,t) -> func(xu[1:3], params)
            B_func = (xu, t) -> TestParticle.dipole(xu[1:3], current_b_dipole_params)
        elseif b_field_type_ref[] == "Current Sheet"
            current_b_cs_params = b_cs_params_ref[]
            # TestParticle.getB_CS_harris(r, B₀, L, Bn)
            # B_func = (xu, t) -> TestParticle.getB_CS_harris(xu[1:3], current_b_cs_params[1], current_b_cs_params[2], current_b_cs_params[3])
            # Simpler:
            B_func = (xu, t) -> TestParticle.getB_CS_harris(xu[1:3], current_b_cs_params...)
        end

        # --- Determine Particle Species for prepare() ---
        species_kwargs = Dict{Symbol, Any}()
        if species_type_ref[] == "Proton"
            species_kwargs[:species] = TestParticle.Proton
        elseif species_type_ref[] == "Electron"
            species_kwargs[:species] = TestParticle.Electron
        elseif species_type_ref[] == "Custom"
            # Assuming TestParticle.prepare can handle q and m directly,
            # potentially with a generic species type if needed.
            # The doc for prepare states `species::Species = Proton, q = 1.0, m = 1.0`
            # So we provide q and m, and potentially TestParticle.User if it exists and is required.
            # For now, let's assume omitting `species` or using a placeholder + q/m works.
            # Based on `prepare` signature, we should be able to just pass q and m.
            # Let's try passing TestParticle.Ion (assuming it's a general type for this)
            # or if TestParticle.User is a thing. If not, just q and m.
            # The example uses Proton/Electron explicitly.
            # The safest is to provide q and m, and let species default if User is not a type.
            # Let's assume TestParticle.Species is an enum and Proton/Electron are values.
            # If custom q/m is provided, the `species` parameter might be ignored or used as a base.
            # For safety, let's use Proton as a base species and override q, m.
            species_kwargs[:species] = TestParticle.Proton # Base type
            species_kwargs[:q] = custom_species_q_ref[]
            species_kwargs[:m] = custom_species_m_ref[]
        end
        
        # --- Call TestParticle.prepare() ---
        # Note: The example showed `param` (singular) from prepare, used in TraceProblem.
        # ODEProblem uses `p` which is often `params` plural. Let's call it `sim_params`.
        sim_params = TestParticle.prepare(E_func, B_func; species_kwargs...)

        # --- Construct Initial State ---
        # stateinit must be a SVector for Boris, and regular Vector for ODEProblem
        # For now, let's use a regular vector, adaptable later if Boris needs SVector specifically here.
        stateinit = vcat(initial_pos_ref[], initial_vel_ref[]) # Regular Vector

        # --- Create Problem and Solve ---
        try
            if solver_type_ref[] == "Boris (fixed-step)"
                # Boris method specific setup from example
                # stateinit_boris = SVector{6, Float64}(stateinit...) # Boris might prefer SVector
                # prob_boris = TestParticle.TraceProblem(stateinit_boris, t_span_ref[], sim_params)
                # sol_ref[] = TestParticle.solve(prob_boris; dt = solver_dt_ref[])[1] # [1] if it returns a tuple
                # Re-checking example: TestParticle.solve(prob_boris; dt) returns sol directly.
                # The example `demo_proton_dipole.jl` uses `sol = TestParticle.solve(prob_boris; dt)[1]`.
                # This implies `TestParticle.solve` might return a tuple (sol, some_other_info).
                # Let's assume it returns the solution object that can be plotted.
                # And `param` (singular) was used in example.
                
                # Reconciling: TestParticle.jl README uses `params` for ODEProblem,
                # and `param` for `TraceProblem`. Let's assume `sim_params` from `prepare` works for both.
                stateinit_svector = SVector{6,Float64}(stateinit)
                prob_boris = TestParticle.TraceProblem(stateinit_svector, t_span_ref[], sim_params)
                sol_ref[] = TestParticle.solve(prob_boris; dt = solver_dt_ref[])

            else # OrdinaryDiffEq solvers
                prob = ODEProblem(TestParticle.trace!, stateinit, t_span_ref[], sim_params)
                current_solver_reltol = solver_reltol_ref[]
                current_solver_dt = solver_dt_ref[]

                if solver_type_ref[] == "Vern9 (adaptive)"
                    sol_ref[] = OrdinaryDiffEq.solve(prob, Vern9(), reltol = current_solver_reltol, abstol=1e-8, save_everystep=false) #abstol and save_everystep are common additions
                elseif solver_type_ref[] == "Tsit5 (adaptive)"
                    sol_ref[] = OrdinaryDiffEq.solve(prob, Tsit5(), reltol = current_solver_reltol, abstol=1e-8, save_everystep=false)
                elseif solver_type_ref[] == "ImplicitMidpoint (fixed-step)"
                    # ImplicitMidpoint requires dt. Check if reltol is also used or ignored.
                    # Typically for fixed step, reltol is not the primary control.
                    sol_ref[] = OrdinaryDiffEq.solve(prob, ImplicitMidpoint(), dt = current_solver_dt, save_everystep=false)
                else
                    println("Error: Unknown solver '", solver_type_ref[], "'. Cannot solve.")
                    return # Exit if solver is not recognized
                end
            end

            if sol_ref[] !== nothing && sol_ref[].retcode == :Success
                error_display_label.text[] = "Simulation successful! Retcode: $(sol_ref[].retcode)"
                error_display_label.color[] = :green
                
                # 3D Trajectory Plot (already being cleared by clear_all_plots)
                current_trajectory_plot[] = lines!(plot_ax, sol_ref[], idxs=(1,2,3), color=:blue, linewidth=2)
                autolimits!(plot_ax)

                # --- 2D Time Series Plots ---
                # Extract solution data
                sol = sol_ref[]
                t_vals = sol.t # Time points
                s_vals = sol.u # Array of state vectors [x,y,z,vx,vy,vz] at each time point

                # Extract individual position and velocity components over time
                x_pos = [state[1] for state in s_vals]
                y_pos = [state[2] for state in s_vals]
                z_pos = [state[3] for state in s_vals]
                
                vx_vel = [state[4] for state in s_vals]
                vy_vel = [state[5] for state in s_vals]
                vz_vel = [state[6] for state in s_vals]

                # Position vs Time plots
                push!(pos_ts_plots_ref[], lines!(pos_ts_ax, t_vals, x_pos, label="x(t)"))
                push!(pos_ts_plots_ref[], lines!(pos_ts_ax, t_vals, y_pos, label="y(t)"))
                push!(pos_ts_plots_ref[], lines!(pos_ts_ax, t_vals, z_pos, label="z(t)"))
                if !isempty(pos_ts_plots_ref[]) axislegend(pos_ts_ax, position=:rt); end
                autolimits!(pos_ts_ax)

                # Velocity vs Time plots
                push!(vel_ts_plots_ref[], lines!(vel_ts_ax, t_vals, vx_vel, label="vx(t)"))
                push!(vel_ts_plots_ref[], lines!(vel_ts_ax, t_vals, vy_vel, label="vy(t)"))
                push!(vel_ts_plots_ref[], lines!(vel_ts_ax, t_vals, vz_vel, label="vz(t)"))
                if !isempty(vel_ts_plots_ref[]) axislegend(vel_ts_ax, position=:rt); end
                autolimits!(vel_ts_ax)

            elseif sol_ref[] !== nothing # Simulation ran but might not be successful
                error_display_label.text[] = "Simulation completed with issues. Retcode: $(sol_ref[].retcode)"
                error_display_label.color[] = :orange
                # Plots already cleared by clear_all_plots() at the start of the callback
            else # sol_ref[] is nothing, meaning simulation didn't run or produce output
                error_display_label.text[] = "Simulation failed to produce a solution object."
                error_display_label.color[] = :red
                # Plots already cleared
            end

        catch e
            error_display_label.text[] = "Simulation failed: " * sprint(showerror, e)
            error_display_label.color[] = :red
            # sol_ref[] is already nothing or will be reset if error happened mid-simulation
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
