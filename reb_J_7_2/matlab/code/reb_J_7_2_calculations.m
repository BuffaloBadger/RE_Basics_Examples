function reb_J_7_2_calculations
%REB_J_7_2_CALCULATIONS Reaction Engineering Basics Example J.7.2

    % results file specification
    results_file_spec = '../results/reb_J_7_2_results.csv';

    % given and known constants
    T_f = 325.; % K
    D = 5.; % cm
    dH = -14000.; % cal /mol
    Cp = 1.3; % cal /cm^3 /K
    nDot_A_in = 1.0; % mol /min
    nDot_Z_in = 0.0; % mol /min
    T_in = 300.; % K
    L = 50.0; % cm
    k_0 = 4.2E15; % cm^3 /mol /min
    E = 18000.; % cal /mol
    Re = 1.987; % cal /mol /K

    % make unknown Vdot available to the derivatives function
    Vdot = nan;

    % derivatives function
    function derivs = eval_derivs(~, cur_dep)
        % extract ind and dep vars for the current integration step
        nDot_A_cur = cur_dep(1);
        nDot_Z_cur = cur_dep(2);
        T_cur = cur_dep(3);

        % calculate rate
        r = k_0*exp(-E/Re/T_cur)*(nDot_A_cur^2)/Vdot^2;

        % evaluate the derivatives
        dnDotAdz = -pi()*D^2/4*r;
        dnDotZdz = pi()*D^2/4*r;
        dTdz = -pi()*D^2/4*r*dH/Vdot/Cp;

        % return the derivatives
        derivs = [dnDotAdz; dnDotZdz; dTdz];
    end
    
    % residuals function
    function residuals = eval_resids(Vdot_guess)
        % make Vdot_guess available to the derivatives function
        Vdot = Vdot_guess;

        % define initial values
        ind_0 = 0.0;
        dep_0 = [nDot_A_in; nDot_Z_in; T_in];

        % define stopping criterion
        f_var = 0;
        f_val = L;

        % solve the IVODEs
        odes_are_stiff = false;
        [~, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @eval_derivs, odes_are_stiff);

        % check tha the solution is valid
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract Tf
        Tf_calc = dep(end,3);

        % evaluate and return the residual
        residuals = T_f - Tf_calc;
    end

    % initial guess
    init_guess = 1150.;

    % solve the ATE
    [Vdot, flag, message] = solve_ates(@eval_resids, init_guess);

    % check that the solution converged
    if flag <= 0
        disp(' ')
        disp(['The ATE solver did not converge: ',message])
    end

    % define initial values
    ind_0 = 0.0;
    dep_0 = [nDot_A_in; nDot_Z_in; T_in];

    % define stopping criterion
    f_var = 0;
    f_val = L;

    % solve the IVODEs
    odes_are_stiff = false;
    [~, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
        , @eval_derivs, odes_are_stiff);

    % check tha the solution is valid
    if flag <= 0
        disp(' ')
        disp('WARNING: The ODE solution may not be accurate!')
    end

    % extract the quantities of interest
    nDot_A_f = dep(end,1);
    nDot_Z_f = dep(end,2);

    % tabulate the results
    item = ["$\dot{V}$"; "$\dot{n}_{A,f}$"; "$\dot{n}_{A,f}$"];
    value = [Vdot; nDot_A_f; nDot_Z_f];
    units = ["cm^3^ min^-1^"; "mol min^-1^"; "mol min^-1^"];
    results_table = table(item,value,units);

    % display the results
    disp(' ')
    disp(['Volumetric Flow Rate: ', num2str(Vdot,3), ' cm^3^ min^-1^'])
    disp(['Final Molar Flow of A: ', num2str(nDot_A_f,3), ' mol min^-1^'])
    disp(['Final Molar Flow of Z: ', num2str(nDot_Z_f,3), ' mol min^-1^'])

    % save the results
    writetable(results_table,results_file_spec);
end