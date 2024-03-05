function reb_K_3_matlab
%REB_K_3_2_CALCULATIONS Reaction Engineering Basics Example K.3

    % given and known constants
    D = 5.; % cm
    dH = -14000.; % cal /mol
    V_flow = 1150.; % cm^3 /min
    Cp = 1.3; % cal /cm^3 /K
    R_R = 1.3; %
    nDot_A_feed = 1.0; % mol /min
    L = 50.; % cm
    T_feed = 300.; % K
    k_0 = 4.2E15; % cm^3 /mol /min
    E = 18000.; % cal /mol
    Re = 1.987; % cal /mol /K

    % derivatives function
    function derivs = eval_derivs(~, cur_dep)
        % extract ind and dep vars for the current integration step
        nDot_A_cur = cur_dep(1);
        nDot_Z_cur = cur_dep(2);
        T_cur = cur_dep(3);

        % calculate rate
        r = k_0*exp(-E/Re/T_cur)*nDot_A_cur*nDot_Z_cur/V_flow^2;

        % evaluate the derivatives
        dnDotAdz = -pi()*D^2/4*r;
        dnDotZdz = pi()*D^2/4*r;
        dTdz = -pi()*D^2/4*r*dH/V_flow/Cp;

        % return the derivatives
        derivs = [dnDotAdz; dnDotZdz; dTdz];
    end
    
    % residuals function
    function residuals = eval_resids(guess)
        % extract the individual guesses
        nAin_guess = guess(1);
        nZin_guess = guess(2);
        Tin_guess = guess(3);

        % define initial values
        ind_0 = 0.0;
        dep_0 = [nAin_guess; nZin_guess; Tin_guess];

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

        % extract the outlet values
        nAf = dep(end,1);
        nZf = dep(end,2);
        Tf = dep(end,3);

        % evaluate the residuals
        residual1 = nDot_A_feed + R_R/(1+R_R)*nAf - nAin_guess;
        residual2 = R_R/(1+R_R)*nZf - nZin_guess;
        residual3 = Tin_guess - T_feed - R_R*(Tf - Tin_guess);

        % return the residuals
        residuals = [residual1; residual2; residual3];
    end

    % initial guess
    init_guess = [nDot_A_feed; nDot_A_feed; T_feed + 5.];

    % solve the ATE
    [soln, flag, message] = solve_ates(@eval_resids, init_guess);

    % check that the solution converged
    if flag <= 0
        disp(' ')
        disp(['The ATE solver did not converge: ',message])
    end

    % extract the results
    n_dot_A_in = soln(1);
    n_dot_Z_in = soln(2);
    T_in = soln(3);

    % initial values
    ind_0 = 0.0;
    dep_0 = [n_dot_A_in; n_dot_Z_in; T_in];

    % stopping criterion
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

    % extract the outlet values
    n_dot_A_out = dep(end,1);
    n_dot_Z_out = dep(end,2);
    T_out = dep(end,3);

    % Display the results
    disp(' ')
    disp(['Inlet molar flow of A: ',num2str(n_dot_A_in,3),' mol/min'])
    disp(['Inlet molar flow of Z: ',num2str(n_dot_Z_in,3),' mol/min'])
    disp(['Inlet temperature: ',num2str(T_in,3),' K'])
    disp(['Outlet molar flow of A: ',num2str(n_dot_A_out,3),' mol/min'])
    disp(['Outlet molar flow of Z: ',num2str(n_dot_Z_out,3),' mol/min'])
    disp(['Outlet temperature: ',num2str(T_out,3),' K'])
end