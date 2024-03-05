function reb_J_7_3_matlab
%REB_J_7_3_MATLAB Reaction Engineering Basics Example J.7.3

    % given and known constants
    T_out = 400.; % K
    D = 1. ; % in
    k_0 = 7.49E9*61.02; % in^3 /mol /min
    E = 15300.; % cal /mol
    P = 4.; % atm
    dH = -14500.; % cal /mol
    Cp_A = 10.9; % cal /mol /K
    Cp_Z = 21.8; % cal /mol /K
    L = 100.; % in
    nDot_A_in = 1.5; % mol /min
    Re = 1.987; % cal /mol /K
    Rw = 0.08206*61.02; % in^3 atm /mol /K

    % derivatives function
    function derivs = eval_derivs(~, cur_dep)
        % extract ind and dep vars for the current integration step
        nDot_A_cur = cur_dep(1);
        nDot_Z_cur = cur_dep(2);
        T_cur = cur_dep(3);

        % calculate rate
        r = k_0*exp(-E/Re/T_cur)...
            *(nDot_A_cur*P/Rw/T_cur/(nDot_A_cur + nDot_Z_cur))^2;

        % evaluate the derivatives
        dnDotAdz = -2*pi()*D^2/4*r;
        dnDotZdz = pi()*D^2/4*r;
        dTdz = -pi()*D^2/4*r*dH/(nDot_A_cur*Cp_A + nDot_Z_cur*Cp_Z);

        % return the derivatives
        derivs = [dnDotAdz; dnDotZdz; dTdz];
    end
    
    % residuals function
    function residuals = eval_resids(T_in_guess)
        % define initial values
        ind_0 = 0.0;
        dep_0 = [nDot_A_in; 0.0; T_in_guess];

        % define stopping criterion
        f_var = 0;
        f_val = L;

        % solve the IVODEs
        odes_are_stiff = false;
        [z, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @eval_derivs, odes_are_stiff);

        % check tha the solution is valid
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract Tf
        Tf = dep(end,3);

        % evaluate and return the residual
        residuals = Tf - T_out;
    end

    % initial guess
    init_guess = T_out - 100.0;

    % solve the ATE for T_in
    [T_in_soln, flag, message] = solve_ates(@eval_resids, init_guess);

    % check that the solution converged
    if flag <= 0
        disp(' ')
        disp(['The ATE solver did not converge: ',message])
    end

    % define initial values
    ind_0 = 0.0;
    dep_0 = [nDot_A_in; 0.0; T_in_soln];

    % define stopping criterion
    f_var = 0;
    f_val = L;

    % solve the IVODEs
    odes_are_stiff = false;
    [~, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
        , @eval_derivs, odes_are_stiff);

    % check that the solution is valid
    if flag <= 0
        disp(' ')
        disp('WARNING: The ODE solution may not be accurate!')
    end

    % extract the results
    nDot_A = dep(:,1);
    nDot_Z = dep(:,2);
    T = dep(:,3);

    % tabulate the results
    profile_results_table = table(z, nDot_A, nDot_Z, T);

    % display the results
    disp(' ')
    disp(['Inlet Temperature: ', num2str(T_in_soln,3), ' K'])
    disp(' ')
    disp('Molar flow and Temperature Profiles')
    disp(profile_results_table)
end