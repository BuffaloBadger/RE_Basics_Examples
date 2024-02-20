function reb_J_7_calculations
%REB_J_7_CALCULATIONS solve IVODEs in Reaction Engineering Basics 
%   Example J.7
    
    % set filepath
    filepath_to_results = '../results/';

    % given and known constants
    V = 2.0; % m^3
    dH_1 = -3470; % cal/mol
    Cp_A = 7.5; % cal/mol/K
    Cp_B = 8.5; % cal/mol/K
    Cp_Y = 12.1; % cal/mol/K
    Cp_Z = 5.7; % cal/mol/K
    k_0_1 = 83; % m^3 /mol /h
    E_1 = 10200; % cal/mol
    % Known
    Re = 1.987; % cal/mol
    Rw = 8.206E-5; % m^3 atm/mol/K

    % derivatives function
    function derivs = eval_derivs(~, cur_dep)
        % extract the dependent variables needed to evaluate the
        % derivatives for the current integration step
        nA_cur = cur_dep(1);
        nB_cur = cur_dep(2);
        nY_cur = cur_dep(3);
        nZ_cur = cur_dep(4);
        T_cur = cur_dep(5);
        
        % create mass matrix, settinng all elements to zero
        mass_matrix = zeros(6,6);

        % add 1 on the diagonal for the first 4 rows
        mass_matrix(1,1) = 1;
        mass_matrix(2,2) = 1;
        mass_matrix(3,3) = 1;
        mass_matrix(4,4) = 1;

        % add the elements for the energy balance
        mass_matrix(5,5) = nA_cur*Cp_A + nB_cur*Cp_B + nY_cur*Cp_Y ...
            + nZ_cur*Cp_Z;
        mass_matrix(5,6) = -V*Re/Rw;

        % add the elements for the ideal gas law equation
        mass_matrix(6,1) = Rw*T_cur;
        mass_matrix(6,2) = Rw*T_cur;
        mass_matrix(6,3) = Rw*T_cur;
        mass_matrix(6,4) = Rw*T_cur;
        mass_matrix(6,5) = Rw*(nA_cur + nB_cur + nY_cur + nZ_cur);
        mass_matrix(6,6) = -V;

        % calculate the rate
        CA = nA_cur/V;
        CB = nB_cur/V;
        r = k_0_1*exp(-E_1/Re/T_cur)*CA*CB;

        % create right side vector
        rhs = [-V*r; -V*r; V*r; V*r; -V*r*dH_1; 0];

        % evaluate the derivatives
        derivs = mass_matrix\rhs;
    end

    % initial values
    ind_0 = 0.0;
    dep_0 = [190.0; 190.0; 0.0; 0.0; 450.0; 7.0];

    % stopping criterion
    f_var = 0;
    f_val = 2.0;

    % solve the IVODEs
    odes_are_stiff = false;
    [t, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
        , @eval_derivs, odes_are_stiff);
    if flag <= 0
        disp(' ')
        disp('WARNING: The ODE solution may not be accurate!')
    end

    % extract the results
    nA = dep(:,1);
    nB = dep(:,2);
    nY = dep(:,3);
    nZ = dep(:,4);
    T = dep(:,5);
    P = dep(:,6);

    % tabulate the results
    results_table = table(t,nA,nB,nY,nZ,T,P);

    % display and save the results
    disp(' ')
    disp(results_table)
    results_file =strcat(filepath_to_results,"reb_J_7_results.csv");
    writetable(results_table,results_file);
end