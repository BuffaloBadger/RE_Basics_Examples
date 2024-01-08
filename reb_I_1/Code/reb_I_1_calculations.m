function reb_I_1_calculations
% solve 3 ATEs in Reaction Engineering Basics Example I.6
    % Given and known constants
    nDotAin = 500.; % mol/h
    V = 500.; % L
    nDotZin = 0.; % mol/h
    Cp = 1170.; % cal/L/K
    Vdot = 250.; % L/h
    Tin = 423.; % K
    deltaH = 18200.; % cal/mol
    k0 = 1.14E9; % L/mol/h
    E = 16200.; % cal/mol
    R = 1.987; % cal/mol/K

    % Initial guess for the solution
    initialGuess = [nDotAin/2; nDotAin/2; Tin - 20.];

    % Solve the ATEs
    [solution, flag, message] = solveATEs(@evalResiduals...
        , initialGuess);

    % Check that the solution is converged
    if flag < 0
        disp(' ')
        disp(['The ATE solver did not converge: ',message])
    end

    % Display and save the results
    nA = solution(1);
    nZ = solution(2);
    T = solution(3);
    disp(' ')
    disp(['Flow Rate of A: ',num2str(nA,3),' mol/h'])
    disp(['Flow Rate of Z: ',num2str(nZ,3),' mol/h'])
    disp(['Temperature: ',num2str(T-273.15,3),' °C'])

    resultsFile ="../Results/reb_I_1_results.csv";
    item = ["Flow Rate of A";"Flow Rate of Z";"Temperature"];
    value = [nA;nZ;T-273.15];
    units = ["mol h^-1^";"mol h^-1^";"°C"];
    resultsTable = table(item,value,units);
    writetable(resultsTable,resultsFile);

    % Function that evaluates the residuals functions.
    function residuals = evalResiduals(guess)
        % Extract the guess values
        nDotAguess = guess(1);
        nDotZguess = guess(2);
        Tguess = guess(3);

        % Calculate the concentration of A
        CA = nDotAguess/Vdot;

        % Calculate the rate
        r = k0*exp(-E/R/Tguess)*CA^2;
        
        % Evaluate the residuals
        residuals(1) = nDotAin - nDotAguess - V*r;
        residuals(2) = nDotZin - nDotZguess + V*r;
        residuals(3) = Cp*Vdot*(Tguess - Tin) + V*r*deltaH;
    end
end