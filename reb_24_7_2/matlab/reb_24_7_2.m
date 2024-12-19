function reb_24_7_2()
% Calculations for Reaction Engineering Basics Example 24.7.2
    % given and known constants
    k0 = 2.65e6; % m^3/kmol/s
    E = 15e6; % cal/kmol
    L = 15.0; % m
    D = 0.1; % m
    VFR = 0.1/1000.0/60.0; % m^3/s
    dH = -0.75E7; % cal/kmol
    Cp = 1700; % cal/kg/K
    rho = 650; % kg/m^3
    CAfeed = 5.0; % kmol/m^3
    CBfeed = 5.0; % kmol/m^3
    Tfeed = 293.15; % K
    R = 1987; % cal/kmol/K
    V = pi*D^2/4*L; % m^3

    nAfeed = CAfeed*VFR;
    nBfeed = CBfeed*VFR;

    % global constants
    g_Dax = nan;
    g_lambda = nan;

    % CSTR residuals function
    function epsilon = residuals(guess)
        % extract the individual guesses
        nA = guess(1);
        nB = guess(2);
        nZ = guess(3);
        T = guess(4);

        % rate
        k = k0*exp(-E/R/T);
        CA = nA/VFR;
        CB = nB/VFR;
        r = k*CA*CB;

        % evaluate and return the residuals
        epsilon_1 = nAfeed - nA - V*r;
        epsilon_2 = nBfeed - nB - V*r;
        epsilon_3 = -nZ + V*r;
        epsilon_4 = -VFR*rho*Cp*(T - Tfeed) - V*r*dH;
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4];
    end

    % CSTR reactor function
    function [nA1, nB1, nZ1, T1] = CSTR_unknowns()
        guess = [0.1*nAfeed; nBfeed - 0.9*nAfeed; 0.9*nAfeed...
            ; Tfeed + 40.];

        % solve the ATEs
        [soln, flag, message] = solve_ates(@residuals, guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end
        
        % extract results
        nA1 = soln(1);
        nB1 = soln(2);
        nZ1 = soln(3);
        T1 = soln(4);
    end

    % derivatives function
    function ddz = pfr_derivatives(~,dep)
        % extract the dependent variables
        nA = dep(1);
        nB = dep(2);
        T = dep(4);

        % calculate the unknowns
        k = k0*exp(-E/R/T);
        CA = nA/VFR;
        CB = nB/VFR;
        r = k*CA*CB;

        % evaluate the derivatives
        dnAdz = -pi()*D^2/4*r;
        dnBdz = -pi()*D^2/4*r;
        dnZdz = pi()*D^2/4*r;
        dTdz = -pi()*D^2/4*r*dH/VFR/rho/Cp;

        % return the derivatives
        ddz = [dnAdz; dnBdz; dnZdz; dTdz];
    end

    % PFR model
    function [z, nA, nB, nZ, T] = pfr_profiles()
        % initial values
        ind_0 = 0;

        dep_0 = [nAfeed; nBfeed; 0; Tfeed];

        % stopping criterion
        f_var = 0;
        f_val = L;

        % solve the design 
        odes_are_stiff = false;
        [z, dep, flag, message] = solve_ivodes(ind_0, dep_0, f_var...
            , f_val, @pfr_derivatives, odes_are_stiff);

        % check that the solver was successful
        if (flag <= 0.0)
            disp(" ")
            disp(["error solving the ivodes: ",message])
        end

        % extract and return the profiles
        nA = dep(:,1);
        nB = dep(:,2);
        nZ = dep(:,3);
        T = dep(:,4);
    end

    % axial dispersion derivatives function
    function ddz = ad_derivatives(~,dep)
        % extract the dependent variables that are needed
        CA = dep(1);
        CB = dep(2);
        T = dep(4);
        wA = dep(5);
        wB = dep(6);
        wZ = dep(7);
        wT = dep(8);

        % set the dispersion coefficients
        Dax = g_Dax;
        lambda = g_lambda;

        % rate
        k = k0*exp(-E/R/T);
        r = k*CA*CB;

        % evaluate the derivatives
        dCAdz = wA;
        dCBdz = wB;
        dCZdz = wZ;
        dTdz = wT;
        dwAdz = -1.0/Dax*(-4*VFR*wA/pi/D^2 - r);
        dwBdz = -1.0/Dax*(-4*VFR*wB/pi/D^2 - r);
        dwZdz = -1.0/Dax*(-4*VFR*wZ/pi/D^2 + r);
        dwTdz = -1.0/lambda*(-4*VFR*rho*Cp*wT/pi/D^2 - r*dH);

        % combine the derivatives as a vector and return
        ddz = [dCAdz; dCBdz; dCZdz; dTdz; dwAdz; dwBdz; dwZdz; dwTdz];
    end

    % boundary value residuals functions
    function epsilon = ad_residuals(dep_lb, dep_ub)
        % extract the boundary values that are needed
        CA0 = dep_lb(1);
        CB0 = dep_lb(2);
        CZ0 = dep_lb(3);
        T0 = dep_lb(4);
        wA0 = dep_lb(5);
        wB0 = dep_lb(6);
        wZ0 = dep_lb(7);
        wT0 = dep_lb(8);
        wA1 = dep_ub(5);
        wB1 = dep_ub(6);
        wZ1 = dep_ub(7);
        wT1 = dep_ub(8);

        % set the dispersion coefficients
        Dax = g_Dax;
        lambda = g_lambda;

        % evaluate the residuals
        epsilon_1 = CA0 - (CAfeed + Dax*pi*D^2*wA0/4/VFR);
        epsilon_2 = CB0 - (CBfeed + Dax*pi*D^2*wB0/4/VFR);
        epsilon_3 = CZ0 - (Dax*pi*D^2*wZ0/4/VFR);
        epsilon_4 = T0 - (Tfeed + lambda*pi*D^2*wT0/4/VFR/rho/Cp);
        epsilon_5 = wA1;
        epsilon_6 = wB1;
        epsilon_7 = wZ1;
        epsilon_8 = wT1;
        
        % combine the residuals as a vector and return
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4...
            ; epsilon_5; epsilon_6; epsilon_7; epsilon_8];
    end

    % axial dispersion reactor function
    function [z, CA, CB, CZ, T, wA, wB, wZ, wT] = ad_profiles(guess)
        % solve the bvodes
        % set the upper and lower bounds of the independent variable
        xl = 0.0;
        xu = L;
        % solve the bvodes
        [z, dep] = solve_bvodes(xl, xu, guess, @ad_derivatives...
            , @ad_residuals);

        % extract the profiles
        CA = dep(1,:);
        CB = dep(2,:);
        CZ = dep(3,:);
        T = dep(4,:);
        wA = dep(5,:);
        wB = dep(6,1);
        wZ = dep(7,1);
        wT = dep(8,1);
    end

    % master function
    function perform_the_calculations()
        % solve the pfr design equations
        [zpfr, nApfr, ~, ~, Tpfr] = pfr_profiles();
        CApfr = nApfr/VFR;
        Tpfr = Tpfr - 273.15;
        lpfr = "0 (PFR)";
        fApfr = 100.0*(nAfeed - nApfr(end))/nAfeed;
        Tf_pfr = Tpfr(end);

        % solve the CSTR design equations
        [nA, ~, nZcstr, Tcstr] = CSTR_unknowns();
        CAcstr = nA/VFR;
        Tcstr = Tcstr - 273.15;
        lcstr = "\infty (CSTR)";
        fAcstr = 100.0*(nAfeed - nA)/nAfeed;
        Tf_cstr = Tcstr;

        % make a guess for the axial dispersion solution
        guess = [CAcstr; CAcstr; VFR*nZcstr; Tcstr + 273.15...
            ; (CAcstr-CAfeed)/L; (CAcstr-CAfeed)/L; VFR*nZcstr/L ...
            ; (Tfeed - Tcstr)/L];

        % solve the axial dispersion equations using CSTR-like values for
        % the dispersion coefficients
        g_Dax = 1.0E0;
        g_lambda = 1.0E6;
        [zad1, CAad1, ~, ~, Tad1, ~, ~, ~, ~] = ad_profiles(guess);
        Tad1 = Tad1 - 273.15;
        lad1 = "1 x 10^0, 1 x 10^6";
        fAad1 = 100.0*(nAfeed - VFR*CAad1(end))/nAfeed;
        Tf_ad1 = Tad1(end);

        % vary lambda from the CSTR-like value
        g_Dax = 1.0E0;
        g_lambda = 1.0E4;
        [zad2, CAad2, ~, ~, Tad2, ~, ~, ~, ~] = ad_profiles(guess);
        Tad2 = Tad2 - 273.15;
        lad2 = "1 x 10^{0}, 1 x 10^4";
        fAad2 = 100.0*(nAfeed - VFR*CAad2(end))/nAfeed;
        Tf_ad2 = Tad2(end);

        g_Dax = 1.0E0;
        g_lambda = 1.0E3;
        [zad3, CAad3, ~, ~, Tad3, ~, ~, ~, ~] = ad_profiles(guess);
        Tad3 = Tad3 - 273.15;
        lad3 = "1 x 10^{0}, 1 x 10^3";
        fAad3 = 100.0*(nAfeed - VFR*CAad3(end))/nAfeed;
        Tf_ad3 = Tad3(end);

        g_Dax = 1.0E0;
        g_lambda = 1.0E1;
        [zad4, CAad4, ~, ~, Tad4, ~, ~, ~, ~] = ad_profiles(guess);
        Tad4 = Tad4 - 273.15;
        lad4 = "1 x 10^{0}, 1 x 10^1";
        fAad4 = 100.0*(nAfeed - VFR*CAad4(end))/nAfeed;
        Tf_ad4 = Tad4(end);

        Dax_and_lambda = [lcstr; lad1; lad2; lad3; lad4; lpfr];
        Pct_conversion = [fAcstr; fAad1; fAad2; fAad3; fAad4; fApfr];
        Tf = [Tf_cstr; Tf_ad1; Tf_ad2; Tf_ad3; Tf_ad4; Tf_pfr];
        comparison_table = table(Dax_and_lambda,Pct_conversion,Tf);
        disp(comparison_table)

        fa = figure;
        hold on
        yline(CAcstr,'k:','LineWidth',3)
        plot(zad1,CAad1,'r','LineWidth',2)
        plot(zad2,CAad2,'g','LineWidth',2)
        plot(zad3,CAad3,'b','LineWidth',2)
        plot(zad4,CAad4,'Color','#EDB120','LineWidth',2)
        plot(zpfr,CApfr,'k--','LineWidth',3)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Axial Position (m)','FontSize', 14)
        ylabel('Concentration of A (kmol m^-^3)','FontSize', 14)
        lgd = legend({lcstr, lad1, lad2, lad3, lad4, lpfr},'Location'...
            ,'east','FontSize',14);
        lgd.Title.String = 'D_{ax},  \lambda_{ax}';
        saveas(gcf,"CSTR_like_CA_profile_varying_lambda.png")
        set(fa,'Units','Inches');
        pos = get(fa,'Position');
        set(fa,'PaperPositionMode','Auto','PaperUnits','Inches'...
            ,'PaperSize',[pos(3), pos(4)])
        print(fa,'CSTR_like_CA_profile_varying_lambda.pdf','-dpdf','-r0')

        ft = figure;
        hold on
        yline(Tcstr,'k:','LineWidth',3)
        plot(zad1,Tad1,'r','LineWidth',2)
        plot(zad2,Tad2,'g','LineWidth',2)
        plot(zad3,Tad3,'b','LineWidth',2)
        plot(zad4,Tad4,'Color','#EDB120','LineWidth',2)
        plot(zpfr,Tpfr,'k--','LineWidth',3)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Axial Position (m)','FontSize', 14)
        ylabel('Temperature (°C)','FontSize', 14)
        lgd = legend({lcstr, lad1, lad2, lad3, lad4, lpfr},'Location'...
            ,'southeast','FontSize',14);
        lgd.Title.String = 'D_{ax},  \lambda_{ax}';
        saveas(gcf,"CSTR_like_T_profile_varying_lambda.png")
        set(ft,'Units','Inches');
        pos = get(ft,'Position');
        set(ft,'PaperPositionMode','Auto','PaperUnits','Inches'...
            ,'PaperSize',[pos(3), pos(4)])
        print(ft,'CSTR_like_T_profile_varying_lambda.pdf','-dpdf','-r0')

        % vary Dax from the CSTR-like value
        g_Dax = 1.0E-2;
        g_lambda = 1.0E6;
        [zad2, CAad2, ~, ~, Tad2, ~, ~, ~, ~] = ad_profiles(guess);
        Tad2 = Tad2 - 273.15;
        lad2 = "1 x 10^{-2}, 1 x 10^6";
        fAad2 = 100.0*(nAfeed - VFR*CAad2(end))/nAfeed;
        Tf_ad2 = Tad2(end);

        g_Dax = 1.0E-3;
        g_lambda = 1.0E6;
        [zad3, CAad3, ~, ~, Tad3, ~, ~, ~, ~] = ad_profiles(guess);
        Tad3 = Tad3 - 273.15;
        lad3 = "1 x 10^{-3}, 1 x 10^6";
        fAad3 = 100.0*(nAfeed - VFR*CAad3(end))/nAfeed;
        Tf_ad3 = Tad3(end);

        g_Dax = 1.0E-5;
        g_lambda = 1.0E6;
        [zad4, CAad4, ~, ~, Tad4, ~, ~, ~, ~] = ad_profiles(guess);
        Tad4 = Tad4 - 273.15;
        lad4 = "1 x 10^{-5}, 1 x 10^6";
        fAad4 = 100.0*(nAfeed - VFR*CAad4(end))/nAfeed;
        Tf_ad4 = Tad4(end);

        Dax_and_lambda = [lcstr; lad1; lad2; lad3; lad4; lpfr];
        Pct_conversion = [fAcstr; fAad1; fAad2; fAad3; fAad4; fApfr];
        Tf = [Tf_cstr; Tf_ad1; Tf_ad2; Tf_ad3; Tf_ad4; Tf_pfr];
        comparison_table = table(Dax_and_lambda,Pct_conversion,Tf);
        disp(comparison_table)

        fa = figure;
        hold on
        yline(CAcstr,'k:','LineWidth',3)
        plot(zad1,CAad1,'r','LineWidth',2)
        plot(zad2,CAad2,'g','LineWidth',2)
        plot(zad3,CAad3,'b','LineWidth',2)
        plot(zad4,CAad4,'Color','#EDB120','LineWidth',2)
        plot(zpfr,CApfr,'k--','LineWidth',3)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Axial Position (m)','FontSize', 14)
        ylabel('Concentration of A (kmol m^-^3)','FontSize', 14)
        lgd = legend({lcstr, lad1, lad2, lad3, lad4, lpfr},'Location'...
            ,'east','FontSize',14);
        lgd.Title.String = 'D_{ax},  \lambda_{ax}';
        saveas(gcf,"CSTR_like_CA_profile_varying_Dax.png")
        set(fa,'Units','Inches');
        pos = get(fa,'Position');
        set(fa,'PaperPositionMode','Auto','PaperUnits','Inches'...
            ,'PaperSize',[pos(3), pos(4)])
        print(fa,'CSTR_like_CA_profile_varying_Dax.pdf','-dpdf','-r0')

        ft = figure;
        hold on
        yline(Tcstr,'k:','LineWidth',3)
        plot(zad1,Tad1,'r','LineWidth',2)
        plot(zad2,Tad2,'g','LineWidth',2)
        plot(zad3,Tad3,'b','LineWidth',2)
        plot(zad4,Tad4,'Color','#EDB120','LineWidth',2)
        plot(zpfr,Tpfr,'k--','LineWidth',3)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Axial Position (m)','FontSize', 14)
        ylabel('Temperature (°C)','FontSize', 14)
        lgd = legend({lcstr, lad1, lad2, lad3, lad4, lpfr},'Location'...
            ,'east','FontSize',14);
        lgd.Title.String = 'D_{ax},  \lambda_{ax}';
        saveas(gcf,"CSTR_like_T_profile_varying_Dax.png")
        set(ft,'Units','Inches');
        pos = get(ft,'Position');
        set(ft,'PaperPositionMode','Auto','PaperUnits','Inches'...
            ,'PaperSize',[pos(3), pos(4)])
        print(ft,'CSTR_like_T_profile_varying_Dax.pdf','-dpdf','-r0')

        % % solve the axial dispersion equations using PFR-like values for
        % the dispersion coefficients
        g_Dax = 1.0E-5;
        g_lambda = 1.0E1;
        [zad1, CAad1, ~, ~, Tad1, ~, ~, ~, ~] = ad_profiles(guess);
        Tad1 = Tad1 - 273.15;
        lad1 = "1 x 10^{-5}, 1 x 10^1";
        fAad1 = 100.0*(nAfeed - VFR*CAad1(end))/nAfeed;
        Tf_ad1 = Tad1(end);

        % vary lambda from the PFR-like value
        g_Dax = 1.0E-5;
        g_lambda = 1.0E3;
        [zad2, CAad2, ~, ~, Tad2, ~, ~, ~, ~] = ad_profiles(guess);
        Tad2 = Tad2 - 273.15;
        lad2 = "1 x 10^{-5}, 1 x 10^3";
        fAad2 = 100.0*(nAfeed - VFR*CAad2(end))/nAfeed;
        Tf_ad2 = Tad2(end);

        g_Dax = 1.0E-5;
        g_lambda = 1.0E4;
        [zad3, CAad3, ~, ~, Tad3, ~, ~, ~, ~] = ad_profiles(guess);
        Tad3 = Tad3 - 273.15;
        lad3 = "1 x 10^{-5}, 1 x 10^4";
        fAad3 = 100.0*(nAfeed - VFR*CAad3(end))/nAfeed;
        Tf_ad3 = Tad3(end);

        g_Dax = 1.0E-5;
        g_lambda = 1.0E6;
        [zad4, CAad4, ~, ~, Tad4, ~, ~, ~, ~] = ad_profiles(guess);
        Tad4 = Tad4 - 273.15;
        lad4 = "1 x 10^{-5}, 1 x 10^6";
        fAad4 = 100.0*(nAfeed - VFR*CAad4(end))/nAfeed;
        Tf_ad4 = Tad4(end);

        Dax_and_lambda = [lcstr; lad1; lad2; lad3; lad4; lpfr];
        Pct_conversion = [fAcstr; fAad1; fAad2; fAad3; fAad4; fApfr];
        Tf = [Tf_cstr; Tf_ad1; Tf_ad2; Tf_ad3; Tf_ad4; Tf_pfr];
        comparison_table = table(Dax_and_lambda,Pct_conversion,Tf);
        disp(comparison_table)

        fa = figure;
        hold on
        yline(CAcstr,'k:','LineWidth',3)
        plot(zad1,CAad1,'r','LineWidth',2)
        plot(zad2,CAad2,'g','LineWidth',2)
        plot(zad3,CAad3,'b','LineWidth',2)
        plot(zad4,CAad4,'Color','#EDB120','LineWidth',2)
        plot(zpfr,CApfr,'k--','LineWidth',3)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Axial Position (m)','FontSize', 14)
        ylabel('Concentration of A (kmol m^-^3)','FontSize', 14)
        lgd = legend({lcstr, lad1, lad2, lad3, lad4, lpfr},'Location'...
            ,'east','FontSize',14);
        lgd.Title.String = 'D_{ax},  \lambda_{ax}';
        saveas(gcf,"PFR_like_CA_profile_varying_lambda.png")
        set(fa,'Units','Inches');
        pos = get(fa,'Position');
        set(fa,'PaperPositionMode','Auto','PaperUnits','Inches'...
            ,'PaperSize',[pos(3), pos(4)])
        print(fa,'PFR_like_CA_profile_varying_lambda.pdf','-dpdf','-r0')

        ft = figure;
        hold on
        yline(Tcstr,'k:','LineWidth',3)
        plot(zad1,Tad1,'r','LineWidth',2)
        plot(zad2,Tad2,'g','LineWidth',2)
        plot(zad3,Tad3,'b','LineWidth',2)
        plot(zad4,Tad4,'Color','#EDB120','LineWidth',2)
        plot(zpfr,Tpfr,'k--','LineWidth',3)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Axial Position (m)','FontSize', 14)
        ylabel('Temperature (°C)','FontSize', 14)
        lgd = legend({lcstr, lad1, lad2, lad3, lad4, lpfr},'Location'...
            ,'east','FontSize',14);
        lgd.Title.String = 'D_{ax},  \lambda_{ax}';
        saveas(gcf,"PFR_like_T_profile_varying_lambda.png")
        set(ft,'Units','Inches');
        pos = get(ft,'Position');
        set(ft,'PaperPositionMode','Auto','PaperUnits','Inches'...
            ,'PaperSize',[pos(3), pos(4)])
        print(ft,'PFR_like_T_profile_varying_lambda.pdf','-dpdf','-r0')

        % vary Dax from the PFR-like value
        g_Dax = 1.0E-3;
        g_lambda = 1.0E1;
        [zad2, CAad2, ~, ~, Tad2, ~, ~, ~, ~] = ad_profiles(guess);
        Tad2 = Tad2 - 273.15;
        lad2 = "1 x 10^{-3}, 1 x 10^1";
        fAad2 = 100.0*(nAfeed - VFR*CAad2(end))/nAfeed;
        Tf_ad2 = Tad2(end);

        g_Dax = 1.0E-2;
        g_lambda = 1.0E1;
        [zad3, CAad3, ~, ~, Tad3, ~, ~, ~, ~] = ad_profiles(guess);
        Tad3 = Tad3 - 273.15;
        lad3 = "1 x 10^{-2}, 1 x 10^1";
        fAad3 = 100.0*(nAfeed - VFR*CAad3(end))/nAfeed;
        Tf_ad3 = Tad3(end);

        g_Dax = 1.0E0;
        g_lambda = 1.0E1;
        [zad4, CAad4, ~, ~, Tad4, ~, ~, ~, ~] = ad_profiles(guess);
        Tad4 = Tad4 - 273.15;
        lad4 = "1 x 10^{0}, 1 x 10^1";
        fAad4 = 100.0*(nAfeed - VFR*CAad4(end))/nAfeed;
        Tf_ad4 = Tad4(end);

        Dax_and_lambda = [lcstr; lad1; lad2; lad3; lad4; lpfr];
        Pct_conversion = [fAcstr; fAad1; fAad2; fAad3; fAad4; fApfr];
        Tf = [Tf_cstr; Tf_ad1; Tf_ad2; Tf_ad3; Tf_ad4; Tf_pfr];
        comparison_table = table(Dax_and_lambda,Pct_conversion,Tf);
        disp(comparison_table)

        fa = figure;
        hold on
        yline(CAcstr,'k:','LineWidth',3)
        plot(zad1,CAad1,'r','LineWidth',2)
        plot(zad2,CAad2,'g','LineWidth',2)
        plot(zad3,CAad3,'b','LineWidth',2)
        plot(zad4,CAad4,'Color','#EDB120','LineWidth',2)
        plot(zpfr,CApfr,'k--','LineWidth',3)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Axial Position (m)','FontSize', 14)
        ylabel('Concentration of A (kmol m^-^3)','FontSize', 14)
        lgd = legend({lcstr, lad1, lad2, lad3, lad4, lpfr},'Location'...
            ,'east','FontSize',14);
        lgd.Title.String = 'D_{ax},  \lambda_{ax}';
        saveas(gcf,"PFR_like_CA_profile_varying_Dax.png")
        set(fa,'Units','Inches');
        pos = get(fa,'Position');
        set(fa,'PaperPositionMode','Auto','PaperUnits','Inches'...
            ,'PaperSize',[pos(3), pos(4)])
        print(fa,'PFR_like_CA_profile_varying_Dax.pdf','-dpdf','-r0')

        ft = figure;
        hold on
        yline(Tcstr,'k:','LineWidth',3)
        plot(zad1,Tad1,'r','LineWidth',2)
        plot(zad2,Tad2,'g','LineWidth',2)
        plot(zad3,Tad3,'b','LineWidth',2)
        plot(zad4,Tad4,'Color','#EDB120','LineWidth',2)
        plot(zpfr,Tpfr,'k--','LineWidth',3)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Axial Position (m)','FontSize', 14)
        ylabel('Temperature (°C)','FontSize', 14)
        lgd = legend({lcstr, lad1, lad2, lad3, lad4, lpfr},'Location'...
            ,'southeast','FontSize',14);
        lgd.Title.String = 'D_{ax},  \lambda_{ax}';
        saveas(gcf,"PFR_like_T_profile_varying_Dax.png")
        set(ft,'Units','Inches');
        pos = get(ft,'Position');
        set(ft,'PaperPositionMode','Auto','PaperUnits','Inches'...
            ,'PaperSize',[pos(3), pos(4)])
        print(ft,'PFR_like_T_profile_varying_Dax.pdf','-dpdf','-r0')
    end

    % perform the calculations
    perform_the_calculations();
end