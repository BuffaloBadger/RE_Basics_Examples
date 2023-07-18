function reb_13_4_calculations

    % set a range of coolant flow rates
    mDotEx = linspace(100.0,250.0,100);

    % allocate storage for the corresponding net rates
    rNet = nan(100,1);

    % calculate the net rates
    for i = 1:100
        [rNet(i),~,~,~] = reb_13_4_response(mDotEx(i));
    end

    % find the maximum rate
    [rMax, iMax] = max(rNet);
    mDotExMax = mDotEx(iMax);

    % report and save the results
    disp(' ')
    disp(['Maximum Net Rate: ',num2str(rMax,3),' mol/min'])
    disp(['Coolant Flow: ',num2str(mDotExMax,3),' g/min'])
    results_file ="../Results/reb_13_4_results.csv";
    item = ["maximum rate";"coolant flow rate"];
    value = [round(rMax,4);round(mDotExMax,0)];
    units = ["mol/min";"g/min"];
    results_table = table(item,value,units);

    % plot rate vs coolant flow
    figure
    plot(mDotEx,rNet,'k','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Coolant Flow Rate (g/min)','FontSize', 14)
    ylabel('Net Rate (mol/min)','FontSize', 14)
    saveas(gcf,"../Results/reb_13_4_net_rate_plot.png")

    % plot conversion and temperature profiles
    [~, t, fA, T] = reb_13_4_response(mDotEx(iMax));

    figure
    plot(t, fA, 'k','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Reaction Time (min)','FontSize', 14)
    ylabel('Conversion of A (%)','FontSize', 14)
    saveas(gcf,"../Results/reb_13_4_conversion_profile.png")

    figure
    plot(t, T, 'k','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Reaction Time (min)','FontSize', 14)
    ylabel('Temperature (°C)','FontSize', 14)
    saveas(gcf,"../Results/reb_13_4_temperature_profile.png")

    % calculate and plot the rate profile
    k01 = 2.59E9; % /min
    E1 = 16500; % cal/mol
    R = 1.987; % cal/mol/K
    CA0 = 2.0; % mol/L
    CA = CA0*(1.0 - fA/100.0);
    r = k01*exp(-E1/R./T).*CA;

    figure
    semilogy(t, r, 'k','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Reaction Time (min)','FontSize', 14)
    ylabel('Reaction Rate (M/min)','FontSize', 14)
    saveas(gcf,"../Results/reb_13_4_rate_profile.png")

    [rMax, iMax] = max(r);
    rdata = {'maximum rate',round(rMax,3,'significant'),'mol/min/L';
        'time of maximum r',round(t(iMax),1),'min.'};
    [Tmax, iMax] = max(T);
    Tdata = {'maximum T',round(Tmax,1),'°C';
        'time of maximum T',round(t(iMax),1),'min.'};
    convData = {'conversion at maximum T',round(fA(iMax),1),'%';
        'final conversion',round(fA(end),1),'%'};
    results_table = [results_table;rdata;Tdata;convData];
    writetable(results_table,results_file);
end
