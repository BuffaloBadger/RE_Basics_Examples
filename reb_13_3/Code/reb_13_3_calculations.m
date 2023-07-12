function reb_13_3_calculations
% performs the calculations for the 2022 CE 329 Evaluation 3
    % set a range of reaction times
    times = linspace(1,100.0,1000);

    % calculate the yields and conversions for those times
    [yield, conversion] = reb_13_3_response(times);

    % find and report the maximum yield
    [yMax, iMax] = max(yield);
    tMax = times(iMax);
    convMax = conversion(iMax);

    % report and save the results
    disp(' ')
    disp(['Maximum Yield: ',num2str(100*yMax,3),'%'])
    disp(['Conversion: ',num2str(100*convMax,3),'%'])
    disp(['Reaction Time: ',num2str(tMax,3),' min'])
    results_file ="../Results/reb_13_3_results.csv";
    item = ["reaction time";"maximum yield";
        "conversion"];
    value = [round(tMax,1);round(100*yMax,1);round(100*convMax,1)];
    units = ["min";"%";"%"];
    results_table = table(item,value,units);
    writetable(results_table,results_file);

    % plot yield vs time
    figure
    plot(times,yield,'k','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Reaction Time (min)','FontSize', 14)
    ylabel('Yield (D/A)','FontSize', 14)
    saveas(gcf,"../Results/reb_13_3_yield_plot.png")
end
