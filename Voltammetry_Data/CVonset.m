%% Enter file name and the name you would like for the exported data
clear all; close all; clc;


file_name = 'PhCN_DiPh_PMDI_50.csv'; %To do

CVdata = importdata(file_name);

save_name = 'PhCN_DiPh_PMDI_50';

%% First run this section to find the approximate max and baseline locations from the plot. 

data = CVdata.data;
data = sortrows(data,[1 2]);
voltage = CVdata.data(:,1);
current = CVdata.data(:,2);
figure(1)
plot(voltage,current,'r-');
volt_max = max(voltage);
volt_min = min(voltage);
volt_range = range(voltage);
volt_partitions = [volt_min:volt_range/100:volt_max];
volt_reductive = [];
current_reductive = [];

%% Enter the approximate locations of the max and baseline
approximate_max_location = -1.75;
baseline_bounds = [-0.806 -0.78];


%% Now run this section to compute the onset potential

%isolates the lower half of the scan (reductive portion)
for i = 1:length(volt_partitions)-1;
    test_range = [volt_partitions(i) volt_partitions(i+1)];
    test_volts = [];
    test_current = [];
    for j = 1:length(voltage);
        if voltage(j) >= test_range(1) & voltage(j) < test_range(2);
            test_volts = [test_volts voltage(j)];
            test_current = [test_current current(j)];
            continue
        end
    end
    test_current_mean = mean(test_current);
    accepted_current_vals = []; 
    accepted_volt_vals = [];
    for k = 1:length(test_current)
        if test_current(k) < test_current_mean;
            accepted_current_vals = [accepted_current_vals test_current(k)];
            accepted_volt_vals = [accepted_volt_vals test_volts(k)];
            continue
        end
    end
    current_reductive = [current_reductive accepted_current_vals];
    volt_reductive = [volt_reductive accepted_volt_vals];
end
plot(voltage, current);
figure(2)
plot(volt_reductive, current_reductive);
hold on 


%baseline region
baseline_volts = [];
baseline_current = [];
for m = 1:length(volt_reductive);
    if volt_reductive(m) >= baseline_bounds(1) & volt_reductive(m) < baseline_bounds(2)
        baseline_volts = [baseline_volts volt_reductive(m)];
        baseline_current = [baseline_current current_reductive(m)];
        continue
    end
end
current_min = max(baseline_current);

%baseline best_fit 
[baseline_params, Set_baseline] = Zlstsq(baseline_volts,baseline_current,1);


baseline_plot_x = volt_reductive;
baseline_plot_y = baseline_params(1).*baseline_plot_x + baseline_params(2); 

current_reductive_corrected = current_reductive - baseline_plot_y;

%find max
max_region = [];
for m = 1:length(volt_reductive);
    if (abs(volt_reductive(m) - approximate_max_location) < 0.1)
        max_region = [max_region m];
        continue
    end
end

[current_max current_max_index] = min(current_reductive_corrected(max_region));
max_region_voltage = volt_reductive(max_region);
peak_voltage = max_region_voltage(current_max_index);
peak_current = current_reductive_corrected(current_max_index);



%halfmax region
current_reductive2 = current_reductive_corrected(volt_reductive>peak_voltage);
volt_reductive2 = volt_reductive(volt_reductive>peak_voltage);

halfmax_volts = [];
halfmax_current = [];
[closest_current closest_index] = min(abs(current_reductive2 - (current_max/2)));
halfmax_voltage = volt_reductive2(closest_index);
halfmax_current_val = current_reductive2(closest_index) + baseline_params(1)*halfmax_voltage + baseline_params(2);
for m = 1:length(volt_reductive);
    if abs(volt_reductive(m) - halfmax_voltage) < 0.03;
        halfmax_volts = [halfmax_volts volt_reductive(m)];
        halfmax_current = [halfmax_current current_reductive(m)];
        continue
    end
end

%halfmax best_fit 
[halfmax_params, Set_halfmax] = Zlstsq(halfmax_volts,halfmax_current,1);

%Onset potential
onset_volt_val = (halfmax_params(2) - baseline_params(2))/(baseline_params(1)-halfmax_params(1));
onset_current_val = baseline_params(1)*onset_volt_val + baseline_params(2);

%bestfit xy data
halfmax_plot_x = linspace(halfmax_voltage-0.1,onset_volt_val+0.05,length(volt_reductive));
halfmax_plot_y = halfmax_params(1).*halfmax_plot_x + halfmax_params(2);


plot(halfmax_plot_x,halfmax_plot_y,'k-')
plot(baseline_plot_x,baseline_plot_y,'r-')
plot(onset_volt_val,onset_current_val,'b*')

plot_table_name = strcat(save_name,'_plot_table.csv');
plot_lines = [baseline_plot_x' baseline_plot_y' halfmax_plot_x' halfmax_plot_y'];
plot_table = array2table(plot_lines,'VariableNames',{'baseline_V' 'baseline_I' 'halfmax_line_V' 'halfmax_line_I'});
writetable(plot_table,plot_table_name)

onset_point_name = strcat(save_name,'_special_points.csv');
onset_point = array2table([onset_volt_val onset_current_val peak_voltage peak_current halfmax_voltage halfmax_current_val],'VariableNames',{'onset_V' 'onset_I' 'peak_V' 'peak_I' 'halfmax_V' 'halfmax_I'});
writetable(onset_point,onset_point_name)
