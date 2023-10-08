%% Run this Section First

%This code is adapted from a .m file sent by the Stoll Lab

close all;
clear, clc, clf;

[B_G,spc] = eprload('JG-2648-NT-forbidden-1p9dB_20scan-actual');
B = B_G/10; % convert field values from Gauss to mT
B_length = length(B);
seventh = floor(B_length/7);

%Plot raw data so you can identify the center, left edge, and right edge
figure(1)
subplot(411)
plot(B,spc);
title('Derivative Lineshape (Raw Data)')

%% Now Enter the Approximate X-Values for Left, Right, and Center of the Curve 

%Designate Left End of Curve
B_left = 144.3
%Designate Center of Curve 
B_center = 161.4 
%Designate Right End of Curve
B_right = 180.2

%Also set max and min x-values (for plotting only)
x_minimum = 120; 
x_maximum = 200; 
%% Now run this section to process the cumulative sum and plot everything else.

[left_point left_index] = min(abs(B-B_left));
[center_point center_index] = min(abs(B-B_center));
[right_point right_index] = min(abs(B-B_right));

%Left Baseline
B_sample1 = [B(1:left_index) ; ones(100,1)*B(center_index)];
spc_sample1 = [spc(1:left_index) ; ones(100,1)*spc(center_index)];
base_line1 = fit(B_sample1,spc_sample1,'poly2');
coeff_a1 = base_line1.p1; 
coeff_b1 = base_line1.p2;
coeff_c1 = base_line1.p3;
spc_baseline1 = coeff_a1.*(B.^2) + coeff_b1*B + coeff_c1;
corrected_spc1 = spc - spc_baseline1;


%Right Baseline
B_sample2 = [ones(100,1)*B(center_index) ; B(right_index:end)];
spc_sample2 = [ones(100,1)*spc(center_index) ; spc(right_index:end)];
base_line2 = fit(B_sample2,spc_sample2,'poly2');
coeff_a2 = base_line2.p1; 
coeff_b2 = base_line2.p2;
coeff_c2 = base_line2.p3;
spc_baseline2 = coeff_a2.*(B.^2) + coeff_b2*B + coeff_c2;
corrected_spc2 = spc - spc_baseline2;

%Combine into one Baseline
new_baseline = zeros(1,length(spc_baseline1));
for i = 1:length(spc_baseline1);
    if spc_baseline1(i) >= spc_baseline2(i);
        new_baseline(i) = spc_baseline2(i);
    else 
        new_baseline(i) = spc_baseline1(i);
    end
end

%Obtain Cumulative Sum
corrected_spc = spc - new_baseline';
cum_sum = cumsum(corrected_spc);

%Obtain Baseline for Cumulative Sum
new_slope = (cum_sum(right_index) - cum_sum(left_index))/(B(right_index) - B(left_index));
new_intercept = cum_sum(left_index) - new_slope*B(left_index);
cumsum_baseline_B = B(left_index:right_index);
cumsum_baseline_spc = new_slope*cumsum_baseline_B + new_intercept;

%Baseline Correction for Cumulative Sum 
corrected_cumsum = cum_sum;
corrected_cumsum(1:left_index-1) = corrected_cumsum(1:left_index-1) - cum_sum(left_index);
corrected_cumsum(right_index+1:end) = corrected_cumsum(right_index+1:end) - cum_sum(right_index);
corrected_cumsum(left_index:right_index) = corrected_cumsum(left_index:right_index) - cumsum_baseline_spc;

%Calculate Integral for Corrected Cumulative Sum 
cumulative_integral = sum(corrected_cumsum(left_index:right_index)');

%Plot Data
figure(1)
subplot(411)
plot(B,spc);
hold on
plot(B,new_baseline);
ylim([5, 10])
xlim([x_minimum, x_maximum])
title('Derivative Lineshape (Raw Data)')
subplot(412)
plot(B,corrected_spc)
ylim([-0.7, 1.2])
xlim([x_minimum, x_maximum])
title('Derivative Lineshape (Corrected)')
subplot(413)
plot(B,cum_sum)
hold on
plot(cumsum_baseline_B, cumsum_baseline_spc)
title('Integrated Lineshape (Cumulative Sum)')
ylim([-10, 90])
xlim([x_minimum, x_maximum])

subplot(414)
plot(B,corrected_cumsum)
hold on
title('Integrated Lineshape (Corrected)')
for i = left_index:10:right_index;
    x_data = B(i)*ones(1,25);
    y_data = linspace(0,corrected_cumsum(i),25);
    plot(x_data,y_data,'k')
end
ylim([-10, 90])
xlim([x_minimum, x_maximum])

dB = B_G(2)-B_G(1);
double_integral = sum(corrected_cumsum(left_index:right_index))*dB
% Now this is positive!


set(gcf,'Position',[100 100 800 1200])

