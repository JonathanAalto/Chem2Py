%% Run This Section First

%This code is adapted from a .m file sent by the Stoll Lab

close all;
clear, clc, clf;

[B_G,spc,Pars] = eprload('JG-2648-NT-allowed-9db-20scan'); %Get data from file
B = B_G/10;              % convert field values from Gauss to mT

figure(1)
subplot(411)
plot(B,spc);
title('Derivative Lineshape (Raw Data)')

%% Enter Left and Right Sides of the EPR Signature
left_edge = 275;
right_edge = 355;
coefficients = polyfit(B(B<left_edge),spc(B<left_edge), 1);
spc = spc - coefficients(1).*(B) - coefficients(2); %Corrects baseline in spectrum

x_max = max(B); 
x_min = min(B);
y_max = max(spc)+10; 
y_min = min(spc)-10; 

%Plots the raw data to make sure it looks as expected. 
figure(1)
subplot(411)
plot(B,spc);
xlim([x_min x_max]);
ylim([y_min y_max])
title('Derivative Lineshape (Raw Data)')

%% In this section, enter the G-value and G-strain of radical peak and the approximate edges of the downward peak.

radical_G = 2.002860333; %These values are obtained from "Radical fit.xlsx" document in Jack's EPR folder
radical_gStrain = 0.002444853;

peak_left = 335; 
peak_right = 350; 

%% This section contains parameter ranges that the code will scan to to obtain the optimal peak to subtract. 

amp_range = [0.001 0.5]; %Range of amplitudes that the code will try
lw1_range = [0.08832384]; %Range of parameters for Gaussian smoothing. 
%Right now this range is the original broadening parameter from the radical
%peak raw data, but if you want to vary it, set lw1_range = [0 4] or
%something like that. 


lw2_range = [0.69617629]; %Range of parameters for Lorentzian smoothing
%Right now this range is the original broadening parameter from the radical
%peak raw data, but if you want to vary it, set lw2_range = [0 4] or
%something like that. 

test_increments = [0.001, 0.05, 0.05]; 
%This dictates the step size when varying the amplitude and broadening
%parameters. 

%% This section uses the data provided in previous sections to find and subtract the optimal radical peak
tic
amp_list1 = amp_range(1):test_increments(1):amp_range(end);

%Isolate the narrow regions of the raw data around the artifact 
B_new = B(B>peak_left);
B_new = B_new(B_new<peak_right);

spc_new = spc(B>peak_left);
spc_new = spc_new(B_new<peak_right);

%In this for loop, a set of test peaks is generated, and each test peak is
%evaluated based on whether, when subtracted, the resulting curve has any
%peak artifacts (identified as places where the derivative changes sign
%where it shouldn't) 

Sys1.g = radical_G; Sys1.gStrain = radical_gStrain;
Exp1.mwFreq = 9.636032; Exp1.nPoints = 10*numel(B_new); Exp1.Range = [min(B_new) max(B_new)];
valid_pars1 = [];
for i1 = 1:length(amp_list1)
    for j1 = lw1_range(1):test_increments(2):lw1_range(end)
        for k1 = lw2_range(1):test_increments(3):lw2_range(end)
            Sys1.lw = [j1 k1];
            [Field1,Spec1] = pepper(Sys1,Exp1);
            spc_sub = spc_new - amp_list1(i1)*Spec1(1:10:end)';
            first_deriv = diff(spc_sub)/(B_new(2)-B_new(1));
            second_deriv = diff(first_deriv)/(B_new(2)-B_new(1));
            second_deriv_sum = sum(abs(second_deriv));
            num_signchange = numel(find(diff(sign(first_deriv))));
            if num_signchange == 1
                valid_pars1 = [valid_pars1 ; amp_list1(i1) j1 k1 second_deriv_sum];
            end
        end
    end
end
toc
%% Find and plot the best set of parameters


%Find Set of Parameters with Smallest Second Derivative Area (smoothest
%curve)
second_derivs = valid_pars1(:,4);
[min_ min_index] = min(second_derivs);
best_parameters = valid_pars1(min_index,:);
    
%% Construct the radical peak with the best set of parameters. 
close all
best_parameters(1)=0.29;
Sys2.lwpp = [0.08832384 1.69617629];
Sys2.g = radical_G; Sys2.gStrain = radical_gStrain; 
Exp2.mwFreq = 9.636032; Exp2.nPoints = 10*numel(B_new); Exp2.Range = [min(B_new) max(B_new)];
[Field2,Spec2] = pepper(Sys2,Exp2);
background = ones(1,length(B));
radical_on_background = [background(B<peak_left)*Spec2(1) Spec2(1:10:end) background(B>peak_right)*Spec2(end)];

spc_subtracted = spc' - best_parameters(1)*radical_on_background; 
figure(2)
%subplot(311)
%plot(B,spc)
hold on 
plot(B,best_parameters(1)*radical_on_background)
hold on 
plot(B,spc_subtracted)




subplot(312)
plot(B,best_parameters(1)*radical_on_background)
%ylim([-100 20])
subplot(313)
plot(B,spc_subtracted)
%ylim([-100 20])