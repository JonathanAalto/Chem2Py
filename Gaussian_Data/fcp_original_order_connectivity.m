%Hi Jack! After copying the cartesion data from the .gjf to a .txt file, this code
%should identify the close phenyl rings and fix the positions of the
%carbons. The output .txt file can then be re-copied into the .gjf file in GaussView

clear all; clc; close all

filename = 'NiCl2 coords.txt';
%For identifying the phenyl rings, use num_nearest = 7
num_nearest = 7;

filename_output = strcat('original_order_output_',filename);
data2 = importdata(filename);

output_matrix = [];
fixed_carbons = [];
second_closest_carbons = [];

%Finds the carbons attached to the carboxyl groups
for i = 1:length(data2.textdata)
    element = data2.textdata{i};
    if element == 'O'
        oxygen_coords = data2.data(i,2:4);
        carbon_list = [];
        for p = 1:length(data2.textdata)
                current_element = data2.textdata{p};
                if current_element == 'C';
                    current_carbon_coords = data2.data(p,2:4);
                    current_distance = norm(current_carbon_coords - oxygen_coords);
                    carbon_list = [carbon_list ; current_distance current_carbon_coords];
                end
        end
        carbon_list = sortrows(carbon_list);
        second_closest_carbons = [second_closest_carbons ; carbon_list(2,:)];
    end
end

carboxylate_carbons = [];

for k = 1:size(second_closest_carbons,1);
    current_carbon = second_closest_carbons(k,2:4);
    for m = 1:size(second_closest_carbons,1);
         if m == k;
             continue
         end
          
         current_carbon2 = second_closest_carbons(m,2:4);
         if sum(current_carbon - current_carbon2) == 0;
             carboxylate_carbons = [carboxylate_carbons ; current_carbon];
         end
    end
end

carboxylate_carbons = unique(carboxylate_carbons,'rows');


%Finds the closest carbons to each of the carboxylate carbons. 

for n = 1:size(carboxylate_carbons,1);
    six_closest = [];
    carbon_list = [];
    for p = 1:length(data2.textdata)
            current_element = data2.textdata{p};
            if current_element == 'C';
                current_carbon_coords = data2.data(p,2:4);
                current_distance = norm(current_carbon_coords - carboxylate_carbons(n,:));
                carbon_list = [carbon_list ; current_distance current_carbon_coords];
            end
    end
    carbon_list = sortrows(carbon_list);
    six_closest = carbon_list(1:num_nearest,:);
    fixed_carbons = [fixed_carbons ; six_closest];
end

%Writes the data to a .txt file in the same atom order as the input, which
%allows for the connectivity data to remain valid. 

fileout = fopen(filename_output, 'w');
for k1 = 1:length(data2.textdata)-1;
    element = data2.textdata{k1};
    if element == 'O';
        fprintf(fileout, ' O               -1    %9.7f   %9.7f    %9.7f\n', data2.data(k1,2:4));
    elseif element == 'C';
        if sum(ismember(fixed_carbons(:,2:4),data2.data(k1,2:4),'rows')) == 1;
            fprintf(fileout, ' C               -1    %9.7f   %9.7f    %9.7f\n', data2.data(k1,2:4));
        else 
            fprintf(fileout, ' C                0    %9.7f   %9.7f    %9.7f\n', data2.data(k1,2:4));
        end
    elseif length(element) == 2;
        fprintf(fileout, ' %s               0    %9.7f   %9.7f    %9.7f\n', element, data2.data(k1,2:4));
    else
        fprintf(fileout, ' %s                0    %9.7f   %9.7f    %9.7f\n', element, data2.data(k1,2:4));
    end
end

k1 = length(data2.textdata);
element = data2.textdata{k1};
if element == 'O';
    fprintf(fileout, ' O               -1    %9.7f   %9.7f    %9.7f', data2.data(k1,2:4));
elseif element == 'C';
    if sum(ismember(fixed_carbons(:,2:4),data2.data(k1,2:4),'rows')) == 1;
        fprintf(fileout, ' C               -1    %9.7f   %9.7f    %9.7f', data2.data(k1,2:4));
    else 
        fprintf(fileout, ' C                0    %9.7f   %9.7f    %9.7f', data2.data(k1,2:4));
    end
elseif length(element) == 2;
    fprintf(fileout, ' %s               0    %9.7f   %9.7f    %9.7f', element, data2.data(k1,2:4));
else
    fprintf(fileout, ' %s                0    %9.7f   %9.7f    %9.7f', element, data2.data(k1,2:4));
end

fclose('all');
