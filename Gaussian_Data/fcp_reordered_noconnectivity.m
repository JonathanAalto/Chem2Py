%Hi Jack! After copying the cartesion data from the .gjf to a .txt file, this code
%should identify the close phenyl rings and fix the positions of the
%carbons. The output .txt file can then be re-copied into the .gjf file in GaussView

clear all; clc; close all

filename = 'test_gjf.txt';
%For identifying the phenyl rings, use num_nearest = 7
num_nearest = 7;

filename_output = strcat('noconnect_output_',filename);
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


%Creates the set of export data with the fixed oxygens and carbons first
export_table = {};

for i = 1:length(data2.textdata)
    element = data2.textdata{i};
    if element == 'O'
        oxygen_coords = data2.data(i,2:4);
        last_line = size(export_table,1);
        export_table{last_line+1,1} = 'O';
        export_table{last_line+1,2} = -1;
        export_table{last_line+1,3} = oxygen_coords;
    end
end

for i = 1:length(data2.textdata)
    element = data2.textdata{i};
    if element == 'C'
        carbon_coords = data2.data(i,2:4);
        last_line = size(export_table,1);
        if sum(ismember(fixed_carbons(:,2:4),carbon_coords)) > 0;
            export_table{last_line+1,1} = 'C';
            export_table{last_line+1,2} = -1;
            export_table{last_line+1,3}= carbon_coords;
        end
    end
end

for i = 1:length(data2.textdata)
    element = data2.textdata{i};
    if element == 'C'
        carbon_coords = data2.data(i,2:4);
        last_line = size(export_table,1);
        if sum(ismember(fixed_carbons(:,2:4),carbon_coords)) == 0;
            export_table{last_line+1,1} = 'C';
            export_table{last_line+1,2} = 0;
            export_table{last_line+1,3} = carbon_coords;
        end
    end
end

for i = 1:length(data2.textdata)
    element = data2.textdata{i};
    %disp(element)
    if max(element ~= 'C') & max(element ~= 'O');
        current_coords = data2.data(i,2:4);
        last_line = size(export_table,1);
        export_table{last_line+1,1} = element;
        export_table{last_line+1,2} = 0;
        export_table{last_line+1,3} = current_coords;      
    end
end

fileout = fopen(filename_output, 'w');
for k1 = 1:size(export_table,1)-1;
    element = export_table{k1,1};
    if length(element) == 1;
        fprintf(fileout, ' %s                %1.0f    %9.7f   %9.7f    %9.7f\n', export_table{k1,:});
    else
        fprintf(fileout, ' %s               %1.0f    %9.7f   %9.7f    %9.7f\n', export_table{k1,:});
    end
end

end_element = export_table{end,1};

if length(end_element) == 1;
    fprintf(fileout, ' %s                %1.0f    %9.7f   %9.7f    %9.7f', export_table{end,:});
else
    fprintf(fileout, ' %s               %1.0f    %9.7f   %9.7f    %9.7f', export_table{end,:});
end

fclose('all');
