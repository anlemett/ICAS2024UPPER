%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICAS24: airspace capacity 
% Anastasia Lemetti
% MATLAB version: MATLAB R2024a
% 
% returns ACC EDUUUTAS polygons
% on all flight levels for day 2023-06-08, for time from 15.00 to 17.30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main_acc_data = icas_function_create_main_acc()

nT = 10;
num_a_bands = 6;
main_acc_data = cell(nT, num_a_bands, 1); % Initialize output variable

% Airspace configuration
upper_sector_filename = fullfile('.', 'code_input', 'airspace_data', 'Upper_airspace', ...
    'fir_EDUU_2023-06-08.json');

upper_sector = jsondecode(fileread(upper_sector_filename));

% For each time, find sector configuration in table 'configuration_upper_20230608_1500_1730.xlsx'

T = readtable(fullfile('.', 'code_input', 'airspace_data',...
        'Upper_airspace',...
        'configuration_upper_20230608_1500_1730.xlsx'), ...
        'FileType', 'spreadsheet', 'ReadVariableNames', true); % Read xlsx file
  
nT = 10; % number of time intervals

config_vec = cell(1, nT);


exp_date = 'x2023_06_08';

acc_struct_arr = [upper_sector.(exp_date)];
acc_names = fieldnames(acc_struct_arr);

num_ACC = numel(acc_names);

flight_levels = [315 325 345 355 365 375 999];
% 315   345   365   999 for ACC EDUUUTAS conf. S6H
num_FL = numel(flight_levels);
num_a_bands = num_FL - 1;

main_airspace_pgons = cell(nT, num_a_bands, 1); 

% 10 time intervals, 13 altitude bands, 1 main ACC

% Iterate through times
for t = 1:nT

    % Iterate through altitude bands
    for h = 1:num_a_bands

        FL_start = flight_levels(h);
        FL_end = flight_levels(h+1);
      
        % Iterate through ACCs
        for a = 1:num_ACC
    
            acc = acc_names{a};

            if ~strcmp(acc, 'EDUUUTAS')
                continue; % skip all except main ACC
            end

            config_vec = T.(acc);

            confs = [upper_sector.(exp_date).(acc).configurations];

            conf_names = fieldnames(confs);
            num_conf = numel(conf_names);

            acc_pgon = polyshape();

            % Iterate through configurations
            for i1 = 1:num_conf

                % Take only configuration for time t
                if ~(strcmp(config_vec{t}, conf_names{i1}))
                    continue;
                end

                conf = confs.(conf_names{i1});

                el_sectors = [conf.elementarySectors];

                el_sectors_names = fieldnames(el_sectors);

                num_sec = numel(el_sectors_names); 

                % Iterate through sectors
                for i2 = 1:num_sec

                    el_sector = el_sectors.(el_sectors_names{i2});

                    airblocks = [el_sector.airblocks];

                    airblocks_names = fieldnames(airblocks);
                    num_blocks = numel(airblocks_names);
        
                    % Iterate through airblocks
                    for i3 = 1:num_blocks
                        airblock = airblocks.(airblocks_names{i3});

                        if (airblock.fl(1) <= FL_start) && (airblock.fl(2) >= FL_end)
              
                            airblock_coord = airblock.polygon;
                            airblock_coord(:, [1, 2]) = airblock_coord(:, [2, 1]);

                            % Remove duplicate vertices
                            airblock_coord = unique(airblock_coord, 'rows', 'stable');
                            airblock_pgon = polyshape(airblock_coord);
                            acc_pgon = union(acc_pgon, airblock_pgon);
                        end
                    end % airblocke
                end % sectors
            end % configurations (used only one)
        
            if ~isempty(acc_pgon.Vertices)
                % Extract vertices from the original polyshape
                vertices = acc_pgon.Vertices;

                % Remove duplicate vertices
                cleanedVertices = unique(vertices, 'rows', 'stable');
                % Create a new polyshape object without duplicate vertices
                acc_pgon = polyshape(cleanedVertices);
            end
            
            main_airspace_pgons{t, h, 1} = acc_pgon;

            main_acc_data{t, h, 1}(1,1) = struct('properties',[], 'geometry', []);

            main_acc_data{t,h,1}(1).properties.DESIGNATOR = acc;
            main_acc_data{t,h,1}(1).properties.LOWER_LIMIT_VALUE = flight_levels(h);
            main_acc_data{t,h,1}(1).properties.UPPER_LIMIT_VALUE = flight_levels(h+1);

            main_acc_data{t,h,1}(1).geometry.coordinates = zeros(1,length(acc_pgon.Vertices),2);
            main_acc_data{t,h,1}(1).geometry.coordinates(:,:,1) = acc_pgon.Vertices(:,1)';
            main_acc_data{t,h,1}(1).geometry.coordinates(:,:,2) = acc_pgon.Vertices(:,2)';


        end % ACCs
    end % altitude bands
end % time intervals
