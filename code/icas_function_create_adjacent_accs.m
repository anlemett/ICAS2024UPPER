function adjacent_sectors_data = icas_function_create_adjacent_sectors

% Create variable with adjacent sector data

% all_adjacent_sectors_data{nT, num_a_band, N}

% adjacent_sectors_data is a nT x num_a_band x N cell
% nT - number of time intervals
% num_a_band - number of altitude bands
% N - number of adjacent sectors
% Each cell contains the information of the M=1 subsectors
% that form each sector, wich is a Mx1 (1x1) structure with fields:
% properties: 1x1 structure with fields DESIGNATOR, UPPER_LIMIT_VALUE,
% LOWER_LIMIT_VALUE
% geometry: 1x1 structure with field coordinates

% N = 3
nT = 10;
num_a_bands = 6;
N = 4;
adjacent_sectors_data = cell(nT, num_a_bands, N); % Initialize output variable

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

adjacent_airspace_pgons = cell(nT, num_a_bands, num_ACC-1); 

% 10 time intervals, 6 altitude bands, 4 ACCs

% Iterate through times
for t = 1:nT

    % Iterate through altitude bands
    for h = 1:num_a_bands

        FL_start = flight_levels(h);
        FL_end = flight_levels(h+1);

        adj_num = 1;
      
        % Iterate through ACCs
        for a = 1:num_ACC
    
            acc = acc_names{a};

            config_vec = T.(acc);

            confs = [upper_sector.(exp_date).(acc).configurations];

            conf_names = fieldnames(confs);
            num_conf = numel(conf_names);

            % Iterate through configurations
            for i1 = 1:num_conf

                % Take only configuration for time t
                if ~(strcmp(config_vec{t}, conf_names{i1}))
                    continue;
                end

                conf = confs.(conf_names{i1});

                acc_pgon = polyshape();

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
            

            if strcmp(acc, 'EDUUUTAS')
                main_airspace_pgons{t, h, 1} = acc_pgon;
            
            else

                adjacent_airspace_pgons{t, h, adj_num} = acc_pgon;

                adjacent_sectors_data{t, h, adj_num}(1,1) = struct('properties',[], 'geometry', []);

                adjacent_sectors_data{t,h,adj_num}(1).properties.DESIGNATOR = acc;
                adjacent_sectors_data{t,h,adj_num}(1).properties.LOWER_LIMIT_VALUE = flight_levels(h);
                adjacent_sectors_data{t,h,adj_num}(1).properties.UPPER_LIMIT_VALUE = flight_levels(h+1);

                adjacent_sectors_data{t,h,adj_num}(1).geometry.coordinates = zeros(1,length(acc_pgon.Vertices),2);
                adjacent_sectors_data{t,h,adj_num}(1).geometry.coordinates(:,:,1) = acc_pgon.Vertices(:,1)';
                adjacent_sectors_data{t,h,adj_num}(1).geometry.coordinates(:,:,2) = acc_pgon.Vertices(:,2)';

                adj_num = adj_num +1;
            end
        end % ACCs
    end % altitude bands
end % time intervals

% Create dummy ACC for all times and altitudes
for t = 1:nT

    % Iterate through altitude bands
    for h = 1:num_a_bands

        FL_start = flight_levels(h);
        FL_end = flight_levels(h+1);

        main_acc = main_airspace_pgons{t, h, 1};
        acc1 = adjacent_airspace_pgons{t, h, 1};
        acc2 = adjacent_airspace_pgons{t, h, 2};
        acc3 = adjacent_airspace_pgons{t, h, 3};

        together = union(main_acc,acc1);
        together = union(together,acc2);
        together = union(together,acc3);

        dummy_pgon = polyshape([6 16 16],[46 46 55]);
        dummy_pgon = subtract(dummy_pgon, together);

        adjacent_sectors_data{t, h, 4}(1,1) = struct('properties',[], 'geometry', []);

        adjacent_sectors_data{t,h,4}(1).properties.DESIGNATOR = "Dummy";
        adjacent_sectors_data{t,h,4}(1).properties.LOWER_LIMIT_VALUE = flight_levels(h);
        adjacent_sectors_data{t,h,4}(1).properties.UPPER_LIMIT_VALUE = flight_levels(h+1);

        adjacent_sectors_data{t,h,4}(1).geometry.coordinates = zeros(1,length(dummy_pgon.Vertices),2);
        adjacent_sectors_data{t,h,4}(1).geometry.coordinates(:,:,1) = dummy_pgon.Vertices(:,1)';
        adjacent_sectors_data{t,h,4}(1).geometry.coordinates(:,:,2) = dummy_pgon.Vertices(:,2)';

    end
end

end % function
