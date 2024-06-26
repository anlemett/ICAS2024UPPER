function [main_acc_data, adjacent_sectors_data] = icas_function_get_accs

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

nT = 10; % number of time intervals
N = 7; % number of adjacent ACCs

flight_levels = [315 325 335 345 355 365 375 385 999];

num_FL = numel(flight_levels);
num_FL_bands = num_FL - 1;

% Initialize output variables
main_acc_data = cell(nT, num_FL_bands, 1); 
adjacent_sectors_data = cell(nT, num_FL_bands, N); 

% Airspace configuration
upper_sector_filename = fullfile('.', 'code_input', 'airspace_data', 'Upper_airspace', ...
    'fir_nextto_EDMMCTAA_upper_2023-06-08.json');

upper_sector = jsondecode(fileread(upper_sector_filename));

% For each time, find sector configuration in table 'configuration_upper_20230608_1500_1730.xlsx'

full_filename = fullfile('.', 'code_input', 'airspace_data',...
        'Upper_airspace',...
        'configuration_upper_20230608_1500_1730.xlsx');

% Create import options for the Excel file
opts = detectImportOptions(full_filename);

opts.VariableNamingRule = 'preserve'; % Preserve the original variable names

% Set the import options to read all columns as text (string)
opts = setvartype(opts, 'string');

T = readtable(full_filename, opts);


nT = 10; % number of time intervals

config_vec = cell(1, nT);


exp_date = 'x2023_06_08';

acc_struct_arr = [upper_sector.(exp_date)];
acc_names = fieldnames(acc_struct_arr);

num_ACC = numel(acc_names);


main_airspace_pgons = cell(nT, num_FL_bands, 1); 

adjacent_airspace_pgons = cell(nT, num_FL_bands, N); 

% 10 time intervals, 8 FL bands, 7 ACCs

% Iterate through times
for t = 1:nT

    % Iterate through flight level bands
    for h = 1:num_FL_bands

        FL_start = flight_levels(h);
        FL_end = flight_levels(h+1);

        adj_num = 1;
      
        % Iterate through all ACCs
        LOVV1CTA_in_use = false;
        LOVVCTA_in_use = false;

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
                    end % airblock
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

                main_acc_data{t, h, 1}(1,1) = struct('properties',[], 'geometry', []);

                main_acc_data{t,h,1}(1).properties.DESIGNATOR = acc;
                main_acc_data{t,h,1}(1).properties.LOWER_LIMIT_VALUE = flight_levels(h);
                main_acc_data{t,h,1}(1).properties.UPPER_LIMIT_VALUE = flight_levels(h+1);

                main_acc_data{t,h,1}(1).geometry.coordinates = zeros(1,length(acc_pgon.Vertices),2);
                main_acc_data{t,h,1}(1).geometry.coordinates(:,:,1) = acc_pgon.Vertices(:,1)';
                main_acc_data{t,h,1}(1).geometry.coordinates(:,:,2) = acc_pgon.Vertices(:,2)';
            
            else
                if ~isempty(acc_pgon.Vertices)
                    if strcmp(acc, 'LOVV1CTA')
                        LOVV1CTA_in_use = true;
                    end
                    if strcmp(acc, 'LOVVCTA')
                        LOVVCTA_in_use = true;
                    end
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
            end
        end % ACCs

        if LOVV1CTA_in_use
            % From EDUUUTAS (1) remove intersection with LOVV1CTA (5)
            acc_pgon = main_airspace_pgons{t, h, 1};
            intersectionPoly = intersect(acc_pgon, adjacent_airspace_pgons{t, h, 5});
            acc_pgon = subtract(acc_pgon, intersectionPoly);
            main_airspace_pgons{t, h, 1} = acc_pgon;
            main_acc_data{t,h,1}(1).geometry.coordinates = zeros(1,length(acc_pgon.Vertices),2);
            main_acc_data{t,h,1}(1).geometry.coordinates(:,:,1) = acc_pgon.Vertices(:,1)';
            main_acc_data{t,h,1}(1).geometry.coordinates(:,:,2) = acc_pgon.Vertices(:,2)';
        end
        if LOVV1CTA_in_use && LOVVCTA_in_use
            % From LOVVCTA (6) remove intersection with LOVV1CTA (5)
            acc_pgon = adjacent_airspace_pgons{t, h, 6};
            intersectionPoly = intersect(acc_pgon, adjacent_airspace_pgons{t, h, 5});
            acc_pgon = subtract(acc_pgon, intersectionPoly);
            adjacent_airspace_pgons{t, h, 6} = acc_pgon;
            adjacent_sectors_data{t,h,6}(1).geometry.coordinates = zeros(1,length(acc_pgon.Vertices),2);
            adjacent_sectors_data{t,h,6}(1).geometry.coordinates(:,:,1) = acc_pgon.Vertices(:,1)';
            adjacent_sectors_data{t,h,6}(1).geometry.coordinates(:,:,2) = acc_pgon.Vertices(:,2)';
            LOVV1CTA_in_use = false;
            LOVVCTA_in_use = false;

        end
    end % altitude bands
end % time intervals


end % function
