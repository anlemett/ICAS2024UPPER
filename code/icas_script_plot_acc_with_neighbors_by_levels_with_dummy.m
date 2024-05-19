%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICAS24: airspace capacity 
% Anastasia Lemetti
% MATLAB version: MATLAB R2024a
% 
% plot ACC EDUUUTAS and its neighbours inside the MUC FIR on all flight levels
% for day 2023-06-08, for time from 15.00 to 17.30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% 10 time intervals, 13 altitude bands, 4 ACCs

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
            struct_acc = acc_struct_arr.(acc);

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
                    end % airblocks
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
                adj_num = adj_num + 1;
            end

        end % ACCs
    end % altitude bands
end % time intervals

% Create dummy sector for all times and altitudes
dummy_airspace_pgons = cell(nT, num_a_bands, 1); 

% Iterate through times
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

        dummy_airspace_pgons{t, h, 1} = dummy_pgon;

    end
end



% Display all ACCs on all band altitudes for all time intervals

for t = 1:nT
    for h=1:num_a_bands

        figure; hold on;
        colors = lines(num_ACC+1);

        min_lon = 5;
        max_lon = 16;
        min_lat = 46;
        max_lat = 55;

        latlim = [min_lat max_lat];
        lonlim = [min_lon max_lon];
        xlim(lonlim);
        ylim(latlim);
        
        daspect([1 cos(mean(latlim)*pi/180) 1]);
        xlabel('Longitude [º]');
        ylabel('Latitude [º]');

        for a=1:num_ACC-1
        
            pgon = adjacent_airspace_pgons{t, h, a};

            if ~isempty(pgon.Vertices)

                fig = plot(pgon, 'FaceColor', 'none', 'EdgeColor', colors(a,:), 'linewidth', 2);
            end
        end

        % Display main acc
        main_pgon = main_airspace_pgons{t, h, 1};
        fig = plot(main_pgon, 'FaceColor', 'none', 'EdgeColor', colors(4,:), 'linewidth', 2);

        % Display dummy acc
        dummy_pgon = dummy_airspace_pgons{t,h,1};
        fig = plot(dummy_pgon, 'FaceColor', 'none', 'EdgeColor', colors(5,:), 'linewidth', 2);
    
        legend_names = acc_names;
        legend_names{length(legend_names)+1} = "dummy";
        legend(legend_names, 'Position', [0.0 0.0 0.2 0.4]);

        FL_start = flight_levels(h);
        FL_end = flight_levels(h+1);

        times = {'15:00', '15:15', '15:30', '15:45', '16:00', '16:15', '16:30', '16:45',...
            '17:00', '17:15', '17:30'};

        fig_title = strcat("Time   ", times{t});
        fig_title = strcat(fig_title, "-");
        fig_title = strcat(fig_title, times{t+1});
        fig_title = strcat(fig_title, " FL ");
        fig_title = strcat(fig_title, string(FL_start));
        fig_title = strcat(fig_title, "-");
        fig_title = strcat(fig_title, string(FL_end));
        title(fig_title);

        times = {'15_00', '15_15', '15_30', '15_45', '16_00', '16_15', '16_30', '16_45',...
            '17_00', '17_15', '17_30'};

        filename = strcat("Time_", times{t});
        filename = strcat(filename, "_");
        filename = strcat(filename, times{t+1});
        filename = strcat(filename, "_FL_");
        filename = strcat(filename, string(FL_start));
        filename = strcat(filename, "_");
        filename = strcat(filename, string(FL_end));

        full_filename = fullfile('.', 'figures', 'accs_by_time_FL_dummy', filename);
    
        saveas(fig, full_filename, 'png');
        clf(fig);
    end
end
close all;
