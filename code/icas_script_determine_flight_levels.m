%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICAS24: airspace capacity 
% Anastasia Lemetti
% MATLAB version: MATLAB R2024a
% 
% determine flight levels in upper airspace, ACC EDUUUTAS,
% and its neighbours  inside the MUC FIR for 2023-06-08
% for time from 15.00 to 17.30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Airspace configuration
upper_sector_filename = fullfile('.', 'code_input', 'airspace_data', 'Upper_airspace', ...
    'fir_EDUU_2023-06-08.json');

upper_sector = jsondecode(fileread(upper_sector_filename));

exp_date = 'x2023_06_08';

acc_struct_arr = [upper_sector.(exp_date)];
acc_arr = fieldnames([upper_sector.(exp_date)]);

% For each time, find sector configuration in table 'configuration_upper_20230608_1500_1730.xlsx'

T = readtable(fullfile('.', 'code_input', 'airspace_data',...
        'Upper_airspace',...
        'configuration_upper_20230608_1500_1730.xlsx'), ...
        'FileType', 'spreadsheet', 'ReadVariableNames', true); % Read xlsx file
  
nT = 10; % number of time intervals

config_vec = cell(1, nT);

flight_levels = {};
flight_levels_from = {};
flight_levels_to = {};

for i = 1:numel(acc_arr)
    
    struct_acc = acc_arr{i};
    acc = char(struct_acc);

    config_vec = T.(acc);

    acc_configs = unique(config_vec);
    
    %struct_acc = acc_struct_arr.EDUUUTAS;
    %acc = 'EDUUUTAS';

    confs = [upper_sector.(exp_date).(acc).configurations];

    conf_names = fieldnames(confs);

    % Loop over the fields
    for ii = 1:numel(conf_names)

        if ~any(strcmp(acc_configs, conf_names{ii}))
            continue;
        end

        %if not(strcmp(conf_names{ii}, 'S6H'))
        %    continue
        %end

        conf = confs.(conf_names{ii});

        el_sectors = [conf.elementarySectors];

        el_sectors_names = fieldnames(el_sectors);

        for j = 1:numel(el_sectors_names)

            el_sector = el_sectors.(el_sectors_names{j});

            airblocks = [el_sector.airblocks];

            airblocks_names = fieldnames(airblocks);

            for jj = 1: numel(airblocks_names)
                airblock = airblocks.(airblocks_names{jj});
                flight_levels{numel(flight_levels)+1} = airblock.fl(1);
                flight_levels{numel(flight_levels)+1} = airblock.fl(2);
                %flight_levels_from{numel(flight_levels_from)+1} = airblock.fl(1);
                %flight_levels_to{numel(flight_levels_to)+1} = airblock.fl(2);
            end
        end
    end
end

flight_levels = unique(cell2mat(flight_levels));
disp(flight_levels); % 315   345   365   999 for ACC EDUUUTAS conf. S6H

%[235 245 255 265 285 295 305 315 325 345 355 365 375 999] for all ACCs

%flight_levels = unique(cell2mat(flight_levels_from));
%disp(flight_levels);
%flight_levels = unique(cell2mat(flight_levels_to));
%disp(flight_levels);
