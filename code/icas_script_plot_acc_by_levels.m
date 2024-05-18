%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICAS24: airspace capacity 
% Anastasia Lemetti
% MATLAB version: MATLAB R2024a
% 
% plot ACC EDUUUTAS and its sectors on all flight levels for 2023-06-08,
% for configuration S6H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Airspace configuration
upper_sector_filename = fullfile('.', 'code_input', 'airspace_data', 'Upper_airspace', ...
    'fir_EDUU_2023-06-08.json');

upper_sector = jsondecode(fileread(upper_sector_filename));

exp_date = 'x2023_06_08';

acc_struct_arr = [upper_sector.(exp_date)];
acc_names = fieldnames(acc_struct_arr);

% Find the coordinates of EDUUUTAS and its sectors, configuration S6H
main_acc = 'EDUUUTAS';
confs = [upper_sector.(exp_date).(main_acc).configurations];

main_conf_name = 'S6H';
main_conf = confs.(main_conf_name);

flight_levels = [315 345 365 999];

num_FL = numel(flight_levels);
num_a_bands = num_FL - 1;

a_band_pgons = cell(1, num_FL-1);

el_sectors = [main_conf.elementarySectors];

sectors_names = fieldnames(el_sectors);

num_sec = numel(el_sectors_names); 


for i1=1:num_a_bands

    FL_start = flight_levels(i1);
    FL_end = flight_levels(i1+1);

    sectors_pgons = cell(1, num_sec);
    
    for i2 = 1:num_sec

        sectors_pgons{i2} = polyshape();

        el_sector = el_sectors.(el_sectors_names{i2});

        airblocks = [el_sector.airblocks];

        airblocks_names = fieldnames(airblocks);
        
        for i3 = 1:numel(airblocks_names)
            airblock = airblocks.(airblocks_names{i3});

            if (airblock.fl(1) <= FL_start) && (airblock.fl(2) >= FL_end)
              
                airblock_coord = airblock.polygon;
                airblock_coord(:, [1, 2]) = airblock_coord(:, [2, 1]);

                % Remove duplicate vertices
                airblock_coord = unique(airblock_coord, 'rows', 'stable');
                airblock_pgon = polyshape(airblock_coord);
                sectors_pgons{i2} = union(sectors_pgons{i2}, airblock_pgon);
            end
        end
        
        if ~isempty(sectors_pgons{i2}.Vertices)
            % Extract vertices from the original polyshape
            vertices = sectors_pgons{i2}.Vertices;

            % Remove duplicate vertices
            cleanedVertices = unique(vertices, 'rows', 'stable');
            % Create a new polyshape object without duplicate vertices
            sectors_pgons{i2} = polyshape(cleanedVertices);
        end

    end

    a_band_pgons{i1} = sectors_pgons;

end

% Display all sectors of EDUUUTAS, configuraion S6H, on all band altitudes

for i1=1:num_a_bands

    figure; hold on;
    colors = lines(num_sec);

    min_lon = 9;
    max_lon = 14;
    min_lat = 46;
    max_lat = 50;

    latlim = [min_lat max_lat];
    lonlim = [min_lon max_lon];
    xlim(lonlim);
    ylim(latlim);
    axis([min_lon max_lon min_lat max_lat]);
    daspect([1 cos(mean(latlim)*pi/180) 1]);
    xlabel('Longitude [º]');
    ylabel('Latitude [º]');

    for i2=1:num_sec
        
        pgon = a_band_pgons{i1}{i2};

        if ~isempty(pgon.Vertices)

            fig = plot(pgon, 'FaceColor', 'none', 'EdgeColor', colors(i2,:), 'linewidth', 2);
        end
    end

    legend(sectors_names, 'Position', [0.0 0.0 0.2 0.4]);

    FL_start = flight_levels(i1);
    FL_end = flight_levels(i1+1);
    plot_title = strcat("EDUUUTAS FL ", string(FL_start));
    plot_title = strcat(plot_title, "-");
    plot_title = strcat(plot_title, string(FL_end));
    title(plot_title);
    axis equal;

    filename = strcat("EDUUUTAS_FL_", string(FL_start));
    filename = strcat(filename, "_");
    filename = strcat(filename, string(FL_end));
    
    saveas(fig, fullfile('.', 'figures', 'sectors', filename), 'png');
    clf(fig)
end
