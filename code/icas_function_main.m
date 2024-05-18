function [ASCR, sector_names, sector_time, sector_data] = FMPMet_WP6_main_function()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICAS2024 Airspace capacity reduction
% Anastasia Lemetti
% MATLAB version: MATLAB R2024a
% 
% Inputs:
%   1) Weather data (xml)
% 
%   2) Flight data: Flight plans (s06)
% 
%   3) Airspace data (json)
%
% Output:
%   Available sector capacity ratio ASCR for EDUUUTAS ACC
%  for day 2023-06-08, time from 15.00 to 17.30


    % Vector of time stamps and airspace configuration

    % Time 
    nT = 10;

    % Time: from 15.00 to 17.15 (end time to 17.30)
    minut_vec = 00:15:15*(nT-1); % Minutes from 15.00
    
    t_string = [repmat('2023-06-08 15:', [nT, 1]), num2str(minut_vec', '%02.0f'),...
        repmat(':00', [nT, 1])];
    t_vec_ini = datenum(t_string, 'yyyy-mm-dd HH:MM:SS');

    %% Read airspace sectors

    % Load the main sector (ACC) and adjacent sectors (ACCs)

    airspace_pgons = icas_function_get_airspace_polygons();

    % TODO: split into main and adjacent ?


    %% Read  weather data

    weather_polygons = icas_function_get_weather_data();

    %% Flight plans

    % TODO from here
    
    %FP_folder = './code_input/flight_plan_data';
    %AC = function_read_FP(FP_folder); 

    %% ASCR Available sector capacity ratio

    nk = numel(sector_names);
    ASCR = cell(nk, 1);

    k_vec_aux = 1:nk;

    for k = 1:nk % For each sector
    %for k = 1 % For one sector
    
        sectors_t = 2; % initialize variable sectors_t
    
        ASCR{k} = nan(nT, nM);

        for t = 1:nT % For each time
        %for t = 1 % Only for the first time (13:00 or 9:00)
    
            if sector_time(t,k) % if sector k exists at time t
        
                if ~isequal(sector_time(t,:), sectors_t) 
                    sectors_t = sector_time(t,:);
                    sector_data_t = sector_data(sectors_t);
                    k_new = find(k_vec_aux(sectors_t)==k);
                    [sector_ab, a_band, flows_j] = function_flows_sector_k(k_new, sector_data_t, adjacent_sectors);
                    [p_in, p_out, AC_in] = function_p_in_out(AC, sector_ab, a_band);
                end
        
                if t == 1 % initialie Wij_m
                    Wij_m1 = function_Wij_ini(flows_j);
                end
        
                [~, Wij, total_ac] = function_Wj(t_vec_ini(t), t_vec_fin(t), p_in, p_out, a_band, flows_j, AC_in);
        
                if total_ac == 0 % If there are no aircraft in the sector, use the previous weights
                    Wij = Wij_m1;
                end
        
                for m = 1:nM % For each member
                %for m = 1 % For one member
                    disp([k,t,m])
                    ASCR_k = function_ASCR_k(sector_ab, flows_j, weather_polygons{t,m}, Wij, a_band);
                    ASCR{k}(t,m) = ASCR_k;
                end
        
                Wij_m1 = Wij;
        
            end
        end

    end

end
