%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICAS24: airspace capacity 
% Anastasia Lemetti
% MATLAB version: MATLAB R2024a
% 
% Read weather data and store it as polygons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function weather_polygons = icas_function_get_weather_data()

T = 10; % Number of weather time intervals from 15.00 to 17.30 (15 minutes intervals)

times = {'1456', '1511', '1526', '1541', '1556', '1611', '1626', '1641', '1656', '1711'};

weather_polygons = cell(T, 1);

for t = 1:T
    
    day = 'DLR_WCB_T_EUR_20230608';
    daytime = strcat(day, times{t});
    filename = strcat(daytime, '.xml');

    weather_filename = fullfile('.', 'code_input', 'meteo_data', filename);

    S = readstruct(weather_filename);
    
    polygons_list = getWeatherPolygons(S, t); 

    num_poly = length(polygons_list); % Number of polygons
    weather_polygons{t,1} = num_poly;
    
    if num_poly>0
        for i = 1:num_poly

            lat_lon = str2num(polygons_list(i));

            lat = lat_lon(1:2:end);
            lon = lat_lon(2:2:end);

            weather_polygons{t,i+1} = polyshape(lon, lat);
        end
    end
end

end


function polygons = getWeatherPolygons(strct, time_interval)

% Define the format of the time strings
timeFormat = 'HH:mm';

start_times = {'15:03', '15:18', '15:33', '15:48', '16:03', '16:18', '16:33', '16:48', '17:03', '17:18'};

polygons = {};

forecast_sets = [strct.wims_forecastSet];

num_of_forecasts = length(forecast_sets);

for i = 1:num_of_forecasts
        
    substrct =  forecast_sets(i);

    start_time_str = substrct.wims_ForecastSet...
                    .wims_status...
                    .wims_StatusWeatherProduct...
                    .wims_metaData...
                    .wims_MetaData...
                    .wims_validityStartTime...
                    .gml_TimeInstant...
                    .gml_timePosition;
    start_time_str = extractAfter(start_time_str, 11);
    start_time_str = extractBefore(start_time_str, 6);

    %disp(start_time);

    if strcmp(start_time_str, start_times{time_interval})

        CBs = [substrct.wims_ForecastSet...
                .wims_weatherProductList...
                .wims_WeatherProductList...
                .gml_featureMembers.wims_CB];

        num_of_CBs = length(CBs);

        for j = 1:num_of_CBs

             CB = CBs(j);

             polygon = getfield(CB.wims_geometry...
                 .gml_Polygon...
                 .gml_exterior...
                 .gml_LinearRing, 'gml_posList');
             polygons = [polygons; polygon];
        end
    end

end
end