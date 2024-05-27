% Function to read flight plan data.

function AC = icas_function_read_FP(flight_plans_folder)

%na = numel(ac_files); % number of AC

%filename = "20230608_m1.so6";
filename = "temp.so6";
full_filename = strcat(flight_plans_folder, filename);
val = readtable(full_filename, 'Delimiter', ' ', 'FileType', 'text',...
    'Format', '%s %s %s %s %s %s %d %d %s %s %s %s %f %f %f %f %s %s %s %s');

val.Properties.VariableNames = {'segmentId', ...
        'origin', 'destination', 'aircraftType', 'beginTime', 'endTime', ...
        'beginFL', 'endFL', 'status', 'callsign', ...
        'beginDate', 'endDate', 'beginLat', 'beginLon', 'endLat', 'endLon', ...
        'segmentLength', 'segmentParityColor', 'beginTimestamp', 'endTimestamp'};

columnsToKeep = [5,6,7,8,10,11,12,13,14,15,16];
flights = val(:, columnsToKeep);

% convert minutes decimals to degrees decimals
% ignore incorrect negative values because it is not our area of interest
flights{:,"beginLat"} = flights{:,"beginLat"} / 60;
flights{:,"beginLon"} = flights{:,"beginLon"} / 60;
flights{:,"endLat"} = flights{:,"endLat"} / 60;
flights{:,"endLon"} = flights{:,"endLon"} / 60;

callsigns_list = flights.callsign;
callsigns_list = unique(callsigns_list);

na = length(callsigns_list);

AC(na) = struct;


max_FL = 0;
for ia = 1:na

    fprintf("Reading %d flight plans, flight plan: %d\n", na, ia);

    callsign = callsigns_list{ia};

    matching_rows = strcmp(flights.callsign, callsign);
    flight_data = flights(matching_rows, :);
           
    nWP = height(flight_data); % Number of segments
    
    AC(ia).WP = zeros(nWP, 4); % [lon, tal, h, time]
   
    for i = 1:nWP
      
        WP_date_time = flight_data.beginDate{i} + " " + flight_data.beginTime{i};

        AC(ia).WP(i,1) = datenum(WP_date_time,'yymmdd HHMMSS');
        AC(ia).WP(i,2) = flight_data.beginLon(i);
        AC(ia).WP(i,3) = flight_data.beginLat(i);
        AC(ia).WP(i,4) = flight_data.beginFL(i);

        if AC(ia).WP(i,4) > max_FL
            max_FL = AC(ia).WP(i,4)
        end
        
    end

    %WP_date_time = flight_data.endDate{nWP} + " " + flight_data.endTime{nWP};

    %AC(ia).WP(nWP+1,1) = datenum(WP_date_time,'yymmdd HHMMSS');
    %AC(ia).WP(nWP+1,2) = flight_data.endLon(nWP);
    %AC(ia).WP(nWP+1,3) = flight_data.endLat(nWP);
    %AC(ia).WP(nWP+1,4) = flight_data.endFL(nWP);

    %fprintf("Maximum FL among all flights: %d\n", max_FL)
    
end

end