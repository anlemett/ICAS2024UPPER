clear; clc;
warning('off');
tic

global nowcast
%nowcast = true; % margins will not be added
nowcast = false; % margins will be added

%[ASCR, sector_names, sector_time, sector_data] = icas_function_main();
ASCR = icas_function_main();
toc