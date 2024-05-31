clear; clc;
warning('off');
tic

global plot_mincut
plot_mincut = true;
%plot_mincut = false;

%[ASCR, sector_names, sector_time, sector_data] = icas_function_main();
ASCR = icas_function_main();
toc