clear; clc;
warning('off');
tic

[ASCR, sector_names, sector_time, sector_data] = icas_function_main();
toc