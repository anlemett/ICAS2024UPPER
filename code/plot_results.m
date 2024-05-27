clear; clc;

x = 1:10;
ASCR = [0.872737511335934;
        1;
        0.945649556742329;
        0.943047094763133;
        0.862010414741266;
        0.719320740765052;
        0.537125994601274;
        0.592788342773918;
        0.855154714058835;
        1]

figure
plot(x, ASCR, 'LineWidth', 2)

xticklabels({'15:00-15:15', '15:15-15:30', '15:30-15:45', '15:45-16:00', '16:00-16:15', ...
    '16:15-16:30', '16:30-16:45', '16:45-17:00', '17:00-17:15', '17:15-17:30'});
xtickangle(90);
ylabel('ASCR');

min_x = 0;
max_x = 11;
min_y = 0.4;
max_y = 1.099;

xlim([min_x max_x]);
ylim([min_y max_y]);
