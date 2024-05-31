function [sector_ab, a_band, flows_j] = icas_function_flows_sector_k(main_acc, adjacent_accs)

disp("icas_function_flows_sector_k")
nSec = 1; % Number of main sectors
nAdj = 7; % Number of adjacent sectors

% Altitude bands within the sector
a_band = 315:10:1000; a_band = a_band';
[nab, ~] = size(a_band); % number of altitude bands

% Polygons, triplets and edges at different altitude bands
sector_ab = cell(nab, 1);
flows_j = cell(nab, 1);

for i = 1:nab

    % Get main and adjacent ACCs for the current altitude
    [main_sector, adjacent_sectors] = get_sectors_at_altitude( ...
        main_acc, adjacent_accs, a_band(i));

    [sectors_pgon, adj_pgon] = icas_function_sector_adjacent_pgon_a_band(main_sector,...
        adjacent_sectors, a_band(i));
    %sector_ab{i} = sectors_pgon{k};
    sector_ab{i} = sectors_pgon{1};
    
    %edge_data = zeros(length(sectors_pgon{k}.Vertices(:,1)), nSec+nAdj);

    edge_data = zeros(length(sectors_pgon{1}.Vertices(:,1)), nSec+nAdj);
    
    % Find edges between adjacent sectors in that altitude band (the edges may
    % change with the altitude) 
    
    for n = 1:nAdj
        in = inpolygon(sectors_pgon{1}.Vertices(:,1),sectors_pgon{1}.Vertices(:,2)...
            ,adj_pgon{n}.Vertices(:,1),adj_pgon{n}.Vertices(:,2));
        edge_data(:,n+nSec) = in;  
    end
    
    % All possible triplets for sector j at altitude band i
    all_edges = 1:(nSec+ nAdj);
    all_edges = all_edges(any(edge_data));

    if length(all_edges)>1
        ntriplets = nchoosek(length(all_edges),2);
        triplets = nchoosek(all_edges,2);
        
        flows_j{i}(2*ntriplets) = struct; 
        
        for n = 1:ntriplets
            
            % Sector triplet
            flows_j{i}(n).triplet = [1, triplets(n,:)];
           
            % Edges
            % - Source S
            flows_j{i}(n).S = create_edge_fun(sectors_pgon{1}, logical(edge_data(:,triplets(n,1))));
            % - Destination
            flows_j{i}(n).D = create_edge_fun(sectors_pgon{1}, logical(edge_data(:,triplets(n,2))));
            % - Top and bottom
            [T_index, B_index] = top_bottom_index_fun(edge_data(:,triplets(n,1)), edge_data(:,triplets(n,2)));
            flows_j{i}(n).T = create_edge_fun(sectors_pgon{1}, T_index);
            flows_j{i}(n).B = create_edge_fun(sectors_pgon{1}, B_index);
            
            flows_j{i}(n).Omincut = function_Omincut(sectors_pgon{1}, flows_j{i}(n).T, flows_j{i}(n).B); % Mincut at altitude band i without weather areas

            
            % Inverse flows
            flows_j{i}(n+ntriplets).triplet = [1, fliplr(triplets(n,:))];
            
            flows_j{i}(n+ntriplets).S = flows_j{i}(n).D;
            flows_j{i}(n+ntriplets).D = flows_j{i}(n).S;
            flows_j{i}(n+ntriplets).T = flows_j{i}(n).B;
            flows_j{i}(n+ntriplets).B = flows_j{i}(n).T;
            
            flows_j{i}(n+ntriplets).Omincut = flows_j{i}(n).Omincut;
                
        end
    end

end
disp("icas_function_flows_sector_k end")
end

function edge = create_edge_fun(sector, index)

    if index(1)&&index(end)
       n = length(index);
       if n>1
           index2 = [find(diff(index)==1):n, 1:find(diff(index)==-1)];
           edge = [sector.Vertices(index2,1), sector.Vertices(index2,2)];
       else 
           edge = [sector.Vertices(index,1), sector.Vertices(index,2)];
       end
    else
        edge = [sector.Vertices(index,1), sector.Vertices(index,2)];
    end

end

function [T_index, B_index] = top_bottom_index_fun(indexS, indexD)

n = length(indexS);

T_index = false(size(indexS));
B_index = false(size(indexS));

ini_T = find(diff(indexS)==-1);
if isempty(ini_T)
    ini_T = n;
end
fin_T = find(diff(indexD)==1);
if isempty(fin_T)
    fin_T = 1;
else
    fin_T = fin_T+1;
end

if ini_T<=fin_T
    T_index(ini_T:fin_T) = true; 
else
    T_index([1:fin_T, ini_T:n]) = true;
end

ini_B = find(diff(indexD)==-1);
if isempty(ini_B)
    ini_B = n;
end
fin_B = find(diff(indexS)==1);
if isempty(fin_B)
    fin_B = 1;
else
    fin_B = fin_B+1;
end

if ini_B<=fin_B
    B_index(ini_B:fin_B) = true; 
else
    B_index([1:fin_B, ini_B:end]) = true;
end

end

function [main_sector, adjacent_sectors] = get_sectors_at_altitude( ...
    main_acc, adjacent_accs, a_band_i)
    
flight_levels = [315 325 335 345 355 365 375 385 1000];
num_FL = length(flight_levels);

for i=1:num_FL-1
    if (a_band_i>=flight_levels(i)) && (a_band_i<flight_levels(i+1))
        main_sector = main_acc{i}; % struct
        adjacent_sectors = adjacent_accs(i,:); % 1x4 cell array containing structs
    end
end

end
