function AFCR_i = function_AFCR_j_onlyweigths(t, weather_polygons, sector_ab, flows_j, Wij, a_band)

disp("function_AFCR_j_onlyweigths")
% [sector_ab, a_band, flows_j] = function_flows_sector_j(k, main_sectors, adjacent_sectors);

% All flows j in sector k:
aux = [flows_j{:}]; aux2 = vertcat(aux.triplet); all_triplets = unique(aux2, 'rows');
% [ntri, ~] = size(all_triplets); % number of triplets

% [nab,ntri] = size(Wij); 

AFCR_i = nan(size(Wij));


vec_ab = find(sum(Wij, 2)>0);
for i = vec_ab'
    
    h = a_band(i);
    vec_trip = find(Wij(i,:)>0);
    trip_j = vertcat(flows_j{i}.triplet);
    
    if numel(weather_polygons)>0
        
        aux = [weather_polygons{:}];
        w_poly = vertcat(aux.pgon);
        index = false(1,numel(weather_polygons));
        for nw = 1:numel(weather_polygons) % If CTH+5000ft<h, weather hazard not considered
                index(nw) = (weather_polygons{nw}.CTH*3.28084+5000)<(h*100);
        end
        w_poly(index) = [];
    
        for j = vec_trip
            [~, j1] = ismember(all_triplets(j,:), trip_j, 'rows');

            Wmincut = function_Wmincut(sector_ab{i}, flows_j{i}(j1).T, flows_j{i}(j1).B, w_poly); % Mincut at altitude band i with weather areas


            if h==315 % draw only for one altitude
                function_Wmincut_draw(t, sector_ab{i}, flows_j{i}(j1).T, flows_j{i}(j1).B, w_poly);
            end
%             Omincut = function_Omincut(sector_ab{i}, flows_j{i}(j1).T, flows_j{i}(j1).B); % Mincut at altitude band i without weather areas
            Omincut = flows_j{i}(j1).Omincut;
            
            AFCR_i(i,j) = Wmincut./Omincut;
        end
    else
        for j = vec_trip
            AFCR_i(i,j) = 1;
        end
    end
end
disp("function_AFCR_j_onlyweigths end")
end