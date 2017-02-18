%*************************************************************************%
% Title: Travelling Salesman Problem Using Genetic ALgorithm
% Description: 
% Version: 1.0
% Author: Group 4
% Year & Section : BSCS 3-3
% [13 22; 20 13; 21 17 ;24 12 ; 25 16;32 19;        33 12 ;35 15; 39 20; 43 8];
%*************************************************************************%

% function TSP :performs genetic algorithm to solve travelling salesman
% problem.
function TSP = tsp_genetic_algo()
    
%*****************************@Initialization*****************************%

    xy_coord           = [7 27; 11 20; 13 29; 17 10; 21 18; 29 21; 33 11; 37 3; 43 26; 43 9]; % City Coordinates
    POP_SIZE           = 60; % population size
    NUM_OF_GENERATIONS = 100 ; % number of generations
    [num_of_cities, ~] = size(xy_coord); % get the number of cities
    global_min         = Inf; % positive integers
    dist_matrix        = []; %point a to b distances/cost
    pop                = zeros(POP_SIZE, 10); % N by N matrix of zeroes
    total_dist         = zeros(1,POP_SIZE); %total distance of the best path
    distHistory = zeros(1,NUM_OF_GENERATIONS);  % initialize an N by N matrix of zeroes, will be used for ploting
    tmp_pop = zeros(4,num_of_cities); % initialize an N by N matrix of zeroes
    new_pop = zeros(POP_SIZE,num_of_cities); % initialize an N by N matrix of zeroes
    current_gen = 0; %number of generations to reach global minima
    gen = 0; %iterator for generations
   
    
%****************************@Code Body***********************************%    
    POP_SIZE = 4*ceil(POP_SIZE/4);
    nPoints = size(xy_coord,1);
    display(nPoints);
    a = meshgrid(1:nPoints);
    dist_matrix = reshape(sqrt(sum((xy_coord(a,:)-xy_coord(a',:)).^2,2))...
    ,nPoints,nPoints);   %compute distances from each points
    display(dist_matrix);
  
    
    % STEP 1: Generate a population using permutation encoding
    
    pop(1,:) = (1:num_of_cities);
    for k = 2:POP_SIZE
        pop(k,:) = randperm(num_of_cities);
    end
    
    % Step 2: FITNESS EVALUATION: After initializing the population, fitness value for each of
    % the individual is calculated. The fitness values are generated using
    % a fitness function. The function provides largest and smallest values
    % for each of the individuals. If the individual has a larger fitness
    % value then result will be a better solution but if a smaller value is
    % obtained then solution obtained will not be better.

    figure('Name','Best Solution');
    for gen = 1:NUM_OF_GENERATIONS
        % Evaluate Each Population Member (Calculate Total Distance)
        for p = 1:POP_SIZE
            tmp_tot_dist = dist_matrix(pop(p,num_of_cities), pop(p, 1)); % temporary variable for total dstance distance
            for k = 2:num_of_cities
                tmp_tot_dist = tmp_tot_dist + dist_matrix(pop(p,k-1),pop(p,k)); % summation of distances of pop(k)
            end
            total_dist(p) = tmp_tot_dist; %select as elite
        end
        
        % Find the least cost if the new computed total distance is less
        % than the latter distance (global_min) then the new global minimum
        % distance is the min_dist
        
        [min_dist,index] = min(total_dist);
        distHistory(gen) = min_dist;
        if min_dist < global_min
            global_min = min_dist;
            optimal_route = pop(index, :);
            path = optimal_route([1 : num_of_cities 1]);
            % Plot the Best Route
            %pause(1);
            plot(xy_coord(path, 1), xy_coord(path, 2), 'r.-'); 
            title(sprintf('Distance = %1.4f',min_dist));     
            drawnow;
            current_gen = gen;
        end
        
        % STEP 4: CROSSOVER AND MUTATION: 
        randomOrder = randperm(POP_SIZE);
        for p = 4:4:POP_SIZE
            paths = pop(randomOrder(p-3:p),:); % select parent 1 from the random population
            dists = total_dist(randomOrder(p-3:p)); %select parent 2 from the matrix of ELITE population
            [~,idx] = min(dists);       % get the index of the path with the min distance from the elite population
            possible_route = paths(idx,:);  %get the values/permutation encoding
            routeInsertionPoints = sort(ceil(num_of_cities*rand(1,2)));
            I = routeInsertionPoints(1); 
            J = routeInsertionPoints(2);
       
            for k = 1:4                         % Mutate and create a 4 new set of population
                tmp_pop(k,:) = possible_route;  
                display(tmp_pop);
                display(size(tmp_pop));
                
                switch k
                    case 2                   
                        tmp_pop(k,I:J) = tmp_pop(k,J:-1:I); % mutate by flipping
                    case 3 
                        tmp_pop(k,[I J]) = tmp_pop(k,[J I]); %mutate by swapping
                    case 4 
                        tmp_pop(k,I:J) = tmp_pop(k,[I+1:J I]);   % mutate by sliding 
                    otherwise % default
                end
            end
            new_pop(p-3:p,:) = tmp_pop;
            %end of generation
        end
        pop = new_pop; % new generation  
        
    end
    %*********************************************************************%  
    
    %**************************@Plotting**********************************%
    figure('Name','RESULTS','Numbertitle','off');
    title('Best Solution History');
    plot(distHistory,'b','LineWidth',3);
    set(gca,'XLim',[0 NUM_OF_GENERATIONS+1],'YLim',[0 1.1*max([1 distHistory])]);
    xlabel('Generation');
    ylabel('Distance');
   
    %*************************@ENDPlotting********************************%
    
    %**************************@Display***********************************%
     display(optimal_route);
     display(min_dist);
     display(current_gen);
    %*********************************************************************%
end    


