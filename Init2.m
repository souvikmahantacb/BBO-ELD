function [OPTIONS, MinCost, AvgCost, InitFunction, CostFunction,FeasibleFunction,Population] = Init2(DisplayFlag, ProblemFunction, RandSeed)

% Initialize population-based optimization software.
if ~exist('DisplayFlag', 'var')
    DisplayFlag = true;
end
% WARNING: some of the optimization routines will not work if population size is odd.
OPTIONS.popsize = 50; % total population size
OPTIONS.Maxgen = 50; % generation count limit
OPTIONS.numVar = 6; % number of genes in each population member
OPTIONS.pmutate = 0.005; % mutation probability
%OPTIONS.bcoefficient = [0.000676 0.0000953 -0.0000507;0.0000953 0.000521 0.0000901;-0.0000507 0.0000901 0.000294];
OPTIONS.pdemand=1263;
OPTIONS.pgmin=[100 50 80 50 50 50];
OPTIONS.pgmax=[500 200 300 150 200 120];
OPTIONS.gama=[240 200 220 200 220 190] ;
OPTIONS.beta=[7 10 8.5 11 10.5 12];
OPTIONS.alpha=[0.00070 0.0095 0.0090 0.0090 0.0080 0.0075];
% OPTIONS.gamas=[13.8593 13.8593 40.2669 40.2669 42.8955 42.8955];
% OPTIONS.betas=[0.32767 0.32767 -0.54551 -0.54551 -0.51116 -0.51116];
% OPTIONS.alphas=[0.00419 0.00419 0.00683 0.00683 0.00461 0.00461];

OPTIONS.bcoefficient=10^-6*[140     17  15  19  26  22;
                            17      60  13  16  15  20;
                            15      13  65  17  24  19;
                            19      16  17  71  30  25;
                            26      15  24  30  69  32;
                            22      20  19  25  32  85;];
%OPTIONS.bcoefficient;

%  if ~exist('RandSeed', 'var')
%     RandSeed = round(sum(100*clock))
% end
% rand('state', RandSeed); % initialize random number generator
% if DisplayFlag
%     disp(['random # seed = ', num2str(RandSeed)]);
% end

% Get the addresses of the initialization, cost, and feasibility functions.
% [InitFunction, CostFunction, FeasibleFunction] = ProblemFunction();
[InitFunction, CostFunction,FeasibleFunction] = eld;


% Initialize the population.
[Population,OPTIONS] = InitFunction(OPTIONS);
% Make sure the population does not have duplicates. 
Population = ClearDups(OPTIONS, Population);
% Compute cost of each individual  
Population = CostFunction(OPTIONS, Population);
% Sort the population from most fit to least fit
Population = PopSort(Population);

% Compute the average cost
AverageCost = ComputeAveCost(Population);
% Display info to screen
MinCost = [Population(1).cost];
AvgCost = [AverageCost];
if DisplayFlag
    disp(['The best and mean of Generation # 0 are ', num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
end

return;