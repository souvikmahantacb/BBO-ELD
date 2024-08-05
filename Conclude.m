function Conclude(DisplayFlag, OPTIONS, Population, nLegal, MinCost)

% Output results of population-based optimization algorithm.

if DisplayFlag
    % Count the number of duplicates
%     NumDups = 0;
%     for i = 1 : OPTIONS.popsize
%          Chrom1 = sort(Population(i).chrom);
% %Chrom1 = (Population(i).chrom);
%         for j = i+1 : OPTIONS.popsize
%              Chrom2 = sort(Population(j).chrom);
% % Chrom2 = sort(Population(j).chrom);
%             if isequal(Chrom1, Chrom2)
%                 NumDups = NumDups + 1;
%             end
%         end
%     end 
%     Population = FeasibleFunction(OPTIONS, Population);
%     disp([num2str(NumDups), ' duplicates in final population.']);
%     disp([num2str(nLegal), ' legal individuals in final population.']);
    % Display the best solution
%     Chrom = sort(Population(1).chrom);
     
    Chrom = (Population(1).chrom);
    XX=sum(Chrom);
    cost1 = 0;
    for i=1:OPTIONS.numVar
         cost1=cost1+OPTIONS.alpha(i)*(Chrom(i))^2+OPTIONS.beta(i)*(Chrom(i))+OPTIONS.gama(i);
       % cost1=cost1+abs(OPTIONS.sin1(i)*sin(OPTIONS.sin2(i)*(-Chrom(i)+OPTIONS.pgmin(i))));          
    end
          
%     pl=Chrom*OPTIONS.bcoefficient;
%     loss=pl*(Chrom)';
    sumpgn=0;
    for j=1:OPTIONS.numVar
        sumpgn=sumpgn+Chrom(j);
    end
%     mismatch=abs(sumpgn-loss-OPTIONS.pdemand);
    misfitness=cost1;
%     misfitness=cost1+(loss)*100+ (mismatch)*1000;
 disp(['Best chromosome = ', num2str(Chrom)]);
 disp(['Best cost = ', num2str(cost1)]);
 disp(['Best generation = ', num2str(XX)]);
%  disp(['Best loss = ', num2str(loss)]);
 disp(['Best misfitness = ', num2str(misfitness)]);
    % Plot some results
%     close all;
figure(1);
    plot([0:OPTIONS.Maxgen], MinCost, 'r');
    title('Misfitness Curve')
    xlabel('Generation');
    ylabel('Misfitness');
    hold on;
end
return;