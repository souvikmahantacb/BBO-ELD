function [Population] = CostFunction(OPTIONS, Population)
      popsize = OPTIONS.popsize;
for popindex = 1 : popsize
 
         Population(popindex).cost1 = 0;
    for i=1:OPTIONS.numVar
            Population(popindex).cost1=Population(popindex).cost1+OPTIONS.alpha(i)*(Population(popindex).chrom(i))^2+OPTIONS.beta(i)*(Population(popindex).chrom(i))+OPTIONS.gama(i);
% Population(popindex).cost1=Population(popindex).cost1+abs(OPTIONS.sin1(i)*sin(OPTIONS.sin2(i)*(-Population(popindex).chrom(i)+OPTIONS.pgmin(i))));
    end
%    Population(popindex).emission = 0;
%     for j=1:OPTIONS.numVar
% 			Population(popindex).emission=Population(popindex).emission+OPTIONS.alphas(j)*(Population(popindex).chrom(j))^2+OPTIONS.betas(j)*(Population(popindex).chrom(j))+OPTIONS.gamas(j);
%             %+OPTIONS.mus(j)*exp(OPTIONS.deltas(j)*Population(popindex).chrom(j));     %Bidding case
%             
%      end
 % cost(i)=0;
  % Population(popindex).cost=Population(popindex).cost1;
   
%Population(popindex).cost=Population(popindex).emission;
 %h= 0.21646;
 %h= 18.8646;
% 
 Population(popindex).cost=Population(popindex).cost1;%+h*Population(popindex).emission;
end

return