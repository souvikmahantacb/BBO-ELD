 function [InitFunction CostFunction,FeasibleFunction] = eld
 InitFunction = @eldInit;
 CostFunction=@eldCost
 FeasibleFunction = @eldFeasible;

 return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population,OPTIONS] = eldInit(OPTIONS)
    popindex=1;
for cy = 1 : 10000000000000000000000000000000000000000000000*OPTIONS.popsize
    
for j=1:OPTIONS.numVar-1
   	r=rand(1);
    pg=OPTIONS.pgmin(j)+r*(OPTIONS.pgmax(j)-OPTIONS.pgmin(j));

  
    Population(popindex).chrom(j) = pg;
end
    %***********************************************
    X=OPTIONS.bcoefficient(OPTIONS.numVar,OPTIONS.numVar);
    Y=0;
    for j=1:OPTIONS.numVar-1
        Y=Y+2*(Population(popindex).chrom(j))*OPTIONS.bcoefficient(j,OPTIONS.numVar);
    end
    Y=Y-1;
    Z=0;
    for j=1:OPTIONS.numVar-1
        Z=Z+(Population(popindex).chrom(j))^2*OPTIONS.bcoefficient(j,j)-Population(popindex).chrom(j);
    end
    losssum=0;
    for ig=1:OPTIONS.numVar-1
        for h=1:OPTIONS.numVar-1
            if ig==h
                losssum=losssum;
            else
                losssum=losssum+Population(popindex).chrom(ig)*Population(popindex).chrom(h)*OPTIONS.bcoefficient(ig,h);
            end
        end
    end
     Z=Z+OPTIONS.pdemand+losssum;
     discriminator=Y*Y-4*X*Z;
     AA=(-Y-sqrt(discriminator))/(2*X);
     if (discriminator>0)
         if AA<=OPTIONS.pgmax(OPTIONS.numVar)
             if AA>=OPTIONS.pgmin(OPTIONS.numVar)
                 Population(popindex).chrom(OPTIONS.numVar)=AA;
                 popindex=popindex+1
                 if popindex>OPTIONS.popsize
                     break;
                 end
             end
         end
     end 
%     OPTIONS.pgmax 
%     Population(popindex).chrom
end
Population.chrom

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPTIONS.OrderDependent = false;
 return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [Population] = eldCost(OPTIONS, Population)
      popsize = OPTIONS.popsize;
for popindex = 1 : popsize
 
         Population(popindex).cost1 = 0;
    for i=1:OPTIONS.numVar
            Population(popindex).cost1=Population(popindex).cost1+OPTIONS.alpha(i)*(Population(popindex).chrom(i))^2+OPTIONS.beta(i)*(Population(popindex).chrom(i))+OPTIONS.gama(i);
% Population(popindex).cost1=Population(popindex).cost1+abs(OPTIONS.sin1(i)*sin(OPTIONS.sin2(i)*(-Population(popindex).chrom(i)+OPTIONS.pgmin(i))));
    end
   
%     for j=1:OPTIONS.numVar
% 			Population(popindex).emission=Population(popindex).emission+OPTIONS.alphas(j)*(Population(popindex).chrom(j))^2+OPTIONS.betas(j)*(Population(popindex).chrom(j))+OPTIONS.gamas(j);
%             %+OPTIONS.mus(j)*exp(OPTIONS.deltas(j)*Population(popindex).chrom(j));     %Bidding case
%             
%      end
  cost(i)=0;
  % Population(popindex).cost=Population(popindex).cost1;
   
%Population(popindex).cost=Population(popindex).emission;
% h= 0.21646;
 h= 18.8646;
% 
 Population(popindex).cost=Population(popindex).cost1;
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = eldFeasible(OPTIONS, Population)
    Population.chrom;

    for popindex=1 :OPTIONS.popsize
    for j=1:OPTIONS.numVar
        AA=Population(popindex).chrom(j);
         Population(popindex).chrom(j)=0;
        Population1(popindex).chrom(j)=AA;
    end
end
    popul=1;
for popindex=1 :OPTIONS.popsize
    for j=1:OPTIONS.numVar-1
        if Population1(popindex).chrom(j)>OPTIONS.pgmax(j)
            Population1(popindex).chrom(j)=OPTIONS.pgmax(j);
        elseif Population1(popindex).chrom(j)<OPTIONS.pgmin(j)
            Population1(popindex).chrom(j)=OPTIONS.pgmin(j);
        end
    end
     X=OPTIONS.bcoefficient(OPTIONS.numVar,OPTIONS.numVar);
    Y=0;
    for j=1:OPTIONS.numVar-1
        Y=Y+2*(Population1(popindex).chrom(j))*OPTIONS.bcoefficient(j,OPTIONS.numVar);
    end
    Y=Y-1;
    Z=0;
    for j=1:OPTIONS.numVar-1
        Z=Z+(Population1(popindex).chrom(j))^2*OPTIONS.bcoefficient(j,j)-Population1(popindex).chrom(j);
    end
    losssum=0;
    for ig=1:OPTIONS.numVar-1
        for h=1:OPTIONS.numVar-1
            if ig==h
                losssum=losssum;
            else
                losssum=losssum+Population1(popindex).chrom(ig)*Population1(popindex).chrom(h)*OPTIONS.bcoefficient(ig,h);
            end
        end
    end
     Z=Z+OPTIONS.pdemand+losssum;
     discriminator=Y*Y-4*X*Z;
     AA=(-Y-sqrt(discriminator))/(2*X);
     if (discriminator>0)
         if AA<=OPTIONS.pgmax(OPTIONS.numVar)&& AA>=OPTIONS.pgmin(OPTIONS.numVar)
             
                 Population1(popindex).chrom(OPTIONS.numVar)=(-Y-sqrt(discriminator))/(2*X);
                for j=1:OPTIONS.numVar
                    Population(popul).chrom(j)=Population1(popindex).chrom(j);
                end
                    
                 popul=popul+1;
         
             
         end
     end
   if popul-1<OPTIONS.popsize
       for i=popul:OPTIONS.popsize
            for j=1:OPTIONS.numVar
            Population(i).chrom(j)=Population(i-popul+1).chrom(j);
            end
       end
   end
   
end

return;
