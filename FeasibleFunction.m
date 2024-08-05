
function [Population] = FeasibleFunction(OPTIONS, Population)
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
