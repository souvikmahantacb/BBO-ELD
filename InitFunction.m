function [Population,OPTIONS] = InitFunction(OPTIONS)
    popindex=1;
    pgaddition=(OPTIONS.pgmax+OPTIONS.pgmin);
    pgmedian=(OPTIONS.pgmax+OPTIONS.pgmin)/2;
for cy = 1 : 10000000000000*OPTIONS.popsize
    
for j=1:OPTIONS.numVar-1
   	r=rand(1);
    pg=OPTIONS.pgmin(j)+r*(OPTIONS.pgmax(j)-OPTIONS.pgmin(j));
    Population(popindex).chrom(j) = pg;
end
%%%%%%%%%%%%%%%%% OBBO %%%%%%%%%%%%%%
% for j=1:OPTIONS.numVar-1
% pgo=(OPTIONS.pgmax(j)+OPTIONS.pgmin(j)-Population(popindex).chrom(j));
% Population2(popindex).chrom(j) =pgo;
% Population(popindex).chrom(j)=Population2(popindex).chrom(j);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QOBBO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for j=1:OPTIONS.numVar-1
%     if Population(popindex).chrom(j)> pgmedian(j)
%     qpgo=pgmedian(j)+rand*(Population(popindex).chrom(j)-pgmedian(j));
%     else
%     qpgo=Population(popindex).chrom(j)+rand*(pgmedian(j)-Population(popindex).chrom(j));
%     Population4(popindex).chrom(j) =qpgo;
%     Population(popindex).chrom(j)=Population4(popindex).chrom(j);
%     end
% end
    Population.chrom;
    %pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%                  Population(popindex).chrom
%  pause
                 popindex=popindex+1;
                 if popindex>OPTIONS.popsize
                     break;
                 end
             end
         end
     end 
%     OPTIONS.pgmax 
%     Population(popindex).chrom
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OPTIONS.OrderDependent = false;
 return;
