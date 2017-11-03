function [ w, mu_p, var_p, time ] = old_school( mu, gamma, graphs_on )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 3
       graphs_on = false; 
    end
    len = size(mu,2);
    
    tic
   
    %some formula for finding the mu corresponding to the min var
    numerator = mu*(gamma\ones(len,1));
    denominator =  ones(1,len)*(gamma\ones(len,1));
    mu_p = numerator/denominator;

    %finding the weights composition for min var
    w = min_var(mu', mu_p, gamma);
    var_p = w'*gamma*w;
    
    time = toc;

    %plot results
    if graphs_on
        figure();
        scatter(diag(gamma),mu'); hold on
        scatter(var_p,mu_p,'r'); hold off
        title('Closed form solution using Matlabs mldivide')
        ylabel('Mean'),xlabel('Variance');
        saveas(gca, 'Figures/oldschool','png') 
        figure,bar(w),xlim([0,length(w)+1]);
        saveas(gca, 'Figures/oldschoolbar','png') 
    end


end

