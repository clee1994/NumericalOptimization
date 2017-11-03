function [  w, mu_p, var_p, time ] = null_space( mu, gamma, graphs_on, null_sep  )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


    if nargin < 3
       graphs_on = false; 
    end
    
    tic
    
    A = ones(1,size(gamma,1));
    b = 1;
    c = zeros(size(gamma,1),1);
    
    x = repmat(1/length(mu),[length(mu) 1]);
    h = A*x - b;
    g = c + gamma*x;


    [Y,Z] = null_range(gamma,null_sep, false);

    %pY and pX
    pY=(A*Y)\-h;
    pZ=(Z'*gamma*Z)\(-Z'*gamma*Y*pY-Z'*g);

    if isnan(pZ); pZ=0; end;
    %solution
    p=(Y*pY)+(Z*pZ);
    w=p+x;
    
    time = toc;

    %portfolio
    mu_p =w'*mu';
    var_p =w'*gamma*w;
    

    if graphs_on
        figure();
        scatter(diag(gamma),mu');
        hold on
        scatter(var_p,mu_p,  'r')
        title('KKT System direct solution Null Space method')
        ylabel('Mean');
        xlabel('Variance');
        saveas(gca, 'Figures/null','png') 
        figure,bar(w),xlim([0,length(w)+1]);
        title('Weights Null Space method')
        saveas(gca, 'Figures/nullbar','png') 
    end
    


end

