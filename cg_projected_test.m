function [  w, mu_p, var_p, time ] = cg_projected_test( mu, gamma, graphs_on, tol )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


    if nargin < 3
       graphs_on = false; 
    end
    
    tic
   
    %initlization
    A = ones(1,size(gamma,1));
    b = 1;
    c = zeros(size(gamma,1),1);

    
    %Null space
    %[~,Z] = null_range(gamma, null_sep, false);

    %initlization
    x = A'*((A*A')\b);
    r = gamma*x +c; 
    %there are various choices of the preconditioner
    P = eye(size(gamma,1)) - A'*((A*A')\A);
    g = P*r;
    d = -g;


    %the algorithm as explained in the pdf
    while (r'*g) > tol 
        alpha = (r' * g) / (d'*gamma*d);
        x = x + alpha*d;
        rp = r + alpha * gamma * d;
        vp = (A*A')\(A*rp);
        gp = rp - A'*vp;
        beta = (rp' * gp)/(r'*g);
        d = -gp + beta*d;
        g = gp;
        r = rp;
    end

    time = toc;
    w = x;

    mu_p =w'*mu';
    var_p =w'*gamma*w;
    
   
    %plot some results
    if graphs_on
        figure();
        scatter(diag(gamma),mu');
        hold on
        scatter(var_p,mu_p,  'r')
        title('Projected CG')
        ylabel('Mean');
        xlabel('Variance');
        saveas(gca, 'Figures/cgproj','png')
        figure,bar(w),xlim([0,length(w)+1]);
        title('Weights Projected CG')
        saveas(gca, 'Figures/cgprojbar','png')
    end


end

