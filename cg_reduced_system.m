function [  w, mu_p, var_p, time ] = cg_reduced_system( mu, gamma, graphs_on, null_sep, tol )
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
    %range and Null space
    [Y,Z] = null_range(gamma, null_sep, false );
    xz = repmat(1/size(Z,2), [size(Z,2),1]);
    Wzz = Z'*eye(size(Z,1))*Z;
    xy = (A*Y)\b;
    cz = (Z'*gamma*Y*xy)+(Z'*c);
    rz = Z' *gamma* Z * xz + cz;
    gz = Wzz\rz;
    dz = -gz;


    %algorithm as described in Nocedal and the pdf
    while (rz'*(Wzz\rz)) > tol 
        alpha = (rz' * gz) / (dz'*Z'*gamma*Z*dz);
        xz = xz + alpha *dz;
        rzp = rz + alpha * Z' * gamma * Z * dz;
        gzp = Wzz\rzp;
        beta = (rzp' * gzp)/(rz'*gz);
        dz = -gzp + beta*dz;
        gz = gzp;
        rz = rzp;
    end

    %final weight vector
    w = Z*xz + Y*xy;
    
    time = toc;
    mu_p =w'*mu';
    var_p =w'*gamma*w;

    %plots of results
    if graphs_on
        figure();
        scatter(diag(gamma),mu');
        hold on
        scatter(var_p,mu_p,  'r')
        title('CG applied to reduced system')
        ylabel('Mean'),xlabel('Variance');
        saveas(gca, 'Figures/cgredmv','png')
        figure,bar(w),xlim([0,length(w)+1]);
        title('Weights CG applied to reduced system')
        saveas(gca, 'Figures/cgredbar','png')
    end


end

