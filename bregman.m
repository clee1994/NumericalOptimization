function [  w_ret, mu_p, var_p, time ] =  bregman( mu, gamma, graphs_on , cdts, tol,rho,lambda)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    %return series
    logr = cdts;

    %in case we didnt specify, if we want to see graphs or not, just dont
    %shows them
    if nargin < 3
       graphs_on = false; 
    end
    
    tic
    
    iter = length(mu);

    %bootstrap parameter
    K = 1000;
    
    %bootstrapping to determine parameter uncertainty
    [bootstat] = bootstrp(K,@mean,logr);
    mu_err = abs(bootstat - repmat(mean(bootstat),[K, 1]));
    beta = mean(mu_err);
    
    %same for variance
    [bootstat] = bootstrp(K,@cov,logr);
    alpha = abs(bootstat - repmat(mean(bootstat),[K, 1]));
    alpha = mean(alpha);
    alpha = reshape(alpha, [],iter);

    Da = alpha;
    B = eye(iter);
    B(logical(eye(size(B)))) = beta;

    %variable initlisation
    k = 2; 
    b(:,k) = zeros(iter,1); 
    w(:,k)= zeros(iter,1);
    w(:,1)= -ones(iter,1); 
    d(:,k) = zeros(iter,1);
    R = rho*gamma + Da; %making gamma easier invertible
    
    %run the bregman
    while norm(w(:,k) - w(:,k-1)) > tol
       w(:,k+1) = (R+B'*B)\(mu'+(lambda*B'*(d(:,k)-b(:,k))));
       for j = 1:iter
           d(j,k+1) = shrinkage_operator((beta(j)*w(j,k+1)+b(j,k)),(1/lambda));
           b(j,k+1) = b(j,k) + beta(j)*w(j,k+1)-d(j,k+1);
       end
       k = k + 1;
    end

    %final weight vector
    w_ret=w(:,k);
    
    time = toc;
    %whats my mean and variance of the portfolio
    mu_p = w(:,k)'*mu';
    var_p = w(:,k)'*gamma*w(:,k);

    %plot some results 
    if graphs_on
        figure();
        scatter(diag(gamma),mu');
        hold on
        scatter(var_p, mu_p, 'r');
        title('Split Bregman');
        ylabel('Mean');
        xlabel('Variance');
        saveas(gca, 'Figures/breg','png')
        figure,bar(w_ret),xlim([0,length(w_ret)+1]);
        title('Weights Split Bregman')
        saveas(gca, 'Figures/bregbar','png')
    end


end

