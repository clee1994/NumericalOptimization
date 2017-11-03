function [ weights ] = min_var( exp_ret, mu, sigma )
    %this funds the min var portfolio for a expected return
    %forumlas come from Gracia 2016
    %includes matlabs inversion, which might be not very accurate
    %two fund spannung theorem is employed here
    S_inv = inv(sigma); %taking inverse
    pi_mu = (sigma\exp_ret)/sum(sigma\exp_ret); 
    pi_1 = (sum(S_inv) / sum(sum(S_inv)))'; 
    lambda_demoninator = (exp_ret'*pi_mu) - (exp_ret'*pi_1); 
    ll = (mu - (exp_ret'*pi_1))/lambda_demoninator; %4.12
    weights = (pi_mu * ll + pi_1 * (1-ll));


end

