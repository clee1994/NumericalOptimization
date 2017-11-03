function [ null_sep, tol_cgred, tol_cgproj, lambda_fin, rho_fin ] = var_calibration(mu, gamma, data_ret)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    %find the best parameters for all the methods given the data

    %see what might be a good nullspace seperator
    sep_test = 0.0001;
    for i = 1:1000
        sep_test = sep_test + 0.0001;
        [w(:,i), mu_p(i), var_p(i), time(i)] = null_space(mu,gamma,false,sep_test);
    end
    figure();
    plot(var_p);
    title('Min. Variance w.r.t. to different Nullspaces')
    saveas(gca, 'Figures/cal_nullsep','png');
    [~,ind] = min(var_p);
    null_sep = ind * 0.0001;


    %tol test cg_reduced 
    tol_test(1) = 0.1;
    for i = 1:30

        [w(:,i), mu_p(i), var_pt(i), time(i)] = cg_reduced_system(mu,gamma,false, null_sep, tol_test(i));
        tol_test(i+1)= tol_test(i)*0.1;
    end
    figure();
    plot(var_pt);
    title('Min. Var. CG reduced w.r.t. to different tolerance');
    saveas(gca, 'Figures/tol_cgred','png');
    [~,ind] = min(round(var_pt,10)); %this rounding assumption....
    tol_cgred = 0.1^ind;
    
    
    %tol test cg_proj 
    tol_test(1) = 0.1;
    for i = 1:30
        [w(:,i), mu_p(i), var_pt(i), time(i)] = cg_projected_test(mu,gamma,false, tol_test(i));
        tol_test(i+1)= tol_test(i)*0.1;
    end
    figure();
    plot(var_pt);
    title('Min. Var. projected CG w.r.t. to different tolerance');
    saveas(gca, 'Figures/tol_cgproj','png');
    [~,ind] = min(round(var_pt,10)); %this rounding assumption....
    tol_cgproj = 0.1^ind;
    
    
    
    %Challenge
    
    %lambda
    lambda = 0.1;
    for i = 1:100
        lambda = lambda + 0.1;
        [w, mu_p, var_pb(i), time] = bregman(mu,gamma,false,data_ret, 10e-5,8,lambda);
    end
    figure();
    plot(var_pb);
    title('Min. Variance Split Bregman w.r.t. to different lambdas')
    saveas(gca, 'Figures/sb_lambda','png');
    [~,ind] = min(var_pb);
    lambda_fin = ind * 0.1;
    
    %rho
    rho = 0.1;
    for i = 1:175
        rho = rho + 0.1;
        [w, mu_p, var_pb2(i), time] = bregman(mu,gamma,false,data_ret, 10e-5,rho,lambda_fin);
    end
    figure();
    plot(var_pb2);
    title('Min. Variance Split Bregman w.r.t. to different rhos')
    saveas(gca, 'Figures/rho_lambda','png');
    [~,ind] = min(round(var_pb2,5)); %rounding...
    rho_fin = ind * 0.1;
    

end

