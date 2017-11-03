clear 
close all

%Final Project Clemens Schaefer
%Numerical Optimisation
%Applied to minimum variance portfolio optimitisation

%This whole programm might take a while to run (4secs, but plotting and 
%saving might make it slower). In order not to be bored, you can go to 
%the folder boot and enjoy my finest collection of memes.

tic

%% Load Data and calibrate variables
%get DAX data to compare final perfromance
filename1 = 'Data/DAX_pure.csv';
data = csvread(filename1,1,1);
data_ret_ind = diff(log(data),1);

%reading csv, stock prices from all DAX constitutes
filename1 = 'Data/DAX.csv';
data = csvread(filename1,1,1);
data_ret = diff(log(data),1); %log returns, thus additive

%remove columns with NaN, e.g. incomplete price history
[~, c]=find(isnan(data_ret));
data_ret(:,c)=[];
dropped_firms = unique(c);

%subsample for illustration, last year (last 250 days)
data_ret_sub = data_ret((end-250):end,:);
mu = mean(data_ret_sub);
gamma = cov(data_ret_sub);


%variable calibration, e.g. find the best parameters for the optimisation
%methdos
[ null_sep, tol_cgred, tol_cgproj, lambda, rho ] = ...
    var_calibration(mu, gamma, data_ret_sub);

%analyse gamma (covariance matrix, e.g. find null and range space)
[ Y,Z ] = null_range( gamma, null_sep, true );



%% Solve the Problem with different methods
% direct solutions
[w1, mu_p1, var_p1, time1] = old_school(mu,gamma,true);

%Null space
[w2, mu_p2, var_p2, time2] = null_space(mu,gamma,true,null_sep);

% iterative solutions
[w3, mu_p3, var_p3, time3] = cg_reduced_system(mu,gamma,true, null_sep, tol_cgred);
[w4, mu_p4, var_p4, time4] = cg_projected_test(mu,gamma,true, tol_cgproj);

% challenge - Split Bregnman
[w5, mu_p5, var_p5, time5] = bregman(mu,gamma,true,data_ret_sub, 10e-15,rho,lambda);




%Put results into latex table
fileID = fopen('Tables/ilres.tex','w');
fprintf(fileID,'\\begin{tabular}{ r|llll }\n');
fprintf(fileID,'& $\\mu \\times 10^{4}$ & $\\sigma^2 \\times 10^{5}$ & used stocks & time  \\\\ \n \\hline\n');
fprintf(fileID,'Closed form& %2.4f & %2.4f & %i & %2.4f s  \\\\ \n',round(mu_p1*10000,4),round(var_p1*100000,4),size(unique(round(w1,6)),1),time1);
fprintf(fileID,'Null space& %2.4f & %2.4f & %i & %2.4f s  \\\\ \n',round(mu_p2*10000,4),round(var_p2*100000,4),size(unique(round(w2,6)),1),time2);
fprintf(fileID,'CG reduced system& %2.4f & %2.4f & %i & %2.4f s  \\\\ \n',round(mu_p3*10000,4),round(var_p3*100000,4),size(unique(round(w3,6)),1),time3);
fprintf(fileID,'CG projected& %2.4f & %2.4f & %i & %2.4f s \\\\ \n',round(mu_p4*10000,4),round(var_p4*100000,4),size(unique(round(w4,6)),1),time4);
fprintf(fileID,'Split Bregman& %2.4f & %2.4f & %i & %2.4f s \\\\ \n',round(mu_p5*10000,4),round(var_p5*100000,4),size(unique(round(w5,6)),1),time5);
fprintf(fileID,'\n \\end{tabular}');
fclose(fileID);


disp('Backtesting performance start');

%% performance test, e.g. backtesting

%monthly rebalancing
stepsize = 20;

%looping over the whole past history
for i = 251:stepsize:(size(data,1)-stepsize)
    %dataset for computations
    data_ret_sub = data_ret((i-250):i,:);
    mu = mean(data_ret_sub);
    gamma = cov(data_ret_sub);

    %compute weigths
    [w1, mu_p, var_p, time] = old_school(mu,gamma,false);
    [w2, mu_p, var_p, time] = null_space(mu,gamma,false,null_sep);
    [w3, mu_p, var_p, time] = cg_reduced_system(mu,gamma,false, null_sep, tol_cgred);
    [w4, mu_p, var_p, time] = cg_projected_test(mu,gamma,false, tol_cgproj);
    [w5, mu_p, var_p, time] = bregman(mu,gamma,false,data_ret_sub, 10e-15,rho,lambda);
    
    %result after a month
    j = ((i-251)/20)+1;
    rets = sum(data_ret(i:(i+20),:));
    mu1(j) = w1'*rets';
    mu2(j) = w2'*rets';
    mu3(j) = w3'*rets';
    mu4(j) = w4'*rets';
    mu5(j) = w5'*rets';
    mu6(j) = sum(data_ret_ind(i:(i+20)));
end

%plot the development of a 100 USD/EUR investment in 2000 until the
%present date
figure();
plot(ret2price(mu1, 100));
hold on
plot(ret2price(mu2, 100));
plot(ret2price(mu3, 100));
plot(ret2price(mu4, 100));
plot(ret2price(mu5, 100));
plot(ret2price(mu6, 100));
title('Performance Comparison')
xlabel('Weeks')
ylabel('Euros')
legend('Closed Form','Null Space','CG Reduced System', 'CG Projected', 'Split Bregman','DAX Index','Location','northwest')
saveas(gca, 'Figures/finalres','png')

toc
