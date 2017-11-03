function [ Y,Z ] = null_range( gamma, null_sep, graphs_on )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    %singular value decomposition
    
    [~, S, V] = svd(gamma);
    sv = diag(S);

    %find the range space
    sv_rel = sv./sv(1);
    ind_rel = find(sv_rel < null_sep);
    %tol= sv_rel(1);
    %Ri=find(sv>tol);%   index vector for range

    %find the null space
    %Ni=find(sv<tol);%   index vector for null

    %range and null
    Y=V(:,setdiff(1:end,ind_rel)');
    Z=V(:,ind_rel');
    
    %when no Nullspace, use zeros
    [m, n]=size(Z); if n<1,Z=zeros(m,1);end;
    
    %checking gamma for singularity
    if graphs_on
        %decay_fig = figure();
        %set(decay_fig, 'Position', [0 0 800 300])
        figure();
        plot(sv);
        title('Singular values')
        saveas(gca, 'Figures/nullsv1','png')
        figure();
        semilogy(sv)
        title('Singular values in Semilog')
        saveas(gca, 'Figures/nullsv2','png')
        if sum(diag(eig(Z'*gamma*Z)))>=0;disp('G is PSD'),end;
    end
    
    %conditions
    %if rank(A)==size(A,1);disp('A is full rank'),end;
    
    
    
    %ill condition
    
end

