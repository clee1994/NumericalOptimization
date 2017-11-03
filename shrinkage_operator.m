function [ retn_val ] = shrinkage_operator( x, gamma )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %definition given in the Ho paper
    retn_val = (x/abs(x)) * max(abs(x)-gamma,0);

end

