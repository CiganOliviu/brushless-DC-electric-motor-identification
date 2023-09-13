function [A, B, C, D] = getStateSpaceForSysOrderOne(K, T)
    A = -1/T;
    B = K/T; 
    C = 1;
    D = 0;
end