function H = poldecsqrt(A)
%POLDECSQRT - Square-root using the Cholesky & polar decomposition
%
%   H = poldecsqrt(A) computes the square root of the input matrix (must 
%   be a Hermitian positive definite) using the polar decomposition
%   of the Cholesky factor (A = R*R,  R = UH,  A = H^2)

    % Compute the Cholesky factor of A
    [R, flag] = chol(A);
    if flag ~= 0
        error("The input matrix is not Hermitian Positive Definite");
    end

    % Compute the Polar Decomposition of R, return H
    [U, H, its] = poldec(R);
    
end