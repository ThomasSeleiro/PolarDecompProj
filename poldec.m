function [U, H, its] = poldec(A)
%POLDEC Polar Decomposition
%   [U, H, ITS] = poldec(A) computes the polar decomposition 
%   A = U*H of the square, nonsingular matrix A. its is the 
%   number of iterations before convergence is achieved.
    n = size(A,1);
    X    = A;
    Xnew = zeros(n);
    its  = 0;
    newtSchulz = false;
    converged = false;
    fprintf("k   \t|X_k-X_{k-1}|/|X_{k}|\t|I - X_k*X_k|\n");
    fprintf("====\t===================\t==============\n");
    while(not(converged) && its < 100)
        if(not(newtSchulz))
            %We use the Newton method until either the 
            %   convergence conditionfor the Newton-Schulz 
            %   iterations is fulfilled, or convergence is
            %   acheived.
            Xnew = (X + inv(X)')/2;
        else
            %We use the Newton-Schulz method, having guaranteed 
            %   it will converge from this point onwards.
            Xnew = X/2 * (3*eye(n) - X' * X);
        end
        
        iterDist = norm(Xnew - X, inf)/norm(Xnew, inf);
        unitDist = norm(eye(n) - Xnew' * Xnew, inf);
        
        newtSchulz = norm(Xnew, 2) < sqrt(3);
        converged = (unitDist <= 1e-16*n) || (iterDist <= 1e-16*n);
        
        X = Xnew;
        its = its + 1;
        fprintf("%4d\t%19.8e\t%13.8e\n", its, iterDist, unitDist);
    end
    U = Xnew;
    Hstar = U' * A;
    H = (Hstar + Hstar')/2;
end
