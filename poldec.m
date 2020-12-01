function [U, H, its] = poldec(A)
%POLDEC Polar Decomposition
%   [U, H, ITS] = poldec(A) computes the polar decomposition A = U*H
%   of the square, nonsingular matrix A. ITS is the number of
%   iterations for convergence.
    n = size(A,1);
    X    = A;
    Xnew = zeros(n);
    its  = 0;
    newtSchulz = false;
    converged = false;
    fprintf("k   \t|X_{k_1}-X_k|/|X_k|\t|I - X_k*X_k|\n");
    fprintf("====\t===================\t==============\n");
    while(not(converged) && its < 1000)
        if(not(newtSchulz))
            %We use the Newton method until either the convergence condition
            %   for the Newton-Schulz iterations is fulfilled, or convergence
            %   is acheived.
            Xnew = (X + inv(X)')/2;
        else
            %We use the Newton-Schulz method, having guaranteed it will
            %   converge from this point onwards.
            Xnew = X/2 * (3*eye(n) - X' * X);
        end
        
        iterDist = norm(Xnew - X, inf)/norm(Xnew, inf);
        unitDist = norm(eye(n) - Xnew' * Xnew, inf);

        
        %Insert Newton-Schulz condition here
        newtSchulz = false;
        converged = unitDist <= 5e-15;
        
        X = Xnew;
        its = its + 1;
        fprintf("%4d\t%19.8e\t%13.8e\n", its, iterDist, unitDist);
    end
    U = Xnew;
    H = U' * A;
end
