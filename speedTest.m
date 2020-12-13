itArray = zeros(2,4);
accuracyArray = zeros(2,4);
matrixArray = {rand(8) rand(20) hilb(6) magic(6) hadamard(8)};
numTries = 100;
timesArray = zeros(2, 5, numTries);
for i = 1:size(matrixArray,2)
    A = matrixArray{i};
    for k = 1:numTries
        newtonTime = tic;
        [U, H, its] = poldecTest(A, "n");
        timesArray(1, i, k) = toc(newtonTime);
        itArray(1, i) = its;
        nsTime = tic;
        [U, H, its] = poldecTest(A, "h", its);
        timesArray(2, i, k) = toc(nsTime);
        itArray(2, i) = its;
    end
    timesArray = timesArray ./ numTries;
end


function [U, H, its] = poldecTest(A, type, iters, conv)
    hybrid = true;
    convCond = [1e-16 1e-16];
    n = size(A,1);
    X    = A;
    Xnew = zeros(n);
    its  = 0;
    newtSchulz = false;
    converged = false;
    
    %We check the input arguments
    if nargin >= 2
        if type == "h"
        elseif type == "n"
            hybrid = false;
        elseif type == "ns"
            if norm(X,2) > sqrt(3)
                fprintf("The supplied matrix will diverge under the Newton-Schulz iteration");
            else
                newtSchulz = true;
            end
        else
            fprintf("The supplied type is not supported")
        end
    end
    if nargin == 4
        convCond = conv;
    end
    
    if nargin == 3
        fprintf("k   \t|X_k-X_{k-1}|/|X_{k}|\t|I - X_k*X_k|\n");
        fprintf("====\t===================\t==============\n");
        while(its < iters)
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

            newtSchulz = (norm(Xnew, 2) < sqrt(3)) && hybrid;
            converged = (unitDist <= convCond(1)) || (iterDist <= convCond(2));

            X = Xnew;
            its = its + 1;
            fprintf("%4d\t%19.8e\t%13.8e\n", its, iterDist, unitDist);
        end
    else

        fprintf("k   \t|X_k-X_{k-1}|/|X_{k}|\t|I - X_k*X_k|\n");
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

            newtSchulz = (norm(Xnew, 2) < sqrt(3)) && hybrid;
            converged = (unitDist <= convCond(1)) || (iterDist <= convCond(2));

            X = Xnew;
            its = its + 1;
            fprintf("%4d\t%19.8e\t%13.8e\n", its, iterDist, unitDist);
        end
    end
    U = Xnew;
    H = U' * A;
end