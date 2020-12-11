x = 5:100;
numTests = 3;
itArray = [];
accuracyArray = [];
for n=x
    acc = 0;
    its = 0;
    for k=1:numTests
        A = rand(n);
        [U, H, tempIts] = poldecTest(A, "h", [1e-16 1e-16]);
        acc = acc + norm(A-U*H, 2);
        its = its + tempIts;
    end
    itArray = [itArray its/numTests];
    accuracyArray = [accuracyArray acc/numTests];
end

plot(x, accuracyArray, "Marker", "s", "MarkerFaceColor", 'b');

function [U, H, its] = poldecTest(A, type, conv)
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
    if nargin == 3
        convCond = conv;
    end
    
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
    U = Xnew;
    H = U' * A;
end