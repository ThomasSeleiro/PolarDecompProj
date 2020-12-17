itArray = zeros(2,4);
accuracyArray = zeros(2,4);
matrixArray = {eye(8) hilb(6) cast(magic(6), "double") hadamard(8)};
namesArray = ["eye(8)" "hilb(6)" "magic(6)" "hadamard(8)"];
for i = 1:size(matrixArray,2)
    A = matrixArray{i};
    [U, H, its] = poldecTest(A, "n", [1e-16 1e-16]);
    acc = norm(A-U*H, 2);
    itArray(1,i) = its;
    accuracyArray(1,i) = acc;
    [U, H, its] = poldecTest(A, "h", [1e-16 1e-16]);
    acc = norm(A-U*H, 2);
    itArray(2,i) = its;
    accuracyArray(2,i) = acc;
end

%clf;
%hold on;
%yyaxis left;
%bar(itArray');
%h1 = ylabel("its", "FontName", "Consolas");
%yyaxis right;
%semilogy(accuracyArray', "Marker", "s", "Color", "r", "MarkerFaceColor", "r");
%ylabel('$\|A - UH\|_2$','Interpreter','latex')
%set(gca, "XTickLabel", namesArray);
%set(gca, "XTick", 1:size(itArray, 2));
%set(gca, 'XTickLabelRotation',45);
%legend("Newton", "N-Schulz", "Newton", "N-Schulz", "Location", "northwest");
%hold off;

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
    while(not(converged) && its < 100)
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
        converged = (unitDist <= convCond(1)*n) || (iterDist <= convCond(2)*n);
        
        X = Xnew;
        its = its + 1;
        fprintf("%4d\t%19.8e\t%13.8e\n", its, iterDist, unitDist);
    end
    U = Xnew;
    H = U' * A;
    H = (H + H') / 2;
end