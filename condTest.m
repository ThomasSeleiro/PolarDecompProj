c = 10 * ones(1, 16);
for i = 2:16
    c(i) = c(i) * c(i-1);
end

itArray = zeros(1, size(c,2));
accuracyArray = zeros(1, size(c,2));
condArray = zeros(1, size(c,2));
skewHermitianArray = zeros(1, size(c,2));
for i = 1:size(c,2)
    [A, calcCond] = condMatrix(10, c(1, i));
    condArray(:, i) = calcCond;
    %fprintf("condNumber: %g\t calcCond:%g\n", 1/condNumber, calcCond);
    [U, H, its] = poldecTest(A);
    itArray(1, i) = its;
    accuracyArray(1, i) = norm(A - U*H);
    skewHermitianArray(1, i) = norm(H - H')/2;
end

clf
box on
%loglog(condArray(1,:), accuracyArray(1, :), "color", "b", "Marker", "s", "MarkerFaceColor", 'b');
%hold on
loglog(condArray(1,:), skewHermitianArray(1,:), "color", "r", "Marker", "o", "MarkerFaceColor", "r");
%hold off
%legend('$\|A - UH\|_2$', '$(\|H - H^*\|_2)/2$', 'Interpreter','latex', "Location", "northwest");
ylabel('$(\|H - H^*\|_2)/2$','Interpreter','latex');
xlabel("$\kappa_2(A)$", 'Interpreter','latex');
grid;
%saveas(gcf, "randnIts", "pdf");


function [A, calcCond] = condMatrix(n, condNumber)
    B = rand(n);
    [P, ~, Q] = svd(B);
    s = (rand(n, 1)*condNumber);
    s = sort(s, "descend");
    s(n) = 1; s(1) = condNumber;
    A = P * diag(s) * Q';
    calcCond = cond(A, 2);
end

function [U, H, its] = poldecTest(A)
%POLDEC Polar Decomposition
%   [U, H, ITS] = poldec(A) computes the polar decomposition 
%   A = U*H of the square, nonsingular matrix A. ITS is the 
%   number ofiterations for convergence.
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
    H = U' * A;
    %H = (H + H')/2;
end