x = 5:100;
numTests = 32;
itArray = zeros(2,size(x,2));
accuracyArray = zeros(2,size(x,2));
for n=x
    i = n - x(1) + 1;
    for k=1:numTests
        A = rand(n);
        [U, H, tempIts] = poldecTest(A, "n", [1e-16 1e-16]);
        itArray(1,i) = itArray(1,i) + tempIts;
        accuracyArray(1,i) = accuracyArray(1,i) + norm(A-U*H,2);
        A = rand(n);
        [U, H, tempIts] = poldecTest(A, "h", [1e-16 1e-16]);
        itArray(2,i) = itArray(2,i) + tempIts;
        accuracyArray(2,i) = accuracyArray(2,i) + norm(A-U*H,2);
    end
    %itArray(:,i) = itArray(:,i) ./ numTests;
    %accuracyArray(:,i) = accuracyArray(:,i) ./ numTests;
end
itArray = itArray ./ numTests;
accuracyArray = accuracyArray ./ numTests;

%Plot the accuracies calculated
clf
hold on
box on
plot(x, accuracyArray(1,:), "color", "b", "Marker", "s", "MarkerFaceColor", 'b');
plot(x, accuracyArray(2,:), "color", "r", "Marker", "o", "MarkerFaceColor", "r");
legend("Newton", "N-Schulz", "Location", "northwest");
ylabel('$\|A - UH\|_2$','Interpreter','latex');
xlabel("n");
grid;
%saveas(gcf, "randnAccuracy", "pdf");
hold off

%fprintf("Press enter to see the number of iterations\n");
%w = waitforbuttonpress;

%Plot the iterations
clf
hold on
box on
plot(x, itArray(1,:), "color", "b", "Marker", "s", "MarkerFaceColor", 'b');
plot(x, itArray(2,:), "color", "r", "Marker", "o", "MarkerFaceColor", "r");
legend("Newton", "N-Schulz", "Location", "northwest");
ylabel("its", "FontName", "Consolas");
xlabel("n");
grid;
saveas(gcf, "randnIts", "pdf");
hold off


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
    %H = (H + H') / 2;
end