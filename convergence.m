p = @(x) (3*x - x.^3)/2;
x = linspace(0, 2.5, 3000);
p_x = p(x);
lims = zeros(size(x));
counter = 1;
for y=x
    temp = y;
    for i=1:100
        temp = p(temp);
    end
    lims(counter) = temp;
    counter = counter+1;
end