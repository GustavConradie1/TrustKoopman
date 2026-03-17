function [coef,breakPt,R2] = fitBogartz(x,y,verbose)
%FITBOGARTZ performs two-segmented linear regression
%   [coef,breakPt,R2] = fitBogartz(x,y,verbose)
%   performs two-segmented linear regression described in Bogartz (1968).
%   R. S. Bogartz (1968): "A least squares method for fitting intercepting
%   line segments to a set of data points." Psychol Bull 70(6), 749-755.
%
%   x      : independent variable
%   y      : dependent variable
%   verbose: outputs calculation process (1) or not (0)
%   coef   : [a1 a0 b1 b0] y=a1x+a0; y=b1x+b0
%   breakPt: [breakPtX breakPtY]
%   R2     : R-squared
%
%(c) S. Okazaki 2023.01

x_order = unique(x,'sorted');

n       = length(x_order);
output  = zeros(n-3,2);
results = zeros(n-3,7);

for k = 2:length(x_order)-2
    x0 = x_order(k);
    i1 = find(x <= x0);
    i2 = find(x  > x0);
    
    x1 = x(i1);
    y1 = y(i1);
    x2 = x(i2);
    y2 = y(i2);
    
    lm1 = polyfit(x1,y1,1);
    lm2 = polyfit(x2,y2,1);
    
    y1fit = polyval(lm1,x1);
    y2fit = polyval(lm2,x2);
    
    y1resid = y1 - y1fit;
    y2resid = y2 - y2fit;
    
    SSresid1 = sum(y1resid.^2);
    SSresid2 = sum(y2resid.^2);
    
    SSresid = SSresid1 + SSresid2;
    SStotal = sum((y - mean(y)).^2);
    
    R2 = 1 - SSresid / SStotal;
    
    a1 = lm1(1);
    a0 = lm1(2);
    b1 = lm2(1);
    b0 = lm2(2);
    
    breakPtX = (b0-a0) / (a1 - b1);
    breakPtY = a1 * breakPtX + a0;
    
    output(k-1,:)  = [x0 SSresid];
    results(k-1,:) = [a1 a0 b1 b0 breakPtX breakPtY R2];
end

output = array2table(output,"VariableNames",["UpTo","SS"]);
if verbose == 1
    output
end

[~,I] = min(output.SS);
coef    = results(I,1:4);
breakPt = results(I,5:6);
R2      = results(I,7);

end

