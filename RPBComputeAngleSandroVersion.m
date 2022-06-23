function [theta,phi] = RPBComputeAngleSandroVersion(NUMBER_PROJECTIONS)
    %ComputeAngleSandroVersion It computes angles for each projection
    %   Detailed explanation goes here
    
    % constant
    C = 3.6;

    
    theta = zeros(NUMBER_PROJECTIONS,1);
    phi = zeros(NUMBER_PROJECTIONS,1);
    h = zeros(NUMBER_PROJECTIONS,1);

    for i=3:(NUMBER_PROJECTIONS)
        h(i) = (-1.0 + (2*(i-2)/NUMBER_PROJECTIONS));
        theta(i) = acos(h(i));
        phi(i) = mod(phi(i-1)+(C/sqrt(NUMBER_PROJECTIONS*(1-h(i)*h(i)))),2*pi); 
    end
    
    % MAKES NO SENSE, WHY HE DOES THAT???
    %theta(NUMBER_PROJECTIONS+1) = pi;
    %phi(NUMBER_PROJECTIONS+1) = pi;
end

