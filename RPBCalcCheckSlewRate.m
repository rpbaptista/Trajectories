function [validate_tag, slew_rate] = RPBCalcCheckSlewRate(gradient,limit, dt)
%CalcCheckSlewRate Calculate and Check if the gradient respect the limits return true if it does and false if it does not.
%   Detailed explanation goes here
%   Author: Renata PORCIUNCULA BAPTISTA
%   e-mail: renata.porciunculabaptista@cea.fr

    
    if ~exist('limit','var') || isempty(limit);  limit = 200e-3;  end;
    if ~exist('dt','var') || isempty(dt); dt = 10e-6;  end;
    
    gradient = gradient/(dt*1000); 
    slew_rate = diff(gradient,1,2); % diff(matrix; order; axis)
    
    fprintf('SlewRateMax: %.02d\n',max(abs(slew_rate(:))));

    if max(abs(slew_rate(:))) > limit
        validate_tag = false;
    else
        validate_tag = true;
    end
end

