function x = clamp(x, max, min)
%CLAMP Clamps the value of x between min and max
%   Can be used for example to model saturation in nonlinear systems.
%   Follows the relationship
%
%           min <= x <= max
%
%   Arguments:
%       x:   Value to saturate
%       max: Upper bound of the saturation
%       min: (optional) Lower bound. Defaults to min = -max

%% Input validation
if length(x) ~= length(max)
    if length(max) == 1
        max = kron(ones(length(x),1), max);
    else
        error('The bounds must be scalar or of the same dimension as x')
    end
end

if nargin <= 2
    min = -max;
else
    if length(x) ~= length(min)
        if length(min) == 1
            min = kron(ones(length(x),1), min);
        else
            error('The bounds must be scalar or of the same dimension as x')
        end
    end
    
    if any(min > max)
        error('The lower bound must be smaller than the upper bound')
    end
end

%% Clamp input to desired range
if any(x < min)
    mask    = x<min;
    x(mask) = min(mask);
elseif x > max
    mask    = x>max;
    x(mask) = max(mask);
end
end
