function [result] = getRelativeError(y, y0)
    result = sqrt(1/length(y)*sum((y-y0).^2));
end