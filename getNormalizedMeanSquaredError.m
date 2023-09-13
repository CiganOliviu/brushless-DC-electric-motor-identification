function [result] = getNormalizedMeanSquaredError(y, y0)
    result = norm(y-y0) / norm(y-mean(y));
end