% x - sources
% y - targets
% check whether y interlaces x
function [flag, fail_idx] = check_interlace(x, y)

    % Ensure column vectors
     x = x(:);
     y = y(:);
   

    n = length(x);

    if length(y) ~= n
        error('x and y must have the same length');
    end

    fail_idx = find( ~(x(1:n-1) < y(1:n-1) & y(1:n-1) < x(2:n)) );
    flag = isempty(fail_idx);

end