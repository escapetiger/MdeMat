function [res] = gaussint(f,a,b,n)
    [x,w] = gausslegd(n,a,b,false);
    res = 0;
    for i = 1:numel(x)
        res = res+w(i)*f(x(i));
    end
%     res = w'*f(x);
end