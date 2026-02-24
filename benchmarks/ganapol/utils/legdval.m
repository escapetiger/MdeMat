%% Function: Legendre values.
function res = legdval(m,x)
    M = max(m); x = x(:);
    P = zeros(numel(x),M+1);
    P(:,1)=1;
    if M>=1
        P(:,2)=x;
    end
    for l=1:M-1
        P(:,l+2)=((2*l+1).*x.*P(:,l+1)-l.*P(:,l))/(l+1);
    end
    res = P(:,m+1);
end

