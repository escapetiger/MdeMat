function S = euler_knopp(a,N,pe)
    S = 0;
    for k = 0:N
        for j = 0:k
            b = binom(k,j);
            u = pe.^j;
            v = (-1).^(k-j).*a(k-j);
            S = S+b*u*v/(1+pe)^(k+1);
        end
    end

%     Version 1
%     S = 0;
%     for k = 0:N
%         b = binom(k,0:k);
%         u = pe.^(0:k);
%         v = (-1).^(k:-1:0).*a(k:-1:0);
%         S = S+sum(b.*u.*v)/(1+pe)^(k+1);
%     end

%     Version 2
%     B = binom(0:N,0:N);
%     U = pe.^(0:N);
%     V = (-1).^(0:N).*a(0:N);
%     C = (1+pe).^(1:N+1);
%     UU = tril(repmat(U,N+1,1));
%     VV = flip(tril(repmat(V,N+1,1)),2);
%     for k=1:N+1
%         VV(k,:) = circshift(VV(k,:),k-N-1,2);
%     end
%     S = sum(sum(B.*UU.*VV,2)./C(:));
end