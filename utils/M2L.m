function [B1, B2] = M2L( r, a, b, dx, dy,scaling)
% centers: a and b diam:  dx and dy
% fun = 1, 1/(x-y)
% fun = 2, 1/(x-y)^2
% B1 translation matrix for 1./(x-y)
% B2 translation matrix for 1./(x-y)^2

% ba = b-a;
ab = a-b;
% S = diag( (-1).^(r-1:-1:0) );
% S = S( 1:r, r:-1:1 );
DD = diag( (-1).^(0:r-1) );
B1 = zeros(r);



switch scaling
    case 0
        % no scaling in B
        cc = -ones(r,1)/ab;
        for k = 1:r-1
            cc(k+1) = k*cc(k)/ab;
        end
        B1 = zeros(r);
        for i = 1:r
            B1(i,1:r-i+1) = cc(i:end);
        end

    case 1
        sr =  (sqrt(pi)*((8*r^3+4*r^2+r+(1/30))^(1/6)))^(1/r)/exp(1);
        B1(1,1) = -1/ab;
        B1(1,2) = -1/ab^2/(sr*2/dy);
        B1(2,1) = -1/ab^2/(sr*2/dx);
        B1(2,2) = 2/ab*B1(2,1)/(sr*2/dy);
        % complete first row of B
        for i = 3:r
            B1( 1,i ) = (i-1)/ab*B1( 1,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
        end
        % complete second row of B
        for i = 3:r-1
            B1( 2,i ) = i/ab*B1( 2,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
        end
        % complete first column of B
        for i = 3:r
            B1( i,1 ) = (i-1)/ab*B1( i-1,1 )/(sr*2/dx)/(i-1)*(1-1/(i-1))^(i-2);
        end
        % complete second column of B
        for i = 3:r-1
            B1( i,2 ) = i/ab*B1( i-1,2 )/(sr*2/dx)/(i-1)*(1-1/(i-1))^(i-2);
        end
        % complete the rest
        for k = 3:r-2
            for i = 3:r-k+1
                B1( k,i ) = (k+i-2)/ab*B1( k,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
            end
        end
end

B1 = DD*B1;

B2 = zeros(r);
switch scaling
    case 0   % without scaling derivative translation matrix
        eta = 1; 
        etaba = 1/( eta*ab );
        cc = 1/ab^2 * ones( r,1 );
        for k = 1:r-1
            cc(k+1) = ( k+1 )*cc(k)*etaba;
        end
         for i = 1:r
            B2(i,1:r-i+1) = cc(i:end);
        end
        % C = triu( toeplitz( cc(end:-1:1) ) );
        B2 = DD * B2;

    case 1  % with scaling derivative translation matrix
        sr =  (sqrt(pi)*((8*r^3+4*r^2+r+(1/30))^(1/6)))^(1/r)/exp(1);
        B2(1,1) = 1/ab^2;
        B2(1,2) = 2/ab^3/(sr*2/dy);
        B2(2,1) = 2/ab^3/(sr*2/dx);
        B2(2,2) = 3/ab*B2(2,1)/(sr*2/dy);
        % complete first row of B
        for i = 3:r
            B2( 1,i ) = i/ab*B2( 1,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
        end
        % complete second row of B
        for i = 3:r-1
            B2( 2,i ) = (i+1)/ab*B2( 2,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
        end
        % complete first column of B
        for i = 3:r
            B2( i,1 ) = i/ab*B2( i-1,1 )/(sr*2/dx)/(i-1)*(1-1/(i-1))^(i-2);
        end
        % complete second column of B
        for i = 3:r-1
            B2( i,2 ) = (i+1)/ab*B2( i-1,2 )/(sr*2/dx)/(i-1)*(1-1/(i-1))^(i-2);
        end
        % complete the rest
        for k = 3:r-2
            for i = 3:r-k+1
                B2(k,i) = (k+i-1)/ab*B2(k,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
            end
        end
        B2 = B2*DD;
end


