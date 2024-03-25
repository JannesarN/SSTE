function [ internalParameters, w_downdate ] = general_qr_opt_func(internalParameters, x, y)
N = internalParameters.n;           %number of filter taps
lambda = internalParameters.alpha;  %regularization parameter

% x is the input vector (m1xN)
% y is the desired output vector (m1x1)
%[n1 m1]=size(x);
x_bar = [x,y,zeros(1,N+1)];
for row = 1:N+1
    [internalParameters, ac, bs] = BCsgr(internalParameters, x_bar(row), row);
    [internalParameters, x_bar]  = ICsgr(internalParameters, x_bar, ac, bs, row);        
end
%final cells
w_downdate = -lambda*internalParameters.r(N+1,N+2:2*N+1)';  %i add lambda via an experimental process
end
%boundary and internal cells functions
function [internalParameters, ac, bs] = BCsgr(internalParameters, x_bar, row)
    r_bar = internalParameters.r(row, row);
    sigma = internalParameters.sigma(row);
    beta = internalParameters.beta;

    r_prime = beta*beta*r_bar + sigma*x_bar*x_bar;
    ac = sigma*x_bar;
    if r_bar == 0
        bs = 0;
    else
        bs = x_bar/r_bar;
    end
    if r_prime == 0
        internalParameters.sigma(row+1) = sigma;
    else
        internalParameters.sigma(row+1) = beta*beta*sigma*r_bar/r_prime;
    end
    internalParameters.r(row, row) = r_prime;    
end

function [internalParameters, x_prime] = ICsgr(internalParameters, x_bar, ac, bs, row)
    beta = internalParameters.beta;
    N = internalParameters.n;
    r_bar = internalParameters.r(row, row+1:row+N+1);
    x_bar_row = x_bar(row+1:row+N+1);

    x_prime_tmp = x_bar_row - bs*r_bar;
    x_bar(row+1:row+N+1) = x_prime_tmp;
    x_prime = x_bar;
    internalParameters.r(row, row+1:row+N+1) = ac*x_bar_row + beta*beta*r_bar;
end
