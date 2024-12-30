% SVEAIR model function
function dp = sveair_model(t, p, par)
    a = par(1);
    w = par(2);
    alpha = par(3);
    bs = par(4);
    ba = par(5);
    zi = par(6);
    e = par(7);
    u = par(8);
    sigma = par(9);
    r = par(10);
    eta = par(11);
    del = par(12);
    phi = par(13);
    
    dp = zeros(6,1);
    N=p(1)+p(2)+p(3)+p(4)+p(5)+p(6);
    dp(1) = a + w * p(2) + alpha * p(6) - ((bs * p(5) + ba * p(4))/N) * p(1) - (u + zi) * p(1);
    dp(2) = zi * p(1) - (1 - e) * ((bs * p(5) + ba * p(4))/N) * p(2) - (u + w) * p(2);
    dp(3) = ((bs * p(5) + ba * p(4)) / N) * p(1) + (1 - e) * ((bs * p(5) + ba * p(4)) / N) * p(2) - (u + sigma) * p(3);
    dp(4) = (1-r)* sigma * p(3) - (u + eta) * p(4);
    dp(5) = (r) * sigma * p(3) - (u + del + phi) * p(5);
    dp(6) = eta * p(4) + phi * p(5) - (u + alpha) * p(6);
end