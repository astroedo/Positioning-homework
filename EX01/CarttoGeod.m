function [geodetic] = CarttoGeod(cartesian, a, e, b)
    n = size(cartesian, 2);
    lambda = atan2(cartesian(2,:), cartesian(1,:));
    r = sqrt(cartesian(1,:).^2 + cartesian(2,:).^2);
    theta = atan2(cartesian(3,:) * a, r * b);
    
    phi = atan2(cartesian(3,:) + e^2 * b * sin(theta).^3, ...
               r - e^2 * a * cos(theta).^3);
    
    N = a ./ sqrt(1 - e^2 * sin(phi).^2);
    h = r ./ cos(phi) - N;
    
    geodetic = [lambda; phi; h];
end