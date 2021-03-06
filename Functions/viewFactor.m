function FE = viewFactor(r,rho)
    % r = R / (R + h), R = Earth radius, h = altitude
    % rho = angle between surface normal and local zenith
    % Angles in degrees
    
    e = - 160.31*(r^6) ...
        + 723.36*(r^5) ...
        - 1380*(r^4) ...
        + 1394.6*(r^3) ...
        - 780.65*(r^2) ...
        + 226.81*r ...
        - 21.232;
    
    FE = (r.^2.1).*(sind(rho/2).^e);
end