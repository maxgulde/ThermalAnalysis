function qA = albedoSolarFlux(Cs,subsolar,r,rho)
    % Returns the albedo's solar flux, multiply by area and absortion
    % coefficient to get power
    qA = Cs*(cosd(0.9*subsolar).^1.5)*viewFactor(r,rho);
end