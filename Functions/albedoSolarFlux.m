function qA = albedoSolarFlux(Cs,subsolar,vF)
    % Returns the albedo's solar flux, multiply by area and absortion
    % coefficient to get power
    % Angles in degrees
    qA = Cs.*(cosd(0.9*subsolar).^1.5).*vF;
end