function f = mass_balance(C, Dz, dz, Dr, dr, uz, ur, EL)
% C: concentration vector
% Dz,Dr: diffusion coefficient vector
% dr,dz: grid spacing
% ur, uz: velocity vector

N = numel(C); % number of nodes (Concentration Values)
Nr = size(C,1)-1;
Nz = size(C,2)-1;
f = zeros(1,size(EL,1)+1);
inlet_flux = zeros(1,size(EL,1)+1); outlet_flux = zeros(1,size(EL,1)+1);
t = 0;

% for 1 = 2:(Nr)
%     for i = 2:(Nz)
for i = 1:(size(EL,1))
        % advection term
        adv_flux1 = uz*(C(1,(EL(i,1))));
        adv_flux2 = uz*(C(1,(EL(i,2))));%+ur*(C(1,i));

        % diffusive term
        diff_flux_in = -Dz*(C(1,(EL(i,2))) - C(1,(EL(i,1))))/dz; %- Dr*(C(1,left_index(EL(i,2))) - C(1-1,left_index(EL(i,1))))/dr;
        diff_flux_out = -Dz*(C(1,(EL(i,2))) - C(1,(EL(i,1))))/dz; % - Dr*(C(1+1,i) - C(1,i))/dr;

        % mass balance
        inlet_flux(EL(i,1)) = adv_flux1 + diff_flux_in;
        outlet_flux(EL(i,2)) = adv_flux2 + diff_flux_out;
        t = t+1;

        f = (inlet_flux-outlet_flux);
end
end
% end
