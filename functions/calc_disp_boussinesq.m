function [ue,un,uz] = calc_disp_boussinesq(x,y,z,G,nu)
% function to calculate displacements due to a unit normal load at the origin in
% a halfspace using Boussinesq solution for point sources
% x,y,z can be vectors or matrices
% Rishav Mallick, EOS, 2018

r = sqrt(x.^2 + y.^2 + z.^2);
iszero = r==0;

ue = 1/(4*pi*G).*((1-2*nu).*x./(r.*(z+r)) - x.*z./r.^3);
un = 1/(4*pi*G).*((1-2*nu).*y./(r.*(z+r)) - y.*z./r.^3);
uz = -1/(4*pi*G).*(2*(1-nu)./(r) + (z.^2)./r.^3);

ue(iszero)=0;
un(iszero)=0;
uz(iszero)=0;
end