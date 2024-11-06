function Cp = calcCp_from_K(K,C_psi,m_target)
% compute Cp from sensitivity kernel, assumed covariance of hidden parameters and target slip distribution 
% Rishav Mallick, 2023 Caltech Seismolab

K_psi = zeros(size(K));
for i = 1:length(m_target)
    K_psi(:,i) = K(:,i).*m_target(i);
end
Cp = (K_psi*C_psi)*K_psi';

end