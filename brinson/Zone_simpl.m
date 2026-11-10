function [ phase ] = Zone_simpl( sigma, sigma_s_cr_p, sigma_f_cr_p, sigma_s_cr_n, sigma_f_cr_n)

phase = 0;

if sigma > sigma_s_cr_p && sigma < sigma_f_cr_p
    phase = 'MD3p'
elseif sigma < sigma_s_cr_n && sigma > sigma_f_cr_n
    phase 'MD3n'
end

end
