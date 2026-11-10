function [ phase ] = Zone( sigma, T, Mf, Ms, As, Af, sigma_s_cr_p, sigma_f_cr_p, sigma_s_cr_n, sigma_f_cr_n, C_M, C_A)

phase = 0;

if T <= Mf && sigma > sigma_s_cr_p && sigma < sigma_f_cr_p
    phase = 'MD3p';
elseif T <= Mf && sigma < sigma_s_cr_n && sigma < sigma_f_cr_n
    phase = 'MD3m'
elseif T > Mf && T<= Ms
    if sigma <= sigma_s_cr_p
        phase = 't';
    elseif sigma > sigma_s_cr_p && sigma < sigma_f_cr_p
        phase = 'ot';
    end
else
    if sigma > C_M * (T - Ms) + sigma_s_cr_p && sigma <= C_M * (T - Ms) + sigma_f_cr_p
        phase = 'M';
    elseif sigma > C_A * (T - Af) && sigma <= C_A * (T - As)
        phase = 'A';
    end
end

end

