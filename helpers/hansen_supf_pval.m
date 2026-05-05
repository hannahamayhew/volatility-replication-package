function pval = hansen_supf_pval(stat, m, pi0)
% =========================================================================
% DESCRIPTION:
% Approximate asymptotic p-value for Hansen (1997) SupF/SupW structural
% break tests using the chi-square approximation:
%
%   p(x) = 1 - chi2cdf(theta0 + theta1*x, eta)
%
% This version codes the coefficient set needed for the Sensier-style
% unconditional-volatility test:
%   m   = 1
%   pi0 = 0.15
%
% INPUTS:
%   stat = SupW (or SupF-family) statistic
%   m    = number of restrictions tested
%   pi0  = trimming parameter
%
% OUTPUT:
%   pval = approximate asymptotic p-value
% =========================================================================

    pval = NaN(size(stat));

    % Sensier UV case: one restriction, 15% trimming
    if ~(m == 1 && abs(pi0 - 0.15) < 1e-10)
        error('This helper currently only implements Hansen Table 2 for m=1, pi0=0.15.');
    end

    % Hansen (1997), Table 2, SupF, m=1, pi0=.15
    theta0 = -0.99;
    theta1 =  1.02;
    eta    =  3.0;

    z = theta0 + theta1 .* stat;

    % Guard against negative transformed values
    z(z < 0) = 0;

    pval = 1 - chi2cdf(z, eta);
end