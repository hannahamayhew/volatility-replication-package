function pval = hansen_supf_pval_general(stat, m, pi0)
% =========================================================================
% DESCRIPTION:
% Hansen (1997) approximate asymptotic p-value for SupF / SupW tests:
%
%   p(x) = 1 - chi2cdf(theta0 + theta1*x, eta)
%
% This helper supports the pi0 = 0.15 case needed for the Sensier-style
% UV and CM work. For CM, m = p + 1 varies by series.
%
% INPUTS:
%   stat = SupW / SupF-family statistic
%   m    = number of restrictions tested
%   pi0  = trimming parameter
%
% OUTPUT:
%   pval = approximate asymptotic p-value
%
% NOTES:
%   Coefficients are from Hansen (1997), Table 2, SupF Distribution,
%   pi0 = 0.15.
% =========================================================================

    if nargin < 3
        error('Usage: pval = hansen_supf_pval_general(stat, m, pi0)');
    end

    % This project currently only needs pi0 = 0.15
    if abs(pi0 - 0.15) > 1e-10
        error('This helper currently supports only pi0 = 0.15.');
    end

    % Hansen (1997), Table 2, SupF distribution, pi0 = 0.15
    % Columns are: [m, theta0, theta1, eta]
    coeffs = [
         1,  -0.99, 1.02,  3.0;
         2,  -1.65, 1.06,  4.7;
         3,  -2.05, 1.13,  6.8;
         4,  -2.52, 1.11,  8.0;
         5,  -3.46, 1.07,  8.3;
         6,  -4.05, 1.08,  9.5;
         7,  -4.42, 1.10, 11.0;
         8,  -5.36, 1.08, 11.3;
         9,  -5.43, 1.10, 13.1;
        10,  -6.47, 1.06, 12.8;
        11,  -6.79, 1.04, 13.5;
        12,  -7.80, 1.02, 13.6;
        13,  -7.93, 1.07, 15.9
    ];

    row = coeffs(coeffs(:,1) == m, :);

    if isempty(row)
        error('No Hansen Table 2 coefficient coded for m = %d at pi0 = 0.15.', m);
    end

    theta0 = row(2);
    theta1 = row(3);
    eta    = row(4);

    z = theta0 + theta1 .* stat;

    % Guard against negative transformed values
    z(z < 0) = 0;

    pval = 1 - chi2cdf(z, eta);
end