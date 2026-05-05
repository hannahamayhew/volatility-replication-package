function vol_obj = build_uv_object(y_window)
% BUILD_UV_OBJECT  Construct the unconditional volatility proxy.
%
%   vol_obj = BUILD_UV_OBJECT(y_window)
%
%   Implements Sensier & van Dijk (2004) definition:
%
%       z_t = sqrt(pi/2) * | y_t - mean(y) |
%
%   The sqrt(pi/2) factor makes z_t an unbiased estimator of sigma_t
%   under the assumption that y_t - mean(y) is conditionally normal
%   (because E[|X|] = sigma * sqrt(2/pi) for X ~ N(0,sigma^2)).
%
%   Demeaning is column-wise with nanmean so that missing observations
%   do not bias the mean estimate within each window.
%
% -------------------------------------------------------------------------
% INPUTS
%   y_window : T x N matrix of transformed series (may contain NaN)
%              Each column is one FRED-MD series within the estimation window.
%
% OUTPUT
%   vol_obj  : T x N matrix of unconditional volatility proxies z_t.
%              NaN entries in y_window propagate as NaN in vol_obj.
%
% REFERENCE
%   Sensier, M. & van Dijk, D. (2004). "Testing for Volatility Changes in
%   U.S. Macroeconomic Time Series." Review of Economics and Statistics,
%   86(3), 833-839.
% =========================================================================

    % Column-wise mean (ignores NaN)
    mu_hat = nanmean(y_window, 1);                               % 1 x N

    % Demean each column
    y_dm = y_window - repmat(mu_hat, size(y_window,1), 1);      % T x N

    % Apply Sensier & van Dijk scaling
    vol_obj = sqrt(pi/2) * abs(y_dm);                           % T x N

end
