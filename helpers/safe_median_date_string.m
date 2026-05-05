function out = safe_median_date_string(x)
% =========================================================================
% DESCRIPTION:
% Safely compute the median date from a string array/cellstr of dates in
% format dd-MMM-yyyy. Returns "" if no usable dates are present.
% =========================================================================

    out = "";

    if isempty(x)
        return;
    end

    try
        dt = datetime(x, 'InputFormat','dd-MMM-yyyy');
        dt = dt(~isnat(dt));
        if isempty(dt)
            return;
        end
        out = string(datetime(median(datenum(dt)), 'ConvertFrom','datenum'), 'dd-MMM-yyyy');
    catch
        out = "";
    end
end