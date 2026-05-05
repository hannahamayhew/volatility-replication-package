function x = get_table_var_or_default(T, varname, default_value)
% =========================================================================
% DESCRIPTION:
% Safely return a table variable if it exists; otherwise return a default.
%
% INPUTS:
%   T             = table
%   varname       = variable name as char/string
%   default_value = scalar or array to use if varname is missing
%
% OUTPUT:
%   x             = T.(varname) if it exists; otherwise a vector/array
%                   filled with default_value and sized to height(T)
% =========================================================================

    if isstring(varname)
        varname = char(varname);
    end

    if ismember(varname, T.Properties.VariableNames)
        x = T.(varname);
    else
        if isscalar(default_value)
            x = repmat(default_value, height(T), 1);
        else
            x = default_value;
        end
    end
end