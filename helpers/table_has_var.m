function tf = table_has_var(T, varname)
% =========================================================================
% DESCRIPTION:
% True if table T has variable varname
% =========================================================================

    if isstring(varname)
        varname = char(varname);
    end

    tf = ismember(varname, T.Properties.VariableNames);
end