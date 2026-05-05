function rec = make_nber_recession_indicator(dates_dt)
% =========================================================================
% DESCRIPTION:
% Create a monthly NBER recession indicator aligned to a datetime vector.
%
% Convention:
% - recession starts in the month after an NBER peak
% - recession ends in the NBER trough month
%
% Based on official NBER chronology through the 2020 recession:
% peaks: 1957-07, 1960-04, 1969-12, 1973-11, 1980-01, 1981-07,
%        1990-07, 2001-03, 2007-12, 2020-02
% troughs: 1958-04, 1961-02, 1970-11, 1975-03, 1980-07, 1982-11,
%          1991-03, 2001-11, 2009-06, 2020-04
% =========================================================================

    dates_dt = dates_dt(:);
    rec = false(size(dates_dt));

    peak_months = [
        datetime(1957,7,1)
        datetime(1960,4,1)
        datetime(1969,12,1)
        datetime(1973,11,1)
        datetime(1980,1,1)
        datetime(1981,7,1)
        datetime(1990,7,1)
        datetime(2001,3,1)
        datetime(2007,12,1)
        datetime(2020,2,1)
    ];

    trough_months = [
        datetime(1958,4,1)
        datetime(1961,2,1)
        datetime(1970,11,1)
        datetime(1975,3,1)
        datetime(1980,7,1)
        datetime(1982,11,1)
        datetime(1991,3,1)
        datetime(2001,11,1)
        datetime(2009,6,1)
        datetime(2020,4,1)
    ];

    for i = 1:length(peak_months)
        rec_start = peak_months(i) + calmonths(1);
        rec_end   = trough_months(i);
        rec = rec | (dates_dt >= rec_start & dates_dt <= rec_end);
    end
end