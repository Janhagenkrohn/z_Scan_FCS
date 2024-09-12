% Author: Jan-Hagen Krohn, MPI for Biochemistry, 2024



function photon_weights = get_blcorr_weights(time_tags)

[trace_x, trace_y] = get_dense_trace(time_tags, 1000);

[~, ~, fitparam, ~] = blcorr_trace(trace_x, trace_y);

trace_y_fit_photons = polyval(fitparam, time_tags);

photon_weights = real(1 ./ sqrt(trace_y_fit_photons./trace_y_fit_photons(1)) + (1 - sqrt(trace_y_fit_photons./trace_y_fit_photons(1))));


end % function





function [trace_y_bleachcorr, trace_y_fit, fitparam, bleachcorr_degree] = blcorr_trace(trace_x, trace_y)


trace_w = sqrt(trace_y);
trace_w(trace_y == 0) = max(trace_w);

bleachcorr_degree = 1;

chi_sq = zeros([10,1]);
fitparam = cell([10, 1]);

while bleachcorr_degree <= 10

    % Fit trace with polynomial of user-defined degree
    fitparam{bleachcorr_degree} = polyfit(trace_x, trace_y, bleachcorr_degree);
    trace_y_fit = polyval(fitparam{bleachcorr_degree}, trace_x);
        
    % Get goodness of fit
    chi_sq(bleachcorr_degree) = sum((trace_y - trace_y_fit ./ trace_w).^2) / (length(trace_y) - bleachcorr_degree);

    % Check improvement
    if bleachcorr_degree > 1
        f_crit = finv(0.95, bleachcorr_degree - 1, bleachcorr_degree);
        f_val = chi_sq(bleachcorr_degree) / chi_sq(bleachcorr_degree - 1);

        if f_val <= f_crit
            % Insignificant improvement: Previous bleachcorr_degree was optimal, stop
            bleachcorr_degree = bleachcorr_degree - 1;
            break

        elseif bleachcorr_degree == 10
            % At 10 we stop either way
            break

        else
            % Still significant improvement: Next round
            bleachcorr_degree = bleachcorr_degree + 1;

        end % if f_val <= f_crit

    elseif bleachcorr_degree == 1
        % We look at least until bleachcorr_degree == 2
        bleachcorr_degree = bleachcorr_degree + 1;

    end % if bleachcorr_degree > 1



end % while bleachcorr_degree <= 10

% Get actual result at optimal bleachcorr_degree
fitparam = fitparam{bleachcorr_degree};
trace_y_fit = polyval(fitparam, trace_x);

% Perform correction and return relevant results
trace_y_bleachcorr = real(trace_y./sqrt(trace_y_fit./trace_y_fit(1)) + trace_y_fit(1).*(1 - sqrt(trace_y_fit./trace_y_fit(1))));

end % subfunction blcorr_trace





function [trace_x, trace_y] = get_dense_trace(time_tags, n_bins)

[trace_y, trace_x_appended] = histcounts(time_tags, n_bins);
trace_x = trace_x_appended(1:end-1);

end % subfunction get_dense_trace



