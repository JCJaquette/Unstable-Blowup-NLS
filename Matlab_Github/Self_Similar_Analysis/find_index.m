function [J_min, J_max ] = find_index(t_min,t_max,tspan)
    % Find the index for start of the time window
    tspan_finding_index =  abs(tspan - t_min);
    [~,J_Index] = sort(tspan_finding_index);
    J_start = J_Index(1);
    J_min = J_start;

    tspan_finding_index =  abs(tspan - t_max);
    [~,J_Index] = sort(tspan_finding_index);
    J_max = J_Index(1);
end