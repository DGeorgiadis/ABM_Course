function [ route_choice ] = Strategies( strategy_tp1 )
    % If new strategy is 1 (Main)
    if strategy_tp1 == 1
        route_choice(player,t) = strategy_tp1;
    % If the new strategy is 2 (Secondary)
    elseif strategy_tp1 == 2
        route_choice(player,t) = strategy_tp1;
    elseif strategy_tp1 == 3
        % If the payoff is lower than before, change route
        if payoff_t < median_all_prev_payoffs
            if route_choice_t == 1
                route_choice(player,t) = 2;
            else
                route_choice(player,t) = 1;
            end
        % If the payoff is higher than before, keep route
        elseif payoff_t > median_all_prev_payoffs
            route_choice(player,t) = route_choice_t;
        % If payoff is the same as before, decision rule depends on
            % nr. of occurances
        else
            if count_smaller_median < count_bigger_median
                route_choice(player,t) = route_choice_t;
            elseif count_smaller_median > count_bigger_median
                if route_choice_t == 1
                    route_choice(player,t) = 2;
                else
                    route_choice(player,t) = 1;
                end
            else
                route_choice(player,t) = randsample([1 2],1,true,[0.5 0.5]);
            end
        end
    % If selected strategy is nr. 4 (Contrarian)
    else
        if payoff_t > median_all_prev_payoffs
            if route_choice_t == 1
                route_choice(player,t) = 2;
            else
                route_choice(player,t) = 1;
            end
        elseif payoff_t < median_all_prev_payoffs
            route_choice(player,t) = route_choice_t;
        else
            if count_smaller_median > count_bigger_median
                route_choice(player,t) = route_choice_t;
            elseif count_smaller_median < count_bigger_median
                if route_choice_t == 1
                    route_choice(player,t) = 2;
                else
                    route_choice(player,t) = 1;
                end
            else
                route_choice(player,t) = randsample([1 2],1,true,[0.5 0.5]);
            end
        end  
    end
end

