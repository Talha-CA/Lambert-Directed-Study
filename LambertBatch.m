clc
clear
tic

% Constants
global mu
mu = 1.327124e11;                                   % km3 s-2
deg = pi/180;
max_iterations = 100;           % 100 (also check sma_max in func)
Bool_List = [false,true];

Rev_Range = 3;                  % 3
Start_Range = 6*365;            % 4*365
End_Range = 6*365;

results = zeros(Start_Range,End_Range);

for N = 0:Rev_Range

    for L = 1:2
        Long = Bool_List(L);
        disp('Solving next Batch')

        for i = 1:Start_Range
            % Earth OEs
            Earth_id = 3;
            Earth_year = 2033;
            Earth_month = 7;
            Earth_day = i;
            Earth_hour = 12;
            Earth_minute = 0;
            Earth_second = 0;
            %...
            %...Algorithm 8.1:
            [coe_Earth, r_Earth, V_Earth, jd_Earth] = planet_elements_and_sv(Earth_id, Earth_year, Earth_month, Earth_day, Earth_hour, Earth_minute, Earth_second);
        
            for j = 1:End_Range
                % Mars OEs
                %...Input data
                Mars_id = 4;            % 4
                Mars_year = 2035;       % 2035
                Mars_month = 2;
                Mars_day = j;
                Mars_hour = 12;
                Mars_minute = 0;
                Mars_second = 0;
                %...
                %...Algorithm 8.1:
                [coe_Mars, r_Mars, V_Mars, jd_Mars] = planet_elements_and_sv(Mars_id, Mars_year, Mars_month, Mars_day, Mars_hour, Mars_minute, Mars_second);
        
                R1_vec = r_Earth;
                R2_vec = r_Mars;
        
        
                ToF_input = (jd_Mars-jd_Earth)*24*60*60; % Calculate off OE dates
        
                [DeltaV1, DeltaV2] = Func_Lambert_Mars(R1_vec, R2_vec, ToF_input, Long, N, mu, max_iterations, V_Earth, V_Mars);
        
                % Store the results for analysis
                results(i, j, :) = DeltaV1;
        
            end
        end
                      
        
        %%
        disp('Drawing Figure')

        figure;
        surf(results,'edgecolor','none'); hold on
        title(append('Porkchop Plot ', num2str(N), ' Revolutions, Long Way ', string(Long)));
        xlabel(append('Arrival Days Since 1/', num2str(Mars_month),'/', num2str(Mars_year)));
        ylabel(append('Departure Days Since 1/', num2str(Earth_month),'/', num2str(Earth_year)));
        zlabel('Departure + Arrival Delta V (km/s)');
        hColorbar = colorbar;
        set(hColorbar, 'Ticks', sort([hColorbar.Limits, hColorbar.Ticks]));
        hColorbar.Label.String = 'Arrival Delta V (km/s)';
        colormap("turbo");
        view(0,90)
        minimum_result = min(min(results));
        maximum_result = max(max(results));

        xlim([0 End_Range])
        ylim([0 Start_Range])
        
        [j,k] = find(results==minimum_result);
        p = plot3(k,j,maximum_result,'o','MarkerSize',6,'color','magenta','LineWidth',1); % Display circle above plot
        % datatip(p);
        
        %%
        % figure;
        % contour(results,1000); hold on
        % xlabel('Arrival Days After May 31, 2035');
        % ylabel('Departure Days After Nov 30, 2034');
        % 
        % hColorbar = colorbar;
        % set(hColorbar, 'Ticks', sort([hColorbar.Limits, hColorbar.Ticks]));
        % hColorbar.Label.String = 'Departure Delta V (km/s)';
        % colormap("turbo");
        % minimum_result = min(min(results));
        % 
        % [j,k] = find(results==minimum_result);
        % p = plot(k,j,'x');
        % datatip(p,k,j);
        
        %%
        
        % Define the output file name
        outputFileName = append('Porkchop_Arrival_', num2str(N), '_', string(Long), '_Results.csv');
        
        % write results
        % writematrix(results, outputFileName);
        
        disp(append('Results have been written to ', outputFileName));
    end
end

disp(['Runtime = ',num2str(toc), ' sec for ', num2str(Rev_Range), ' revolutions'])