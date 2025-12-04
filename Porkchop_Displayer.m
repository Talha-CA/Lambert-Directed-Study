% Porkchop Displayer


Rev_Range = 3;
Bool_List = [false true];
Earth_month = 7;
Earth_year = 2033;
Mars_month = 2;
Mars_year = 2035;

DelV_String = 'Departure Delta V (km/s)';
Title_String = 'Porkchop Plot ';




%% For Individual Figures
% for N = 0:Rev_Range
% 
%     for L = 1:2
%         Long = Bool_List(L);
% 
% 
%         data1 = readmatrix(append('Porkchop_', num2str(N), '_', string(Long), '_Results.csv'));
% 
%         sz = size(data1);
% 
%         %%
% 
% 
%         figure
%         surf(data1,'edgecolor','none'); hold on
%         title(append(Title_String, num2str(N), ' Revolutions, Long Way ', string(Long)));
%         xlabel(append('Arrival Days Since 1/', num2str(Mars_month),'/', num2str(Mars_year)));
%         ylabel(append('Departure Days Since 1/', num2str(Earth_month),'/', num2str(Earth_year)));
%         zlabel(DelV_String);
%         hColorbar = colorbar;
%         set(hColorbar, 'Ticks', sort([hColorbar.Limits, hColorbar.Ticks]));
%         hColorbar.Label.String = DelV_String;
%         colormap("turbo");
%         view(0,90)
%         minimum_result = min(min(data1));
%         maximum_result = max(max(data1));
% 
%         xlim([0 sz(2)])
%         ylim([0 sz(1)])
% 
%         [j,k] = find(data1==minimum_result);
%         p = plot3(k,j,maximum_result,'o','MarkerSize',10,'color','magenta','LineWidth',2);
%         % datatip(p,k,j);
% 
%         exportgraphics(gcf, append('Porkchop_', num2str(N), '_', string(Long), '_Results.png'), 'Resolution', 500);
%     end
% end


%% Subplots
figure;
sgtitle('Porkchop Subplots for Departure Delta V (km/s)') 

for N = 0:Rev_Range

    for L = 1:2
        Long = Bool_List(L);


        data1 = readmatrix(append('Porkchop_', num2str(N), '_', string(Long), '_Results.csv'));

        sz = size(data1);

        %%

        subplot(2,4,(N+1)+(Rev_Range+1)*(L-1))


        surf(data1,'edgecolor','none'); hold on
        title(append(num2str(N), ' Revolutions, Long Way ', string(Long)));
        xlabel(append('Arrival Days Since 1/', num2str(Mars_month),'/', num2str(Mars_year)));
        ylabel(append('Departure Days Since 1/', num2str(Earth_month),'/', num2str(Earth_year)));
        %zlabel(DelV_String);
        hColorbar = colorbar;
        clim([0 60])
        % set(hColorbar, 'Ticks', sort([hColorbar.Limits, hColorbar.Ticks]));
        %hColorbar.Label.String = DelV_String;
        colormap("turbo");
        view(0,90)
        minimum_result = min(min(data1));
        maximum_result = max(max(data1));

        xlim([0 sz(2)])
        ylim([0 sz(1)])

        [j,k] = find(data1==minimum_result);
        p = plot3(k,j,maximum_result,'o','MarkerSize',6,'color','magenta','LineWidth',1); % Display circle above plot
        % datatip(p,k,j);

    end
end
% exportgraphics(gcf,'Porkchop_Subplots_Departure_Common_Colorbar.png', 'Resolution', 500);