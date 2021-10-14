clear
%before run this routine, run interpolate_data_1st_response.m
import plasma.Dip.LK2alpha
import plasma.pc2en

sdate_range = {'03-Sep-2015'};
edate_range = {'05-Sep-2015'};

mu_range = [300 3200 4500];
K_target = 0.172; % must be a number!
mag_field = 'TS07Dmid15';

satellite = 'rbspb';
input_folder = 'Interp_Data/';

print_figure = true;
figure_folder = '../Figure/';
figure_name = ['PSD_', mag_field, '_Inv'];

% Prepare data for plotting
nd = length(sdate_range);
nm = length(mu_range);

%% Choose in- and outbound data
PSD_info(nd,nm) = struct();
for id = 1:nd
    sdate = datenum(sdate_range{id});
    edate = datenum(edate_range{id});
    for im = 1:nm
        mu_target = mu_range(im);
        if mu_target <= 700
            instrument = 'mageis';
        else
            instrument = 'rept';
        end
        if strcmpi(instrument, 'rept') % NOTE: temporary duck tape for MagEIS HIGH
            filename = sprintf('%s_%s_%sto%s_%s_mu=%g_K=%g.mat', satellite, instrument, '20150901', '20150910', mag_field, mu_target, K_target);
        elseif strcmpi(instrument, 'mageis')
            filename = sprintf('%s_%s_%sto%s_%s_mu=%g_K=%g.mat', satellite, instrument, '20150901', '20150910', mag_field, mu_target, K_target);
        end
            
        load([input_folder, filename]);
        ind = find_bounds(xGEO.arr, time, sdate, edate);
        ni = length(ind);
        PSD_info(id,im).arr = cell(ni, 1);
        PSD_info(id,im).time = cell(ni, 1);
        PSD_info(id,im).Lstar = cell(ni, 1);
        for ii = 1:ni
            PSD_info(id,im).arr{ii} = PSD(ind{ii});
            PSD_info(id,im).time{ii} = time(ind{ii});
            PSD_info(id,im).Lstar{ii} = Lstar(ind{ii});
        end
    end
end

%% Calculate energies
en_info(nm) = struct();
k0=0.311;
mc2 = 0.511;
for im = 1:nm
    en_info(im).L_g = 1:0.1:6;
    mu_target = mu_range(im);
    pa_target = LK2alpha(en_info(im).L_g, K_target);
    pc = sqrt(2*mc2*k0./(en_info(im).L_g.^3)*mu_target) ./ sin(pa_target);
    en_info(im).arr = pc2en(pc);
end

%% Plot data
hf = figure('Units','inches','Position',[0 0 9.5 6],'PaperType','usletter',...
    'defaultAxesFontWeight','bold','Color','w');

% colormap for psd profiles
cmap = jet;
cmap = cmap([1:32,40:end], :);
cmap = cat(1, cmap, [0 0 0]);

lw = 1.2; % line width
rjc = 2; % row joining coefficient

nrows = rjc*nd + 1;
ncols = nm;

sp = gobjects(nd + 1, ncols);
% Plot PSD profiles
for id = 1:nd
    for im = 1:nm
        nn = length(PSD_info(id, im).arr);
        sp(id,im) = subplot(nrows, ncols, [rjc*(id-1)*ncols + im, (rjc*(id-1)*ncols + im + (rjc-1)*ncols)]);
        for in = 2:nn-1
            pp = plot(PSD_info(id, im).Lstar{in}, (PSD_info(id, im).arr{in}), 'o-', 'LineWidth', lw, 'MarkerSize', 3);
            % Choose color for plot
            ci = 1 + floor(size(cmap,1)/(nn-2)) * (in - 2);
            if ci <= 0
                ci = 1;
            end
            set(pp, 'Color', cmap(ci, :), 'MarkerFaceColor', cmap(ci,:));
            hold on
            % Fill legend string
            legend_str{in-1} = datestr(PSD_info(id, im).time{in}(end), ' mm/dd HH:MM');
        end
        if im == 1
            lg{id}= legend(legend_str, 'Orientation', 'Vertical');
        end
    end
end

% Plot energy values
for im = 1:nm
    sp(nd + 1, im) = subplot(nrows, ncols, (nrows-1)*ncols + im);
    pp = plot(en_info(im).L_g, en_info(im).arr, 'r-', 'LineWidth', lw);
end

%% Setup parameters for PSD profile subplots
xmin = 3.5;
xmax = 5.5;
set(sp(1:nd,:), 'XLim', [xmin, xmax], 'YScale', 'log', ...
    'XGrid', 'on', 'YGrid', 'on', 'Box', 'on', 'YMinorGrid', 'on', 'Units', 'Centimeter', ...
    'XTick', [xmin:0.5:xmax], 'YMinorTick', 'on');
ymax = [1e-4, 2e-7, 1e-8];
ymin = [6e-7, 2e-9, 1e-10];
for id = 1:nd
    for im = 1:nm
        sp(id,im).YLim = [ymin(im), ymax(im)];        
        sp(id,im).XAxis.LineWidth = lw;
        sp(id,im).YAxis.LineWidth = lw;
        if im == 1
            sp(id,im).YLabel.String = 'PSD, (c/cm/MeV)^3';
        end
        if id == 1
            title(sp(id,im), ['\mu = ', num2str(mu_range(im)), ' MeV/G']);
        end
        uistack(sp(id,im), 'bottom');
    end
end

%% Setup legend location and properties
for id = 1:nd
    set(lg{id}, 'FontSize', 7, 'Location', 'southeast', 'Units', 'Centimeter');
    uistack(lg{id}, 'top');
end

%% Setup parameters for energy subplots
set(sp(nd+1,:), 'XLim', [xmin, xmax], 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'FontWeight', 'bold', ...
    'XTick', [xmin:0.5:xmax], 'Units', 'Centimeter');
 for im = 1:nm
     sp(nd+1,im).XAxis.LineWidth = lw;
     sp(nd+1,im).YAxis.LineWidth = lw;
     sp(nd+1,im).XLabel.String = 'L^*';
     sp(nd+1,im).YLabel.String = 'E, MeV';
     sp(nd+1,im).XLabel.Units = 'Centimeter';
 end

 drawnow
 %% Added labels to plots
 txt = gobjects(nd + 1, nm);
 for id = 1:nd+1
     for im = 1:nm
         x = sp(id, im).Position(3) * 0 + 0.1;
         y = sp(id, im).Position(4) * 0;
         txt(id, im) = text(x, y, [char('a' + (id-1)*nm + im - 1), ')'], 'Units', 'Centimeter', ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'k', 'Layer', 'front', 'Parent', sp(id, im), ...
            'FontSize', 12, 'FontWeight', 'bold');
     end
 end
 
drawnow;
hf.PaperPositionMode = 'auto';
hf.InvertHardcopy = 'off';
hf.PaperSize = hf.Position(3:4);
hf.Renderer = 'opengl';
if print_figure
    print(hf, [figure_folder, figure_name, '.png'], '-dpng', '-r600');
    print(hf, [figure_folder, figure_name, '.pdf'], '-dpdf', '-r600');
end