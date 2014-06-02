function haxis = multiplot(x, yz, xrange, yzrange, color, logscale, xlabeltext, yzlabeltext, yznames, dataname, haxis, yy, yyrange, yylabeltext)

% Check input arguments
ny = size(yz, ndims(yz));

if (nargin < 3) || isempty(xrange)
    xrange = [x(1) x(end)];
end

if (nargin < 5) || isempty(color)
    color = 'b';
end

if (nargin < 6) || isempty(logscale)
    logscale = false;
end

if (nargin < 4) || isempty(yzrange)
    if logscale
        yzrange = [max(min(min(yz)), 1e-6) max(max(yz))];
    else
        yzrange = [min(min(yz)) max(max(yz))];
    end
end

if (nargin < 7) || isempty(xlabeltext)
    xlabeltext = '';
end

if (nargin < 8) || isempty(yzlabeltext)
    yzlabeltext = 'a.u.';
end

if (nargin < 9)
    yznames = [];
end

if (nargin < 10)
    dataname = [];
end

% Plot
for i = 1:ny
    if (nargin < 11) || isempty(haxis)
        figure;
        haxis(i) = gca;
    else
        axes(haxis(i));
    end
    switch ndims(yz)
        case 2
            hnew = plot(x, yz(:,i), 'Color', color);
        case 3
            hnew = surf(x, yy, yz(:,:,i)', 'EdgeColor', 'none', 'LineStyle', 'none', 'FaceLighting', 'phong');
            view(0,90);
    end
    hold on;
    grid on;
    if logscale
        set(haxis(i), 'YScale', 'log');
    else
        set(haxis(i), 'YScale', 'linear');
    end
    axis([xrange yzrange])
    xlabel(xlabeltext, 'FontSize', 14, 'FontWeight', 'bold');
    ylabel(yzlabeltext, 'FontSize', 14, 'FontWeight', 'bold');
    set(haxis(i), 'FontSize', 14);
    if ~isempty(yz)
        [~,~,outh,outm] = legend;
        n = length(outm);
        outm{n+1} = '';
        if ~isempty(yznames)
            outm{n+1} = yznames{i};
        end
        if ~isempty(dataname)
            outm{n+1} = [outm{n+1} ' - ' dataname];
        end
        legend([outh;hnew], outm, 'FontSize', 10);
    end
end
