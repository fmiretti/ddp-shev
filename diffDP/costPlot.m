function costPlot(V0_iter, V0_aug_iter, endpointError, opts)
%costPlot
%
% Cost and endpoint error vs iteration plot.

% Set up layout
npsi = size(endpointError, 2);
nrows = fix((npsi+1)/2) + 1;
colSpans = ones(npsi,1);
if mod(npsi,2)  ~= 0
    colSpans(end) = 2;
end

figure
tiledlayout(nrows, 2)

% Cost vs iteration number
nexttile([1, 1])
hold on
grid on

[tags2, tags1] = ndgrid(0:size(V0_iter,1)-1, 0:size(V0_iter,2)-1);
tags = tags1 + "." + tags2;
tags = tags(:);
tags(V0_iter==0) = [];

iterMajor = cumsum(sum(V0_iter~=0, 1)) + 1;
iterMajor = [1 iterMajor(1:end-1)];
V0_major = V0_iter(1,:);
V0_minor = V0_iter(:);
V0_minor( V0_minor == 0 ) = [];

plot(V0_minor, '--o', 'LineWidth', 1, 'Color', [1	0.498039215686275	0.0549019607843137])
plot(iterMajor, V0_major, '-o', 'LineWidth', 1.5, 'Color', [1	0.498039215686275	0.0549019607843137])

ylabelString = "Total cost";
if ~isempty(opts.CostUnit) && ~strcmp(opts.CostUnit, "")
    ylabelString = ylabelString + ", " + opts.CostUnit;
end
ylabel(ylabelString)
xlabel("Iteration")

ax = gca;
ax.XLim = [0 length(tags)+1];
ax.XTick = 1:length(tags);
ax.XTickLabel = tags;

% Augmented cost vs iteration number
nexttile([1, 1])
hold on
grid on

V0_major = V0_aug_iter(1,:);
V0_minor = V0_aug_iter(:);
V0_minor( V0_minor == 0 ) = [];

plot(V0_minor, '--o', 'LineWidth', 1, 'Color', [1	0.498039215686275	0.0549019607843137])
plot(iterMajor, V0_major, '-o', 'LineWidth', 1.5, 'Color', [1	0.498039215686275	0.0549019607843137])
legend("Minor iter", "Major iter")
ylabelString = "Augmented cost";
if ~isempty(opts.CostUnit) && ~strcmp(opts.CostUnit, "")
    ylabelString = ylabelString + ", " + opts.CostUnit;
end
ylabel(ylabelString)
xlabel("Iteration")

ax = gca;
ax.XLim = [0 length(tags)+1];
ax.XTick = 1:length(tags);
ax.XTickLabel = tags;

% FXEP error
for n = 1:npsi
    nexttile([1 colSpans(npsi)])
    plot(endpointError(:,n), '-o', 'LineWidth', 1, 'Color', '#2ca02c')
    grid on
%     title("Fixed endpoint error")
    ylabelString = opts.StateName(n)+ "(t_f) - " + opts.StateName(n) + "_f";
    if ~isempty(opts.StateUnit(n)) && ~strcmp(opts.StateUnit(n), "")
        ylabelString = ylabelString + ", " + opts.StateUnit(n);
    end
    ylabel(ylabelString)
    xlabel("Iteration")
    ax = gca;
    ax.XLim = [0 length(tags)+1];
    ax.XTick = 1:length(tags);
    ax.XTickLabel = tags;
end

end