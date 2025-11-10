data = readtable("sorted_data.xlsx");

cmask = endsWith(data.Properties.VariableNames,'Concentration');
emask = endsWith(data.Properties.VariableNames,'Error1s');

words = data.Properties.VariableNames(cmask);
for i = 1:length(words)
    words{i} = erase(words{i},'Concentration');
end

concs = data(:,cmask);
concs = table2array(concs) ./ 1e4;

errors = data(:,emask);
errors = table2array(errors) ./ 1e4;

errors = errors(:,any(concs>0.7));
words = words(:,any(concs>0.7));
concs = concs(:,any(concs>0.7));
error_low = errors;
neg_check = concs - errors;
error_low(neg_check<0) = concs(neg_check<0);
error_high = errors;
figure
b = bar(concs','grouped');
b(1).FaceColor = [0.6 0.6 0.6];
b(2).FaceColor = [0 0.4470 0.7410];
b(3).FaceColor = [0 0.4470 0.7410];
b(4).FaceColor = [0.8500 0.3250 0.0980];
b(5).FaceColor = [0.8500 0.3250 0.0980];
b(6).FaceColor = [0.8500 0.3250 0.0980];
hold on
[ngroups, nbars] = size(concs');
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, concs(i,:), error_low(i,:), error_high(i,:), 'k', 'linestyle', 'none');
end
axis([-inf inf 0 100])
xticklabels(words)
legend([b(1), b(2), b(4)],'Container','MGS1 pre testing', 'MGS1 post testing')
ylabel('Percentage Composition $[%]$')
xlabel('Element')