temp=[70.755,   74.228,   74.935,   75.801,   76.031,   76.954,   77.470,   78.698,   80.314,   81.492];
figure
cdfplot(temp);
xlabel('Value');
ylabel('Cumulative Distribution Function');
title('Empirical Cumulative Distribution Function');
grid on;
