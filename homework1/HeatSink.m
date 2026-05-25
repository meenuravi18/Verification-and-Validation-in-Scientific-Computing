clear all
clc


sampling_number=1000;
[Tbase_list] = monte_carlo_sampling(sampling_number);
disp(Tbase_list);

[Tbase_list] = latin_hypercube_sampling(sampling_number);
disp(Tbase_list);

figure
cdfplot(Tbase_list);
xlabel('Value of base temperature');
ylabel('Cumulative Distribution Function');
title('Empirical Cumulative Distribution Function');
grid on;


% Finding the probability that the CPU base temperature
%  will exceed the peak allowable temperature of 85 degrees
i=0;
for element = Tbase_list
    if element>85
      i=i+1
    end
end
disp(i/sampling_number)