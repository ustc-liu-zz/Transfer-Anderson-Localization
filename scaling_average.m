%% compute the slope fo t=1 and T belongs to [1,1.3]
k_slope = zeros(31,10);
for step = 1:10
    scaling_exponent;
    k_slope(:,step) = -p_slope;
end
T = 1:0.01:1.3;
k_mean = mean(k_slope,2);
k_err = std(k_slope,1,2);
errorbar(T,k_mean,k_err,'k-','linewidth',1.5)
xlim([1 1.3])
xlabel('T')
ylabel('slope')