addpath('/usr/share/itpp/');

itload('P_min_far.it');
itload('P_min_near.it');
itload('rate_pairs_due_to_spc.it');

rates=[0.50 0.67 0.75 0.83 1.00 1.33 1.50 1.67 2.00 2.67 3.00 3.33];
highest_p_min_near=P_min_near(end);
highest_p_min_far=P_min_far(max(find(P_min_far<highest_p_min_near)));

x_intercept=rates(end);
y_intercept=rates(find(P_min_far==highest_p_min_far));

slope=-(y_intercept/x_intercept);
x=0:0.01:rates(end);
y=slope*x+y_intercept;
plot(x,y);
xlabel('Near User Rates (in bits/sec/Hz)');
ylabel('Far User Rates (in bits/sec/Hz)');
grid on

hold on;

plot(rate_pairs_due_to_spc(1,:),rate_pairs_due_to_spc(2,:),'r+');

[index_x index_y]=find(rate_pairs_due_to_spc>0);

for i=1:length(index_x)
    rate_pairs_due_to_spc_convhull(index_x(i),index_y(i))=rate_pairs_due_to_spc(index_x(i),index_y(i));
end

rate_pairs_due_to_spc_convhull=[rate_pairs_due_to_spc_convhull [x_intercept;0] [0;y_intercept] [0;0]];

convex_hull=convhull(rate_pairs_due_to_spc_convhull(1,:),rate_pairs_due_to_spc_convhull(2,:));

plot(rate_pairs_due_to_spc_convhull(1,convex_hull),rate_pairs_due_to_spc_convhull(2,convex_hull),'r');

axis([0,4,0,4]);
legend('Line obtainable by just time-sharing','Rate-pairs due to spc','Convex hull of rate-pairs of spc');
hold off;
print -color rate_region.eps
%print(1, "rate_region.png");
