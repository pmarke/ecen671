%% Finding filter with 5 coefficients

[e_5,h_5] = myfilter(5);
h_5

%% Finding filter coefficients with 6 coefficients

[e_6,h_6] = myfilter(6);
h_6

%% Finding filter coefficients with 4 coefficients

[e_4,h_4] = myfilter(4);
h_4


figure(1),clf

hold on
plot(e_5)
plot(e_6)
plot(e_4)

legend('5 coeff','6 coeff', '4 coeff');

%%
% The filters that try to match the impulse response h = [1,2,3,4,5] can do
% a good job at estimating it if they have at least 5 coefficients.
% Otherwise they fail to do it. 
