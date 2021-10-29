% main_amgrams.m

clc

% Load data into workspace
analyze_setup;


%% Prepare
figure(1);
t = data.spec_st.t;
f = data.spec_st.f;
k = 3;
S = data.spec_st.Sft{k};
spec.plot_spectrogram(t, f, S, 'fignum', 1);
title(sprintf('Label: %s', data.spec_st.labels{k}));



%% Create the AM-grams










