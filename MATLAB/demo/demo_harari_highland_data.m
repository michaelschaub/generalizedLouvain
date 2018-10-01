

% load data and labels and assign name to store outputs
load('./gama.mat');
labels = {    
    'GAVEV'
    'KOTUN'
    'OVE'
    'ALIKA'
    'NAGAM'
    'GAHUK'
    'MASIL'
    'UKUDZ'
    'NOTOH'
    'KOHIK'
    'GEHAM'
    'ASARO'
    'UHETO'
    'SEUVE'
    'NAGAD'
    'GAMA'};
outname = 'harari_gama';


%assemble joint signed adjacency matrix and time parameters
A_tot = A_pos - A_neg;
time = logspace(-3,2,120);

% run louvain optimization with signed Laplacian quality function
partition_stability(A_tot,time,'Laplacian','louvain_signedLap',...
    'v','out',outname,'t','p');
   
% postprocess files, write out results and plot
stability_postprocess(['Stability_' outname '.mat'],A_tot)
write_stab_results(['Stability_' outname '_PP.mat'],outname);
load(['Stability_' outname '_PP.mat']);
script_plot_stability_one_panel

% print partition at Markov time 0.165

pvec = C_new(:,50);
for i =1:max(pvec)
    labels(pvec==i)
end