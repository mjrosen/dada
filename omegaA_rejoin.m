function out = omegaA_rejoin(folder)

% produces a list of omegaA thresholds at which each cluster would be
% reabsorbed into some other cluster. to compute this, p-values for
% membership of the reads in the error-free family to each other cluster
% are computed and the largest one is kept.

%takes as its input a folder of .mat files as produced by dada

files = dir([folder '/*.mat']);
files = files(~strcmp({files.name},'ERR.mat'));
A = cell(1,length(files));

for f = 1:length(files)
    load([folder '/' files(f).name]);
    %A is a cell array of the \Omega_a threshold that would be required for 
    %each cluster to be reabsorbed into its nearest neighbor
    A{f} = -inf(1,length(bin));
    for i = 1:length(bin)
        %find the error-free family in this cluster
        for j = 1:length(bin(i).fam)
            if isempty(bin(i).fam(j).raw(1).subPos{i}) %error-free fam
                r = bin(i).fam(j).r;
                %find largest p-value for membership to another cluster
                for k = setdiff(1:length(bin),i) %check each other cluster
                    lambda = bin(i).fam(j).raw(1).lambda(k);
                    %first compute normalization of abundance p-value
                    norm = 1 - poisspdf(0,bin(k).R*lambda);
                    if norm == 0 %lambda*R too small for poisspdf
                        norm = lambda*bin(k).R; %better approx
                    end
                    %compute numerator: the tail of the poisson
                    %distribution
                    tail = 1 - poisscdf(r-1,bin(k).R*lambda);
                    if tail == 0
                        tail = poisspdf(r,bin(k).R*lambda); %use first term
                    end
                    p = tail / norm;
                    if p > A{f}(i)
                        A{f}(i) = p;
                    end
                end
            end
        end        
    end
    %Bonferroni correct the p-values in this cluster by the total number of
    %clusters
    A{f} = A{f} * length([bin.fam]); 
end
out = [A{:}];
end