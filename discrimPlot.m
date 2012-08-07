function discrimPlot(omegaA,omegaR,T,rho,bc,Nfam)

%parameters: significance levels omegaA and omegaR, 4x4 context-independent
%transition matrix T, the number of reads of a genotype rho, the basecount
%of the genotype template (either a 4x1 or a 1x4 vector) bc, and the number
%of families of the genotype Nfam. These can all be acquired easily from
%data after clustering by dada.

%make sure bc is a 4x1:
bc = reshape(bc,4,1);

T1 = diag(T);
T2 = [[T(1,2) T(1,3) T(1,4)] / T(1,1); ...
    [T(2,1) T(2,3) T(2,4)] / T(2,2); ...
    [T(3,1) T(3,2) T(3,4)] / T(3,3); ...
    [T(4,1) T(4,2) T(4,3)] / T(4,4)];
self = prod(T1.^bc);
maxT = max(setdiff(T(:),diag(T))); %the largest error rate

%first entry of lambda should be the maximum possible actual error rate.
%log10 is taken because logspace's arguments are of this form
lambda = logspace(log10(max(T2(:))*self),-20,2000);

% construct omegaA discrimination line: A(lambda) is the minimum number of
% reads of error of probability lambda that are declared significant
A = zeros(size(lambda));
for i = 1:length(lambda)
    r = max(1,floor(lambda(i)*rho));
    while(1)
        tail = 1 - poisscdf(r-1,lambda(i) * rho);
        if tail == 0
            tail = poisspdf(r,lambda(i) * rho);
        end
        norm = 1 - poisspdf(0,rho*lambda(i));
        if norm == 0
            norm = lambda(i) * rho;
        end
        p = tail / norm;
        if p * Nfam < omegaA %Bonferroni correct
            A(i) = r;
            break;
        else
            r = r + 1;
        end
    end
end

% construct omegaR discrimination line: lambdaSig is the upper bound on
% lambda that are declared significant. this is all taken from the main
% dada source code, so see there for better commenting

dmax = 10; %maximum number of errors to consider
LV = 1;
mV = 1;
multi = cell(1,dmax);
bin = zeros(dmax,4);
for d = 1:dmax
    multi{d} = zeros(d+1,d+1,d+1);
    for i = ceil(d/3):d
        for j = ceil((d-i)/2):min(i,d-i)
            k = d - i - j;
            n = factorial(d) / prod(factorial([i j k]));
            multi{d}(i+1,j+1,k+1) = n;
            multi{d}(i+1,k+1,j+1) = n;
            multi{d}(j+1,i+1,k+1) = n;
            multi{d}(j+1,k+1,i+1) = n;
            multi{d}(k+1,i+1,j+1) = n;
            multi{d}(k+1,j+1,i+1) = n;
        end
    end
    for i = 1:4
        if d <= bc(i)
            bin(d,i) = nchoosek(bc(i),d);
        else
            bin(d,i) = 0;
        end
    end
    L = zeros(nchoosek(d+11,11),1);
    m = ones(nchoosek(d+11,11),1);
    M = zeros(4,3);
    M(1) = d;
    i = 1;
    while(1)
        bases = sum(M,2); %bases involved in errors
        L(i) = prod(prod(T2.^M)); %relative error rate of A
        for j = 1:4
            if bases(j) > 0
                m(i) = m(i) * bin(bases(j),j);
                m(i) = m(i) * multi{bases(j)}(...
                    M(j,1)+1,M(j,2)+1,M(j,3)+1);
            end
        end
        if M(end) == d %d T->G errors is termination condition
            break;
        elseif M(end) == 0 %no T->G errs: advance final error
            x = find(M,1,'last');
            M(x) = M(x) - 1;
            M(x+1) = M(x+1) + 1;
        else
            n = M(end);
            M(end) = 0;
            x = find(M,1,'last');
            M(x) = M(x) - 1;
            M(x+1) = n + 1;
        end
        i = i + 1; %advance index
    end
    LV = [LV;L]; %#ok<AGROW>
    mV = [mV;m]; %#ok<AGROW>
end
[LV,I]=sort(LV,'descend');
mV = mV(I);
p = 1 - cumsum(self * LV.*mV);
%find first p-value that would be significant
lambdaSig = lambda(find(lambda<LV(find(p<(omegaR/rho),1)),1));

%plot
plot(log(lambda)/log(maxT),A,'color','k','linewidth',2);hold on;
line([log(lambdaSig)/log(maxT) log(lambdaSig)/log(maxT)],...
    [0 A(1)],'color','k','linewidth',2);
set(gca,'fontsize',16,'fontname','helvetica','fontweight','bold');
xlabel('effective Hamming distance','fontsize',16,...
    'fontname','helvetica','fontweight','bold')
ylabel('reads','fontsize',16,'fontname','helvetica','fontweight','bold')
xlim([0 log(lambdaSig)/log(maxT)+2]);
ylim([0 A(1)]);