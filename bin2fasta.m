function bin2fasta(folder,varargin)
p = inputParser;
p.addRequired('folder',@ischar);
p.addOptional('fasta',true,@(x)islogical(x) && isscalar(x));
p.parse(folder,varargin{:});
folder = p.Results.folder;
fasta = p.Results.fasta;
%if 'fasta' is true, output a single fasta file with seqs from all files 
%called 'all.fasta'. otherwise output a fasta for each file.
files = dir([folder '/*.mat']);
files(strcmp('ERR.mat',{files.name})) = [];
reads = [];
seqs = {};
i = 1;
for f = 1:length(files)
    if ~fasta %need to refresh for every file
        reads = [];
        seqs = {};
        i = 1;
    end
    s = load([folder '/' files(f).name]);
    for r = 1:length(s.reals)
        reads(i) = s.reads(r); %#ok<AGROW>
        seqs{i} = s.reals{r}; %#ok<AGROW>
        i = i + 1;
    end
    if ~fasta %output sorted files in "reads /t seq" format (a .uniques)
        [reads,I] = sort(reads,'descend');
        seqs = seqs(I);
        fid = fopen([folder '/' ...
            regexprep(files(f).name,'mat','uniques')],'w+');
        for i = 1:length(reads)
            fprintf(fid,'%d\t%s\n',reads(i),seqs{i});
        end
        fclose(fid);
    end
end
if fasta %output a single sorted .fasta file
    [reads,I] = sort(reads,'descend');
    seqs = seqs(I);
    out = struct('Header',{},'Sequence',{});
    for i = 1:length(reads)
        out(i).Header = [num2str(i) '_' num2str(reads(i))];
        out(i).Sequence = seqs{i};
    end
    fastawrite([folder '/all.fasta'], out);
end
