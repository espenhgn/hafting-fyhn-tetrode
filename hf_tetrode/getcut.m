function [clust cut_for] = getcut(cutfile)
fid = fopen(cutfile, 'rt');
clust = [];
while ~feof(fid)
    string = fgetl(fid);
    if ~isempty(string)
        if (string(1) == 'E')
            tmp = regexp(regexp(string, ' ', 'split'), ',', 'split');
            cut_for = tmp{2};
            break;
        end
    end
end
while ~feof(fid)
  string = fgetl(fid);
  if ~isempty(string)
     content = sscanf(string,'%u')';
     clust = [clust content]; %#ok
  end
end
fclose(fid);
clust = clust';