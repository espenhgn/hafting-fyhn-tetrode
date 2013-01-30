function out = is_my_logical(n)
%IS_MY_LOGICAL    checks if value is logical. 
%
% IS_MY_LOGICAL(n) returns true if *n* is boolean or integer.
%
% Born in version 0.2, Jan 2013, Eivind Skjønsberg Norheim

out = 0;
if islogical(n)
    out=1;
% integers are doubles if not handled explicitly, so check if the double is
% integer
elseif (n-int8(n))==0
    out=1;
end
