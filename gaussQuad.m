function varargout = gaussQuad(quadfunc, minval1, maxval1, rank)
if nargin < 4
    rank = 6;
end
[XG, WGT] = gauleg(minval1, maxval1, rank);
m = nargout;
temp = cell(1, m);
varargout = cell(1, m);

% firstly for ii = 1
[temp{:}] = quadfunc(XG(1));
for ll = 1:m
    varargout{ll} = WGT(1).*temp{ll};
end
for ii = 2:rank
    [temp{:}] = quadfunc(XG(ii));
    for ll = 1:m
        varargout{ll} = varargout{ll} + WGT(ii).*temp{ll};
    end
end