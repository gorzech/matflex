function varargout = tripleGaussQuad(quadfunc, minval1, maxval1, ...
    minval2, maxval2, minval3, maxval3, rank, rank2, rank3)
if nargin < 8
    rank = 6;
    rank2 = 6;
    rank3 = 6;
elseif nargin < 9
    rank2 = rank;
    rank3 = rank;
end
[XG, WGT] = gauleg(minval1, maxval1, rank);
[XG2, WGT2] = gauleg(minval2, maxval2, rank2);
[XG3, WGT3] = gauleg(minval3, maxval3, rank3);
m = nargout;
temp = cell(1, m);
varargout = cell(1, m);
varargout(1, :) = {0};

for ii = 1:rank
    for jj = 1:rank2
        for kk = 1:rank3
            [temp{:}] = quadfunc(XG(ii), XG2(jj), XG3(kk));
            coeff = WGT(ii).*WGT2(jj).*WGT3(kk);
            for ll = 1:m
                varargout{ll} = varargout{ll} + coeff.*temp{ll};
            end
        end
    end
end