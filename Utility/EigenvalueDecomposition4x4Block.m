% The function for eigenvalue decomposition of 4x4 Hermitian Matrix
function varargout = EigenvalueDecomposition4x4Block(RBlock, frqLen, timeLen)
rowLen = size(RBlock,1);
colLen = size(RBlock,2);
%% Find eigenvalues
indxij = cell(4);
offsetVect = 0:16:16*(frqLen*timeLen-1);
for i = 1:4
    for j = 1:4
      indxij{i,j} = i + (j-1)*4 + offsetVect;  
    end
end

coef4 = ones(size(indxij{1,1}));
coef3 = (-RBlock(indxij{1,1}) - RBlock(indxij{2,2}) - RBlock(indxij{3,3}) - RBlock(indxij{4,4}));
coef2 = (RBlock(indxij{1,1}).*(RBlock(indxij{2,2}) + RBlock(indxij{3,3}) + RBlock(indxij{4,4})) + RBlock(indxij{2,2}).*(RBlock(indxij{3,3}) + RBlock(indxij{4,4})) + RBlock(indxij{3,3}).*RBlock(indxij{4,4}) - abs(RBlock(indxij{1,2})).^2 - abs(RBlock(indxij{1,3})).^2 - abs(RBlock(indxij{1,4})).^2 - abs(RBlock(indxij{2,3})).^2 - abs(RBlock(indxij{2,4})).^2 - abs(RBlock(indxij{3,4})).^2);
coef1 = ((abs(RBlock(indxij{1,3})).^2).*RBlock(indxij{2,2}) - RBlock(indxij{1,4}).*RBlock(indxij{2,1}).*RBlock(indxij{4,2}) - RBlock(indxij{1,4}).*RBlock(indxij{3,1}).*RBlock(indxij{4,3}) - RBlock(indxij{2,4}).*RBlock(indxij{3,2}).*RBlock(indxij{4,3}) - RBlock(indxij{1,1}).*RBlock(indxij{2,2}).*RBlock(indxij{3,3}) - RBlock(indxij{1,1}).*RBlock(indxij{2,2}).*RBlock(indxij{4,4}) - RBlock(indxij{1,1}).*RBlock(indxij{3,3}).*RBlock(indxij{4,4}) - RBlock(indxij{2,2}).*RBlock(indxij{3,3}).*RBlock(indxij{4,4}) - RBlock(indxij{1,2}).*RBlock(indxij{2,3}).*RBlock(indxij{3,1}) - RBlock(indxij{1,3}).*RBlock(indxij{2,1}).*RBlock(indxij{3,2}) - RBlock(indxij{1,2}).*RBlock(indxij{2,4}).*RBlock(indxij{4,1}) + RBlock(indxij{1,4}).*RBlock(indxij{2,2}).*RBlock(indxij{4,1}) + RBlock(indxij{1,1}).*RBlock(indxij{2,3}).*RBlock(indxij{3,2}) + RBlock(indxij{1,2}).*RBlock(indxij{3,3}).*RBlock(indxij{2,1}) + RBlock(indxij{1,1}).*RBlock(indxij{2,4}).*RBlock(indxij{4,2}) - RBlock(indxij{1,3}).*RBlock(indxij{3,4}).*RBlock(indxij{4,1}) + RBlock(indxij{1,4}).*RBlock(indxij{3,3}).*RBlock(indxij{4,1}) + RBlock(indxij{1,2}).*RBlock(indxij{4,4}).*RBlock(indxij{2,1}) + RBlock(indxij{1,3}).*RBlock(indxij{4,4}).*RBlock(indxij{3,1}) + RBlock(indxij{1,1}).*RBlock(indxij{3,4}).*RBlock(indxij{4,3}) - RBlock(indxij{2,3}).*RBlock(indxij{3,4}).*RBlock(indxij{4,2}) + RBlock(indxij{2,4}).*RBlock(indxij{3,3}).*RBlock(indxij{4,2}) + RBlock(indxij{2,2}).*RBlock(indxij{3,4}).*RBlock(indxij{4,3}) + RBlock(indxij{2,3}).*RBlock(indxij{4,4}).*RBlock(indxij{3,2}));
coef0 = - RBlock(indxij{1,4}).*RBlock(indxij{2,1}).*RBlock(indxij{3,2}).*RBlock(indxij{4,3}) - RBlock(indxij{1,2}).*RBlock(indxij{2,3}).*RBlock(indxij{3,4}).*RBlock(indxij{4,1}) + RBlock(indxij{1,2}).*RBlock(indxij{2,4}).*RBlock(indxij{3,3}).*RBlock(indxij{4,1}) + RBlock(indxij{1,3}).*RBlock(indxij{2,2}).*RBlock(indxij{3,4}).*RBlock(indxij{4,1}) - RBlock(indxij{1,4}).*RBlock(indxij{2,2}).*RBlock(indxij{3,3}).*RBlock(indxij{4,1}) + RBlock(indxij{1,1}).*RBlock(indxij{2,3}).*RBlock(indxij{3,4}).*RBlock(indxij{4,2}) - RBlock(indxij{1,1}).*RBlock(indxij{2,4}).*RBlock(indxij{3,3}).*RBlock(indxij{4,2}) + RBlock(indxij{1,2}).*RBlock(indxij{2,3}).*RBlock(indxij{4,4}).*RBlock(indxij{3,1}) - RBlock(indxij{1,3}).*RBlock(indxij{2,2}).*RBlock(indxij{4,4}).*RBlock(indxij{3,1}) - RBlock(indxij{1,1}).*RBlock(indxij{2,2}).*RBlock(indxij{3,4}).*RBlock(indxij{4,3}) - RBlock(indxij{1,1}).*RBlock(indxij{2,3}).*RBlock(indxij{4,4}).*RBlock(indxij{3,2}) - RBlock(indxij{1,2}).*RBlock(indxij{3,3}).*RBlock(indxij{4,4}).*RBlock(indxij{2,1}) + RBlock(indxij{1,3}).*RBlock(indxij{2,4}).*RBlock(indxij{3,1}).*RBlock(indxij{4,2}) - RBlock(indxij{1,3}).*RBlock(indxij{2,4}).*RBlock(indxij{4,1}).*RBlock(indxij{3,2}) - RBlock(indxij{1,4}).*RBlock(indxij{2,3}).*RBlock(indxij{3,1}).*RBlock(indxij{4,2}) + RBlock(indxij{1,4}).*RBlock(indxij{2,3}).*RBlock(indxij{4,1}).*RBlock(indxij{3,2}) - RBlock(indxij{1,2}).*RBlock(indxij{2,4}).*RBlock(indxij{3,1}).*RBlock(indxij{4,3}) - RBlock(indxij{1,3}).*RBlock(indxij{3,4}).*RBlock(indxij{2,1}).*RBlock(indxij{4,2}) + RBlock(indxij{1,4}).*RBlock(indxij{2,2}).*RBlock(indxij{3,1}).*RBlock(indxij{4,3}) + RBlock(indxij{1,4}).*RBlock(indxij{3,3}).*RBlock(indxij{2,1}).*RBlock(indxij{4,2}) + RBlock(indxij{1,1}).*RBlock(indxij{2,4}).*RBlock(indxij{3,2}).*RBlock(indxij{4,3}) + RBlock(indxij{1,2}).*RBlock(indxij{3,4}).*RBlock(indxij{2,1}).*RBlock(indxij{4,3}) + RBlock(indxij{1,3}).*RBlock(indxij{4,4}).*RBlock(indxij{2,1}).*RBlock(indxij{3,2}) + RBlock(indxij{1,1}).*RBlock(indxij{2,2}).*RBlock(indxij{3,3}).*RBlock(indxij{4,4});

indxQuadratic = ones(4,1)*(abs(coef0) > 1e-10);
indxCubic = ones(4,1)*(abs(coef0) <= 1e-10);
%% Quadratic Equation Solver
D0 = coef2.^2 - 3*coef3.*coef1 + 12*coef4.*coef0;
D1 = 2*coef2.^3 - 9*coef3.*coef2.*coef1 + 27*(coef3.^2).*coef0 + 27*coef4.*(coef1.^2) - 72*coef4.*coef2.*coef0;
D = (D1.^2 - 4*D0.^3) / (-27);
p = (8*coef4.*coef2 - 3*coef3.^2) ./ (8*coef4.^2);
q = (coef3.^3 - 4*coef4.*coef3.*coef2 + 8*(coef4.^2).*coef1) ./ (8*coef4.^3);

Q = D1.^(1/3);
S1 = 0.5*sqrt(-2/3*p + 1./(3*coef4).*(Q + D0./Q));
phi = acos(D1./(2*sqrt(D0.^3)));
S2 = 0.5*sqrt(-2/3*p + 2./(3*coef4).*sqrt(D0).*cos(phi/3));
Q = ((D1.^2 + sqrt(D1.^2 - 4*D0.^3))/2).^(1/3);
S3 = 0.5*sqrt(-2/3*p + 1./(3*coef4).*(Q + D0./Q));

indx1 = abs(D) ~= 0 & abs(D0) == 0;
indx2 = abs(D) > 0 & not(indx1);
indx3 = not(indx1) & not(indx2);
S = S1.*indx1 + S2.*indx2 + S3.*indx3;

eigValQuadratic = [-coef3./(4*coef4) - S + 0.5*sqrt(-4*S.^2 - 2*p + q./S);...
    -coef3./(4*coef4) - S - 0.5*sqrt(-4*S.^2 - 2*p + q./S);...
    -coef3./(4*coef4) + S + 0.5*sqrt(-4*S.^2 - 2*p - q./S);...
    -coef3./(4*coef4) + S - 0.5*sqrt(-4*S.^2 - 2*p - q./S)];
%% Cubic Equation Solver
zeroVect = zeros(size(indxij{1,1}));
D0 = coef3.^2 - 3*coef4.*coef2;
D1 = 2*coef3.^3 - 9*coef4.*coef3.*coef2 + 27*coef4.^2.*coef1;
D = (D1.^2 - 4*D0.^3) ./ (27*coef4.^2);
indx1 = abs(D) ~= 0 & abs(D0) == 0;
C = D1.^(1/3);
u1 = 1;
u2 = (-1+sqrt(-3))/2;
u3 = (-1-sqrt(-3))/2;
temp1 = [-1./(3*coef4).*(coef3 + u1*C);-1./(3*coef4).*(coef3 + u2*C);-1./(3*coef4).*(coef3 + u3*C);zeroVect].*[indx1;indx1;indx1;indx1];
indx2 = abs(D) == 0 & abs(D0) == 0 & not(indx1);
x123 = -coef3./(3*coef4);
temp2 = [x123;x123;x123;zeroVect].*[indx2;indx2;indx2;indx2];
indx3 = abs(D) == 0 & abs(D0) ~= 0 & not(indx1) & not(indx2);
x12 = (9*coef4.*coef1 - coef3.*coef2) ./ (2*D0);
x3 = (4*coef4.*coef3.*coef2 - 9*coef4.^2.*coef1 - coef3.^3) ./ (coef4.*D0);
temp3 = [x12;x12;x3;zeroVect].*[indx3;indx3;indx3;indx3];
indx4 = not(indx1) & not(indx2) & not(indx3);
C = ((D1 + sqrt(D1.^2 - 4*D0.^3))/2).^(1/3);
u1 = 1;
u2 = (-1+sqrt(-3))/2;
u3 = (-1-sqrt(-3))/2;
temp4 = [-1./(3*coef4).*(coef3 + u1*C + D0./(u1*C));-1./(3*coef4).*(coef3 + u2*C + D0./(u2*C));-1./(3*coef4).*(coef3 + u3*C + D0./(u3*C));zeroVect].*[indx4;indx4;indx4;indx4];
eigValCubic = temp1 + temp2 + temp3 + temp4;
%% Final Solution
eigVal = sort(eigValQuadratic.*indxQuadratic + eigValCubic.*indxCubic,'descend');
%% Find Eigen Vectors
eigVect = zeros(4,4*size(eigVal,2));
for i = 1:4
    RI = RBlock - kron(eigVal(i,:),eye(4));
    alpha2 = RI(indxij{1,2}).*RI(indxij{2,1}) - RI(indxij{2,2}).*RI(indxij{1,1});
    alpha3 = RI(indxij{1,3}).*RI(indxij{2,1}) - RI(indxij{2,3}).*RI(indxij{1,1});
    alpha4 = RI(indxij{1,4}).*RI(indxij{2,1}) - RI(indxij{2,4}).*RI(indxij{1,1});
    beta2 = RI(indxij{3,2}).*RI(indxij{4,1}) - RI(indxij{4,2}).*RI(indxij{3,1});
    beta3 = RI(indxij{3,3}).*RI(indxij{4,1}) - RI(indxij{4,3}).*RI(indxij{3,1});
    beta4 = RI(indxij{3,4}).*RI(indxij{4,1}) - RI(indxij{4,4}).*RI(indxij{3,1});
    b4 = ones(size(indxij{1,1}));
    b3 = (beta4.*alpha2 - alpha4.*beta2) ./ (alpha3.*beta2 - beta3.*alpha2);
    b3(isinf(b3) | isnan(b3)) = 1;
    b2 = -(alpha3.*b3 + alpha4 + beta3.*b3 + beta4) ./ (alpha2 + beta2);
    b2(isinf(b2) | isnan(b2)) = 1;
    b1 = -((RI(indxij{1,2}) + RI(indxij{2,2}) + RI(indxij{3,2}) + RI(indxij{4,2})).*b2 + ...
           (RI(indxij{1,3}) + RI(indxij{2,3}) + RI(indxij{3,3}) + RI(indxij{4,3})).*b3 + ...
           (RI(indxij{1,4}) + RI(indxij{2,4}) + RI(indxij{3,4}) + RI(indxij{4,4})).*b4) ./ (RI(indxij{1,1}) + RI(indxij{2,1}) + RI(indxij{3,1}) + RI(indxij{4,1}));
    vectLength = sqrt(abs(b1).^2 + abs(b2).^2 + abs(b3).^2 + abs(b4).^2);
    eigVect(:,i:4:end) = [b1;b2;b3;b4] ./ [vectLength;vectLength;vectLength;vectLength];
end
%% OUTPUTS
varargout{1} = 0;
varargout{2} = eigVect;
varargout{3} = eigVal;