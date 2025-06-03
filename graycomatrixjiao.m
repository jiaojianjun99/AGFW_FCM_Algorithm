function [GLCMS,SI] = graycomatrixjiao(varargin)%%%%进来就调用一个函数
%GRAYCOMATRIX Create gray-level co-occurrence matrix.
%   GLCMS = GRAYCOMATRIX(I) analyzes pairs of horizontally adjacent pixels
%   in a scaled version of I.  If I is a binary image, it is scaled to 2
%   levels. If I is an intensity image, it is scaled to 8 levels. In this
%   case, there are 8 x 8 = 64 possible ordered combinations of values for
%   each pixel pair. GRAYCOMATRIX accumulates the total occurrence of each
%   such combination, producing a 8-by-8 output array, GLCMS. The row and
%   column subscripts in GLCMS correspond respectively to the first and
%   second (scaled) pixel-pair values.
%
%   GLCMS = GRAYCOMATRIX(I,PARAM1,VALUE1,PARAM2,VALUE2,...) returns one or
%   more gray-level co-occurrence matrices, depending on the values of the
%   optional parameter/value pairs. Parameter names can be abbreviated, and
%   case does not matter.
%   
%   Parameters include:
%  
%   'Offset'         A p-by-2 array of offsets specifying the distance 
%                    between the pixel-of-interest and its neighbor. Each 
%                    row in the array is a two-element vector, 
%                    [ROW_OFFSET COL_OFFSET], that specifies the 
%                    relationship, or 'Offset', between a pair of pixels. 
%                    ROW_OFFSET is the number of rows between the 
%                    pixel-of-interest and its neighbor.  COL_OFFSET is the
%                    number of columns between the pixel-of-interest and 
%                    its neighbor. For example, if you want the number of
%                    occurrences where the pixel of interest is one pixel 
%                    to the left of its neighbor, then 
%                    [ROW_OFFSET COL_OFFSET] is [0 1].
%  
%                    Because this offset is often expressed as an angle, 
%                    the following table lists the offset values that 
%                    specify common angles, given the pixel distance D.
%                            
%                    Angle     OFFSET
%                    -----     ------  
%                    0         [0 D]   
%                    45        [-D D]
%                    90        [-D 0]
%                    135       [-D -D]  
%  
%                    ROW_OFFSET and COL_OFFSET must be integers. 
%
%                    Default: [0 1]
%            
%   'NumLevels'      An integer specifying the number of gray levels to use
%                    when scaling the grayscale values in I. For example,
%                    if 'NumLevels' is 8, GRAYCOMATRIX scales the values in
%                    I so they are integers between 1 and 8.  The number of
%                    gray levels determines the size of the gray-level
%                    co-occurrence matrix (GLCM).
%
%                    'NumLevels' must be an integer. 'NumLevels' must be 2
%                    if I is logical.
%  
%                    Default: 8 for numeric
%                             2 for logical
%   
%   'GrayLimits'     A two-element vector, [LOW HIGH], that specifies how 
%                    the values in I are scaled into gray levels. If N is
%                    the number of gray levels (see parameter 'NumLevels')
%                    to use for scaling, the range [LOW HIGH] is divided
%                    into N equal width bins and values in a bin get mapped
%                    to a single gray level. Grayscale values less than or
%                    equal to LOW are scaled to 1. Grayscale values greater
%                    than or equal to HIGH are scaled to NumLevels. If
%                    'GrayLimits' is set to [], GRAYCOMATRIX uses the
%                    minimum and maximum grayscale values in I as limits,
%                    [min(I(:)) max(I(:))].
%  
%                    Default: the LOW and HIGH values specified by the
%                    class, e.g., [LOW HIGH] is [0 1] if I is double and
%                    [-32768 32767] if I is int16.
%
%   'Symmetric'      A Boolean that creates a GLCM where the ordering of
%                    values in the pixel pairs is not considered. For
%                    example, when calculating the number of times the
%                    value 1 is adjacent to the value 2, GRAYCOMATRIX
%                    counts both 1,2 and 2,1 pairings, if 'Symmetric' is
%                    set to true. When 'Symmetric' is set to false,
%                    GRAYCOMATRIX only counts 1,2 or 2,1, depending on the
%                    value of 'offset'. The GLCM created in this way is
%                    symmetric across its diagonal, and is equivalent to
%                    the GLCM described by Haralick (1973).
%
%                    The GLCM produced by the following syntax, 
%
%                    graycomatrix(I, 'offset', [0 1], 'Symmetric', true)
%
%                    is equivalent to the sum of the two GLCMs produced by
%                    these statements.
%
%                    graycomatrix(I, 'offset', [0 1], 'Symmetric', false) 
%                    graycomatrix(I, 'offset', [0 -1], 'Symmetric', false) 
%
%                    Default: false
%  
%  
%   [GLCMS,SI] = GRAYCOMATRIX(...) returns the scaled image used to
%   calculate GLCM. The values in SI are between 1 and 'NumLevels'.
%
%   Class Support
%   -------------             
%   I can be numeric or logical.  I must be 2D, real, and nonsparse. SI is
%   a double matrix having the same size as I.  GLCMS is an
%   'NumLevels'-by-'NumLevels'-by-P double array where P is the number of
%   offsets in OFFSET.
%  
%   Notes
%   -----
%   Another name for a gray-level co-occurrence matrix is a gray-level
%   spatial dependence matrix.
%
%   GRAYCOMATRIX ignores pixels pairs if either of their values is NaN. It
%   also replaces Inf with the value 'NumLevels' and -Inf with the value 1.
%
%   GRAYCOMATRIX ignores border pixels, if the corresponding neighbors
%   defined by 'Offset' fall outside the image boundaries.
%
%   References
%   ----------
%   Haralick, R.M., K. Shanmugan, and I. Dinstein, "Textural Features for
%   Image Classification", IEEE Transactions on Systems, Man, and
%   Cybernetics, Vol. SMC-3, 1973, pp. 610-621.
%
%   Haralick, R.M., and L.G. Shapiro. Computer and Robot Vision: Vol. 1,
%   Addison-Wesley, 1992, p. 459.
%
%   Example 1
%   ---------      
%   Calculate the gray-level co-occurrence matrix (GLCM) and return the
%   scaled version of the image, SI, used by GRAYCOMATRIX to generate the
%   GLCM.
%
%        I = [1 1 5 6 8 8;2 3 5 7 0 2; 0 2 3 5 6 7];
%       [GLCMS,SI] = graycomatrix(I,'NumLevels',9,'G',[])
%     
%   Example 2
%   ---------  
%   Calculate the gray-level co-occurrence matrix for a grayscale image.
%  
%       I = imread('circuit.tif');
%       GLCMS = graycomatrix(I,'Offset',[2 0])
%
%   Example 3
%   ---------
%   Calculate gray-level co-occurrences matrices for a grayscale image
%   using four different offsets.
%  
%       I = imread('cell.tif');
%       offsets = [0 1;-1 1;-1 0;-1 -1];
%       [GLCMS,SI] = graycomatrix(I,'Of',offsets); 
%
%   Example 4
%   ---------  
%   Calculate the symmetric gray-level co-occurrence matrix (the Haralick
%   definition) for a grayscale image.
%  
%       I = imread('circuit.tif');
%       GLCMS = graycomatrix(I,'Offset',[2 0],'Symmetric', true)
%  
%   See also GRAYCOPROPS.
  
%   Copyright 1993-2012 The MathWorks, Inc.

[I, Offset, NL, GL, makeSymmetric] = ParseInputs(varargin{:});
%%NL灰度级数'NumLevels'
%%GL，A two-element vector, [LOW HIGH], that specifies how 
%  the values in I are scaled into gray levels.

% Scale I so that it contains integers between 1 and NL.
% I=I';
if GL(2) == GL(1)%%如果两个相等，或为空'GrayLimits' is set to []，SI按原始图像的灰度级，不缩放
  SI = ones(size(I));%SI按原始图像的灰度级，不缩放
else 
    slope = NL / (GL(2) - GL(1));
    intercept = 1 - (slope*(GL(1)));
    SI = floor(imlincomb(slope,I,intercept,'double'));%%线性缩放I到SI
end

% Clip values if user had a value that is outside of the range, e.g.,
% double image = [0 .5 2;0 1 1]; 2 is outside of [0,1]. The order of the
% following lines matters in the event that NL = 0.
SI(SI > NL) = NL;
% SI(SI < 1) = 0;
SI(SI < GL(1)) = 0;

numOffsets = size(Offset,1);%%￥￥￥￥￥￥￥￥￥￥%%Offset[n,m],取其测纹理数量，即行数￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥

if NL ~= 0%%%%%非零

  % Create vectors of row and column subscripts for every pixel and its
  % neighbor.
  s = size(I);%%%%%%%%%%%%%%%要变！！！！！！！！！！！！！！！！！！！！


% [r,c] = meshgrid(1:s(1),1:s(2));%%s(1)元素索引。好函数,如果是不规则图形，取坐标时在这里进行变化？？？？？？？？？？？？？！！！！！！！！！！！！！！！！！！！！！！！

[r,c]=find(I);%取出所有非零点的坐标,非边界
  %%%取好每个像素的坐标点
%   r = r(:);%将矩阵变成列向量，按照列的方式索引
%   c = c(:);

  % Compute GLCMS
  GLCMS = zeros(NL,NL,numOffsets);
  for k = 1 : numOffsets
    GLCMS(:,:,k) = computeGLCM(r,c,Offset(k,:),SI,NL);
    
    if makeSymmetric %%%%%%%%%%%%%%%%%像素对(如1 2或2 1出现的概率)的排序问题，缺省值时由offset决定，不用看
        % Reflect glcm across the diagonal
        glcmTranspose = GLCMS(:,:,k).';
        GLCMS(:,:,k) = GLCMS(:,:,k) + glcmTranspose;
    end
  end

else
  GLCMS = zeros(0,0,numOffsets);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算灰度共生矩阵
% -----------------------------------------------------------------
function oneGLCM = computeGLCM(r,c,offset,si,nl)%原图像的行坐标，列坐标，offset，灰度缩放后图像，灰度级数
% computes GLCM given one Offset
% （r2，c2）与（r，c）两点形成一定角度和距离（移动产生）
r2 = r + offset(1);%offset(1)取第一个元素￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥
c2 = c + offset(2);
E=[r,c];
A=[r2,c2];
F=[E,A];
%%%%%%%%%%查找A中的点那些在边界外，边界就是E组成的图形
TF=[];
for i=1:size(A,1)
    tf=all(E==A(i,:),2);
    TF=[TF,tf];
end

D=all(TF==0,1);
Index=find(D);%找出A中那一行坐标是边界外点
%将所有边界外点放到矩阵里
co=[];
for n=1:size(Index,2)
    
  c=[A(Index(n),1),A(Index(n),2)];  
    
   co=[co;c]; 
end

outsize=co-offset;%还原回去，这些点属于边界
v=si;
%删除这些边界
for i=1:size(outsize,1)
v(outsize(i,1),outsize(i,2))=0;
end
%按行排序
v=v';
In=find(v);
v1=v(In);%得到除去边界的原始图


v0=A;%移动后的坐标
%找出边界外点在A中对应的行
 B=[ ];
for i=1:size(co,1)
  b=find(all(A==co(i,:),2));
   B=[B;b]; 

end
v0(B,:)=[];%删除这些行

si1=zeros(size(si));

for n=1:size(v0,1)
  
 si1(v0(n,1)-offset(1),v0(n,2)-offset(2))=si(v0(n,1),v0(n,2));
    
end

si1=si1';
In2=find(si1);
v2=si1(In2);

% Remove pixel and its neighbor if their value is NaN.
bad = isnan(v1) | isnan(v2);%%B=isnan(A)，返回和A同样大小得数组，A中元素为NAN是该位置返回1，否则该位置返回0
if any(bad)%bad中出现非零元素或1
    warning(message('images:graycomatrix:scaledImageContainsNan'));
end

Ind = [v1 v2];
Ind(bad,:) = [];%%%%%%%%%%如果有，剔除这一行

if isempty(Ind)%如果A是空数组，TF = isempty(A)返回逻辑1 (true)，否则返回逻辑0 (false)。
    oneGLCM = zeros(nl);
else
%     Tabulate the occurrences of pixel pairs having v1 and v2.%列出包含v1和v2的像素对的出现次数。
%     Aq=max(Ind);
    
    oneGLCM = accumarray(Ind, 1, [nl nl]);%%%%%%Ind中数据为输出定义了一个矩阵，其中数据为矩阵的下标，下标相同的加1，矩阵的大小由[nl nl]决定。
end










%%-----------------------------------------------------------------------------解析输入，不用看
function [I, offset, nl, gl, sym] = ParseInputs(varargin)

narginchk(1,9);%验证输入参数个数的正确与否，对了什么都不做

% Check I
I = varargin{1};
validateattributes(I,{'logical','numeric'},{'2d','real','nonsparse'}, ...
              mfilename,'I',1);%%检查数组的有效性

% Assign Defaults
offset = [0 1];
if islogical(I)
  nl = 2;
else
  nl = 8;
end
gl = getrangefromclass(I);
sym = false;

% Parse Input Arguments
if nargin ~= 1
 
  paramStrings = {'Offset','NumLevels','GrayLimits','Symmetric'};
  
  for k = 2:2:nargin

    param = lower(varargin{k});
    inputStr = validatestring(param, paramStrings, mfilename, 'PARAM', k);
    idx = k + 1;  %Advance index to the VALUE portion of the input.
    if idx > nargin
      error(message('images:graycomatrix:missingParameterValue', inputStr));        
    end
    
    switch (inputStr)
     
     case 'Offset'
      
      offset = varargin{idx};
      validateattributes(offset,{'logical','numeric'},...
                    {'2d','nonempty','integer','real'},...
                    mfilename, 'OFFSET', idx);
      if size(offset,2) ~= 2
        error(message('images:graycomatrix:invalidOffsetSize'));
      end
      offset = double(offset);

     case 'NumLevels'
      
      nl = varargin{idx};
      validateattributes(nl,{'logical','numeric'},...
                    {'real','integer','nonnegative','nonempty','nonsparse'},...
                    mfilename, 'NL', idx);
      if numel(nl) > 1
        error(message('images:graycomatrix:invalidNumLevels'));
      elseif islogical(I) && nl ~= 2
        error(message('images:graycomatrix:invalidNumLevelsForBinary'));
      end
      nl = double(nl);      
     
     case 'GrayLimits'
      
      gl = varargin{idx};
      % step 1: checking for classes
      validateattributes(gl,{'logical','numeric'},{},mfilename, 'GL', idx);
      if isempty(gl)
        gl = [min(I(:)) max(I(:))];
      end
      
      % step 2: checking for attributes
      validateattributes(gl,{'logical','numeric'},{'vector','real'},mfilename, 'GL', idx);
      
      if numel(gl) ~= 2
        error(message('images:graycomatrix:invalidGrayLimitsSize'));
      end
      gl = double(gl);
    
      case 'Symmetric'
        sym = varargin{idx};
        validateattributes(sym,{'logical'}, {'scalar'}, mfilename, 'Symmetric', idx);
        
    end
  end
end
