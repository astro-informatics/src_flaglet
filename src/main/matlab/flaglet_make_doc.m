function flaglet_make_doc(flagletpath)

% flaglet_make_doc
% Generate Matlab documentation
%
% Default usage :
%
%   flaglet_make_doc(flagletpath)
%
% flagletpath is the path for the FLAGLET package (root)
%
% FLAGLET package to perform exact wavelets on the ball.
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

cd(flagletpath)
m2html('mfiles', 'src/main/matlab', 'htmldir', 'doc/matlab');

end