%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = rosenbr( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The ever famous 2 variables Rosenbrock "banana valley" problem
%
%   Source: problem 1 in
%      J.J. More, B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%
%   For dimension higher than 2, this version of problem (with the starting
%   point of all ones) appears to be much harder than extrosnb.
%   
%   Ph. Toint 26 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'rosenbr';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in rosenbr: n = ', int2str(n), ' must be > 1!' ] )
      end
   else
      n = 2;
   end
   if ( n == 2 )
      varargout{1} = [ -1.2; ones( n-1,1) ];    % x0
   else
      varargout{1} = [ -ones( n,1) ];           % x0
   end
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AY-V-0';                % class

case 'cpsstr'

   n = varargin{1};
   for iel = 1:n-1
      eldom{ iel } = [ iel iel+1 ];
   end
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = {};
   varargout{1} = cpsstr;
   
case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) );
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', length( x ) );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, cpsstr.param{:} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   r1   = 10*(x(2)-x(1)^2);
   r2   = 1-x(1);
   varargout{1} = r1^2 + r2^2;
   if ( nargout > 1 )
      J1 = 10*[-2*x(1); 1 ];
      J2 = [ -1; 0 ];
      varargout{2} = 2* ( J1*r1 + J2*r2 );
      if ( nargout > 2 )
         H1 = [ -20, 0; 0, 0  ];
	 varargout{3} = 2 * ( J1*J1.' + r1*H1 + J2*J2.' );
      end
   end
end

return

end