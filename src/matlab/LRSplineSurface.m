classdef LRSplineSurface < handle
% LRSplineSurface Matlab wrapper class for c++ LR-spline object
%     Locally Refined (LR) B-splines is a technique to achive local adaptivity while using smooth spline 
%     functions. This is a sample library which implements these techniques and can be used for geometric
%     representations or isogeometric analysis.
%     
% LRSplineSurface Properties: 
%     p        - polynomial degree
%     knots    - knot vectors
%     cp       - control points
%     w        - weights
%     lines    - mesh lines, (u0,v0, u1,v1, m), where m is the multiplicity
%     elements - fintite elements (u0, v0, u1, v1)
%     support  - element to basis function support list
%    
% LRSplineSurface Methods:
%     copy                 - Performs a deep copy of the spline object
%     refine               - Performs local refinements
%     raiseOrder           - Performs global degree elevation
%     getEdge              - Extracts functions with support on one of the four parametric edges
%     getElementContaining - Get element index at parametric point (u,v)
%     point                - Evaluates the physical coordinates (x,y) corresponding to a parametric point (u,v)
%     computeBasis         - Compute all basis functions (and their derivatives)
%     getBezierExtraction  - Get the bezier extraction matrix for one element
%     setContinuity        - Performs global continutiy reduction
%     L2project            - L2-project results onto the spline basis 
%     surf                 - Plot scalar results in a surface plot (per element or per controlpoint)
%     plot                 - Plot the mesh structure 
%     print                - Prints raw c++ lr data structure

	properties(SetAccess = private, Hidden = false)
		p        % polynomial degree
		knots    % knot vectors
		cp       % control points
		w        % weights
		lines    % mesh lines, (u0,v0, u1,v1, m), where m is the multiplicity
		elements % fintite elements (u0, v0, u1, v1)
		support  % element to basis function support list
	end
	properties(SetAccess = private, Hidden = true)
		objectHandle;
	end

	methods
		function this = LRSplineSurface(varargin)
		% LRSplineSurface  Constructor, initialize a tensor product LRSplinSurface object
		% LRSplineSurface(n,p)
		% LRSplineSurface(n,p, knotU, knotV)
		% LRSplineSurface(n,p, knotU, knotV, controlpoint)
		% 
		%   parameters 
		%     n            - number of basis functions in each direction (2 components)
		%     p            - polynomial degree in each direction (2 components)
		%     knotU        - global open knot vector in u-direction (n(1)+p(1)+1 components)
		%     knotV        - global open knot vector in v-direction (n(2)+p(2)+1 components)
		%     controlpoint - list of control points (matrix of size dim x n(1)*n(2)), where dim is dimension in physical space


			% error check input
			if(nargin == 0)
				objectHandle = 0;
				return;
			end
			n = varargin{1};
			p = varargin{2};
			if(nargin ~= 2 && nargin ~=4 && nargin ~= 5)
				throw(MException('LRSplineSurface:constructor',  'Error: Invalid number of arguments to LRSplineSurface constructor'));
			end
			if(length(p) ~=2 || length(n) ~=2)
				throw(MException('LRSplineSurface:constructor', 'Error: p and n should have 2 components'));
			end
			if(nargin > 3)
				for i=1:2
					if(~(size(varargin{i+2}) == [1, p(i)+n(i)+1]) )
						throw(MException('LRSplineSurface:constructor', 'Error: Knot vector should be a row vector of length p+n+1'));
					end
				end
			end
			if(nargin > 4)
				if(size(varargin{5},2) ~= n(1)*n(2))
					throw(MException('LRSplineSurface:constructor', 'Error: Control points should have n(1)*n(2) columns'));
				end
			end
			
			this.objectHandle = lrsplinesurface_interface('new', varargin{:});
			this.updatePrimitives();
		end

		function delete(this)
		% LRSplineSurface  Destructor clears object from memory
			lrsplinesurface_interface('delete', this.objectHandle);
		end


		function copyObject = copy(this)
		% COPY  peforms a deep copy of the spline object
		% LRSplineSurface.copy()
		%
		%   returns:
		%     new LRSpline object
			newHandle  = lrsplinesurface_interface('copy', this.objectHandle);
			copyObject = LRSplineSurface();
			copyObject.setHandle(newHandle);
		end


		function print(this)
		% PRINT  Dumps the backend c++ representation of this LR-spline object to screen
		% LRSplineSurface.print()
		%
		%   parameters:
		%     none
			lrsplinesurface_interface('print', this.objectHandle);
		end


		function refine(this, indices, varargin)
		% REFINE  Performs local refinement of elements or basis functions
		% LRSplineSurface.refine(indices)
		% LRSplineSurface.refine(indices, 'elements')
		% LRSplineSurface.refine(indices, 'basis')
		% LRSplineSurface.refine(indices, 'continuity', n)
		%
		%   parameters:
		%     indices      - index of the basis function or elements to refine
		%     'elements'   - perform full span refinement on the input elements
		%     'basis'      - perform structure mesh refinement on the input functions
		%     'continuity' - set the refinement continuity to n (less than polynomial degree)
		%   returns
		%     none
			mult     = 1;
			elements = false;
			i        = 0;
			% read input parameters
			while(i<nargin-2)
				i=i+1;
				if     strcmp(varargin{i}, 'elements')
					elements = true;
				elseif strcmp(varargin{i}, 'basis')
					elements = false;
				elseif strcmp(varargin{i}, 'continuity')
					mult = max(this.p)-varargin{i+1};
					i=i+1;
				else
					throw(MException('LRSplineSurface:refine',  'Error: Unknown refine parameter'));
				end
			end

			if(elements)
				% error control
				if(min(indices)<0 || max(indices)>size(this.elements, 1))
					throw(MException('LRSplineSurface:refine',  'Error: Invalid refinement index'));
				end
				% perform refinement
				lrsplinesurface_interface('refine_elements', this.objectHandle, indices, mult);
			else
				% error control
				if(min(indices)<0 || max(indices)>size(this.knots, 1))
					throw(MException('LRSplineSurface:refine',  'Error: Invalid refinement index'));
				end
				% perform refinement
				lrsplinesurface_interface('refine_basis', this.objectHandle, indices, mult);
			end

			% new LR-mesh... update static variables
			this.updatePrimitives();
		end

		function setContinuity(this, newCont, newCont2)
		% SETCONTINUITY  Lowers the global continuity to max C^{newCont}
		%
		%   parameters:
		%     newCont  - new continuity for the global solution space
		%     newCont2 - new continuity in second parameter direction
			c = newCont;
			if nargin > 2
				c(2) = newCont2;
			elseif numel(c)==1
				c(2) = newCont;
			end
			lrsplinesurface_interface('set_continuity', this.objectHandle, c);
			this.updatePrimitives();
		end


		function cp = L2project(this, u,v,z, w)
		% L2project  Performs a global L2 projection into LR spline space
		% LRSplineSurface.raiseOrder(u,v,z)
		% LRSplineSurface.raiseOrder(u,v,z,w)
		%
		%   parameters:
		%     u - vector of first parametric point
		%     v - vector of second parametric point
		%     z - vector of L2-projection points
		%     w - [optional] list of weights
		%   returns
		%     cp - list of control points corresponding to this projection
			nCP = size(this.cp,2);
			if(nCP > length(u))
				throw(MException('LRSplineSurface:L2project', 'Error: too few evaluation points to do global L2 projection'));
			end
			A = sparse(nCP, nCP);
			b = sparse(nCP, size(z,2));
			for i=1:length(u)
				el  = this.getElementContaining(u(i), v(i));
				ind = this.support{el};
				N   = this.computeBasis(u(i), v(i));
				if(nargin > 4) % continuous L2 projection, include weights
					A(ind, ind) = A(ind, ind) + N'*N      * w(i);
					b(ind,:)    = b(ind,:)    + N'*z(i,:) * w(i);
				else
					A(ind, ind) = A(ind, ind) + N'*N;
					b(ind,:)    = b(ind,:)    + N'*z(i,:);
				end
			end
			cp = A \ b;
		end


		function raiseOrder(this, dp, dq)
		% RAISEORDER  Performs global degree elevation
		% LRSplineSurface.raiseOrder(dp)
		% LRSplineSurface.raiseOrder(dp, dq)
		%
		%   parameters:
		%     dp - amount to increase in the first parametric direction
		%     dq - amount to increase in the second parametric direction
		%   returns
		%     none
			oldGuy = this.copy();
			newHandle = lrsplinesurface_interface('raise_order', this.objectHandle, dp, dq);
			lrsplinesurface_interface('delete', this.objectHandle);
			this.objectHandle = newHandle;
			this.updatePrimitives();

			nElms  = size(this.elements,1);
			nBasis = size(this.knots,1);
			newCP  = zeros(size(this.cp,1), nBasis);
			% ideally we would like to do an greville interpolation, or even quasi interpolation would
			% work, but sometimes the greville points seem to stack on top of each other. We'll do nxn
			% evaluation points for each element and hope this suffices for an L2-projection
			nPts  = ceil(sqrt(nBasis / nElms));
			uAll  = zeros(nPts*nPts*nElms,1);
			vAll  = zeros(nPts*nPts*nElms,1);
			cpAll = zeros(nPts*nPts*nElms,2);

			k = 1;
			for iEl=1:nElms,
				% make a tensor grid of evaluation points on this element
				u = linspace(this.elements(iEl,1), this.elements(iEl,3), nPts+2);
				v = linspace(this.elements(iEl,2), this.elements(iEl,4), nPts+2);
				u = u(2:end-1);
				v = v(2:end-1);
				for i=1:nPts
					for j=1:nPts
						uAll(k)    = u(i);
						vAll(k)    = v(j);
						cpAll(k,:) = oldGuy.point(u(i), v(j));
						k = k+1;
					end
				end
			end
			newCP = this.L2project(uAll, vAll, cpAll);

			lrsplinesurface_interface('set_control_points', this.objectHandle, newCP');
			clear oldGuy;
			this.updatePrimitives();
		end

		function x = point(this, u, v)
		% POINT  Evaluates the mapping from parametric to physical space
		% x = LRSplineSurface.point(u,v)
		%
		%   parameters:
		%     u - first parametric coordinate
		%     v - second parametric coordinate
		%   returns
		%     the parametric point mapped to physical space
			x = lrsplinesurface_interface('point', this.objectHandle, [u,v]);
		end


		function N = computeBasis(this, u, v, varargin)
		% COMPUTEBASIS  Evaluates all basis functions at a given parametric point, as well as their derivatives
		% N = LRSplineSurface.computeBasis(u, v)
		% N = LRSplineSurface.computeBasis(u, v, derivs)
		%
		%   parameters:
		%     u      - first parametric coordinate
		%     v      - second parametric coordinate
		%     derivs - number of derivatives (greater or equal to 0)
		%   returns
		%     the value of all nonzero basis functions at a given point
		%     in case of derivatives, a cell is returned with all derivatives requested
			N = lrsplinesurface_interface('compute_basis', this.objectHandle, [u, v], varargin{:});
		end


		function C = getBezierExtraction(this, element)
		% GETBEZIEREXTRACTION  Returns the bezier extraction matrix for this element
		% C = LRSplineSurface.getBezierExtraction(element)
		%
		%   parameters:
		%     element - global index to the element 
		%   returns
		%     a matrix with as many rows as there is active basis functions and (p(1)+1)*(p(2)+1) columns
			C = lrsplinesurface_interface('get_bezier_extraction', this.objectHandle, element);
		end


		function iel = getElementContaining(this, u,v)
		% GETELEMENTCONTAINING  Returns the index of the element containing the parametric point (u,v)
		% iel = getElementContaining(u,v)
		%
		%   parameters:
		%     u - first parametric coordinate
		%     v - second parametric coordinate
		%   returns
		%     index to the element containint this parametric point
			iel = lrsplinesurface_interface('get_element_containing', this.objectHandle, [u,v]);
		end


		function index = getEdge(this, edge, varargin)
		% GETEDGE  Returns a list of all basis functions with nonzero value at one of the four parametric edges
		% index = LRSplineSurface.getEdge()
		% index = LRSplineSurface.getEdge(n)
		% index = LRSplineSurface.getEdge(n, 'elements')
		%
		%   parameters:
		%     n - the local edge number (all=0, umin=1, umax=2, vmin=3, vmax=4)
		%   returns
		%     list of all elements or basis function with support on this edge
			elements = false;
			index = [];
			if(nargin < 2)
				edge = 0;
			end
			if(nargin > 2)
				if(strcmp(varargin{1}, 'elements'))
					elements = true;
				else 
					throw(MException('LRSplineSurface:getEdge', 'Error: Unkown parameters'));
				end
			end
			%%% error test input
			if(edge~=0 && edge~=1 && edge~=2 && edge~=3 && edge~=4)
				throw(MException('LRSplineSurface:getEdge', 'Error: Invalid edge enumeration'));
			end
			if(edge == 1 || edge == 0)
				umin = min(this.elements(:,1));
				if(elements)
					index = [index; find(this.elements(:,1) == umin)];
				else
					index = [index; find(this.knots(:, this.p(1)+1) == umin)];
				end
			end
			if(edge == 2 || edge == 0)
				umax = max(this.elements(:,3));
				if(elements)
					index = [index; find(this.elements(:,3) == umax)];
				else
					index = [index; find(this.knots(:, 2) == umax)];
				end
			end
			if(edge == 3 || edge == 0)
				vmin = min(this.elements(:,2));
				if(elements)
					index = [index; find(this.elements(:,2) == vmin)];
				else
					index = [index; find(this.knots(:, end-1) == vmin)];
				end
			end
			if(edge == 4 || edge == 0)
				vmax = max(this.elements(:,4));
				if(elements)
					index = [index; find(this.elements(:,4) == vmax)];
				else
					index = [index; find(this.knots(:, this.p(1)+4) == vmax)];
				end
			end
		end

		function H = plot(this, varargin)
		% PLOT  Creates a plot of the LRSplineSurface as mapped into the physical coordinate space
		% H = LRSplineSurface.plot()
		% H = LRSplineSurface.plot('enumeration')
		% H = LRSplineSurface.plot('nviz', n)
		% H = LRSplineSurface.plot('parametric')
		%
		%   parameters:
		%     'enumeration' - tags all elements with their corresponding enumeration index
		%     'basis'       - plots control points as dots (greville points if 'parametric' is specified)
		%     'parametric'  - prints the elements in the parametric space instead of the physical 
		%     'nviz'        - sets the line resolution for plots to use n points for drawing each line
		%   returns
		%     handle to the figure
			
			nPtsPrLine  = 41;
			enumeration = false;
			parametric  = false;
			basis       = false;

			i = 1;
			while i<nargin
				if strcmp(varargin{i}, 'enumeration')
					enumeration = true;
				elseif strcmp(varargin{i}, 'nviz')
					i = i+1;
					nPtsPrLine = varargin{i};
				elseif strcmp(varargin{i}, 'parametric')
					parametric = true;
				elseif strcmp(varargin{i}, 'basis')
					basis = true;
				else
					throw(MException('LRSplineSurface:plot',  'Error: Unknown input parameter'));
				end
				i = i+1;
			end

			if(parametric)
				nPtsPrLine = 2;
			end
			nLines     = size(this.lines, 1);
			x = zeros(nPtsPrLine, nLines);
			y = zeros(nPtsPrLine, nLines);
			for i=1:nLines
				u = linspace(this.lines(i,1), this.lines(i,3), nPtsPrLine);
				v = linspace(this.lines(i,2), this.lines(i,4), nPtsPrLine);
				for j=1:nPtsPrLine
					if(parametric)
						x(j,i) = u(j);
						y(j,i) = v(j);
					else
						res = this.point(u(j), v(j));
						x(j,i) = res(1);
						y(j,i) = res(2);
					end
				end
			end
			holdOnReturn = ishold;
			H = plot(x,y, 'k-');

			% enumerate elements if specified
			if(enumeration && ~basis)
				hold on;
				for i=1:size(this.elements, 1),
					if(parametric)
						x = [sum(this.elements(i, [1,3]))/2, sum(this.elements(i,[2,4]))/2];
					else 
						x = this.point(sum(this.elements(i, [1,3]))/2, sum(this.elements(i,[2,4]))/2);
					end
					text(x(1), x(2), num2str(i));
				end
			end

			% plot basis if specified
			if basis
				hold on;
				for i=1:size(this.knots,1)
					if parametric
						x = [sum(this.knots(i,2:(this.p(1)+1))); sum(this.knots(i,(this.p(1)+4):(end-1)))] ./ this.p;
					else
						x = this.cp(:,i);
					end
					plot(x(1), x(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', [0.9843137, 0.8384314, 0.40117648], 'MarkerEdgeColor', 'black');
					if enumeration
						text(x(1), x(2), num2str(i));
					end
				end
			end
			if ~holdOnReturn
				hold off;
			end
		end


		function H = surf(this, u, varargin)
		% SURF  Creates a surface plot of scalar results u given by control point values OR per element values
		% H = LRSplineSurface.surf(u)
		% H = LRSplineSurface.surf(u, 'nviz', n)
		% H = LRSplineSurface.surf(u, 'secondary', f)
		%
		% If the number of components passed is equal to the number of elements, this is interpreted as per-element
		% results (i.e. error norms). Else, it is treated as scalar control-point variables (i.e. primary solution field)
		%
		%   parameters:
		%     u            - control point results
		%     'nviz'       - sets the plotting resolution to n points per element
		%     'diffX'      - plots the derivative with respect to X 
		%     'diffY'      - plots the derivative with respect to Y
		%     'secondary'  - plots secondary solutions such as functions of u and dudx
		%     'parametric' - displays results in parametric space (and parametric derivatives)
		%   returns
		%     handle to the figure
			nviz               = 6;
			diffX              = false;
			diffY              = false;
			parametric         = false;
			mononomial         = false;
			per_element_result = false;
			function_result    = false;
			secondary          = false;
			sec_function       = 0;

			i = 1;
			while i<nargin-1
				if strcmp(varargin{i}, 'diffX')
					diffX = true;
				elseif strcmp(varargin{i}, 'diffY')
					diffY = true;
				elseif strcmp(varargin{i}, 'mononomial')
					mononomial = true;
				elseif strcmp(varargin{i}, 'secondary')
					secondary = true;
					i = i+1;
					sec_function = varargin{i};
				elseif strcmp(varargin{i}, 'nviz')
					i = i+1;
					nviz = varargin{i};
				elseif strcmp(varargin{i}, 'parametric')
					parametric = true;
				else
					throw(MException('LRSplineSurface:surf',  'Error: Unknown input parameter'));
				end
				i = i+1;
			end
			xg = linspace(-1,1,nviz);

			if strcmp(class(u), 'function_handle')
				function_result = true;
			elseif ~mononomial
				u = u(:)'; % make u a row vector
				if(numel(u) == size(this.elements,1))
					per_element_result = true;
				end
			end

			if mononomial,
				nDOF = size(this.knots,1);
				grevU = zeros(2, nDOF);
				grevX = zeros(2, nDOF);
				for i=1:nDOF
					grevU(:,i) = [sum(this.knots(i,2:(this.p(1)+1))); sum(this.knots(i,(this.p(1)+4):(end-1)))] ./ this.p;
					grevX(:,i) = this.point(grevU(1,i), grevU(2,i));
				end
			end

			holdOnReturn = ishold;
			H = gcf;
			hold on;

			Xlines = zeros(size(this.elements, 1)*4, nviz);
			Ylines = zeros(size(this.elements, 1)*4, nviz);
			Zlines = zeros(size(this.elements, 1)*4, nviz);

			bezierKnot1 = [ones(1, this.p(1)+1)*-1, ones(1, this.p(1)+1)];
			bezierKnot2 = [ones(1, this.p(2)+1)*-1, ones(1, this.p(2)+1)];
			[bezNu, bezNu_diff] = getBSplineBasisAndDerivative(this.p(1), xg, bezierKnot1); 
			[bezNv, bezNv_diff] = getBSplineBasisAndDerivative(this.p(2), xg, bezierKnot2); 
			for iel=1:size(this.elements, 1)
				umin = this.elements(iel,1);
				vmin = this.elements(iel,2);
				umax = this.elements(iel,3);
				vmax = this.elements(iel,4);
				hu = umax-umin;
				hv = vmax-vmin;
				ind  = this.support{iel}; % indices to nonzero basis functions
				C  = this.getBezierExtraction(iel);
				X  = zeros(nviz);
				Y  = zeros(nviz);
				U  = zeros(nviz);
				Ux = zeros(nviz);
				Uy = zeros(nviz);
				% for all gauss points
				for i=1:nviz
					for j=1:nviz
						xi  = (.5*xg(i)+.5)*(umax-umin)+umin;
						eta = (.5*xg(j)+.5)*(vmax-vmin)+vmin;

						% compute all basis functions
						N     = bezNu(:,i)       * bezNv(:,j)';
						dNdu  = bezNu_diff(:,i)  * bezNv(:,j)';
						dNdv  = bezNu(:,i)       * bezNv_diff(:,j)';
						N     = N(:); % and make results colum vector
						dN    = [dNdu(:)*2/hu, dNdv(:)*2/hv];

						% evaluates physical mapping and jacobian
						x  = this.cp(:,ind) * C * N;
						Jt = this.cp(:,ind) * C * dN; % transpose jacobian matrix [dx/du,dy/du; dx/dv, dy/dv]

						% physical derivatives
						dNdx = dN * inv(Jt'); 

						% write results depending on type of plot
						if(parametric)
							X(i,j) = xi;
							Y(i,j) = eta;
						else
							X(i,j) = x(1);
							Y(i,j) = x(2);
						end
						if function_result || secondary
							if secondary
								if nargin(sec_function)==2 % input parameters x and u
									U(i,j) = sec_function(x, u(ind) * C * N);
								elseif nargin(sec_function)==3 % input parameters x, u and dudx
									U(i,j) = sec_function(x, u(ind) * C * N, (u(ind) * C * dNdx)');
								end
							else
								U(i,j) = u(x);
							end
						elseif per_element_result
							U(i,j) = u(iel);
						elseif diffX && parametric
							U(i,j)  = u(ind) * C * dN(:,1);
						elseif diffX 
							U(i,j)  = u(ind) * C * dNdx(:,1);
						elseif diffY && parametric
							U(i,j)  = u(ind) * C * dN(:,1);
						elseif diffY
							U(i,j)  = u(ind) * C * dNdx(:,2);
						elseif mononomial
							nFun = 0;
							for iBasis=ind
								if ~isnan( u(1,iBasis) )
									nFun = nFun + 1;
									p = sqrt(size(u,1))-1;
									mon = ones(2,p+1);
									for k=1:p
										mon(:,k+1) = mon(:,k) .* (x-grevX(:,iBasis));
									end
									monAll = mon(1,:)' * mon(2,:);
									U(i,j) = U(i,j) + u(:,iBasis)'*monAll(:);
								end
							end
							U(i,j) = U(i,j) / nFun;
						else
							U(i,j) = u(ind) * C * N;
						end
					end
				end
				if sum(sum(isnan(U)))>0,
					iel
					U
					pause;
				end
				surf(X,Y,U, 'EdgeColor', 'none');
				Xlines((iel-1)*4+1,:) = X(1,:);
				Ylines((iel-1)*4+1,:) = Y(1,:);
				Zlines((iel-1)*4+1,:) = U(1,:);

				Xlines((iel-1)*4+2,:) = X(end,:);
				Ylines((iel-1)*4+2,:) = Y(end,:);
				Zlines((iel-1)*4+2,:) = U(end,:);

				Xlines((iel-1)*4+3,:) = X(:,1);
				Ylines((iel-1)*4+3,:) = Y(:,1);
				Zlines((iel-1)*4+3,:) = U(:,1);

				Xlines((iel-1)*4+4,:) = X(:,end);
				Ylines((iel-1)*4+4,:) = Y(:,end);
				Zlines((iel-1)*4+4,:) = U(:,end);
				% plot3(X(1,:),   Y(1,:),   U(1,:),   'k-');
				% plot3(X(end,:), Y(end,:), U(end,:), 'k-');
				% plot3(X(:,1),   Y(:,1),   U(:,1),   'k-');
				% plot3(X(:,end), Y(:,end), U(:,end), 'k-');
			end
			plot3(Xlines', Ylines', Zlines', 'k-');
			if(per_element_result)
				view(2);
			else	
				view(3);
			end
		end
	end

	methods (Hidden = true)
		function insertLine(this, u,v,m)
			if(numel(u) ~=2 || numel(v) ~=2)
				throw(MException('LRSplineSurface:insertLine',  'Error: Invalid arguments'));
			end
			if nargin<4
				m = 1;
			end
				
			lrsplinesurface_interface('insert_line', this.objectHandle, u,v,m);
			this.updatePrimitives();
		end
	end

	methods (Access = private, Hidden = true)
		function updatePrimitives(this)
			[this.knots, this.cp, this.w, ...
			 this.lines, this.elements,   ...
			 this.support, this.p] = lrsplinesurface_interface('get_primitives', this.objectHandle);
		end

		function setHandle(this, handle)
			if this.objectHandle ~= 0
				lrsplinesurface_interface('delete', this.objectHandle);
			end
			this.objectHandle = handle;
			this.updatePrimitives();
		end
	end
end
