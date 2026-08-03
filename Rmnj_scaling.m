

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 1D FMM with adaptive domain partition  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [u1, u2, post_ord,ch, Px, Py,I, nflops, neighbor] = Rmnj_scaling(r, x, y, q, scaling)

%%%%%%%%%%% Input %%%%%%%%%%%%

% r: truncation order, number of terms in Taylor expansion

% x: source points, real vectors
% y: target points, real vectors

% q: charges, vectors

% scaling = 1 : including diagonal scaling for stability
% scaling = 0 : NOT including diagonal scaling for stability


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% set up %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert x, y to column vectors
if size( x,1 ) < size( x,2 )
    x = x.';
end
if size( y,1 ) < size( y,2 )
    y = y.';
end

% % sort x, y
% x = sort( x );
% y = sort( y );

N0 = 64;  % best leaf size for max speedup factors

% separation ratio
tau = 0.6;

% computation cell
z2 = max([x;y]);
z1 = min([x;y]);
z2 = z2 + 0.1*abs(z2);
z1 = z1 - 0.1*abs(z1);

nflops = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  adaptive partition %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [center, radius, left_end_point, right_end_point, tree_index, parent, level]
I = [(z1+z2)/2, (z2-z1)/2, 1, 1, 1, 0, 1];
S = {};
S = push(S, I);
Px = {[1:length(x)]};
Py = {[1:length(y)]};

while ~isempty(S)
    [S, i] = pop(S);
    c = i(1); d = i(2); el = i(3); er = i(4); idx = i(5); lvl = i(7);
    n = size(I, 1);
    px = Px{idx}; Px{idx} = [];
    py = Py{idx}; Py{idx} = [];
    
    %%% bisect the interval
    i1 = [c-d/2, d/2, el, 1, n+1, idx, lvl+1];
    i2 = [c+d/2, d/2, 0, er, n+2, idx, lvl+1];
    I = [I; i1; i2];
    if el == 1
        Px{n+1} = px(c-d <= x(px) & x(px) <= c);
        Py{n+1} = py(c-d <= y(py) & y(py) <= c);
    elseif el == 0
        Px{n+1} = px(c-d < x(px) & x(px) <= c);
        Py{n+1} = py(c-d < y(py) & y(py) <= c);
    end
     
    if er == 1
        Px{n+2} = px(c < x(px) & x(px) <= c+d);
        Py{n+2} = py(c < y(py) & y(py) <= c+d);
    elseif er == 0
        Px{n+2} = px(c < x(px) & x(px) < c+d);
        Py{n+2} = py(c < y(py) & y(py) < c+d);
    end
    
    if length(Px{n+1}) > N0 || length(Py{n+1}) > N0
        S = push(S, i1);
    end
    
    if length(Px{n+2}) > N0 || length(Py{n+2}) > N0
        S = push(S, i2);
    end 
end


% tree information
tr = I(:, 6);
ch = child1(tr);
n = length(tr);
post_ord = postorder(1, tr, ch);
pre_ord = preorder(1, tr, ch);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% neighbors and interaction list %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neighbor = cell(n,1);
interlist = cell(n,1);
for i = pre_ord
    if tr(i) == 0
        neighbor{i} = [];
        interlist{i} = [];
    else
        neighbor{i} = [sib(tr,ch,i)];
        
        cousins = [];
        for j = neighbor{tr(i)}
            if isempty(ch{j})
                cousins = [cousins, j];
            else
                cousins = [cousins, ch{j}(1), ch{j}(2)];
            end
        end
        
        ci = I(i, 1); 
        ri = I(i, 2);
        for j = cousins
            cj = I(j, 1); 
            rj = I(j, 2);
            if ri + rj <= tau * abs(ci - cj)
                interlist{i} = [interlist{i}, j];
            else
                neighbor{i} = [neighbor{i}, j];
            end
        end
    end
end


for i = pre_ord
    if isempty(ch{i})
        for j = neighbor{i}
            neighbor{j} = unique([neighbor{j}, i]);
        end
    end 
    
    for j = interlist{i}
        interlist{j} = unique([interlist{j}, i]);
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% save fmm coefficients %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u1 = cell(n, 1);    % local expansion
u2 = cell(n, 1);
v1 = cell(n, 1);    % multipole expansion


%%%% post order (bottom up) traversal
for i = post_ord
    lvl = I(i, 7);
    if lvl >= 2
        if isempty(ch{i})
            xi = x(Px{i});
            % Vi1 = particle_to_local(xi, r, I(i,1), 2*I(i,2), scaling);
            Vi = P2L(xi, r, I(i,1), 2*I(i,2), scaling);
            v1{i} = Vi' * q(Px{i}, :);
            nflops = nflops + size(Vi,2)*(2*size(Vi,1)-1)*size(q,2);
        else
            c1 = ch{i}(1); 
            c2 = ch{i}(2);
             Wc1 = M2M( r, I(c1,1), I(i,1), 2*I(c1,2), 2*I(i,2), scaling);  % upward pass
             Wc2 = M2M( r, I(c2,1), I(i,1), 2*I(c2,2), 2*I(i,2), scaling);  % upward pass

            v1{i} = Wc1' * v1{c1} + Wc2' * v1{c2};
            nflops = nflops + flops('prod', Wc1, 't', v1{c1}, 'n')...
                                         + flops('prod', Wc2, 't', v1{c2}, 'n') + numel(v1{i});
        end
    end
end

%%%% pre order (top down) traversal
for i = pre_ord
    lvl = I(i, 7);
    if lvl >= 2
        list = interlist{i};
        for j = list
       % [BB1_S, BB2_S] = multipole_to_local(r, I(i,1), I(j,1),2*I(i,2),2*I(j,2), scaling); % Stirling scaling 
       [Bij, B2] = M2L(r, I(i,1), I(j,1), 2*I(i,2), 2*I(j,2), scaling); % Ramanujan scaling
            if isempty(u1{i})
                u1{i} = Bij * v1{j};
                u2{i} = B2 * v1{j};
                nflops = nflops + flops('prod', Bij, 'n', v1{j}, 'n');
                nflops = nflops + flops('prod', B2, 'n', v1{j}, 'n');
            else
                u1{i} = u1{i} + Bij * v1{j};
                u2{i} = u2{i} + B2 *  v1{j};
                nflops = nflops + flops('prod', Bij, 'n', v1{j}, 'n') + numel(u1{i});
                nflops = nflops + flops('prod', B2, 'n', v1{j}, 'n') + numel(u1{i});
            end
        end
        
        if lvl > 3
            p = tr(i);
            if ~isempty(u1{p})
                % Ri_1 = multipole_to_multipole(r, I(i,1), I(p,1),  2*I(i,2), 2*I(p,2), scaling); % downward pass
                Ri = M2M(r, I(i,1), I(p,1),  2*I(i,2), 2*I(p,2), scaling); % downward pass
                if isempty(u1{i})
                    u1{i} = Ri * u1{p};
                    u2{i} = Ri * u2{p};
                    nflops = nflops + flops('prod', Ri, 'n', u1{p}, 'n');
                    nflops = nflops + flops('prod', Ri, 'n', u2{p}, 'n');
                else
                    u1{i} = u1{i} + Ri * u1{p};
                    u2{i} = u2{i} + Ri * u2{p};
                    nflops = nflops + flops('prod', Ri, 'n', u1{p}, 'n') + numel(u1{i});
                    nflops = nflops + flops('prod', Ri, 'n', u2{p}, 'n') + numel(u2{i});
                end
            end
        end
    end
end







