

clc
clear

addpath(genpath('./'))

load('data2.mat')  % random starting points between each interval

setup_time = zeros(10,1);
elapsed_time = cell(10,2);
fmm_fun_val = cell(10,1);
fmm_y = cell(10,1);  % iterates at the end of the loop
flop_counts = cell(10,1);

scaling = 1;    
n_iter = 10;  % Newton iterations
for r = 20 %5:5:30

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% fmm evaluation %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj =  1 %1:length(data2)
jj
x = data2{jj,1};  % x -- sources
y = data2{jj,2};  % y -- targets 
q = data2{jj,3};

q = q.^2;

% [interlace_flag, interlace_fail_idx] = check_interlace(x, y);  % check whether y interlaces x

tic
[u1, u2, post_ord,ch, Px, Py, I, nflops, neighbor] = test_scaling(r, x, y, q, scaling);
% [u1, u2, post_ord,ch, Px, Py, I, nflops, neighbor] = Stirling_scaling(r, x, y, q, scaling);
% [u1, u2, post_ord,ch, Px, Py, I, nflops, neighbor] = Rmnj_scaling(r, x, y, q, scaling);
setup_time(jj) = toc;

fmm_timer = zeros(n_iter,1);
for iter = 1:n_iter
ftt = tic;
z1 = zeros(length(q),1);  % final product for function
z2 = zeros(length(q),1);  % final product for derivative

for i = post_ord
    if isempty(ch{i})
        py = Py{i}; 
        if ~isempty(py)
            yj = y(py);
            %%% evaluate local expansion
            if ~isempty(u1{i})
                Ui = P2L(yj, r, I(i,1), 2*I(i,2), scaling);
                z1(py, :) = Ui * u1{i};
                z2(py, :) = Ui * u2{i};
                nflops = nflops + flops('prod', Ui, 'n', u1{i}, 'n');
                nflops = nflops + flops('prod', Ui, 'n', u2{i}, 'n');
            end


            %%% near field evaluation
            for j = union(i, neighbor{i})
                px = Px{j};
                if ~isempty(px)
                    xi = x(px);
                    D1 = 1 ./ (xi.' - yj);
                    D2 = 1 ./ (xi.' - yj).^2;
                    D1(isinf(D1)) = 0;
                    D2(isinf(D2)) = 0;
                    z1(py, :) = z1(py, :) + D1 * q(px, :);
                    z2(py, :) = z2(py, :) + D2 * q(px, :);
                    nflops = nflops + size(D1,1)*(2*size(D1,2)-1)*size(q,2) + numel(py)*size(z1,2);
                    nflops = nflops + size(D2,1)*(2*size(D2,2)-1)*size(q,2) + numel(py)*size(z2,2);
                end
            end
        end
    end
end


% Newton iteration
w1 = 1 + z1;
ynew = y - (w1./z2);
y = ynew;

fmm_timer(iter) = toc(ftt);
end  % end of Newton iterations 
fmm_fun_val{jj} = w1;  % fmm function values at the end of the loop
fmm_y{jj} = ynew;           % fmm target values at the end of the loop
elapsed_time{jj,1} = fmm_timer;
flop_counts{jj} = nflops;

end % end of data loop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% direct evaluation %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% exact_fun_val = cell(10,1);
% exact_y = cell(10,1);
% 
% for jj = 1 %1:length(data2)
% jj
% x = data2{jj,1};  % x -- sources
% y = data2{jj,2};  % y -- targets 
% q = data2{jj,3};
% 
% q = q.^2;
% 
% n = size(y,1);
% 
% exact_timer = zeros(n_iter,1);
% 
% for iter = 1:n_iter
% dtt = tic;
% 
% ne = size(x,1);
% exact_fun = zeros(ne,1);
% exact_der = zeros(ne,1);
% 
% for k = 1:ne
% dn = x - y(k);  % (x-y)
% der_dn = (x - y(k)).^2;  % (x-y)^2
% exact_fun(k) = sum(q./dn); % function estimation
% exact_der(k) = sum(q./der_dn); % derivative estimation
% end
% 
% v1 = 1 + exact_fun;
% ynew = y - (v1 ./exact_der);
% y = ynew;
% exact_timer(iter) = toc(dtt);
% end % end of newton iterations
% exact_fun_val{jj} = v1;
% exact_y{jj} = ynew; 
% elapsed_time{jj,2} = exact_timer;
% 
% end % end of data loop 
% 
% for i1 = 1:10
%      % 5000*i1
%  % [elapsed_time{i1,1}, elapsed_time{i1,2}]
%  % mean(elapsed_time{i1,2}./elapsed_time{i1,1})
% flop_counts{i1};
% end


b_fmm = z1;
A = 1./(x.'-y);
b_exact = A*q;

norm(b_fmm-b_exact,1)./norm(b_exact,1)

end  % end of r loop

% flop_counts
% x1 = 5000:5000:50000;
% y1 = zeros(10,1);
% for i2 = 1:10
% y1(i2) = flop_counts{i2};
% end
% plot(x1,y1)




