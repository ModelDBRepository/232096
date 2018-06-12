%% Figure: LTP Hippocampal DPhil data vs statLTSP predictions 

clear all;
close all;
loadColors


%% Init
titles = 'Data: Hippocampal LTP';
file = 'Larkman1992_LTP.xls';

[data, txt] = xlsread(file);

Ni = 1;
p1i = 2;
p2i = 3;
q1i = 4;
q2i = 5;


q1_i = data(:,q1i);
p1_i = data(:,p1i);
q2_i = data(:,q2i);
p2_i = data(:,p2i);
Ns = data(:,Ni);


dKL_dp = @(p,q,t)  ((1-2.*p)./(2.*p-2.*p.^2)  +  ((p.*q-t).*(p.*(q-2.*t)+t))./(2.*(1-p).^2 .* p.^2 .* q.^2));
dKL_dq = @(p,q,t)  ((1./q)  +  (t.*(t-p.*q)) ./ ((p-1).*p.*q.^3));
KL = @(p,q,t) log(sqrt((q.^2).*p.*(1-p))) + ((t-p.*q).^2)./(2.*((q.^2).*p.*(1-p)));


LTPth = 1.05;
idxLTD = find(((q2_i./q1_i) .* (p2_i./p1_i))<=1);
idxLTP = find(((q2_i./q1_i) .* (p2_i./p1_i))>LTPth);

qscale = mean(q1_i(idxLTP));
pscale = 1;

q1_i(idxLTP) = q1_i(idxLTP)./qscale;
q2_i(idxLTP) = q2_i(idxLTP)./qscale;
p1_i(idxLTP) = p1_i(idxLTP)./pscale;
p2_i(idxLTP) = p2_i(idxLTP)./pscale;


lr = -0.0001;
qmin = 0.01;
pmin = 0.01;

t = 3.4697;


%%
nrows = 1;
ncols = 3;

% Matlab figure representation
fig_ind = reshape(1:nrows*ncols, ncols, nrows)';
A = 1;
B = A+1;
C = B+1;
fig_subind = [A B C];

A4_width = 600;
figure('Name', titles, 'Position',  [1046 287 A4_width 250-100]);
colors = [red1; blue1];



%% Panel A: Quiver plot: pred vs obs
fig_a = subplot(nrows, ncols, fig_ind(fig_subind==A), 'replace'); %Weight



%Normalize data

% Estimate all p&q
estpF = zeros(size(idxLTP));
estqF = zeros(size(idxLTP));

t = ones(size(idxLTP)).*t';
tq = t;


    eN = 5.5;

    tt = t(1);

    minp = 0.01;
    minq = 0.01;
    maxp = 1;
    %maxq = 1;
    maxq = max(q2_i);


    npartsq = 1;
    npartsp = 1;
    incp = 0.01;
    incq = 0.02;

    nq = maxq/npartsq;
    np = maxp/npartsp;

    for i=1:npartsq

        for j=1:npartsp

            Ps = [minp + np*(j-1): incp : maxp - np*(npartsp-j)];
            Qs = [minq + nq*(i-1): incq : maxq - nq*(npartsq-i)];

            KL_M = zeros(length(Ps), length(Qs));
            for x=1:length(Ps)
                for y=1:length(Qs)
                    KL_M(x,y) = -KL(Ps(x), Qs(y), tt);
                end
            end

            Y = repmat(Ps, size(Qs,2),1)';
            X = repmat(Qs, size(Ps,2),1);

            [U,V] = gradient(KL_M,0.1);

            hold on
            hsl = streamslice(X,Y,U,V, 0.5, 'arrows');
        end    
    end
    
    for i=1:size(hsl,1)
        hsl(i).Color = [0.9 0.9 0.9];
        hsl(i).LineWidth = 0.75;
    end
    

tmpd_q = 0;
tmpd_p = 0;
q1_i = q1_i - tmpd_q;
p1_i = max(p1_i - tmpd_p,0.02);
q2_i = q2_i - tmpd_q;
p2_i = p2_i - tmpd_p;

for i=1:size(idxLTP,1)
    
    estp = [];
    estq = [];
    estKL = [];
    estp(1) = p1_i(idxLTP(i));
    estq(1) = q1_i(idxLTP(i));
    estKL(1) = KL(estp(1), estq(1), tq(i));
    
    obsmean2 = (p2_i(idxLTP(i))*q2_i(idxLTP(i)));
    KL1 = KL(p1_i(idxLTP(i)), q1_i(idxLTP(i)), tq(i));
    KL2 = KL(p2_i(idxLTP(i)), q2_i(idxLTP(i)), tq(i));

    z = 1;
    while((estp(z)*estq(z))<obsmean2)
        z = z+1;
        
        estp(z) = estp(z-1) + lr.*dKL_dp(estp(z-1), estq(z-1), t(i));
        estq(z) = estq(z-1) + lr.*dKL_dq(estp(z-1), estq(z-1), tq(i));
        
        estKL(z) = KL(estp(z), estq(z), tq(i));
    end
    estpF(i) = estp(end);
    estqF(i) = estq(end);

end

scale = 0.5;
h1 = quiver(q1_i(idxLTP), p1_i(idxLTP), (q2_i(idxLTP) - q1_i(idxLTP)).*scale, (p2_i(idxLTP) - p1_i(idxLTP)).*scale, 0, 'Color', (red1+blue2)./2, 'LineWidth', 1.5);
h1.MaxHeadSize = 0.25;
hold on;
h2 = quiver(q1_i(idxLTP), p1_i(idxLTP), (estqF - q1_i(idxLTP)).*scale, (estpF - p1_i(idxLTP)).*scale, 0, 'Color', black1, 'LineWidth', 1.5);
h2.MaxHeadSize = 0.25;

KL2_HCLTP = KL(p2_i(idxLTP), (q2_i(idxLTP)), t);
KL1_HCLTP = KL(p1_i(idxLTP), (q1_i(idxLTP)), t);

title(fig_a, 'hippocampal LTP', 'fontsize', fontsize);
xlabel('q, quantal amp. (mV)', 'fontsize', minifontsize);
ylabel('P_{rel}, release prob.', 'fontsize', minifontsize);
box off
hla=legend([h1 h2], {'data', 'pred'}, 'Location', 'SouthEast');
legend('boxoff');
hla.Position = [0.2617    0.2673    0.1000    0.1767];

xlim(fig_a, [0.5 2.2]);

xticks = get(fig_a, 'XTick');
aux = num2cell(round(xticks.*qscale./ 1e3, 2));
set(fig_a, 'XTickLabel', aux);

ylim(fig_a, [0 1]);
set(fig_a, 'TickDir', 'out');





%% Panel B: q2-q1 & p2-p1
fig_b = subplot(nrows, ncols, fig_ind(fig_subind==B), 'replace'); %Weight

set(fig_b, 'TickDir', 'out');


qdif = q2_i(idxLTP)./q1_i(idxLTP).*100;
qdifp = estqF./q1_i(idxLTP).*100;
mine = min([qdif; qdifp]) - 0.05;
maxe = max([qdif; qdifp]) + 0.05;
h=plot([mine-5 maxe], [mine-5 maxe], ':k');
h.HandleVisibility= 'off';
hold on

h1 = scatter(qdif, qdifp, ssize, red1);
h1.MarkerEdgeColor = red1;

hold on
xlim([mine maxe])
ylim([mine maxe])
[rq,pq] = corrcoef(qdif, qdifp);

qstars = 'n.s.';
if(pq(1,2)<0.001)
    qstars = '***';
elseif(pq(1,2)<0.01)
    qstars = '**';
elseif(pq(1,2)<0.05)
    qstars = '*';
end

pdif = p2_i(idxLTP)./p1_i(idxLTP).*100;
pdifp = estpF./p1_i(idxLTP).*100;

h2 = scatter(pdif, pdifp, ssize);
h2.MarkerEdgeColor  = blue1;

h=lsline;
h(1).Color = blue2;
h(1).LineWidth = 1.5;

h(2).Color = red2;
h(2).LineWidth = 1.5;

xlabel('parameter_{data} (%)', 'fontsize', minifontsize)
ylabel('parameter_{pred} (%)', 'fontsize', minifontsize)
xlim([mine-10 maxe])
ylim([mine-10 maxe])
[rp,pp] = corrcoef(pdif, pdifp);


pstars = 'n.s.';
if(pp(1,2)<0.001)
    pstars = '***';
elseif(pp(1,2)<0.01)
    pstars = '**';
elseif(pp(1,2)<0.05)
    pstars = '*';
end      

%Assess significance: obs vs pred
ppc = ranksum(pdif, pdifp);
qpc = ranksum(qdif, qdifp);
disp(['HC LTP; p=' num2str(ppc) ', q=' num2str(qpc)]);

lh = legend([h1 h2], {[''], ['']}, 'Location', 'NorthEast');
lh.Position(2) = lh.Position(2)+0.065;
lh.Position(1) = lh.Position(1)-0.105;

th = text(96, 197, ['P_{rel}']);
th.FontSize = minifontsize-2;
th = text(170, 100, ['r=' num2str(rp(1,2),2) ' ' pstars], 'Color', blue1);

th = text(96, 222, ['q']);
th.FontSize = minifontsize-2;
th = text(170, 125, ['r=' num2str(rq(1,2),2) ' ' qstars], 'Color', red1);


legend('boxoff')

box off

set(fig_b, 'TickDir', 'out');




    
    

%% Panel C: obs-pred theta density
fig_c = subplot(nrows, ncols, fig_ind(fig_subind==C), 'replace'); %Weight

before_v_obs = [q1_i(idxLTP), p1_i(idxLTP)];
after_v_obs = [q2_i(idxLTP), p2_i(idxLTP)];
after_v_pred = [estqF, estpF];

deg_dif = zeros(size(before_v_obs,1), 1);
for i=1:size(before_v_obs,1)    
    u = [after_v_obs(i,:) - before_v_obs(i,:)];
    v = [after_v_pred(i,:) - before_v_obs(i,:)];
    
    CosTheta = dot(u,v)/(norm(u)*norm(v));
    deg_dif(i) = acos(CosTheta)*180/pi;
end




delta = [after_v_obs(:,1).*after_v_obs(:,2) - before_v_obs(:,1).*before_v_obs(:,2)];
after_v_predNaive = [delta./2 delta./2]; 
deg_difNaive = zeros(size(before_v_obs,1), 1);
for i=1:size(before_v_obs,1)    
    u = [after_v_obs(i,:) - before_v_obs(i,:)];
    v = [after_v_predNaive(i,:)];
    
    CosTheta = dot(u,v)/(norm(u)*norm(v));
    deg_difNaive(i) = acos(CosTheta)*180/pi;
end

[f,xi] = ksdensity(deg_difNaive);



after_v_pred = [];
% Rel change control
deg_difControl1 = zeros(size(before_v_obs,1), 1);
    p = [0.01:0.01:0.9 1];
    newp = zeros(size(idxLTP));
    newq = zeros(size(idxLTP));
    for i=1:size(idxLTP,1)
        change = (q2_i(idxLTP(i)).*p2_i(idxLTP(i)))-(q1_i(idxLTP(i)).*p1_i(idxLTP(i)));
        
        We = q1_i(idxLTP(i)).*p1_i(idxLTP(i));
        We = We + change;
        
        qtmp = We./p;
        d = sqrt(((qtmp - q1_i(idxLTP(i))).^2) + (p - p1_i(idxLTP(i))).^2); %Abs
        id = find(d==min(d));
        newq(i) = qtmp(id);
        newp(i) = p(id);    
        
        u = [after_v_obs(i,:) - before_v_obs(i,:)];
        after_v_pred(i,:) = [newq(i), newp(i)];
        v = [after_v_pred(i,:) - before_v_obs(i,:)];

        CosTheta = dot(u,v)/(norm(u)*norm(v));
        deg_difControl1(i) = real(acos(CosTheta)*180/pi);
    end




nruns = 1000;
deg_difRandn = zeros(size(before_v_obs,1), nruns);
deg_difRandn_norm = zeros(size(before_v_obs,1), nruns);
deg_difRandn_normCtrl = zeros(size(before_v_obs,1), nruns);
for r=1:nruns
    after_v_predRandn = [abs(randn(size(idxLTP))) abs(randn(size(idxLTP)))];
    for i=1:size(before_v_obs,1)    
        u = [after_v_obs(i,:) - before_v_obs(i,:)];
        v = [after_v_predRandn(i,:)];

        CosTheta = dot(u,v)/(norm(u)*norm(v));
        deg_difRandn(i,r) = acos(CosTheta)*180/pi;
    end
    deg_difRandn_norm(:,r) = deg_difRandn(:,r)./deg_dif;
    deg_difRandn_normCtrl(:,r) = deg_difRandn(:,r)./deg_difControl1;
end    

[f,xi] = ksdensity(reshape(deg_difRandn, size(deg_difRandn, 1).*size(deg_difRandn, 2), 1), 'Bandwidth', 10);
hRand = plot(xi,f, 'linewidth', 2, 'Color', orange1);
hold on;

[f,xi] = ksdensity(reshape(deg_difControl1, size(deg_difControl1, 1).*size(deg_difControl1, 2), 1));
hCtrl = plot(xi,f, ':k', 'linewidth', 2, 'Color', (orange2.*0.5+orange1)./2);   

[f,xi] = ksdensity(deg_dif);
hModel = plot(xi,f, 'linewidth', 2, 'Color', black1);


ylabel('freq.');
xlabel('angle_{data, pred}', 'fontsize', minifontsize)
xlim(fig_c, [0 125]);

    [ppr,h] = ranksum(deg_dif, reshape(deg_difRandn, size(deg_difRandn, 1).*size(deg_difRandn, 2), 1));
    aux = reshape(deg_difRandn_norm, size(deg_difRandn_norm, 1).*size(deg_difRandn_norm, 2), 1);
    [ppr_norm,h] = ranksum(aux./aux, aux);
    [ppr_norm2,h] = ranksum(deg_dif./deg_dif, reshape(deg_difRandn_normCtrl, size(deg_difRandn, 1).*size(deg_difRandn, 2), 1));
    [pps,h] = ranksum(deg_dif./deg_dif, deg_difControl1./deg_dif);
    
    prstars = 'n.s.';
    if(ppr_norm<0.001)
        prstars = '***';
    elseif(ppr_norm<0.01)
        prstars = '**';
    elseif(ppr_norm<0.05)
        prstars = '*';
    end
    
    psstars = 'n.s.';
    if(pps<0.001)
        psstars = '***';
    elseif(pps<0.01)
        psstars = '**';
    elseif(pps<0.05)
        psstars = '*';
    end
    

lh=legend([hModel, hRand, hCtrl], {'statLTSP', ['random (' prstars ')'], ['shortest (' psstars ')']}, 'Location', 'NorthEast');
lh.Position(1) = lh.Position(1)+0.015;
lh.Position(2) = lh.Position(2)+0.04;
legend('boxoff');
box off;

set(fig_c, 'TickDir', 'out');

set(gcf, 'Color', 'w')