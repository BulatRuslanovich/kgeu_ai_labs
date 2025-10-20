clc;
clear;

function mu = trapmf_val(x, p)
	a = p(1); b = p(2); c = p(3); d = p(4);
	mu = zeros(x);
	for i = 1:length(x)
		xi = x(i);
		if xi <= a then
			mu(i) = 0;
		elseif xi >= d then
			mu(i) = 0;
		elseif xi >= b & xi <= c then
			mu(i) = 1;
		elseif xi > a & xi < b then
			mu(i) = (xi - a) / (b - a);
		else // xi > c & xi < d
			mu(i) = (d - xi) / (d - c);
		end
	end
endfunction

function mu = trimf_val(x, p)
	a = p(1); b = p(2); c = p(3);
	mu = zeros(x);
	for i = 1:length(x)
		xi = x(i);
		if xi <= a | xi >= c then
			mu(i) = 0;
		elseif xi == b then
			mu(i) = 1;
		elseif xi > a & xi < b then
			mu(i) = (xi - a) / (b - a);
		else
			mu(i) = (c - xi) / (c - b);
		end
	end
endfunction


TEMP_C  = [0 0 10 30];      
TEMP_SW = [20 35 50];       
TEMP_W  = [40 50 60];       
TEMP_SH = [50 60 70];       
TEMP_H  = [60 70 100 100];  


CORN_BL = [-90 -90 -72 -36]; 
CORN_SL = [-54 -27 0];       
CORN_Z  = [-18 0 18];        
CORN_SR = [0 27 54];         
CORN_BR = [36 72 90 90];     


corn_x = -90:0.5:90;

function mu = corn_mf(termName)
	select termName
		case "BL" then mu = trapmf_val(corn_x, CORN_BL);
		case "SL" then mu = trimf_val(corn_x, CORN_SL);
		case "Z"  then mu = trimf_val(corn_x, CORN_Z);
		case "SR" then mu = trimf_val(corn_x, CORN_SR);
		case "BR" then mu = trapmf_val(corn_x, CORN_BR);
		otherwise mu = zeros(corn_x);
	end
endfunction

ruleTempIdx = [5 4 3 2 1];
ruleTerm    = ["BR" "SR" "Z" "SL" "BL"];

function mus = fuzzify_temp(t)
	mus = zeros(1,5);
	mus(1) = trapmf_val([t], TEMP_C)(1);
	mus(2) = trimf_val([t], TEMP_SW)(1);
	mus(3) = trimf_val([t], TEMP_W)(1);
	mus(4) = trimf_val([t], TEMP_SH)(1);
	mus(5) = trapmf_val([t], TEMP_H)(1);
endfunction

function agg = aggregated_output(t)
	alpha = fuzzify_temp(t); // [C, SW, W, SH, H]
	agg = zeros(corn_x);
	for r=1:length(ruleTempIdx)
		tempIdx = ruleTempIdx(r);
		term = ruleTerm(r);
		w = alpha(tempIdx);
		if w > 0 then
			mu_out = corn_mf(term);
			agg = max(agg, min(w, mu_out));
		end
	end
endfunction

function y = eval_corn_mamdani(t)
	agg = aggregated_output(t);
	num = sum(corn_x .* agg);
	den = sum(agg);
	if den == 0 then
		y = 0;
	else
		y = num / den;
	end
endfunction

T_test = 55;
Corn_test = eval_corn_mamdani(T_test);

mC  = trapmf_val([T_test], TEMP_C)(1);
mSW = trimf_val([T_test], TEMP_SW)(1);
mW  = trimf_val([T_test], TEMP_W)(1);
mSH = trimf_val([T_test], TEMP_SH)(1);
mH  = trapmf_val([T_test], TEMP_H)(1);

active_rules = [];
if mH  > 0 then active_rules = [active_rules, 1]; end
if mSH > 0 then active_rules = [active_rules, 2]; end
if mW  > 0 then active_rules = [active_rules, 3]; end
if mSW > 0 then active_rules = [active_rules, 4]; end
if mC  > 0 then active_rules = [active_rules, 5]; end

mprintf("Temp=%.3f -> Corn=%.3f\n", T_test, Corn_test);
mprintf("Memberships at T=55: C=%.3f, SW=%.3f, W=%.3f, SH=%.3f, H=%.3f\n", mC, mSW, mW, mSH, mH);
mprintf("Active rules at T=55: ");
if size(active_rules, "*") > 0 then
	for k=1:size(active_rules, "*")
		if k>1 then mprintf(", "); end
		mprintf("%d", active_rules(k));
	end
else
	mprintf("none");
end
mprintf("\n");

Corn_target = 40;

function err = target_error(t)
	y = eval_corn_mamdani(t);
	err = abs(y - Corn_target);
endfunction

Ts = 0:0.5:100;
errs = zeros(Ts);
for i=1:length(Ts)
	errs(i) = target_error(Ts(i));
end
[emin, idx] = min(errs);
T0 = Ts(idx);

step = 0.05;
bestT = T0;
bestE = target_error(bestT);
for t = max(0, T0-5):step:min(100, T0+5)
	et = target_error(t);
	if et < bestE then
		bestE = et;
		bestT = t;
	end
end

Corn_best = eval_corn_mamdani(bestT);

bC  = trapmf_val([bestT], TEMP_C)(1);
bSW = trimf_val([bestT], TEMP_SW)(1);
bW  = trimf_val([bestT], TEMP_W)(1);
bSH = trimf_val([bestT], TEMP_SH)(1);
bH  = trapmf_val([bestT], TEMP_H)(1);

active_rules_best = [];
if bH  > 0 then active_rules_best = [active_rules_best, 1]; end
if bSH > 0 then active_rules_best = [active_rules_best, 2]; end
if bW  > 0 then active_rules_best = [active_rules_best, 3]; end
if bSW > 0 then active_rules_best = [active_rules_best, 4]; end
if bC  > 0 then active_rules_best = [active_rules_best, 5]; end

mprintf("Variant 7 target Corn=+40 -> best Temp=%.3f gives Corn=%.3f (|err|=%.3f)\n", bestT, Corn_best, bestE);
mprintf("Active rules at best Temp: ");
if size(active_rules_best, "*") > 0 then
	for k=1:size(active_rules_best, "*")
		if k>1 then mprintf(", "); end
		mprintf("%d", active_rules_best(k));
	end
else
	mprintf("none");
end
mprintf("\n");

try
	temp_x = 0:0.5:100;
	muC  = trapmf_val(temp_x, TEMP_C);
	muSW = trimf_val(temp_x, TEMP_SW);
	muW  = trimf_val(temp_x, TEMP_W);
	muSH = trimf_val(temp_x, TEMP_SH);
	muH  = trapmf_val(temp_x, TEMP_H);
	scf(); clf();
	plot(temp_x, muC, 'r-');
	hold on;
	plot(temp_x, muSW, 'g-');
	plot(temp_x, muW, 'b-');
	plot(temp_x, muSH, 'm-');
	plot(temp_x, muH, 'k-');
	title('Функции принадлежности Temp');
	xlabel('Temp, °C'); ylabel('μ');
	legend(['C','SW','W','SH','H'], -1);
	xs2png(gcf(), 'temp_mf.png');
catch
end

try
	muBL = corn_mf('BL');
	muSL = corn_mf('SL');
	muZ  = corn_mf('Z');
	muSR = corn_mf('SR');
	muBR = corn_mf('BR');
	scf(); clf();
	plot(corn_x, muBL, 'r-');
	hold on;
	plot(corn_x, muSL, 'g-');
	plot(corn_x, muZ,  'b-');
	plot(corn_x, muSR, 'm-');
	plot(corn_x, muBR, 'k-');
	title('Функции принадлежности Corn');
	xlabel('Corn, градусы'); ylabel('μ');
	legend(['BL','SL','Z','SR','BR'], -1);
	xs2png(gcf(), 'corn_mf.png');
catch
end


try
	agg55 = aggregated_output(55);
	scf(); clf();
	plot(corn_x, agg55, 'b-');
	hold on;
	plot([Corn_test Corn_test],[0 1],'r--');
	title('Агрегированная МФ Corn при Temp=55');
	xlabel('Corn, градусы'); ylabel('μ');
	xs2png(gcf(), 'agg_55.png');
catch
end

try
	aggBest = aggregated_output(bestT);
	scf(); clf();
	plot(corn_x, aggBest, 'b-');
	hold on;
	plot([Corn_best Corn_best],[0 1],'r--');
	title('Агрегированная МФ Corn при Temp=best');
	xlabel('Corn, градусы'); ylabel('μ');
	xs2png(gcf(), 'agg_var7.png');
catch
end
