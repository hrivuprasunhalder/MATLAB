%      initial final  average time        
data= [0       89       72;
       89      139      124;
       139     192      75;
       192     229      122.30;
       229     288      49;
       288     339      56;
       339     380      73.20;
       380     416      113.50;
       416     443      249.80;
       443     481      103.60;
       481     519      126.40;
       519     582      81.50;
       582     656      88.00;
       656     729      143.10;
       729     999      123.50;];

l= .000032;
beta_eff= 0.007;
beta= [.000231 .00153 .001372 .00276 .000805 .000294 ];
lambda = [0.0124, 0.0305, 0.1115, 0.301, 1.138, 3.01];

rho_dollar= zeros(15,1);
diff_worth= zeros(15,1);
int_worth= zeros(15,1);
pcm_conv= beta_eff*1e5;

for i= 1:15
    x1= data(i,1);
    x2= data(i,2);
    t= data(i,3);
    % as known P= P0e^(t/T),where T is the reactor period. P0=200w and
    % P=800W. Thus e^(t/T)=4. Taking ln we get, t/T= ln4
    T= t/log(4);
    sum_term= 0;
    for j= 1:6
        sum_term= beta(j)/(1+(lambda(j)*T));
    end
    rho_dollar(i)=  (1 / beta_eff) * ((l/T + sum_term) / (1 + l/T));
    % Differential worth/step 
    dx= x2-x1;
    diff_worth(i)= rho_dollar(i)/dx;
    %Integral worth
    if i==1
        int_worth(i)= rho_dollar(i);
    else
        int_worth(i)= rho_dollar(i) + rho_dollar(i-1);
    end
end    


diff_worth_pcm = diff_worth * pcm_conv;
int_worth_pcm = int_worth * pcm_conv;

disp(' Stroke  dx     avg t(s)      rho($)    differential worth($/step)   integral worth($)');
for i = 1:15
    fprintf('%3d    %4.0f    %7.2f    %8.5f    %8.5f    %8.5f\n', i, data(i, 2)-data(i, 1), data(i, 3), rho_dollar(i), diff_worth(i), int_worth(i));
end

fprintf('\nTotal Integral Worth = %.5f $ = %.2f pcm\n', int_worth(end), int_worth(end)*pcm_conv);