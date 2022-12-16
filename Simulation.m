clc, clear, close all
sim_param=1;                      % sim_param = 1: masse della scheda, sim_param = n: masse * n
fprintf(">> Simulazione con sim_param = %d\n\n", sim_param);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inizializzazione dati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pm = 20.5*10^6;                   % Pressione mandata [Pa]
Ps = 0.5*10^6;                    % Pressione di scarico [Pa]
rho = 850;                        % Densita' del fluido [kg/m^3]
ks = 2;                           % Parametro di sintesi
m = 10*sim_param;                 % Massa del martinetto sul portello [kg]
M = 100*sim_param;                % Massa del martinetto sul carrello [kg]
S = 500*10^(-6);                  % Superficie attuatori (p e c) [m^2]
c = 0.1;                          % Corsa del carrello e del portello [m]
F0 = 5000;                        % Coefficiente dell'espressione delle forze agenti sul carrello (cond. eq.) [N]
As = 0.1*10^-6;                   % Area dello strozzamento fisso [m^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametri di simulazione
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 10^-5;                       % Passo di integrazione [s]
Tsim = 10;                        % Tempo di simulazione [s]
t = linspace(0,Tsim,1+Tsim/dt);   % Vettore contenente la scansione temporale della simulazione.
                                  % Tale vettore viene usato come riferimento temporale 
                                  % durante tutta la simulazione. Essso è un vettore colonna contenente
                                  % (1+Tsim/dt) elementi. Ogni elemento è incrementato rispetto 
                                  % al precedente del passo base di simulazione (dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inizializzazione variabili
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Per comodità vengono inizializzate tutte a 0 e poi successivamente
% modificate nel corso della simulazione. Si noti che tutti i vettori hanno
% la medesima dimensione: 1+Tsim/dt, ovvero la 2° componente del vettore 't'

Pa = zeros(size(t, 2),1);         % Pressioni a monte della valvola (sx)
Pb = zeros(size(t, 2),1);         % Pressioni a monte della valvola (dx)
Fc = zeros(size(t, 2),1);         % Forza aerodinamica agente sul carrello
x2p = zeros(size(t, 2),1);        % Accelerazioni del portello (derivata delle velocità)
x4p = zeros(size(t, 2),1);        % Accelerazioni del carrello (derivata delle velocità)
x1p = zeros(size(t, 2),1);        % Velocità del portello (derivata delle posizioni)
x3p = zeros(size(t, 2),1);        % Velocità del carrello (derivata delle posizioni)
errimp_x1 = zeros(size(t,2),1);   % Errore di implementazione sulla posizione del portello
errimp_x2 = zeros(size(t,2),1);   % Errore di implementazione sulla velocità del portello
errimp_x3 = zeros(size(t,2),1);   % Errore di implementazione sulla posizione del carrello
errimp_x4 = zeros(size(t,2),1);   % Errore di implementazione sulla velocità del carrello
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = zeros(size(t, 2),1);         % Posizioni del portello
x2 = zeros(size(t, 2),1);         % Velocità del portello
x3 = zeros(size(t, 2),1);         % Posizioni del carrello
x4 = zeros(size(t, 2),1);         % Velocità del carrello
Pp1 = zeros(size(t, 2),1);        % Pressioni nella camera inferiore dell'attuatore del portello
Pp2 = zeros(size(t, 2),1);        % Pressioni nella camera superiore dell'attuatore del portello
Pc1 = zeros(size(t, 2),1);        % Pressioni nella camera inferiore dell'attuatore del carrello
Pc2 = zeros(size(t, 2),1);        % Pressioni nella camera superiore dell'attuatore del carrello
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input e condizioni iniziali
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xv = t>=1 & t<8;                  % Posizioni della valvola 
                                  % (0: valvola in posizione centrale, 1: mandata a sinistra, scarico a destra)
                                  % xv è un vettore booleano che vale 1 quando t è compreso tra 1 e 8 (secondi)

x1(1) = 0;                        % Condizione iniziale della posizione del portello
x2(1) = 0;                        % Condizione iniziale della velocità del portello
x3(1) = 0;                        % Condizione iniziale della posizione del carrello
x4(1) = 0;                        % Condizione iniziale della velocità del carrello
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulazione 0<=t<1 (fase 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(">> Simulazione fase iniziale ...\n\n");
faseIniziale = ~xv & t<=1; % è un vettore booleano che vale 1 quando xv=0 e t<=1
                           % e 0 altrimenti, in questo modo nei
                           % comandi successivi, un certo valore viene assegnato
                           % ad un certo vettore solo nei "punti" in cui faseIniziale vale 1
x1(faseIniziale) = 0;
x1p(faseIniziale) = 0;
x2(faseIniziale) = 0;
x2p(faseIniziale) = 0;
x3(faseIniziale) = 0;
x3p(faseIniziale) = 0;
x4(faseIniziale) = 0;
x4p(faseIniziale) = 0;
Fc(faseIniziale) = F0.*(c-x1(faseIniziale))./(c+x1(faseIniziale));
Pa(faseIniziale) = (Pm+Ps)/2;
Pb(faseIniziale) = (Pm+Ps)/2;
Pp1(faseIniziale) = Pa(faseIniziale);
Pp2(faseIniziale) = Pb(faseIniziale);
Pc2(faseIniziale) = Pb(faseIniziale);
Pc1(faseIniziale) = Pc2(faseIniziale)+(Fc(faseIniziale)/S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulazione 1<=t<8 (estrazione del portello e del carrello, fase 2a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pa(xv) = Pm;
Pb(xv) = Ps;

fprintf(">> Simulazione estrazione portello ...\n\n");

%% Estrazione del portello
i2a=find(xv==1, 1); % minimo valore per cui xv=1, cioè indica quando inizia la fase 2a
i=i2a; % i2a è tenuto fisso come valore di riferimento
       % nel successivo while si usa un'altro indice denominato 'i'
while x1(i)<c
    x2p(i) = (Pa(i)-Pb(i))*(S/m)-((rho*ks*(S^3)*x2(i)^2)/(2*m*(As^2))); % equazione dell'accelerazione
    x1p(i) = x2(i);
    x2(i+1) = x2(i) + x2p(i)*dt;
    x1(i+1) = x1(i) + x1p(i)*dt;
    errimp_x1(i) = 0.5*x2p(i)*dt^2; % errore di implementazione sulla posizione
    errimp_x2(i) = -0.5*((rho*ks*(S^3)*x2p(i)*x2(i))/(m*(As^2)))*dt^2; % errore di implementazione sulla velocità
    %fprintf(">> t=%f, xp=%f, xpp=%f, xppp=%f \n", 1+dt*i, x1(i), x2(i), x2p(i));
    %pause;
    i=i+1;
end

% Una volta calcolati x1, x2 (=x1p) e x2p si determinano i funzione di questi
% anche gli altri output (y=h(x,u))

estrazionePortello = i2a:i-1; % l'estrazione del portello va dall'indice i2a all'indice i-1 al quale si ferma il while.
                              % Nei successivi comandi solo i valori dei vettori corrispondenti agli indici i2a,
                              % i2a+1, i2a+2, ..., i-1 saranno modificati e/o presi in considerazione                    
x3(estrazionePortello) = 0;
x3p(estrazionePortello) = 0;
x4(estrazionePortello) = 0;
x4p(estrazionePortello) = 0;
Fc(estrazionePortello) = F0.*(c-x3(estrazionePortello))./(c+x3(estrazionePortello));
Pp1(estrazionePortello) = Pa(estrazionePortello);
Pp2(estrazionePortello) = Pb(estrazionePortello)+(rho*ks*S^2*x2(estrazionePortello).^2)/(2*As^2);
Pc1(estrazionePortello) = Pb(estrazionePortello)+Fc(estrazionePortello)/S;
Pc2(estrazionePortello) = Pb(estrazionePortello);

fprintf(">> Simulazione estrazione carrello ...\n\n");

%% Estrazione del carrello
i2b=i; % indice per il quale inizia l'estrazione del carrello
while x3(i)<c % qui è stata apportata una leggera modifica rispetto alla scheda del progetto, in quanto
              % il carrello arriva a fine corsa circa 0.07 secondi prima di quanto previsto nella suddetta
              % dunque si è optato per terminare la fase di estrazione in tale istante,
              % analogamente con la fase precedente (2a)

    Fc(i) = F0.*(c-x3(i))./(c+x3(i)); % equazione delle forze agenti sul carrello
    x4p(i) = (Pa(i)-Pb(i))*(S/M)-((rho*ks*S^3*x4(i)^2)/(2*M*As^2))-(Fc(i)/M); % equazione dell'accelerazione
    x3p(i) = x4(i);
    x4(i+1) = x4(i) + x4p(i)*dt;
    x3(i+1) = x3(i) + x3p(i)*dt;
    errimp_x3(i) = 0.5*x4p(i)*dt^2; % errore di implementazione sulla posizione
    errimp_x4(i) = -0.5*(((rho*ks*(S^3)*x4p(i)*x4(i))/(m*(As^2)))+(-2*F0*c*x4(i))/(c+x3(i))^2)*dt^2; % errore di implementazione sulla velocità
    i=i+1;
end

% Una volta calcolati x3, x4 (=x3p) e x4p si determinano i funzione di questi
% anche gli altri output (y=h(x,u))

estrazioneCarrello = i2b:i-1;  % l'estrazione del carrello va dall'indice i2b all'indice i-1 al quale si ferma il while.
                               % Nei successivi comandi solo i valori dei vettori corrispondenti agli indici i2b,
                               % i2b+1, i2b+2, ..., i-1 saranno modificati e/o presi in considerazione    
x1(estrazioneCarrello,1) = c;
x1p(estrazioneCarrello,1) = 0;
x2(estrazioneCarrello,1) = 0;
x2p(estrazioneCarrello,1) = 0;
Pp1(estrazioneCarrello) = Pa(estrazioneCarrello);
Pp2(estrazioneCarrello) = Pa(estrazioneCarrello);
Pc1(estrazioneCarrello) = Pa(estrazioneCarrello);
Pc2(estrazioneCarrello) = Pb(estrazioneCarrello)+(rho*ks*S^2*x4(estrazioneCarrello).^2)/(2*As^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulazione t>=8 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(">> Simulazione fase finale ...\n\n");
faseFinale = i:size(t,2); % la fase finale parte dall'indice i al quale termina il while precedente
                          % e si protrae fino alle fine della simulazione (cioè alla dimensione del vettore t)
                          % Tali indici sono gli unici che vengono presi in considerazione per le operazioni di
                          % assegnazione tra vettori nella righe successive
x1(faseFinale) = c;
x1p(faseFinale) = 0;
x2(faseFinale) = 0;
x2p(faseFinale) = 0;
x3(faseFinale) = c;
x3p(faseFinale) = 0;
x4(faseFinale) = 0;
x4p(faseFinale) = 0;
Pa(faseFinale) = (Pm+Ps)/2;
Pb(faseFinale) = (Pm+Ps)/2;
Fc(faseFinale) = F0.*(c-x3(faseFinale))./(c+x3(faseFinale));
Pa(faseIniziale) = (Pm+Ps)/2;
Pb(faseIniziale) = (Pm+Ps)/2;
Pp1(faseFinale) = Pa(faseFinale);
Pp2(faseFinale) = Pb(faseFinale);
Pc2(faseFinale) = Pb(faseFinale);
Pc1(faseFinale) = Pc2(faseFinale)+(Fc(faseFinale)/S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot dei risultati 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(">> Plot dei risultati ...\n\n");
figure Name 'Posizioni e velocità'
subplot(211)
plot(t, x1, t, x3, 'LineWidth', 1.5)
legend('Portello','Carrello')
xlabel('t [s]')
ylabel('Spostamento [m]')
subplot(212)
plot(t, x2, t, x4, 'LineWidth', 1.5)
legend('Portello', 'Carrello')
xlabel('t [s]')
ylabel('Velocità [m/s]')

figure Name 'Pressioni attuatori'
subplot(211)
plot(t, Pp2, t, Pp1, 'LineWidth', 1.5)
set(gca,'YScale','log') % per le pressioni il plot con ordinate logaritmiche  
                        % è risultato descrivere meglio i risultati ottenuti
                        % rispetto ad un plot con ordinate lineari
axis([0 10 2*10^5 10^9])
legend('Camera superiore attuatore portello','Camera inferiore attuatore portello')
xlabel('t [s]')
ylabel('Pressione [Pa]')
subplot(212)
plot(t, Pc2, t, Pc1, 'LineWidth', 1.5)
set(gca,'YScale','log') % per le pressioni il plot con ordinate logaritmiche  
                        % è risultato descrivere meglio i risultati ottenuti 
                        % rispetto ad un plot con ordinate lineari
legend('Camera superiore attuatore carrello','Camera inferiore attuatore carrello', 'Location','northwest')
xlabel('t [s]')
ylabel('Pressione [Pa]')
axis([0 10 2*10^5 10^9]) 

figure Name 'Errore x1'
subplot(211)
plot(t, errimp_x1, 'LineWidth', 1.5)
title('Errore di implementazione sulla posizione del portello')
xlabel('t [s]')
ylabel('Errore [m]')
subplot(212)
plot(t(i2a-5:i2a+21), errimp_x1(i2a-5:i2a+21), 'LineWidth', 1.5)
title('Ingrandimento sul tratto significativo')
xlabel('t [s]')
ylabel('Errore [m]')

figure Name 'Errore x2'
subplot(211)
plot(t, errimp_x2, 'LineWidth', 1.5)
title('Errore di implementazione sulla velocità del portello')
xlabel('t [s]')
ylabel('Errore [m/s]')
subplot(212)
plot(t(i2a-5:i2a+20), errimp_x2(i2a-5:i2a+20), 'LineWidth', 1.5)
title('Ingrandimento sul tratto significativo')
xlabel('t [s]')
ylabel('Errore [m/s]')

figure Name 'Errore x3'
subplot(211)
plot(t, errimp_x3, 'LineWidth', 1.5)
title('Errore di implementazione sulla posizione del carrello')
xlabel('t [s]')
ylabel('Errore [m]')
subplot(212)
plot(t(i2b-200:i2b+500), errimp_x3(i2b-200:i2b+500), 'LineWidth', 1.5)
title('Ingrandimento sul tratto significativo')
xlabel('t [s]')
ylabel('Errore [m]')

figure Name 'Errore x4'
subplot(211)
plot(t, errimp_x4, 'LineWidth', 1.5)
title('Errore di implementazione sulla velocità del carrello')
xlabel('t [s]')
ylabel('Errore [m/s]')
subplot(212)
plot(t(i2b-200:i2b+500), errimp_x4(i2b-200:i2b+500), 'LineWidth', 1.5)
title('Ingrandimento sul tratto significativo')
xlabel('t [s]')
ylabel('Errore [m/s]')
