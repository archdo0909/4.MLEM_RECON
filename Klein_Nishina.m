function [Pro_scat]=Klein_Nishina(theta, E_o)
    a = E_o/0.511;
    Pro_scat = (1/(1+a*(1-theta)))^2*((1 + theta^2)/2)*( 1+ (a^2*(1-theta)^2/((1+theta^2)*(1+a*(1-theta)))));
end