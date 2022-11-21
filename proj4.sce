/*  Pedro Henrique de Queiroz Ramos
    20/11/2022
*/

clear
clc
s=%s;

// Entradas do programa

R = .05;
b = 2;

A=[-.625 -.375; 
    .003 -.2];
    
B=[12.5; 
    0];
    
C=[0 1];

rho = 50;

t = linspace( 0 , 30 , 3001 );

// a) Determinar o mínimo do critério quadrático para

Q = [R  0; 
     0  b];

// Determinação do compensador

P = riccati( A , inv( rho ) * B * B' , Q , 'c' );
Ka = -inv( rho ) * B' * P;

// Determinação da matriz de ganho de referência

Ma = -1 / ( (C * inv( A + B * Ka ) * B) * norm(Ka')^2 ) * Ka';

// Montagem do sistema em malha fechada

sys = syslin( 'c' , A + B * Ka , B , eye( A ));

// Simulação

ya = csim( zeros( t ) , t , sys , [ .5 3 ]' );

// b) Determinar o mínimo do critério quadrático para

Q = [R  0; 
     0  0];

// Determinação de V baseado em Q

V = [sqrt(.05) 0]; // v1^2 = R

// Determinação do LR simétrico

phi = V * inv(s*eye(A) - A ) * B;       // função phi(s)
phi_ = V * inv((-s)*eye(A) - A ) * B;   // função phi(-s)

Pol = 1 + rho^(-1) * phi * phi_;        // Polinômio do LR simétrico
rPol = roots( Pol.den );                // Raízes do denominador de Pol

// Lugar das raízes simétrico

//evans( Pol );

// Obtém todas as raízes do SPE

Sd=[];
for i = 1 : size(rPol)(1),                  // Percorre todas as raízes de Pol
    if ( rPol(i, 1) < 0 ) then              // Se a raíz atual estiver no SPE
        Sd( 1, size(Sd)(1, 2) + 1 ) = rPol(i, 1); // Armazena em Sd
    end
end

// Determinação do compensador

Kb = -ppol( A, B, Sd );

// Determinação da matriz de ganho de referência

Mb = -1 / ( (C * inv( A + B * Kb ) * B) * norm(Kb')^2 ) * Kb';

// Montagem do sistema em malha fechada

sys = syslin( 'c' , A + B * Kb , B , eye( A ));

// Simulação

yb = csim( zeros( t ) , t , sys , [ .5 3 ]' );

// Plotagem dos gráficos

subplot( 121 );
plot( t, ya( 1, : ), 'b', t, ya( 2, : ), 'r'); xgrid;
legend(['x1'; 'x2'])
title('a) Minimização com ARE');
subplot( 122 );
plot( t, yb( 1, : ), 'b', t, yb( 2, : ), 'r'); xgrid;
legend(['x1'; 'x2'])
title('b) Minimização com LR simétrico');
