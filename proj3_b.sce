/*  Pedro Henrique de Queiroz Ramos
    14/11/2022
*/

clear
clc
s=%s;

// Entradas do programa
A=[-.625 -.375; .003 -.2];
B=[12.5; 0];
C=[0 1];
D=[0];
p1=complex(-.3, .4);
p2=complex(-.3, -.4);
vrpm=15; // Velocidade em RPM
c_inicial_x1=2;
c_inicial_x2=3;


tsim=60;                            // Tempo de simulação em segundos
t=linspace(0, tsim, tsim*1000+1);   // Tempo de simulação
vel=vrpm*%pi/30;                    // Velocidade em rad/s
n=size(A)(1);                       // Ordem da planta

// Simulação em realimentação dinâmica de saída
disp('Simulando o sistema em realimentação dinâmica de saída');
// Sistema sem compensação em FTMF
G=syslin('c', A, B, C);

//Passo 1: Verificar se o sistema é controlável
if contr(G.A, G.B) == size(G.A) then
    disp('Passo 1: Sistema é controlável')
else
    disp('Passo 1: Sistema não é controlável')
    abort
end

// Passo 2: Determinar o polinômio alocador desejado
Pa=(s-p1)*(s-p2);
// Aloca os demais polos 10x mais a esquerda
for i=size(coeff(Pa))(2)-1:n-1,
    Pa = Pa *(s+10);
end
disp('Passo 2: Polinômio alocador desejado');
disp('Pa(s) = '+string(Pa));

// Passo 3: Determinar a matriz phi
function P = phi(A, Pa)
    P=0;
    for i=0:size(coeff(Pa))(2)-1,
        P = P + A^i*coeff(Pa)(1, i+1);
    end
endfunction
disp('Passo 3: Matriz de Ackerman');
disp('phi(A) = ', phi(G.A, Pa));

// Passo 4: Determinar a matriz K
K=(-eye(n, n)(n, :))*inv(cont_mat(G.A, G.B))*phi(G.A, Pa);
disp('Passo 4: Matriz K');
disp('K = ', K);

// Passo 5: Determinar a matriz M
M=vel/((-(C)*inv(G.A+G.B*K)*G.B)*norm(K')^2)*K';
disp('Passo 5: Matriz M');
disp('M = ', M);

// Passo 6: Projeto do observador transposto
O_t=cont_mat(A', C');
disp('Passo 6: Matriz O''');
disp('O'' = ', O_t);

// Determinar o polinômio alocador para o controlador dual
// Coloquemos os polos do controlador dual 2x mais a esquerda dos polos do controlador primal
polo_dual=min(real(roots(Pa)'))*n;
Pa_ob=(s-polo_dual)^n;
disp('Polinômio do observador = ', Pa_ob);

// Passo 9: Determinar a matriz L com Ackerman
L=((-eye(n, n)(n, :))*inv(O_t)*phi(G.A', Pa_ob))';
disp('Passo 9: Matriz L');
disp('L = ', L);

// Matrizes do sistema compensado
Ac=[G.A+G.B*K -G.B*K; zeros(n, n) G.A+L*G.C];
Bc=[B*K*M; zeros(n, 1)];
Cc=[G.C];

// Preparando simulação para o Xcos
// Carrega biblioteca de blocos e configuração de simulação
loadXcosLibs(); loadScicos();
// Importa o arquivo Xcos que está no mesmo caminho do script .sce
importXcosDiagram(get_absolute_file_path()+'proj3_b.zcos')

typeof(scs_m)
scs_m.props.context;
// Passagem das variáveis para o SetContext do Xcos
// Valores das condições iniciais para o erro em x1 e x2
Context.c_ini_x1=c_inicial_x1;
Context.c_ini_x2=c_inicial_x2;

// Sistema sem compensação
Context.a11=G.A(1, 1);
Context.a12=G.A(1, 2);
Context.a21=G.A(2, 1);
Context.a22=G.A(2, 2);
Context.b11=G.B(1, 1);

// Sistema com compensação
Context.ac11=Ac(1, 1);
Context.ac12=Ac(1, 2);
Context.ac13=Ac(1, 3);
Context.ac14=Ac(1, 4);
Context.ac21=Ac(2, 1);
Context.ac22=Ac(2, 2);
Context.ac23=Ac(2, 3);
Context.ac24=Ac(2, 4);
Context.ac31=Ac(3, 1);
Context.ac32=Ac(3, 2);
Context.ac33=Ac(3, 3);
Context.ac34=Ac(3, 4);
Context.ac41=Ac(4, 1);
Context.ac42=Ac(4, 2);
Context.ac43=Ac(4, 3);
Context.ac44=Ac(4, 4);
Context.bc11=Bc(1, 1);
Context.bc21=Bc(2, 1);
Context.bc31=Bc(3, 1);
Context.bc41=Bc(4, 1);
Context.cc11=Cc(1, 1);
Context.cc12=Cc(1, 2);

scicos_simulate(scs_m,list(),Context,'nw'); // Simula sem display

// Plotagem das simulações com 2 variáveis de estado
clf;
subplot(3, 2, 1);
plot(out.time, out.values(:,1), 'r', out.time, out.values(:,2), 'b'); xgrid
title('x1 sem compensação')
ylabel('Corrente(Ampére)')
xlabel('Tempo(segundos)')
subplot(3, 2, 3);
plot(out.time, out.values(:,1), 'r', out.time, out.values(:,3), 'b'); xgrid
title('x1 com compensação')
ylabel('Corrente(Ampére)')
xlabel('Tempo(segundos)')
subplot(3, 2, 5);
plot(out.time, out.values(:,1), 'r', out.time, out.values(:,6), 'b'); xgrid
title('Erro em x1')
ylabel('Corrente(Ampére)')
xlabel('Tempo(segundos)')
subplot(3, 2, 2);
plot(out.time, out.values(:,1), 'r', out.time, out.values(:,4), 'b'); xgrid
title('x2 sem compensação')
ylabel('Vel. Angular(rad/segundo)')
xlabel('Tempo(segundos)')
subplot(3, 2, 4);
plot(out.time, out.values(:,1), 'r', out.time, out.values(:,5), 'b'); xgrid
title('x2 com compensação')
ylabel('Vel. Angular(rad/segundo)')
xlabel('Tempo(segundos)')
subplot(3, 2, 6);
plot(out.time, out.values(:,1), 'r', out.time, out.values(:,7), 'b'); xgrid
title('Erro em x2')
ylabel('Vel. Angular(rad/segundo)')
xlabel('Tempo(segundos)')
