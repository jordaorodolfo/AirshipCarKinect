O arquivo inicial a ser chamado é o "tslip1.m"

A rotina "ftrack3.m" contém o modelo dinâmico do veículo.

As missões (estipuladas em trajetórias por pontos de passagem)
são definidas na matriz "mis" que tá na na linha 308 da tslip1.m:

mis=         [0  0    0;  % ==> primeiro ponto
             10  0 pi/2;  % ===> segundo ponto
             15  5 pi/2;  % etc...
             10 10    0;
              0 10 pi/2;
             -5  5 pi/2;
              0  0    0];

Cada linha indica dessa matriz indica um ponto de passagem.
A primeira coluna é o Y (Norte), a segunda é o X( Leste) 
e a terceira é o ângulo do segmento entre o ponto atual (linha atual)
e o proximo ponto (linha seguinte).
Se esse terceiro elemento for "0" será feito um segmento de reta,
e se for "pi/2" será feito um arco circulo de 90 graus.
Com isso voce pode criar diferentes tipos de trajetorias, inclusive
uma pista com um "8", etc.
