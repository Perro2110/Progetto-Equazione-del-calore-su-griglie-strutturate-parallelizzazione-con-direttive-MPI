 % Read Heat2D.f90 output file 
 clear  all 
 close  all 
 clc 
 % INITIAL CONDITION 
 % Open file 
 FileID  =  fopen  (  'Heat2D_output-0000-0000.dat'  ); 
 % Read array dimension 
 % Read data 
 N  =  fscanf  (  FileID  ,  '%d \n'  ,  1  ); 
 TempoI  =  fscanf  (  FileID  ,  '%d \n'  ,  1  ); 
 nCPU  =  fscanf  (  FileID  ,  '%d \n'  ,  1  ); 
 x  =  fscanf  (  FileID  ,  '%f \n'  ,  N  ); 
 y  =  fscanf  (  FileID  ,  '%f \n'  ,  N  ); 
 T0  =  zeros  (  N  ,  N  ); 
 TR  =  50 
 TL  =  100 
 t  =  0.05  %  tempo  finale 
 for  k  =  1  :  N 
 T0  (  k  ,:) =  fscanf  (  FileID  ,  '%f \n'  ,  N  )'; 
 end 
 % FINAL SOLUTION 
 % Open file 
 FileID  =  fopen  (  'Heat2D_output-0003-0000.dat'  );  %inserire  nome file a t(end) 
 % Read array dimension 
 N  =  fscanf  (  FileID  ,  '%d \n'  ,  1  ); 
 TempoF  =  fscanf  (  FileID  ,  '%d \n'  ,  1  ); 
 nCPU  =  fscanf  (  FileID  ,  '%d \n'  ,  1  ); 
 x  =  fscanf  (  FileID  ,  '%f \n'  ,  N  ); 
 y  =  fscanf  (  FileID  ,  '%f \n'  ,  N  ); 
 T1  =  zeros  (  N  ,  N  ); 
 Te  =  0.5  .* (  TR  +  TL  )+  0.5  .*  erf  ((  x  -  1  )/(  2  .*  sqrt  (  t  )))  .* (  TR  -  TL  ); 
 Tt  =  repmat  (  Te' ,N,1); 
 for  k  =  1  :  N 
 
 T1  (  k  ,:) =  fscanf  (  FileID  ,  '%f \n'  ,  N  )' ; 
 end 
 Tdiff  =  abs  (  Tt  -  T1  )./  abs  (  Tt  ); 
 fprintf  (  "Numero di Cpu: %g  \n"  ,  nCPU  ) 
 fprintf  (  "Tempo di esecuzione: %g \n"  ,  TempoF  ) 
 fprintf  (  "Valore max abs: %g \n"  ,  max  (  max  (  abs  (  Tdiff  )))) 
 % Plot data 
 fg  =  figure  (  1  ); 
 fg  .  Position  = [  100  100  1300  650  ]; 
 tiledlayout  (  2  ,  2  ) 
 nexttile 
 surf  (  x  ,  y  ,  T0  ) 
 view  (  45  ,  22  ) 
 ylim  ([  0  2  ]) 
 xlim  ([  0  2  ]) 
 title  (  'Initial condition'  ) 
 nexttile 
 surf  (  x  ,  y  ,  T1  ) 
 ylim  ([  0  2  ]) 
 xlim  ([  0  2  ]) 
 view  (  45  ,  22  ) 
 title  (  'Final condition '  ) 
 nexttile 
 surf  (  x  ,  y  ,  Tt  ) 
 ylim  ([  0  2  ]) 
 xlim  ([  0  2  ]) 
 view  (  45  ,  22  ) 
 title  (  'Final condition teorica'  ) 
 nexttile 
 surf  (  x  ,  y  ,  Tdiff  ) 
 ylim  ([  0  2  ]) 
 xlim  ([  0  2  ]) 
 view  (  45  ,  22  ) 
 title  (  'Differenza'  ) 