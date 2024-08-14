# Precondiciones: matriz cuadrada, diagonal distinto de cero, diagonal dominante

clear all

#--------------------------------------------------------------#

function [res] = grilla (num, t)
 if (num == 0)
    res = 0;
 else
    num_signo = sign(num);
    num_abs = abs(num);
    cor = log10(num_abs);
    cor = fix(cor) + 1;
    res = num_abs + 0.5 * 10^(cor-t);
    res = fix(res * 10^(t-cor));
    res = num_signo * res * 10^(cor-t);
 endif
endfunction
#--------------------------------------------------------------#

function [res] = matriz_aplicar_grilla (A, t)
 for i=1 :rows(A)
   for j=1 :columns(A)
     res(i,j) = grilla(A(i,j), t);
   endfor
 endfor
endfunction
#--------------------------------------------------------------#

function [res] = verificar_convergencia_GS (A)
  D = diag(diag(A));
  L = (-1) * tril(A, -1);
  U = (-1) * triu(A, 1);
  T = (D - L)^-1 * U;
  autovalores = eig(T);
  valores_abs_autovalores = autovalores;
  for i=1 : length(valores_abs_autovalores)
    valores_abs_autovalores(i) = abs(valores_abs_autovalores(i));
  endfor
  radio_espectral = max(valores_abs_autovalores);
  if (radio_espectral < 1)
    display("La matriz converge");
    res = true;
  else
    display("La matriz no converge");
    res = false;
  endif
endfunction
#--------------------------------------------------------------#

function solucion = gauss_seidel (A, b, semilla, cant_max_iter, delta_min, t)

  x = semilla;
  cant_iter = 1;
  continuar = true;
  n = rows(A);

  A = matriz_aplicar_grilla(A, t);
  b = matriz_aplicar_grilla(b, t);
  x = matriz_aplicar_grilla(x, t);

  while (continuar && cant_iter < cant_max_iter)
    x_anterior = x;
    for i=1 : n
      sigma = 0;
      for j=1 : i-1
        sigma = grilla(sigma + grilla(A(i,j) * x(j), t), t); # sigma = sigma + A(i,j) * x(j);
      endfor
      for j=i+1 : n
        sigma = grilla(sigma + grilla(A(i,j) * x_anterior(j), t), t); # sigma = sigma + A(i,j) * x_anterior(j);
      endfor
      x(i) = grilla(grilla(b(i)- sigma, t) / A(i,i), t); # x(i) = (b(i)- sigma) / A(i,i);
    endfor

    delta = norm(x - x_anterior);
    if (delta < delta_min)
      continuar = false;
    endif
    cant_iter = cant_iter + 1;
  endwhile

  solucion = matriz_aplicar_grilla(x, t);

  if (cant_iter >= cant_max_iter)
    disp("Se dejó de iterar porque se superó la cantidad máxima de iteraciones")
  endif
  if (!continuar)
    disp("Se dejó de iterar porque la diferencia entre iteraciones es menor que la mínima buscada")
  endif
  fprintf("Cantidad de iteraciones: %d\n", cant_iter);
  fprintf("La solucion del sistema es:\n");
  fprintf("%f\n", solucion);

endfunction
#--------------------------------------------------------------#

function comparar_solucion (A, b, solucion)
  x = A \ b;
  diferencia = x - solucion;
  fprintf("La diferencia entre mi solucion y la que me da Octave es:\n");
  fprintf("%f\n", diferencia);
endfunction
#--------------------------------------------------------------#

A = [3, 0, 2, 0; 0, 6, 0, 5; 0, 0, 2, 0; 7, 0, 0, 1]
b = [1; 1; 1; 1]
semilla = [0; 0; 0; 0];
cant_max_iter = 200;
delta_min = 0.5e-5;
t = 4;

if (!verificar_convergencia_GS(A))
  return
endif
solucion = gauss_seidel(A, b, semilla, cant_max_iter, delta_min, t);
comparar_solucion(matriz_aplicar_grilla(A, t), matriz_aplicar_grilla(b, t), solucion);
