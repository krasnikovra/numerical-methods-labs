hilbert_vectors_csv = readmatrix("../csv/hilbert_vectors.csv");
hilbert_solutions_csv = readmatrix("../csv/hilbert_solutions.csv");

n = hilbert_vectors_csv(:, 1)';

error = zeros([1, length(n)]);
error_Ax_b = zeros([1, length(n)]);

for i = 1:length(n)
    A = hilb(n(i));
    b = hilbert_vectors_csv(i, 2:(n(i)+1))';
    my_sol = hilbert_solutions_csv(i, 2:(n(i)+1))';
    act_sol = A\b;
    error(i) = norm(act_sol - my_sol);
    error_Ax_b(i) = norm(A * my_sol - b);
end

figure
semilogy(n, error)
grid on
xlabel("Размерность матрицы Гильберта")
ylabel("Ошибка ||x-x^*||")
title("||x-x^*||")

figure
semilogy(n, error_Ax_b)
grid on
xlabel("Размерность матрицы Гильберта")
ylabel("Невязка ||Ax-b||")
title("||Ax-b||")