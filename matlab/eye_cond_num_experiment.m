my_solutions = readmatrix('../csv/eye_cond_num_solutions.csv');
matrices = readmatrix('../csv/eye_cond_num_matrices.csv');
vector = readmatrix('../csv/eye_cond_num_vector.csv');

cond_nums = zeros([1, 10]);
errors = zeros([1, 10]);
errors_Ax_b = zeros([1, 10]);

for i=1:10
    matrix = matrices(10*(i-1)+1:10*i, :)
    my_sol = my_solutions(i, :)';
    actual_sol = matrix\vector;
    cond_nums(i) = cond(matrix);
    errors(i) = norm(my_sol - actual_sol);
    errors_Ax_b(i) = norm(matrix * my_sol - vector);
end

figure
loglog(cond_nums, errors)
title('Зависимость фактической ошибки от числа обусловленности')
xlabel('Число обусловленности')
ylabel('Ошибка ||x-x^*||')
grid on

figure
loglog(cond_nums, errors_Ax_b)
title('Зависимость невязки от числа обусловленности')
xlabel('Число обусловленности')
ylabel('Невязка ||Ax-b||')
grid on