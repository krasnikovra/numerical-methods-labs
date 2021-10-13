actual_solution = readmatrix('../csv/eye_cond_actual_solution.csv');
my_solutions = readmatrix('../csv/eye_cond_solutions.csv');
matrices = readmatrix('../csv/eye_cond_matrices.csv');
vectors = readmatrix('../csv/eye_cond_vectors.csv');

cond_nums = zeros([1, 10]);
errors = zeros([1, 10]);
errors_Ax_b = zeros([1, 10]);

for i=1:length(actual_solution)
    matrix = matrices(10*(i-1)+1:10*i, :);
    vector = vectors(i, :)';
    my_sol = my_solutions(i, :)';
    cond_nums(i) = cond(matrix);
    errors(i) = norm(my_sol - actual_solution);
    errors_Ax_b(i) = norm(matrix * my_sol - vector);
end

figure
loglog(cond_nums, errors, cond_nums, errors_Ax_b)
xlabel('Число обусловленности')
ylabel('Норма')
grid on
legend('факт. ошибка', 'невязка', 'Location', 'northwest')