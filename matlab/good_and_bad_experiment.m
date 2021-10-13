good_matrix = readmatrix('../csv/good_matrix.csv')
bad_matrix = readmatrix('../csv/bad_matrix.csv')
file_res = readmatrix('../csv/good_and_bad_res.csv')
db_on_b = file_res(:, 1)'
dx_on_x_good = file_res(:, 2)'
dx_on_x_bad = file_res(:, 3)'

figure
hold on
plot(db_on_b, dx_on_x_good, '*')
line([0, db_on_b(length(db_on_b))], [0, cond(good_matrix) * db_on_b(length(db_on_b))],...
     'LineStyle', '--')
%line([0, db_on_b(length(db_on_b))], [0, cond(bad_matrix) * db_on_b(length(db_on_b))],...
      %'Color', 'red', 'LineStyle', '--')
plot(db_on_b, dx_on_x_bad, 'r*')

grid on