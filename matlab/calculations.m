A = readmatrix("../csv/matrix.csv");
b = readmatrix("../csv/vector.csv");
x = readmatrix("../csv/solution.csv");
db = readmatrix("../csv/delta_vector.csv");
dx = readmatrix("../csv/delta_solution.csv");

minors = zeros([1, length(A)])
for i = 1:length(A)
    minors(i) = det(A(1:i, 1:i))
end

det(A)

real_x = A\b;
error_x = norm(x - real_x)
error_Ax_b = norm(A*x - b)

db_on_b = norm(db)/norm(b)
dx_on_x = norm(dx)/norm(x)
db_on_b_multiplied_by_condA = db_on_b * cond(A)
dx_on_x <= db_on_b_multiplied_by_condA
