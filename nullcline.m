% dK/dt = 0.004 + 0.07*K^2/(0.04 + K^2) - K/(1 + K + S)
% dS/dt = 0.82/(1 + (K/0.222)^5) - S/(1 + K + S)

n = 45
array_k = zeros(1, n)
array_s = zeros(1, n)
k = 0.01
for i = 1:n
    K = k * i
    syms S
    array_k(i) = solve(0.004 + 0.07*(K^2)/(0.04 + K^2) - K/(1 + K + S) == 0, S)
    array_s(i) = solve (0.82/(1 + (K/0.222)^5) - S/(1 + K + S) == 0, S)
end
x = k * (1 : n)
hold on
plot(x, array_s, 'g')
plot(x, array_k, 'b')
hold off

    
