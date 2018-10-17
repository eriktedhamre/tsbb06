function A = make_A_matrix(N, y1, y2)
A = [];
for cnt=1:N,
    coord1 = y1(:, cnt);
    coord2 = y2(:, cnt);
    A = [A;
        [coord1(1) 0 -coord1(1)*coord2(1)  coord1(2) 0 -coord1(2)*coord2(1) 1 0 -coord2(1)];
        [0  coord1(1)  -coord1(1)*coord2(2) 0 coord1(2) -coord1(2)*coord2(2) 0 1 -coord2(2)]];
end
end

