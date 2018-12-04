function e1 = get_geometric_error(H, y1, y2)
y2b = vgg_get_nonhomg(H*y1);
y1b = vgg_get_nonhomg(inv(H)*y2);
figure(1);plot(y1b(1,:),y1b(2,:),'rx');
figure(2);plot(y2b(1,:),y2b(2,:),'rx');

e1 = 0;
for k = 1:length(y1),
    e1 = e1 + norm(vgg_get_nonhomg(y2(:,k))-y2b(:,k))^2 + ...
    + norm(vgg_get_nonhomg(y1(:,k))-y1b(:,k))^2;
end
e1 = sqrt(e1);

end

