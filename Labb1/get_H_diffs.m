function e1 = get_H_diffs(Ha, Hb, y1, y2)

y2c = Ha*y1;
y1c = inv(Ha)*y2;

y2b = vgg_get_nonhomg(Hb*y1);
y1b = vgg_get_nonhomg(inv(Hb)*y2);
%figure(1);plot(y1b(1,:),y1b(2,:),'rx');
%figure(2);plot(y2b(1,:),y2b(2,:),'rx');

figure(1);plot(y1b(1,:),y1b(2,:),'rx');
figure(1);plot(y1c(1,:),y1c(2,:),'rx');


e1 = 0;
for k = 1:length(y1),
    e1 = e1 + norm(vgg_get_nonhomg(y2c(:,k))-y2b(:,k))^2 + ...
    + norm(vgg_get_nonhomg(y1c(:,k))-y1b(:,k))^2;
end
e1 = sqrt(e1);


end

