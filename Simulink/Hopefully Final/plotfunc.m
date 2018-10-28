function plot(Y)
    t2 = Y(1);
    t3 = Y(5);
    t4 = cos(t2);
    t5 = sin(t2);
    t6 = Y(2);
    t7 = t2+t6;
    t8 = cos(t7);
    t9 = sin(t7);
    t10 = Y(3);
    t11 = t2+t6+t10;
    t12 = cos(t11);
    t13 = Y(4);
    t14 = t2+t6+t10+t13;
    t15 = sin(t11);
    tester = reshape([0.0,t4,t4+t8,t4+t8+t12,t4+t8+t12+cos(t14),t3,t3+t5,t3+t5+t9,t3+t5+t9+t15,t3+t5+t9+t15+sin(t14)],[5,2]);
    P = plot(tester(:,1),tester(:,2));
end

