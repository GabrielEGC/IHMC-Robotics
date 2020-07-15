function R = eul2rotm1(eul)
ct = cos(eul);
st = sin(eul);
R = [  ct(2)*ct(1)   st(2)*st(3)*ct(1)-st(1)*ct(3)    st(2)*ct(3)*ct(1)+st(1)*st(3)
       ct(2)*st(1)   st(2)*st(3)*st(1)+ct(1)*ct(3)    st(2)*ct(3)*st(1)-ct(1)*st(3)
         -st(2)            ct(2)*st(3)             ct(2)*ct(3)];
end

