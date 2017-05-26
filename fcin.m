function dPdt=fcin(inp)

uv=inp(1:2);r=inp(3);psi=inp(4);
dXYdt=[cos(psi) -sin(psi);sin(psi) cos(psi)]*uv;
dpsidt=r;
dPdt=[dXYdt;dpsidt];
