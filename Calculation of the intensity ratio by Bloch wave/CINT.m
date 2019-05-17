function faR=CINT(X,XX,F)%cubic interpolation

X1=XX-X(1);
X2=XX-X(2);
X3=XX-X(3);
X4=XX-X(4);
X12=X(1)-X(2);
X13=X(1)-X(3);
X14=X(1)-X(4);
X23=X(2)-X(3);
X24=X(2)-X(4);
X34=X(3)-X(4);

faR=F(1)*X2*X3*X4/(X12*X13*X14)-F(2)*X1*X3*X4/(X12*X23*X24)+F(3)*X1*X2*X4/(X13*X23*X34)-F(4)*X1*X2*X3/(X14*X24*X34);

end