function QQINT=QQINT(X,XX,F)
X1=XX-X(1);
X2=XX-X(2);
X3=XX-X(3);
X12=X(1)-X(2);
X13=X(1)-X(3);
X23=X(2)-X(3);

QQINT=F(1)*X2*X3/(X12*X13)-F(2)*X1*X3/(X12*X23)+F(3)*X1*X2/(X13*X23); % F size mismatch

end
