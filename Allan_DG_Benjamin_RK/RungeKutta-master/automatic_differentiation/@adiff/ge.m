function c = ge(a,b)
% GE implements greater or equal than, where either a or bor both is an adiff object.

switch [class(a),class(b)]
   
case 'adiffdouble'
   c = a.x>=b;
   
case 'doubleadiff'
   c = a>=b.x;
   
case 'adiffadiff'
   c = a.x>=b.x;
   
otherwise
   error(['Can''t do ',class(a),'>=',class(b)]);
   
end
