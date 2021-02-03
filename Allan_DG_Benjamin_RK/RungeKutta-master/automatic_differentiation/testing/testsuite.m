% This script tests the ado objects against a set of objective function


[x,dx] = autodiff(1:20,'ackley');
[x,dx] = autodiff(1:2,'boh1');
[x,dx] = autodiff(1:4,'colville');
[x,dx] = autodiff(1:2,'easom');
[x,dx] = autodiff(1:3,'hart3');
[x,dx] = autodiff(1:12,'mich');
[x,dx] = autodiff_sparse(1:10,'powell');
[x,dx] = autodiff(1:4,'rosen');
[x,dx] = autodiff(1:4,'schw');
