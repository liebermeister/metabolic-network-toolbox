function s = cs_print(kinetics,n,actions,network)

s = [sprintf('u: %f\n',kinetics.u); ...
     sprintf('c: %f\n',kinetics.c)];