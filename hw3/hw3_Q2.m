% question 2 assistance
clc
clear all
syms xs xn ys yn yss ysn ynn xss xsn xnn

yns=ysn
xns=xsn

J=xs*yn-xn*ys

sx=yn/J
nx=-ys/J
sy=-xn/J
ny=xs/J

Js=xss*yn+xs*ysn-xsn*ys-xn*yss
Jn=xsn*yn+xs*ynn-xnn*ys-xn*ysn



sxs=(J*yns-Js*yn)/(J^2)
sxn=(J*ynn-Jn*yn)/(J^2)
sys=-(J*xns-Js*xn)/(J^2)
syn=-(J*xnn-Jn*xn)/(J^2)

nxs=-(J*yss-Js*ys)/(J^2)
nxn=-(J*ysn-Jn*ys)/(J^2)
nys=(J*xss-Js*xs)/(J^2)
nyn=(J*xsn-Jn*xs)/(J^2)


sxx=sx*sxs+nx*sxn
syy=sy*sys+ny*syn
nxx=sx*nxs+nx*nxn
nyy=sy*nys+ny*nyn

     RHS1=((sxx+syy)*xs+(nxx+nyy)*xn)*(-J^2)
simpleRHS1=simplify(RHS1)

RHS2=((sxx+syy)*ys+(nxx+nyy)*yn)*(-J^2)
simpleRHS2=simplify(RHS2)
