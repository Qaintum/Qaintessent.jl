// quantum ripple-carry adder from Cuccaro et al, quant-ph/0410184
OPENQASM 2.0;
include "qelib1.inc";
gate majority a,b,c
{
  cx a,b;
  cx a,c;
  ccx a,c,b;
}
gate unmaj a,b,c
{
  h a;
  h b;
  h c;
}
gate test(var1,var2) a,e,d
{
  U(0, var1, var2) e;
  cx d,e;
  ccx a,d,e;
}
qreg cin[1];
qreg a[4];
qreg b[4];
qreg cout[1];
creg ans[4];
creg out[1];
// set input states
x a[0];
x b;
// add a to b, storing result in b
majority cin[0],b[0],a[0];
cx a[3],cout[0];
unmaj a[2],b[3],a[3];
test(0.1,0.3) a[3], b[0], a[1];

measure b -> ans;
measure a[0] -> out;
measure cin -> out[0];
