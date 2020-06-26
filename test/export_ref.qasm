OPENQASM 2.0;
include "qelib1.inc";

//Defining registers
qreg q[4];

x q[2];
swap q[2],q[1];
ry(0.9424777960769379) q[0];
rz(2.4) q[1];
t q[1];
cx q[0],q[1];
crz(0.3) q[2],q[0];
