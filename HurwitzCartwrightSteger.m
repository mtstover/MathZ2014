/*

This program enumerates the smooth arithmetic ball quotients S with |Aut(S)| = 288e(S) and (topological) Euler characteristic e(S) <= e as a list called torsion_free.

One can run this program by typing `magma HurwitzCartwrightSteger.m' at the command line.

*/

G<j, u, v, b> := Group<j, u, v, b | u^4, v^8, (u,j), (v,j), j^-3*v^2, u*v*u*v^-1*u*v^-1, (b*j)^2*(v*u^2)^-1, (b,v*u^2), b^3, (b*v*u^3)^3>;

H := sub<G | v*u*b*j*u^-1, u^-1*j^-1*b*j^2, u^2*v*b*u*j^-2>;

e := 252;
N := 288*e;
norms := LowIndexNormalSubgroups(G, N); 
torsion_free := [* *]; 
K1:=u; 
K2:=u^2; 
K3:=u^3; 
K4:=v*u^2*v*j^9; 
K5:=u*v*u^2*v*j^9; 
K6:=u^2*v*u^2*v*j^9; 
K7:=u^3*v*u^2*v*j^9; 
K8:=v; 
K9:=u*v; 
K10:=u^2*v; 
K11:=u^3*v; 
K12:=u^2*v*u; 
K13:=u*v*u; 
K14:=u^3*v*u; 
K15:=v*u; 
K16:=u*v*u^2; 
K17:=u^2*v*u^2; 
K18:=u^3*v*u^2; 
K19:=v*u^2; 
K20:=u^2*v*u^3; 
K21:=u*v*u^3; 
K22:=u^3*v*u^3; 
K23:=v*u^3; 
KList:=[Id(G),K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11, 
	  K12,K13,K14,K15,K16,K17,K18,K19,K20,K21,K22,K23]; 
KTors := [* *]; 
for i1 in [1..24] do 
for i2 in [0..11] do 
	 Append(~KTors, KList[i1]*j^i2); 
end for; 
end for; 
TorsList := [K20*j^5*b^-1*K17, K20*j^5*b^-1*K5*j^4, 
		K20*j^5*b^-1*K21*j^2, K20*j^6*b^-1*K23*j^5, 
		K20*j^6*b^-1*K4*j^7, K20*j^6*b^-1*K16*j^3, 
		K18*j^7*b^-1*K12*j^2, K18*j^7*b^-1*K2*j^4,
		K18*j^7*b^-1*K9*j^3, K17*j^7*b^-1*K2*j^3, 
		K17*j^7*b^-1*K12*j, K17*j^7*b^-1*K9*j^2, 
		K18*j^6*b^-1*K8*j^3, K18*j^6*b^-1*K1*j^4, 
		K18*j^6*b^-1*K13*j^2, K17*j^8*b^-1*K17*j, 
		K17*j^8*b^-1*K5*j^5, K17*j^8*b^-1*K21*j^3, 
		K17*j^7*b^-1*K20*j^8*b^-1*K13*j^2, 
		K17*j^7*b^-1*K20*j^8*b^-1*K8*j^3, 
		K17*j^7*b^-1*K20*j^8*b^-1*K1*j^4, 
		K20*j^5*b^-1*K12*j, K18*j^7*b^-1*K21*j^3, 
		K20*j^5*b^-1*K21*j^2, K20*j^5*b^-1*K16*j^2, 
		K20*j^5*b^-1*K9*j^2, K20*j^5*b^-1*K13*j^2, 
		K20*j^6*b^-1*K16*j^3, K20*j^6*b^-1*K21*j^3, 
		K20*j^6*b^-1*K13*j^3, K20*j^6*b^-1*K9*j^3, 
		K20*j^5*b^-1*K21*j^3, K20*j^5*b^-1*K16*j^3, 
		K20*j^5*b^-1*K9*j^3, K20*j^5*b^-1*K13*j^3]; 
MasterTors := [* *]; 
for t in TorsList do Append(~MasterTors, t); end for; 
for i in [2..288] do Append(~MasterTors, KTors[i]); end for; 
for m in norms do
      M := m`Group; 
for g in MasterTors do 
if g in M then 
torsion_exists := true; 
break g; 
 else 
 torsion_exists := false; 
end if; 
end for;
if torsion_exists eq false then 
Append(~torsion_free, M); 
end if; 
end for;
