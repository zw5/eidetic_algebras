# Reconstruct the rational-class examples in the Families section.
SizeScreen([1000000,1000000]);;
LoadPackage("ctbllib");;
table:=fail;;
Print("{\"groups\":[");;
qs:=[3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32];;
for j in [1..Length(qs)] do
 q:=qs[j];;
 table:=CharacterTable(Concatenation("L2(",String(q),")"));;
 if table=fail then table:=CharacterTable(PSL(2,q)); fi;
 if j>1 then Print(","); fi;
 orders:=OrdersClassRepresentatives(table);;
 Print("{\"q\":",q,",\"class_orders\":",orders,",\"power_maps\":[");
 for n in [1..Maximum(orders)] do
  if n>1 then Print(","); fi;
  Print(List(PowerMap(table,n),x->x));
 od;
 Print("]}");
od;
Print("]}\n");;
QUIT_GAP(0);
