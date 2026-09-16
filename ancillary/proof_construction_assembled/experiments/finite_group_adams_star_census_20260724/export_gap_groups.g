# Exact conjugacy-class power-map export for the generalized Adams-star census.

SizeScreen([1000000, 1000000]);
LoadPackage("ctbllib");

PrimePowerDirections := function(exponent)
    local primes, output, prime, value;
    primes := Set(FactorsInt(exponent));
    output := [];
    for prime in primes do
        value := prime;
        while exponent mod value = 0 do
            Add(output, value);
            value := value * prime;
        od;
    od;
    return output;
end;

GroupRows := [
    rec(label := "C2", family := "cyclic", group := CyclicGroup(2)),
    rec(label := "C3", family := "cyclic", group := CyclicGroup(3)),
    rec(label := "C4", family := "cyclic", group := CyclicGroup(4)),
    rec(label := "C5", family := "cyclic", group := CyclicGroup(5)),
    rec(label := "C6", family := "cyclic", group := CyclicGroup(6)),
    rec(label := "C8", family := "cyclic", group := CyclicGroup(8)),
    rec(label := "C12", family := "cyclic", group := CyclicGroup(12)),
    rec(label := "C2xC2", family := "elementary_abelian",
        group := DirectProduct(CyclicGroup(2), CyclicGroup(2))),
    rec(label := "C2^3", family := "elementary_abelian",
        group := DirectProduct(CyclicGroup(2), CyclicGroup(2), CyclicGroup(2))),
    rec(label := "C3xC3", family := "elementary_abelian",
        group := DirectProduct(CyclicGroup(3), CyclicGroup(3))),
    rec(label := "D8", family := "dihedral", group := DihedralGroup(8)),
    rec(label := "D10", family := "dihedral", group := DihedralGroup(10)),
    rec(label := "D12", family := "dihedral", group := DihedralGroup(12)),
    rec(label := "D16", family := "dihedral", group := DihedralGroup(16)),
    rec(label := "Q8", family := "quaternion", group := QuaternionGroup(8)),
    rec(label := "S3", family := "symmetric", group := SymmetricGroup(3)),
    rec(label := "S4", family := "symmetric", group := SymmetricGroup(4)),
    rec(label := "S5", family := "symmetric", group := SymmetricGroup(5)),
    rec(label := "S6", family := "symmetric", group := SymmetricGroup(6)),
    rec(label := "S7", family := "symmetric", group := SymmetricGroup(7)),
    rec(label := "S8", family := "symmetric", group := SymmetricGroup(8)),
    rec(label := "A4", family := "alternating", group := AlternatingGroup(4)),
    rec(label := "A5", family := "alternating", group := AlternatingGroup(5)),
    rec(label := "A6", family := "alternating", group := AlternatingGroup(6)),
    rec(label := "A7", family := "alternating", group := AlternatingGroup(7)),
    rec(label := "A8", family := "alternating", group := AlternatingGroup(8)),
    rec(label := "PSL(2,7)", family := "lie_type_simple",
        group := PSL(2,7)),
    rec(label := "PSL(2,8)", family := "lie_type_simple",
        group := PSL(2,8)),
    rec(label := "PSL(2,11)", family := "lie_type_simple",
        group := PSL(2,11)),
    rec(label := "Heis27", family := "nilpotent_nonabelian",
        group := SmallGroup(27,3)),
    rec(label := "Frob20", family := "frobenius",
        group := SmallGroup(20,3))
];

Print("{\"schema\":\"gap-finite-group-power-maps-v2\",");
Print("\"software\":{\"gap\":\"", GAPInfo.Version, "\",");
Print("\"ctbllib\":\"", PackageInfo("ctbllib")[1].Version, "\"},");
Print("\"class_ordering\":\"GAP ConjugacyClasses order\",");
Print("\"groups\":[");
classes := [];; group := fail;; direction := 1;;
for row_index in [1..Length(GroupRows)] do
    row := GroupRows[row_index];
    group := row.group;
    classes := ConjugacyClasses(group);
    exponent := Exponent(group);
    directions := PrimePowerDirections(exponent);
    if row_index > 1 then
        Print(",");
    fi;
    Print("{\"label\":\"", row.label, "\",");
    Print("\"family\":\"", row.family, "\",");
    Print("\"order\":", Size(group), ",");
    Print("\"class_count\":", Length(classes), ",");
    Print("\"exponent\":", exponent, ",");
    Print("\"directions\":", directions, ",");
    Print("\"class_sizes\":", List(classes, Size), ",");
    Print("\"class_orders\":",
        List(classes, class -> Order(Representative(class))), ",");
    Print("\"power_maps\":[");
    for direction_index in [1..Length(directions)] do
        direction := directions[direction_index];
        if direction_index > 1 then
            Print(",");
        fi;
        Print(List(classes, class ->
            Position(classes,
                ConjugacyClass(group, Representative(class)^direction))));
    od;
    Print("]}");
od;
Print("]}\n");
QUIT_GAP(0);
