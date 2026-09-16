# Export exact class-power data for all 26 sporadic simple groups.

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
    ["M11", "M11"], ["M12", "M12"], ["M22", "M22"],
    ["M23", "M23"], ["M24", "M24"], ["J1", "J1"],
    ["J2", "J2"], ["J3", "J3"], ["J4", "J4"],
    ["HS", "HS"], ["McL", "McL"], ["He", "He"],
    ["Ru", "Ru"], ["Suz", "Suz"], ["O'N", "ON"],
    ["Co1", "Co1"], ["Co2", "Co2"], ["Co3", "Co3"],
    ["Fi22", "Fi22"], ["Fi23", "Fi23"], ["Fi24'", "Fi24'"],
    ["HN", "HN"], ["Ly", "Ly"], ["Th", "Th"],
    ["Baby Monster", "B"], ["Monster", "M"]
];

Print("{\"schema\":\"gap-sporadic-power-maps-v2\",");
Print("\"software\":{\"gap\":\"", GAPInfo.Version, "\",");
Print("\"ctbllib\":\"", PackageInfo("ctbllib")[1].Version, "\"},");
Print("\"class_ordering\":\"CTblLib character-table class order\",");
Print("\"groups\":[");
for rowIndex in [1..Length(GroupRows)] do
    label := GroupRows[rowIndex][1];
    tableName := GroupRows[rowIndex][2];
    table := CharacterTable(tableName);
    if table = fail then
        Error("missing character table ", tableName);
    fi;
    exponent := Exponent(table);
    directions := PrimePowerDirections(exponent);
    if rowIndex > 1 then
        Print(",");
    fi;
    Print("{\"label\":\"", label, "\",");
    Print("\"table_name\":\"", tableName, "\",");
    Print("\"order\":", Size(table), ",");
    Print("\"class_count\":", NrConjugacyClasses(table), ",");
    Print("\"exponent\":", exponent, ",");
    Print("\"directions\":", directions, ",");
    Print("\"class_sizes\":", SizesConjugacyClasses(table), ",");
    Print("\"class_orders\":", OrdersClassRepresentatives(table), ",");
    Print("\"power_maps\":[");
    for directionIndex in [1..Length(directions)] do
        if directionIndex > 1 then
            Print(",");
        fi;
        Print(PowerMap(table, directions[directionIndex]));
    od;
    Print("]}");
od;
Print("]}\n");
QUIT_GAP(0);
