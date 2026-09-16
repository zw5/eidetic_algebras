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

UnitRows := [[991, 661, 881, 1057, 1201], [991, 661, 881, 1057, 1201], [2311, 4621, 6161, 3697, 5281, 2521], [159391, 106261, 70841, 170017, 60721, 57961, 64681], [159391, 106261, 70841, 170017, 60721, 57961, 64681], [29261, 17557, 37621, 27931, 16171], [631, 421, 281, 337, 241], [29071, 58141, 25841, 23257, 41041, 61201], [531946220191, 151984634341, 202646179121, 364763122417, 521090174881, 386869978321, 581506427041, 251560774081, 549105775681, 394338510721, 226209688321], [2311, 4621, 6161, 3697, 5281, 2521], [20791, 13861, 15401, 22177, 23761, 2521], [10711, 7141, 9521, 2857, 6121, 8401], [554191, 158341, 211121, 380017, 542881, 146161, 218401], [270271, 180181, 320321, 216217, 51481, 196561, 277201], [1360591, 8163541, 3628241, 8707777, 3109921, 2968561, 4010161, 3160081], [10360351, 12432421, 14734721, 9945937, 4736161, 4520881, 8925841, 7927921], [159391, 956341, 708401, 1020097, 910801, 695521, 277201], [159391, 318781, 70841, 382537, 273241, 57961, 277201], [270271, 180181, 320321, 576577, 411841, 196561, 277201], [739728991, 211351141, 31311281, 507242737, 120772081, 153709921, 390186721, 696215521, 588107521], [9193774591, 6129183061, 12712379681, 9806692897, 3502390321, 11143969201, 3771804961, 5768642881, 19187007841, 6763236481], [658351, 1316701, 292601, 1685377, 1881001, 1915201, 2079001], [7988453551, 5325635701, 4733898401, 8947067977, 9129661201, 2904892201, 687178801, 3166594201, 635896801], [43415191, 28943461, 27871481, 46309537, 49617361, 31169881, 21326761, 29877121], [131643756193951, 204779176301701, 225365442702401, 28084001321377, 33433334906401, 127654551460801, 72010259798401, 96366671200801, 209398255466401, 81402902380801, 45296776324801, 224074478628001], [655680820994086837951, 1019947943768579525701, 863447994724723408001, 839271450872431152577, 333044226536679028801, 847748940275182982401, 537994519790019969601, 1097086863885530918401, 552152270310809968801, 202722572674500278401, 522534907152375717601, 376017675122056968001, 739195722313116868801, 471222150312535221601, 395137217924873424001, 755212964400074980801]];;
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
table := fail;;
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
    Print("\"class_names\":[");
    names := ClassNames(table);
    for j in [1..Length(names)] do
      if j>1 then Print(","); fi;
      Print("\"",names[j],"\"");
    od;
    Print("],\"unit_generators\":",UnitRows[rowIndex],",");
    Print("\"unit_power_maps\":[");
    for j in [1..Length(UnitRows[rowIndex])] do
      if j>1 then Print(","); fi;
      Print(List([1..NrConjugacyClasses(table)], i -> PowerMap(table,UnitRows[rowIndex][j] mod OrdersClassRepresentatives(table)[i])[i]));
    od;
    Print("],\"local_power_maps\":[");
    for j in [1..Maximum(OrdersClassRepresentatives(table))] do
      if j>1 then Print(","); fi;
      Print(List(PowerMap(table,j), x -> x));
    od;
    Print("],\"table_automorphisms\":",
      List(GeneratorsOfGroup(AutomorphismsOfTable(table)),
        a -> List([1..NrConjugacyClasses(table)], i -> i^a)),",");
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
