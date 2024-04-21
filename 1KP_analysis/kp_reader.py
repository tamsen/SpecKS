def get_species_code(Kp_out_file, sample_lookup):
    splat = Kp_out_file.split("_")
    if (len(splat)) < 4:
        return False

    species_code = splat[3][0:4]

    if species_code in sample_lookup:
        species = sample_lookup[species_code]
    else:
        species = "Unknown"

    return (species_code, species)

def get_Ks_from_file(Kp_out_file):

    Ks_values = []
    print("reading file.."+ Kp_out_file)
    with open(Kp_out_file, "r") as f:
        lines=f.readlines()

    for l in lines[1:len(lines)]:
        data=l.split('\t')
        if len(data) < 4:
            break
        if len(data[2]) < 2:
            print("strange line!")
            print("line: " + l)
            continue
        try:
            Ks=float(data[2])
            #print(str(Ks))
        except:
            print("parsing issue!")
            print("line: " + l)
            print(data[2])
            continue

        Ks_values.append(Ks)

    return Ks_values

def get_sample_code_lookup(lookup_file):

    look_up_dict = {}
    print("reading file.." + lookup_file)
    with open(lookup_file, "r") as f:
        lines = f.readlines()

    for l in lines[1:len(lines)]:
        #print(l)
        data = l.split(",")
        #print(str(data))
        #print(str('\n------\n'))
        if len(data) < 5:
            continue
        code=data[0]
        species_name=data[3]
        #sanitize.
        species_name=species_name.replace("/","_")
        look_up_dict[code] = species_name

    return look_up_dict
