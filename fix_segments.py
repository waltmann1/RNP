name = "all.pdb"

f = open(name, "r")
new = "all_fixed.pdb"


def group_id(number):
    ids = ["A", "B", "C", "D", "E", "F", 'G', 'H', 'I', 'J', 'K', 'L' 'M',
           'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    if number < len(ids):
        return ids[number]
    else:
        first = ids[int(number % len(ids))]
        second = ids[int(number / len(ids))]
        return first + second


data = f.readlines()
for line in data:
    s = line.split()
    if len(s) > 5:
        if s[0] == "ATOM":
            chain_number = int(int(s[1]) / 3457)
            line=list(line)
            line[21] = group_id(chain_number)
            line = ''.join([x for x in line])
            print(line)

