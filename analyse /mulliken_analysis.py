from itertools import islice

number = int(input("Enter the atom number:"))
lookup = "Atom number    " + str(number)
lookup_end = "Atom number    " + str(number + 1)

homo = int(input("Enter HOMO state:"))
lumo = int(input("Enter LUMO state:"))
homo_lookup = str(homo - 1) + "         "
lumo_lookup = str(lumo) + "         "

with open("Mulliken.out") as f:
    for num, line in enumerate(f, 1):
        if lookup in line:
         val_1 = num
        if lookup_end in line:
         val_2 = num

f = open('section.txt', 'w')

with open("Mulliken.out") as fin:
           lines = islice(fin, val_1, val_2)
           for line in lines:
               print(line)
               f.write(line + '\n')

with open("section.txt") as f:
    for num, line in enumerate(f, 1):
        if homo_lookup in line:
            val_3 = num
        if lumo_lookup in line:
            val_4 = num

with open("section.txt") as fin:
           lines = islice(fin, val_3, val_4)
           for line in lines:
                print(line)

