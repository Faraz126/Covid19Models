import csv

def csv_to_txt(csv_file):
    with open(csv_file, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        numbers = {"total": []}
        prov_names = []

        num = 0
        for row in spamreader:
            if num == 0:
                prov_names = row[1:]
                for i in prov_names:
                    numbers[i] = []
            else:
                data = row[2:] #first two index have date
                for i in range(len(prov_names)):
                    numbers[prov_names[i]].append(data[i])
                numbers["total"].append(sum([int(i) for i in data]))
            num += 1

    for i in numbers:
        file = open(i+'.txt', "w")
        data = numbers[i]
        lines = [str(i + 1) + "\t" + str(data[i]) + '\n' for i in range(len(data))]
        file.write("day\tcases\n")
        file.writelines(lines)


csv_to_txt("data.csv")




            #print(', '.join(row))