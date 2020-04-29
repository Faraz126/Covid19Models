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
                    numbers[prov_names[i]].append(- sum(numbers[prov_names[i]]) + int(data[i]))
                numbers["total"].append(- sum(numbers["total"]) + sum([int(i) for i in data]))
            num += 1



    for i in numbers:
        file = open(i+'.txt', "w")
        data = numbers[i]
        lines = [str(i + 1) + "\t" + str(data[i]) + '\n' for i in range(len(data))]
        file.write("day\tcases\n")
        file.writelines(lines)


def world_data_to_txt(country_name):
    countires = ['China', 'US', 'United Kingdom', 'Italy', 'France', 'Germany', 'Spain', 'Iran']
    country_index=  countires.index(country_name)
    with open("worldData.csv", newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        spamreader.__next__() #skipping first line
        numbers = {country_name: []}
        epidemic_start = False
        start_date = str()

        for row in spamreader:
            if int(row[country_index + 1]) != 0 and not epidemic_start:
                start_date = row[0]
                num = 1
                epidemic_start = True
            if epidemic_start:
                #print(num, row[country_index + 1])
                numbers[country_name].append(- sum(numbers[country_name]) + int(row[country_index + 1]))
                num += 1

    #print(numbers)
    for i in numbers:
        file = open(i+'.txt', "w")
        data = numbers[i]
        lines = [str(i + 1) + "\t" + str(data[i]) + '\n' for i in range(len(data))]
        file.write("day\tcases started on " + start_date + "\n")
        file.writelines(lines)

#csv_to_txt("data2.csv")
#world_data_to_txt("China")
#world_data_to_txt("Spain")



            #print(', '.join(row))