import csv

def csv_to_txt(csv_file):
    with open(csv_file, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        numbers = {"Pakistan": []}
        prov_names = []
        start_date, end_date = str(), str()

        num = 0
        for row in spamreader:
            if num == 0:
                prov_names = row[1:]
                for i in prov_names:
                    numbers[i] = []
            else:
                data = row[1:] #first two index have date
                if not start_date:
                    start_date = row[0]
                if not data:
                    print(f"File {csv_file} ended")
                    return
                for i in range(len(prov_names)):
                    numbers[prov_names[i]].append(- sum(numbers[prov_names[i]]) + int(data[i]))
                numbers["Pakistan"].append(- sum(numbers["Pakistan"]) + sum([int(i) for i in data]))

            end_date = row[0]
            num += 1



    for i in numbers:
        file = open(i+'.txt', "w")
        data = numbers[i]
        lines = [str(i + 1) + "\t" + str(data[i]) + '\n' for i in range(len(data))]
        file.write(f"{start_date}, {end_date}\n")
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
        end_date = str()

        for row in spamreader:
            if int(row[country_index + 1]) != 0 and not epidemic_start:
                start_date = row[0]
                num = 1
                epidemic_start = True
            if epidemic_start:
                #print(num, row[country_index + 1])
                numbers[country_name].append(- sum(numbers[country_name]) + int(row[country_index + 1]))
                num += 1
            end_date = row[0]

    #print(numbers)
    for i in numbers:
        file = open(i+'.txt', "w")
        data = numbers[i]
        lines = [str(i + 1) + "\t" + str(data[i]) + '\n' for i in range(len(data))]
        file.write(f"{start_date}, {end_date}\n")
        file.writelines(lines)

if __name__ == "__main__":
    csv_to_txt("pakistan_data.csv")
    #world_data_to_txt("China")
    #world_data_to_txt("Spain")
    #world_data_to_txt("United Kingdom")
    #world_data_to_txt("Italy")


